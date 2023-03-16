import networkx as nx #networkx library is neccesary for managing the networks
import networkxPlus as nxP #This library improves some networkx functions

import numpy as np

def nodes_list_from_dict(comm_dict,comms_join=(-1,-1)):
    list_out=[]
    
    if comms_join[0]==-1 or comms_join[1]==-1:
        for _,feats in comm_dict.items():
            list_out.append(feats[0])
    else:
        for comm_id,feats in comm_dict.items():
            if comm_id in comms_join:
                continue
            list_out.append((feats[0])[:])
        list_aux=(comm_dict[comms_join[0]][0])[:]
        list_aux.extend(comm_dict[comms_join[1]][0])
        list_out.append(list_aux)   
    
    return list_out

def deep_len(lst):
    return sum(deep_len(el) if isinstance(el, list) else 1 for el in lst)
     
    
        


comms_orig={}# In this dictionary is saved all the information of the communities

G=nx.Graph()

mod_growth_lim=0.0 #\Delta Q min

GN_evolution_filename="GN_evolution.txt" #Name of the filename with the GN algorithm evolution
#Each line of this file corresponds to a community of the GN final partition.
#Each community have 4 fields in the file separated by ;
    #The first one is the index of the community
    #The second one is the index of the parent community
    #The third one is the \Delta Q correponding to the "birth" of the community
    #The fourth are all the nodes that form the community, separated by ,
    
Network_filename="network_for_input.txt" #Name of the filename where are all the edges of the network.
#The network is required to calculate modularity during the TCM performance
#The format of this file is: node_origin,node_destination,edge_weight

out_modularity_filename="out_modularity.txt"
out_communities_filename="out_community.txt"



in_file=open(GN_evolution_filename,'r')
for line in in_file:
    line=line[:-1]
    line=line.split(";")
    comms_orig[int(line[0])]=[[int(aux) for aux in line[3].split(",")],int(line[1]),float(line[2]),False,False]
    
in_file.close()


nodo_for_index=[]
in_file=open(Network_filename,'r')
for line in in_file:
    line=line[:-1]
    line=line.split(",")
    
    node_or=int(line[0])
    node_dest=int(line[1])
    G.add_weighted_edges_from( [( node_or,node_dest,float(line[2]) )] )
    nodo_for_index.append(node_or)
    nodo_for_index.append(node_dest)
    
in_file.close()

N_nodos=len(nodo_for_index)
    
first_soft=True
second_soft=True



out_file=open(out_modularity_filename,'w')

comms = {key: [value[0][:],value[1],value[2],value[3],value[4]] for key, value in comms_orig.items()}
nodes_list=nodes_list_from_dict(comms)
N_comms_old=len(nodes_list)



comms = {key: [value[0][:],value[1],value[2],value[3],value[4]] for key, value in comms_orig.items()}


#### First filtering ####

ordered_keys=sorted(list(comms.keys()),reverse=True)


for keys in ordered_keys:
    feat=comms[keys]
    if feat[2]>mod_growth_lim :
        comms[keys][3]=True
        comms[ comms[keys][1] ][4]=True
    elif comms[keys][4] :
        comms[ comms[keys][1] ][4]=True
    elif not comms[ comms[keys][1] ][4] :
        comms[ comms[keys][1] ][0].extend(comms[keys][0])
        comms.pop(keys)       

nodes_list=nodes_list_from_dict(comms)
suma_Bools=[0,0,0]
not_relevant_comms=[]
rep_test=set()
for comm_id,nodes_feat in comms.items():
    rep_test.update(nodes_feat[0])
    
    if nodes_feat[3]:
        suma_Bools[0]+=1
    elif nodes_feat[4]:
        suma_Bools[1]+=1
    else:
        suma_Bools[2]+=1
        not_relevant_comms.append(comm_id)
        

#### second filtering ####

ordered_keys=sorted(list(comms.keys()),reverse=True)   

for keys in ordered_keys:
    feat=comms[keys]
    if not feat[3] and not feat[4] :
        ordered_keys_2=sorted(list(comms.keys()),reverse=True)
        
        ref_comm=feat[1]
        possible_union={ref_comm}
        for possible_key in ordered_keys_2:
            act_comm=possible_key
            
            update_list=[possible_key]
            flag=0
            while True:
                if comms[act_comm][1] in possible_union:
                    flag=1
                    break
                
                if comms[act_comm][1]==act_comm:
                    break
                    
                    
                update_list.append(act_comm)
                act_comm=comms[act_comm][1]
            possible_union.update(update_list)
        
        max_delta_mod=[0.0,-1]
        
        if keys in possible_union:
            possible_union.remove(keys)
            
        for possible_key in possible_union:
            nodes_list=nodes_list_from_dict(comms,comms_join=(keys,possible_key))

            mod_var=nxP.modularity(G,nodes_list)

            
            acum_len=sum([len(aux_nod_lis) for aux_nod_lis in nodes_list])
            if mod_var>max_delta_mod[0]:
                max_delta_mod=[mod_var,possible_key]
        
        if not max_delta_mod[1]==-1:
            comms[ max_delta_mod[1] ][0].extend(comms[keys][0])
            comms.pop(keys)

        

nodes_list=nodes_list_from_dict(comms)
suma_Bools=[0,0,0]
rep_test=set()
for comm_id,nodes_feat in comms.items():
    rep_test.update(nodes_feat[0])
    
    if nodes_feat[3]:
        suma_Bools[0]+=1
    elif nodes_feat[4]:
        suma_Bools[1]+=1
    else:
        suma_Bools[2]+=1

modula=nxP.modularity(G,nodes_list)
     

Nodes_list_Mut=np.zeros(N_nodos,dtype=int)
for label,aux in enumerate(nodes_list):
    for auxNod in aux:
        ind=nodo_for_index.index(auxNod)
        Nodes_list_Mut[ind]=label


print(mod_growth_lim,len(nodes_list),modula)
print(suma_Bools)
print(" ")    

out_file.write(str(mod_growth_lim)+","+str(modula)+","+str(len(nodes_list)))
out_file.close()




out_file=open(out_communities_filename,'w')
for i in range(len(nodes_list)):
    out_file.write(str(i)+";"+str(nodes_list[i][0]))
    for aux in nodes_list[i][1:]:
        out_file.write(","+str(aux))
    out_file.write("\n")
    

out_file.close()



