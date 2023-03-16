# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 12:30:29 2019
"""
import time
import networkx as nx
import math
import numpy as np
import operator
import general_functions as GF
from networkx.algorithms.assortativity.mixing import degree_mixing_dict
from networkx.algorithms.shortest_paths.weighted import _weight_function
from heapq import heappush, heappop
from itertools import count
import random
__author__ = ' '.join(['Francho bauza <fbm.prof@hotmail.com>'])
__all__ = ['dist_link',
           'unweighted_func',
           'decodi_edge_extrem',
           'giant_component_size',
           'modularity',
           'premilimar_comm_partition',
           'dijkstra_withPaths',
           'all_pairs_dijkstra',
           'efficiency_weighted',
           'global_efficiency_weighted',
           'local_efficiency_weighted',
           'degree_assortativity_coefficient_weighted',
           'numeric_ac_ValueList',
           'degree_mixing_matrix_floatKeys',
           'wallace_index',
           'wallace_matrix',
           'joint_distribution_discrete']

def dist_link(u,v,att):
    wei=att['weight']
    if wei<1E-10:
        return 1E10
    else:
        return 1/att['weight']
    
def unweighted_func(u,v,att):
    return 1.0

def decodi_edge_extrem(code,n,type="triang_min"):
    coord=[]
    if type=="square":
        coord.append(code/n)
        coord.append(code%n)
        return coord
    elif type=="triang":
        if code>=n*(n+1)/2:
            return [n-1,n-1]
        a=int(math.ceil((-1+math.sqrt(1+4*(n*(n+1)-2*code)))/2))
        row=n-a
        col=n-(n*(n+1)/2-a*(a-1)/2-code)
        return [row,col]
    elif type=="triang_min":
        if code>=n*(n+-1)/2:
            return [n-2,n-1]
        a=int(math.floor((1+math.sqrt(1+4*(n*(n-1)-1-2*code)))/2))
        row=n-a-1
        col=n-(n*(n-1)/2-a*(a-1)/2-code)
        return [row,col]
    else:
        return [-1,-1]
    
def toy_graph_for_testing():
    Gaux=nx.Graph()
    Gaux.add_weighted_edges_from([(0,1,1.),(1,2,1.),(1,3,1.),(2,4,1.0),(3,4,1.0)])
    Gaux.add_weighted_edges_from([(4,5,1.0),(4,6,1.0),(6,7,1.),(6,8,1.)])
    Gaux.add_weighted_edges_from([(7,9,1.0),(8,9,1.0),(9,10,1.)])
    return Gaux


def giant_component_size(G):
    """
    Returns a 2-tuple list (N node, list of nodes) of net components
    
    """
    
    G_copy=G.copy()
    components=[]
    node_list=list(G_copy.nodes(data=False))
    if node_list==[]:
        return[(0,[])]
    
    prev_len=len(G_copy.nodes[node_list[0]])
    
    while(len(node_list)!=0):
        root_node=node_list[0]
        component_list=[]
        component_list.append(root_node)
        queue=[]
        queue.append(root_node)
        G_copy.nodes[root_node]["visited"]=True
        while(len(queue)>0):
            working_node=queue.pop(0)
            for n in G_copy.neighbors(working_node):
                if len(G_copy.nodes[n])==prev_len:
                    queue.append(n)
                    component_list.append(n)
                    G_copy.nodes[n]["visited"]=True
        components.append((len(component_list),component_list))
        for i in component_list:
            node_list.remove(i)
    components.sort(reverse=True)
    return components

def modularity(G,list_nodes,weight='weight',autoloop=False):
    
    KK=0.0
    KK_=0.0
    W=0.0
    sum_strn=0.0
    
    if len(list_nodes)<2 or len(list_nodes)>len(G)-1:
        return 0.0
    for nodes_line in list_nodes:
        KKaux=0.0
        check_nodes_line=set(nodes_line)
        for node in nodes_line:
            deg1=G.degree(node,weight=weight)
            KKaux+=deg1
            sum_strn+=deg1
            KK_+=deg1*deg1
            for node2,nodDict in G[node].items():
                if not node2 in check_nodes_line:
                    continue
                try:
                    aux=nodDict[weight]
                except KeyError:
                    aux=1.0
                W+=aux
        KK+=KKaux*KKaux
    
    if autoloop:
        Q=W/sum_strn-KK/(sum_strn*sum_strn)
    else:
        Q=W/sum_strn-KK/(sum_strn*sum_strn)+KK_/(sum_strn*sum_strn)
    
    return Q

def premilimar_comm_partition(G,control_inc=0.01,com_min=3,max_size=0,max_aver=4.0):
    mod_difs=[]
    mod_difs_dict={}
    sum_strength=0.0
    degreeDict={}
    
    if max_size==0:
        max_size=(int(max_aver*2))
        
    if max_size<max_aver:
        max_size=max_aver
    
    for node in list(G.nodes(data=False)):
        aux=G.degree(node,weight='weight')
        degreeDict[node]=aux
        sum_strength+=aux
        
    for node in list(G.nodes(data=False)):
        deg1=degreeDict[node]
        for node2,nodeDict in G[node].items():
            deg2=degreeDict[node2]
            W=nodeDict['weight']
            mod_difs.append((node,node2,W-(deg1*deg2/sum_strength)))
            if not node in mod_difs_dict:
                mod_difs_dict[node]={}
            mod_difs_dict[node][node2]=W-(deg1*deg2/sum_strength)
    
    mod_difs.sort(reverse=True,key=operator.itemgetter(2))
    
    time_0=time.time()
    G_aux=nx.Graph()
    
    control_value=0.0
    
    ## Pensar si la condicion del numero de nodos por comunidad la hacemos
    #   cogiendo el maximo o la media
    counter_aux=0
    for linkInf in mod_difs:
        if linkInf[2]<=0.0:
            break
        G_aux.add_edge(linkInf[0],linkInf[1])
        counter_aux+=1
        if float(counter_aux)/len(mod_difs) > control_value:
            control_value=float(counter_aux)/len(mod_difs)
            Gcc=giant_component_size(G_aux)
            N_com=len(Gcc)
            nodos_max=Gcc[0][0]
            nodos_aver=G_aux.number_of_nodes()/float(N_com)
            if (N_com>com_min and nodos_aver>max_aver) or nodos_max>max_size:
                break
            control_value+=control_inc
    
    print(N_com,G_aux.number_of_nodes(),nodos_aver,nodos_max,control_value,time.time()-time_0)
    print(" ")


def dijkstra_withPaths(G, source, weight, paths={},perc_nodes_forPaths=0.5, pred=None, cutoff=None):
    
    """
    Check networkx.algorithms.shortest_paths.weighted._dijkstra_multisource()

    """
    
    value_for_float_comparison=0.00000001
    ## Value for float comparison sirve para comparar longitudes de caminos
    # deberia ser un valor en relacion con las longitudes de los caminos
    
    if perc_nodes_forPaths<1.0:
        len_paths_max=int(len(G)*perc_nodes_forPaths)
    else:
        len_paths_max=int(perc_nodes_forPaths)
    
    G_succ = G._succ if G.is_directed() else G._adj

    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {}
    # fringe is heapq with 3-tuples (distance,c,node)
    # use the count c to avoid comparing nodes (may not be able to)
    c = count()
    fringe = []
    
    if source not in G:
        raise nx.NodeNotFound("Source {} not in G".format(source))
    
    ## Primer nodo y sus vecinos
    paths_source={source: [[[source]], 0.0] }
    seen[source] = 0
    dist[source]= 0.0
    
    for u, e in G_succ[source].items():
        cost = weight(source, u, e)
        if cost is None:
            continue
        vu_dist = cost
        if cutoff is not None:
            if vu_dist > cutoff:
                continue
        if u in dist:
            if vu_dist+value_for_float_comparison < dist[u]:
                raise ValueError('Contradictory paths found:',
                                 'negative weights?')
        elif u not in seen or vu_dist < seen[u]:
            seen[u] = vu_dist
            push(fringe, (vu_dist, next(c), u))
            paths_source[u] = [ [paths_source[source][0][0] + [u]], vu_dist]
            if pred is not None:
                pred[u] = [source]
        elif vu_dist == seen[u]:
            if pred is not None:
                pred[u].append(source)
                
    ## Let's work with the input paths to pre-update fringe, seen and dist
    # arrays.

    paths_len_list=[]
    for u, auxdat in paths.items():
        paths_len_list.append((u,-max([len(aux) for aux in auxdat[source][0]])))
    paths_len_list.sort(key=GF.get_element(1))


    main_count=0
    for u, _ in paths_len_list:
        if u in dist:
            continue
        
        main_count+=1
        if main_count>len_paths_max:
            break
        
        ## En esta primera comprobacion no hace falta comprobar multiple paths
        # suponemos que el dict paths, ya tiene todos los caminos multiples
        
        try:
            whole_paths_orig,whole_len=paths[u][source]
        except KeyError:
            continue
        
        whole_paths=[aux2[:] for aux2 in whole_paths_orig]
        ##Aqui deberia ir un bucle a todos los caminos multiples y repetir
        # el mismo procedimiento para todos los caminos
        
        for whole_path in whole_paths:
            
            
            whole_path_aux=whole_path[:-1]
            
            nodo_end_pre=''
            whole_len_aux=whole_len
            while whole_path_aux:
                
                path_to_add=[source] + list(reversed(whole_path_aux))
                nodo_end=whole_path_aux.pop(0)
                if not nodo_end in dist:
                    dist[nodo_end]=whole_len_aux
                    paths_source[nodo_end]=[[path_to_add],whole_len_aux]
                    
                else:
                    path_source_data=paths_source[nodo_end]
                    if path_to_add in path_source_data[0]:
                        break
                    elif not abs(whole_len_aux-path_source_data[1])<value_for_float_comparison:
                        raise ValueError('Contradictory paths found:',
                                             'negative weights?',
                                             whole_len_aux,
                                             path_source_data)
                    else:
                        path_source_data[0].append(path_to_add)                   
                
                cost_subs=0.0
                for nodo_end_neig, e in G_succ[nodo_end].items():
                    if nodo_end_neig==nodo_end_pre:
                        continue
                    
                    cost = weight(nodo_end, nodo_end_neig, e)
                    
                    
                    ## AQUI SE TIENE QUE PROBAR A PONER UN DEEP IN DE TODOS LOS
                    # CAMINOS MULTIPLES, NO SOLO DEL ACTUAL
                    
                    if nodo_end_neig in whole_path_aux:
                        try:
                            if nodo_end_neig==whole_path_aux[0]:
                                cost_subs=cost
                        except IndexError:
                            pass
                        continue
                    
                    if cost is None:
                        continue
                    vu_dist = whole_len_aux + cost
                    if cutoff is not None:
                        if vu_dist > cutoff:
                            continue
                    if nodo_end_neig in dist:
                        if vu_dist+value_for_float_comparison < dist[nodo_end_neig]:
                            raise ValueError('Contradictory paths found:',
                                             'negative weights?',
                                             vu_dist,
                                             dist[nodo_end_neig])
                    elif nodo_end_neig not in seen or vu_dist < seen[nodo_end_neig]:
                        seen[nodo_end_neig] = vu_dist
                        push(fringe, (vu_dist, next(c), nodo_end_neig))
                        
                        paths_source[nodo_end_neig] = [ [path_to_add + [nodo_end_neig]], vu_dist]
                        if pred is not None:
                            pred[nodo_end_neig] = [nodo_end]
                    elif vu_dist == seen[nodo_end_neig]:
                        paths_source[nodo_end_neig][0].append(path_to_add + [nodo_end_neig])

                        if pred is not None:
                            pred[nodo_end_neig].append(nodo_end)
                
                nodo_end_pre=nodo_end
                whole_len_aux-=cost_subs
                
                
    ## Esta es la parte antigua del algoritmo Dijkstra        
    
    while fringe:
        (d, _, v) = pop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        for u, e in G_succ[v].items():
            cost = weight(v, u, e)
            if cost is None:
                continue
            vu_dist = dist[v] + cost
            if cutoff is not None:
                if vu_dist > cutoff:
                    continue
            if u in dist:
                if vu_dist+value_for_float_comparison < dist[u]:
                    raise ValueError('Contradictory paths found:',
                                     'negative weights?',
                                     vu_dist,
                                     dist[u])
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if paths is not None:
                    paths_source[u] = [ [ aux_path + [u] for aux_path in paths_source[v][0] ] , vu_dist]
                if pred is not None:
                    pred[u] = [v]
            elif vu_dist == seen[u]:
                paths_source[u][0].extend([ aux_path + [u] for aux_path in paths_source[v][0] ])
                if pred is not None:
                    pred[u].append(v)


    ## AQUI TIENE QUE IR LA COMPROBACION DE SI NOS PASAMOS DE MEMORIA
    # PARA DECIDIR SI AÑADIMOS EL SOURCE_DICT AL PATHS_DICT
    
    paths[source]=paths_source

def all_pairs_dijkstra(G, cutoff=None, weight='weight', perc_nodes_forPaths=0.5):
    """
    Check networkx.algorithms.shortest_paths.weighted._dijkstra_multisource()

    """
    paths={}
    weight_func = _weight_function(G, weight)
    for n in G:
        dijkstra_withPaths(G, n, weight_func,paths,perc_nodes_forPaths=perc_nodes_forPaths ,cutoff=None)

    return paths

def efficiency_weighted(G,u,v=None,weight=None,path_length=False):
    if v==None:
        eff_dict={}
        path_dict=nx.shortest_path_length(G,u,v,weight)
        for k,i in path_dict.items():
            if u==k: continue
            eff_dict[k]=1.0/i
        if path_length:
            return eff_dict,path_dict
        else:
            return eff_dict
    else:
        try:
            path=nx.shortest_path_length(G,u,v,weight)
        except nx.NetworkXNoPath:
            eff=0.0
            path_length=False
        else:
            eff=1.0/path
    
        if path_length:
            return eff,path
        else:
            return eff

def global_efficiency_weighted(G, weight=None, proportion=1.0,path_length=False):
    
    nodes=list(G.nodes(data=False))
    N=G.number_of_nodes()
    control=N*(N-1)
    aver=[0.0,0.0]
    if control!=0:
        n_selected=int(N*proportion)
        selected_codes=random.sample(range(N), n_selected)
        for code in selected_codes:
            incremento=efficiency_weighted(G,nodes[code],weight=weight,path_length=path_length)
            if isinstance(incremento,tuple):
                incre_eff=0.0
                incre_path=0.0
                for _,eff in incremento[0].items():
                    incre_eff+=eff
                for _,path in incremento[1].items():
                    incre_path+=path
                aver[0]+=incre_eff/float(N-1)
                aver[1]+=incre_path/float(N-1)
            else:
                incre_eff=0.0
                for _,eff in incremento.items():
                    incre_eff+=eff
                aver[0]+=incre_eff/float(N-1)
        aver[0]/=float(n_selected)
        aver[1]/=float(n_selected)    
        if aver[1]==0.0:
            return aver[0]
        else:
            return aver
    else:
        if path_length:
            return aver
        else:
            return aver[0]
    # TODO This can be made more efficient by computing all pairs shortest
    # path lengths in parallel.
    #
    # TODO This summation can be trivially parallelized.
        
def global_efficiency_weighted_2(G, weight=None, proportion=1.0,path_length=False):
    nodes=list(G.nodes(data=False))
    N=G.number_of_nodes()
    aver=0.0
    for i in range(len(nodes)):
        for j in range(i):
            aver+=efficiency_weighted(G,nodes[i],nodes[j],weight=weight,path_length=path_length)
    aver/=float(N*(N-1))
    return aver


def local_efficiency_weighted(G, weight=None, proportion=1.0):
    if proportion<1.0:
        nodes=list(G.nodes(data=False))
        N=G.number_of_nodes()
        n_selected=int(N*proportion)
        if n_selected>0:
            selected_nodes=random.sample(range(N), n_selected)
            efficiency_list = (global_efficiency_weighted(G.subgraph(G[nodes[i]]),weight) for i in selected_nodes)
            return sum(efficiency_list) / float(n_selected)
        else:
            return 0.0
    else:
        efficiency_list = (global_efficiency_weighted(G.subgraph(G[v]),weight) for v in G)
        return sum(efficiency_list) / len(G)
    # TODO This can be made more efficient by computing all pairs shortest
    # path lengths in parallel.
    #
    # TODO This summation can be trivially parallelized.

def degree_assortativity_coefficient_weighted(G, weight, x='out', y='in',
                                     nodes=None):
    """
    For documentation check degree_assortativity_coefficient function from networkx module
    """
    M,ValVec = degree_mixing_matrix_floatKeys(G, x=x, y=y, nodes=nodes, weight=weight)
    return numeric_ac_ValueList(M,ValVec)

def numeric_ac_ValueList(M,ValVec):
    # M is a numpy matrix or array
    # numeric assortativity coefficient, pearsonr

    if M.sum() != 1.0:
        M = M / float(M.sum())
    x = np.array(ValVec)
    y = np.array(ValVec)
    a = M.sum(axis=0)
    b = M.sum(axis=1)
    vara = (a * x**2).sum() - ((a * x).sum())**2
    varb = (b * x**2).sum() - ((b * x).sum())**2
    xy = np.outer(x, y)
    ab = np.outer(a, b)
    return (xy * (M - ab)).sum() / np.sqrt(vara * varb)

def degree_mixing_matrix_floatKeys(G, x='out', y='in', weight=None,
                         nodes=None, normalized=True):
    """Return mixing matrix for attribute.

    Parameters
    ----------
    G : graph
       NetworkX graph object.

    x: string ('in','out')
       The degree type for source node (directed graphs only).

    y: string ('in','out')
       The degree type for target node (directed graphs only).

    nodes: list or iterable (optional)
        Build the matrix using only nodes in container.
        The default is all nodes.

    weight: string or None, optional (default=None)
       The edge attribute that holds the numerical value used
       as a weight.  If None, then each edge has weight 1.
       The degree is the sum of the edge weights adjacent to the node.

    normalized : bool (default=True)
       Return counts if False or probabilities if True.

    Returns
    -------
    m: numpy array
       Counts, or joint probability, of occurrence of node degree.
    """
    ValVec=[]
    d = degree_mixing_dict(G, x=x, y=y, nodes=nodes, weight=weight)
    a,ValVec = GF.dict_to_numpy_array_listBased(d)
    if normalized:
        a = a / float(a.sum())        
    return a,ValVec

def wallace_index(matrix,vector_1,vector_2):
    N=0
    P=0
    T=0
    Q=0
    
    
    
    for i in vector_1:
        N+=i
        Q+=i*(i-1)
    
    for i in vector_2:
        P+=i*(i-1)
        
    for i in matrix:
        for j in i:
            if j>0:
                T+=j*(j-1)
    
    U=N*(N-1)-P-Q+T
    
    a_ARI=T/2.
    b_ARI=(Q-T)/2.
    c_ARI=(P-T)/2.
    d_ARI=U/2.
    
    b_ARI2=Q/2.
    c_ARI2=P/2.
    
    
    R=float(U+T)/(N*(N-1))
    J=float(T)/(P+Q-T)
    k=len(vector_1)
    h=len(vector_2)
    minimo_T=0
    B=T/math.sqrt(P*Q)
    B_1=T/(Q+np.finfo('float64').eps)
    B_2=T/(P+np.finfo('float64').eps)
    if N>(len(vector_1)*len(vector_2)):
        minimo_T=k*h*((N/(k*h))*((N/(k*h))-1))+(N%(k*h))*2*(N/(k*h))
    try:
        B_norm=(T-minimo_T)/math.sqrt((P-minimo_T)*(Q-minimo_T))
    except ValueError:
        B_norm=0.0
    
    comb_2_ARI=N*(N-1)/2.
    
    exp_ARI=b_ARI2*c_ARI2/comb_2_ARI
    ARI=( a_ARI - exp_ARI )/( 0.5*(b_ARI2+c_ARI2) - exp_ARI )
    index_vector=[J,R,B,B_norm,B_1,B_2,ARI]
    return index_vector
    
def wallace_matrix(list_nodes_2,list_nodes):
    vector_row=[]
    vector_column=[]
    line2=[]
    copy_list_nodes=[]
    copy_list_nodes_2=[]
    #realizamos copias de las listas para no machacarlas
    for line in list_nodes:
        copy_list_nodes.append(set(line))
    for line in list_nodes_2:
        copy_list_nodes_2.append(set(line))
        
    for i in list_nodes:
        vector_column.append(len(i))
    matrix_M=np.zeros((len(list_nodes_2), len(list_nodes)),dtype=int)
    i=0
    for cluster_row in copy_list_nodes_2:
        vector_row.append(len(cluster_row))
        while(len(cluster_row)>0):
            node_row=cluster_row.pop()
            j=0
            for cluster_column in copy_list_nodes:
                if (node_row in cluster_column):
                    matrix_M[i][j]+=1
                    cluster_column.remove(node_row)
                    break
                j+=1
        i+=1
    return wallace_index(matrix_M,vector_row,vector_column)

def joint_distribution_discrete(vect1,vect2,normalized=True,returnDict=False):
    aux_dict={}
    
    norm=1.0
    if normalized:
        norm=float(len(vect1))
        
    for auxind in range(len(vect1)):
        aux_key=(vect1[auxind],vect2[auxind])
        if aux_key in aux_dict:
            aux_dict[aux_key]+=1/norm
        else:
            aux_dict[aux_key]=1/norm
    
    if returnDict:
     return aux_dict
    else:
        return list(aux_dict.values())