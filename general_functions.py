import numpy as np
import sys

__author__ = ' '.join(['Francho bauza <fbm.prof@hotmail.com>'])
__all__ = ['deep_in',
           'dict_to_numpy_array_listBased',
           'list_overlapp',
           'get_size']

def get_element(index):
    return lambda elem: elem[index]

def deep_in(obj,lst):
    flag=False
    for el in lst:
        if isinstance(el, list):
            flag=deep_in(obj,el)
        else:
            flag=obj in lst
        if flag:
            break
    return flag

def dict_to_numpy_array_listBased(d):    
    """Convert a dictionary of dictionaries to a numpy array
    with list."""
    import numpy
    s=[]
    s.extend(d.keys())
    for k, v in d.items():
        for key2 in v.keys():
            if key2 not in s: s.append(key2)
    n = len(s)
    a = numpy.zeros((n, n))
    for i in range(n):
        for j in range(n):
            try:
                a[i, j] = d[s[i]][s[j]]
            except KeyError:
                pass                
    return a,s

def list_overlapp(lista1,lista2):
    """ Comprueba si en dos lista hay algun termino comun
    """
    
    listaux=[aux in lista2 for aux in lista1]
    
    if True in listaux:
        return True
    else:
        return False

def get_size(obj, seen=None):
    """Recursively finds size of objects"""

    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()

    obj_id = id(obj)
    if obj_id in seen:
        return 0

    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)

    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, ( str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])

    return size