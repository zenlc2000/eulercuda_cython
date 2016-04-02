# debruijn.pyx

# Import low-level C declarations


cimport pydebruijn
from libc.stdlib cimport malloc, free

cdef extern from "common.h":
    #define KEY_TYPE_VALUE_TYPE
    ctypedef int KEY_T 
    ctypedef int * KEY_PTR 
    ctypedef int VALUE_T
    ctypedef int * VALUE_PTR 
    
    
cdef extern from "graph.h":

    ctypedef struct Vertex:
        pass
 #       int vid
 #       int n1
 #       int n2
        
    
    
    ctypedef struct EulerVertex:
        pass
#        KEY_T	vid
#        int  ep
#        int  ecount
#        int  lp
#        int  lcount
        
    
    ctypedef struct EulerEdge:
        pass
#        KEY_T eid
#        int v1
#        int v2
#        int s
#        int pad
    
    ctypedef struct CircuitEdge:
        pass
#        int ceid
#        int e1
#        int e2
#        int c1
#        int c2

cdef:
    unsigned int *ecount = <unsigned int *>malloc(sizeof(int))
    unsigned int bucketCount = 0
    EulerVertex ** d_ev = <EulerVertex **>malloc(sizeof(EulerVertex))
    unsigned int ** d_l = <unsigned int**>malloc(sizeof(unsigned int))
    unsigned int ** d_e	= <unsigned int**>malloc(sizeof(unsigned int))								
    EulerEdge ** d_ee = <EulerEdge **>malloc(sizeof(EulerEdge))

def constructDebruijnGraphDevice(int[:] d_lmerKeys, int[:] d_lmerValues, lmerCount, int[:] d_kmerKeys, kmerCount, l,int[:] d_TK,int[:] d_TV, int[:] d_bucketSeed):

    global ecount 
    global bucketCount
    global d_ev 
    global d_l 
    global d_e								
    global d_ee 
    
    return pydebruijn.constructDebruijnGraphDevice(
    <int *> &ecount[0],
    <KEY_PTR> &d_lmerKeys[0],
    <VALUE_PTR> &d_lmerValues[0],
    <int>lmerCount,
    <KEY_PTR> &d_kmerKeys[0],
    <int>kmerCount,
    <int>l,
    <KEY_PTR> &d_TK[0], 
    <VALUE_PTR> &d_TV[0],
    <int *> &d_bucketSeed[0],
    <int>bucketCount,
    <EulerVertex **> &d_ev[0],
    <int **> &d_l[0],
    <int **> &d_e[0], 
    <EulerEdge **> &d_ee[0])


#int * ecount,KEY_PTR d_lmerKeys,VALUE_PTR d_lmerValues, int lmerCount,KEY_PTR d_kmerKeys,  int kmerCount,  int l,
#KEY_PTR   d_TK,VALUE_PTR d_TV, int * d_bucketSeed, int bucketCount,EulerVertex ** d_ev, int ** d_l, int ** d_e,EulerEdge ** d_ee)





#d_lmerKeys,d_lmerValues,lmerCount,d_kmerKeys,kmerCount,l,d_TK, d_TV,* d_bucketSeed,bucketCount,** d_ev,** d_l,** d_e, ** d_ee,* ecount