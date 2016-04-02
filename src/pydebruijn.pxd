# cdebruijn.pxd 
#
# Declarations of "external" C functions and structures
#

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

cdef extern from "debruijn.h":
    void constructDebruijnGraphDevice(
            unsigned int *  ecount,
            KEY_PTR         d_lmerKeys,
            VALUE_PTR       d_lmerValues, 
            unsigned int    lmerCount,
            KEY_PTR         d_kmerKeys,  
            unsigned long   kmerCount,  
            unsigned int    l,
            KEY_PTR         d_TK,
            VALUE_PTR       d_TV, 
            unsigned int *  d_bucketSeed, 
            unsigned int    bucketCount,
            EulerVertex **  d_ev, 
            unsigned int ** d_l, 
            unsigned int ** d_e,
            EulerEdge **    d_ee)



