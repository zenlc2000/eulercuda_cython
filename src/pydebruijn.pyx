# debruijn.pyx

# Import low-level C declarations


cimport pydebruijn
from libc.stdlib cimport malloc, free
from cpython.pycapsule cimport *

    

def constructDebruijnGraphDevice(int[:] d_lmerKeys, int[:] d_lmerValues, lmerCount, int[:] d_kmerKeys, kmerCount, l,int[:] d_TK,int[:] d_TV, int[:] d_bucketSeed):

    cdef:
      int[:] ecount = None
      int bucketCount = 0
     #   EulerVertex ** d_ev = <EulerVertex **>malloc(sizeof(EulerVertex))
      EulerVertex[:,:] d_ev = None
     #   unsigned int ** d_l = <unsigned int**>malloc(sizeof(unsigned int
      int[:,:] d_l = None
     #   unsigned int ** d_e	= <unsigned int**>malloc(sizeof(unsigned int))			
      int[:,:] d_e = None
      #  EulerEdge ** d_ee = <EulerEdge **>malloc(sizeof(EulerEdge))
      EulerEdge[:,:] d_ee = None
    
    return pydebruijn.constructDebruijnGraphDevice(
    ecount,
    d_lmerKeys,
    d_lmerValues,
    lmerCount,
    d_kmerKeys,
    kmerCount,
    l,
    d_TK, 
    d_TV,
    d_bucketSeed,
    bucketCount,
    d_ev,
    d_l,
    d_e, 
    d_ee)


