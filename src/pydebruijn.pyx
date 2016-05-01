# debruijn.pyx

# Import low-level C declarations


cimport pydebruijn
from libc.stdlib cimport malloc, free
from cpython.pycapsule cimport *
from pydebruijn cimport constructDebruijnGraphDevice


# void constructDebruijnGraphDevice(
#       unsigned int * ecount,
# 		KEY_PTR d_lmerKeys,         //in lmer keys
# 		VALUE_PTR d_lmerValues,     //in lmer values
# 		unsigned int lmerCount,     //in total lmers
# 		KEY_PTR d_kmerKeys,         //in
# 		unsigned long kmerCount,    //in  total kmers
# 		unsigned int l,             //in k
# 		KEY_PTR d_TK,
# 		VALUE_PTR d_TV,
# 		unsigned int * d_bucketSeed,
# 		unsigned int bucketCount,
# 		EulerVertex ** d_ev,        //out
# 		unsigned int ** d_l,        //out
# 		unsigned int ** d_e,        //out
# 		EulerEdge ** d_ee           //out
# 		)

def construct_Debruijn_Graph_Device(
    KEY_T[:] d_lmerKeys,
    VALUE_T[:] d_lmerValues,
    unsigned int lmerCount,
    KEY_T[:] d_kmerKeys,
    unsigned long kmerCount,
    unsigned int l,
    KEY_T[:] d_TK,
    VALUE_T[:] d_TV,
    unsigned int[:] d_bucketSeed,
    unsigned int bucketCount):

    cdef:
      unsigned int ecount
      # int bucketCount = 0
     #   EulerVertex ** d_ev = <EulerVertex **>malloc(sizeof(EulerVertex))
      EulerVertex *d_ev
     #   unsigned int ** d_l = <unsigned int**>malloc(sizeof(unsigned int
      unsigned int *d_l
     #   unsigned int ** d_e	= <unsigned int**>malloc(sizeof(unsigned int))			
      unsigned int *d_e
      #  EulerEdge ** d_ee = <EulerEdge **>malloc(sizeof(EulerEdge))
      EulerEdge *d_ee
    
    return pydebruijn.constructDebruijnGraphDevice(
    & ecount,
    &d_lmerKeys[0],
    &d_lmerValues[0],
    lmerCount,
    &d_kmerKeys[0],
    kmerCount,
    l,
    &d_TK[0],
    &d_TV[0],
    &d_bucketSeed[0],
    bucketCount,
    &d_ev,
    &d_l,
    &d_e,
    &d_ee)


