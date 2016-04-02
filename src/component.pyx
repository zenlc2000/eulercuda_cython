# component.pyx
#
# Import low-level C declarations
#
# Inputs: Vertex * v, unsigned int length
# Outputs: unsigned int * D / d_D
from libc.stdlib cimport malloc, free
from libc.stdint cimport uintptr_t
cimport component

def findComponent(Vertex[:] v, unsigned int length ):
    cdef:
        unsigned int *D = <unsigned int *>malloc(sizeof(unsigned int))
    result = component.findComponent(<Vertex *> &v[0], <unsigned int *> D, length)
    return  <uintptr_t>result


def findComponentDevice(Vertex[:] d_v, unsigned int length):
    cdef:
        unsigned int ** d_D = <unsigned int **>malloc(sizeof(unsigned int))
    return <uintptr_t>component.findComponentDevice(<Vertex *> &d_v[0], <unsigned int **> d_D, length)
