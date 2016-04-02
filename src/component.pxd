# Vertex * v,  /* in */
# unsigned int * D, d_D /*out*/
# unsigned int length /*in*/

cdef extern from "graph.h":

	cdef struct Vertex:
		unsigned int vid
		unsigned int n1
		unsigned int n2


cdef unsigned int * findComponent(Vertex * v, 
							 unsigned int * D,
							 unsigned int length )


cdef unsigned int * findComponentDevice(Vertex * d_v,	
								   unsigned int ** d_D,
								   unsigned int length	)

#endif //COMPONENT_H