#ifndef COMPONENT_H
#define COMPONENT_H

extern "C" 
unsigned int * findComponent(Vertex * v,  /* in */
							 unsigned int * D, /*out*/
							 unsigned int length /*in*/);

extern "C"
unsigned int * findComponentDevice(Vertex * d_v,	/*in */
								   unsigned int ** d_D, /* out*/
								   unsigned int length	/* in */);

#endif //COMPONENT_H