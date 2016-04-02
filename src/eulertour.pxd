#ifndef EULERTOUR_H
#define EULERTOUR_H

cdef extern from "common.h":
	#define KEY_TYPE_VALUE_TYPE
	ctypedef int KEY_T 
	ctypedef int * KEY_PTR 
	ctypedef int VALUE_T
	ctypedef int * VALUE_PTR 
	

cdef extern from "graph.h":
	
	cdef struct EulerEdge:
		KEY_T eid
		unsigned int v1
		unsigned int v2
		unsigned int s
		unsigned int pad
	
	cdef struct EulerVertex:
		KEY_T	vid
		int  ep
		int  ecount
		int  lp
		int  lcount
		
	cdef struct CircuitEdge:
		unsigned int ceid
		unsigned e1
		unsigned e2
		unsigned c1
		unsigned c2	
		
		

cdef void findEulerDevice(EulerVertex * d_ev,
					 unsigned int * d_l, 
					 unsigned int * d_e, 
					 unsigned int vcount,
					 EulerEdge * d_ee,
					 unsigned int ecount,
					 CircuitEdge ** d_cg_edge, 
					 unsigned int * cg_edgeCount,
					 unsigned int * cg_vertexCount,
					 unsigned int kmerLength)
cdef executeSwipeDevice(EulerVertex * d_ev,
						unsigned int * d_e, 
						unsigned int vcount, 
						EulerEdge * d_ee, 
						unsigned int ecount, 
						CircuitEdge * d_cg_edge,
						unsigned int cg_edgeCount , 
						unsigned int * d_tree,
						unsigned int treeCount)

cdef markContigStart(EulerEdge * d_ee, 
					 unsigned char * d_contigStart, 
					 unsigned int ecount)


