#ifndef EULERTOUR_H
#define EULERTOUR_H


//extern "C"
void findEulerDevice(EulerVertex * d_ev,
					 unsigned int * d_l, 
					 unsigned int * d_e, 
					 unsigned int vcount,
					 EulerEdge * d_ee,
					 unsigned int ecount,
					 CircuitEdge ** d_cg_edge, 
					 unsigned int * cg_edgeCount,
					 unsigned int * cg_vertexCount,
					 unsigned int kmerLength);
//extern "C"
void executeSwipeDevice(EulerVertex * d_ev,
						unsigned int * d_e, 
						unsigned int vcount, 
						EulerEdge * d_ee, 
						unsigned int ecount, 
						CircuitEdge * d_cg_edge,
						unsigned int cg_edgeCount , 
						unsigned int * d_tree,
						unsigned int treeCount);

//extern "C"
void markContigStart(EulerEdge * d_ee, 
					 unsigned char * d_contigStart, 
					 unsigned int ecount);
//extern "C"
void findEulerGold(EulerVertex * h_ev,
		unsigned int * h_l,
		unsigned int * h_e,
		unsigned int vcount,
		EulerEdge * h_ee,
		unsigned int ecount,
		unsigned int kmerLength);


#endif
