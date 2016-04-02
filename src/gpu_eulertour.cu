#define CUDPP_STATIC_LIB
#include <algorithm>
#include "graph.h"
//#include <cutil_inline.h>
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_cuda.h" 
#include <cuda.h>
#include <cuda_runtime.h>
#include "common.h"
#include "utils.h"
#include "cudpp.h"
#include "component.h"
#ifdef EULER_NDEBUG
#define DEBUG_EULER_CU(x)
#else
#define DEBUG_EULER_CU(x) x
#endif
//#define DEBUG_EULER_CU(x) x
#define DEBUG_CALL(x)  DEBUG_EULER_CU(x)

void printSuccessorGraph(Vertex * d_v , unsigned int length){
	
	Vertex * h_v =NULL;
	h_v=(Vertex * ) malloc(length* sizeof(Vertex));
	checkCudaErrors(cudaMemcpy(h_v,d_v,length*sizeof(Vertex),cudaMemcpyDeviceToHost));
	printf("$graph G {\n");
	for (unsigned int i =0; i< length; i++){	
		if(h_v[i].n1 < length)	printf("$\t%u -- %u\n",h_v[i].vid, h_v[i].n1);
		if(h_v[i].n2 < length)  printf("$\t%u -- %u\n",h_v[i].vid, h_v[i].n2);
	}
	printf("$}\n");
	free(h_v);
}
void printCircuitGraph(CircuitEdge * d_ce , unsigned int length){
	
	CircuitEdge * h_ce =NULL;
	h_ce=(CircuitEdge * ) malloc(length* sizeof(CircuitEdge));
	checkCudaErrors(cudaMemcpy(h_ce,d_ce,length*sizeof(CircuitEdge),cudaMemcpyDeviceToHost));
	printf("$graph G {\n");
	for (unsigned int i =0; i< length; i++){	
		printf("$\t%u -- %u [ label= e1:%u-e2:%u ]\n",h_ce[i].c1,h_ce[i].c2,h_ce[i].e1,h_ce[i].e2);
	}
	printf("$}\n");
	free(h_ce);
}
/*** Assig Successor**/

__global__  void assignSuccessor(EulerVertex * ev,unsigned int * l, unsigned int * e, unsigned vcount, EulerEdge * ee ,unsigned int ecount){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	unsigned int eidx=0;
	if(tid<vcount){		
		while(eidx<ev[tid].ecount && eidx<ev[tid].lcount){
			ee[e[ev[tid].ep+eidx]].s=l[ev[tid].lp+eidx] ;
			eidx++;
		}
	}
}
void validateSuccessors(EulerEdge * d_ee, unsigned int ecount) {
	EulerEdge * h_ee;

	h_ee= (EulerEdge * ) malloc( sizeof(EulerEdge) *ecount);
	checkCudaErrors(cudaMemcpy(h_ee,d_ee,ecount * sizeof(EulerEdge),cudaMemcpyDeviceToHost));
	
	unsigned int snot=0;
	for(unsigned int i =0;i< ecount; i++){
		if( h_ee[i].s==ecount) {snot++;}
	}
	printf("total edges with succesors not set :%u\n",snot);
	free(h_ee);
}

/** Constuct Succesor Graph**/ //Redundant
__global__ void constructSuccessorGraphP1(EulerEdge* e, Vertex * v, unsigned int ecount){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	
	if(tid<ecount){
		v[tid].n1=ecount;v[tid].n2=ecount;//v[tid].n3=ecount;v[tid].n4=ecount;
		v[tid].vid=e[tid].eid;
		v[tid].n1=e[tid].s;
	}
}

__global__ void constructSuccessorGraphP2(EulerEdge* e, Vertex * v, unsigned int ecount){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	
	if(tid<ecount){
	/*	if(v[v[tid].n1].n1 < ecount){
			v[v[tid].n1].n2=v[tid].vid;
		}else{
			v[v[tid].n1].n1=v[tid].vid;
		}*/
		if(v[tid].n1 <ecount ){
			v[v[tid].n1].n2=v[tid].vid;
		}
	}
}
/***   Calculate Circuit Graph Vertex  ***/
__global__ void calculateCircuitGraphVertexData( unsigned int * D,unsigned int * C,unsigned int ecount){
	
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid <ecount)
	{
		unsigned int c=D[tid];
		atomicExch(C+c,1);
	}	
}
/*** construct circuit graph vertex **/
__global__ void constructCircuitGraphVertex(unsigned int * C,unsigned int * offset,unsigned int ecount, unsigned int * cv, unsigned int cvCount){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if(tid < ecount){
		if(C[tid]!=0){
			cv[offset[tid]]=tid;
		}
	}
}

/*** Calculate Circuit Graph Edges***/
__global__ void calculateCircuitGraphEdgeData(EulerVertex* v,unsigned int * e,unsigned vCount,unsigned int * D,unsigned int * map,unsigned int ecount, unsigned int * cedgeCount/*, unsigned int cvCount*/){
	
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	unsigned int index=0;
	unsigned int maxIndex=0;
	index=0;
	maxIndex=0;
	if(tid<vCount && v[tid].ecount>0 ){
		index=v[tid].ep;
		maxIndex=index+v[tid].ecount-1;
		while (index < maxIndex ){
			unsigned int c1=map[D[e[index]]];
			unsigned int c2=map[D[e[index+1]]];
			if( c1 !=c2){
				unsigned int c=min(c1,c2);
				atomicInc(cedgeCount+c,ecount);
			}
	
			index++;
		}
	}

}
__global__ void assignCircuitGraphEdgeData(EulerVertex* v,
					   unsigned int * e,
					   unsigned vCount,
					   unsigned int * D,
					   unsigned int * map,
					   unsigned int ecount, 
					   unsigned int * cedgeOffset,
					   unsigned int * cedgeCount, 
					   unsigned int cvCount, 
					   CircuitEdge * cedge,  
					   unsigned int cecount){

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	unsigned int index=0;
	unsigned int maxIndex=0;
	if(tid<vCount && v[tid].ecount>0){
		index=v[tid].ep;
		maxIndex=index+v[tid].ecount-1;
		while (index<maxIndex   ){			
			unsigned int c1=map[D[e[index]]];
			unsigned int c2=map[D[e[index+1]]];
			if( c1 !=c2){
				unsigned int c=min(c1,c2);
				unsigned int t=max(c1,c2);
				unsigned int i=atomicDec(cedgeCount+c,ecount);
				i=i-1;
				cedge[cedgeOffset[c]+i].c1=c;
				cedge[cedgeOffset[c]+i].c2=t;
				cedge[cedgeOffset[c]+i].e1=e[index];
				cedge[cedgeOffset[c]+i].e2=e[index+1];
			}				
			index++;
		}
	}
}

/*
__global__ void markSegments(	unsigned short * d_mark,
				unsigned int 	circuitGraphEdgeCount,
				unsigned int * 	d_cg_edge_start,
				unsigned int *	d_cedgeCount,
				unsigned int 	circuitVertexSize){

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if(tid<circuitVertexSize){
		d_mark[ d_cg_edge_start[tid]]=d_cedgeCount[tid];
	}
	
}
*/
/*
__global__ void sortCircuit(	unsigned int cedgeOffset,
				unsigned int cedgeCount,
				unsigned int circuitCount,
				CircuitEdge * cedge ){

	unsigned int bid=0;
	unsigned int tid=0;
	
	unsigned int keyCount=0;
	unsigned int offset=cedgeOffset[bid];
	unsigned int itemCount=circuitCount[bid];
	unsigned int chunks=itemCount/256;
	unsigned int chunkSize=0;
	unsigned int chunkIdx=0;
	__shared__ unsigned int keys[256];

	if(bid<circuitCount){
		while(chunkIdx<Chunks){
			if(tid<	itemCount)
				keys[tid]=edge[offset+tid].e2;
			__syncthreads();
			if(tid<itemCount){
				for(int i=0;i<256;i++){
					if(keys[tid]>keys[i]) keyCount++;
				}
			}
			__syncthreads();
			CircuitEdge temp;
			if(tid<itemCount){
				temp=cedge[tid];			
			}
			__syncthreads();
		
		}
	
	}
}*/


__device__ unsigned int getValue(CircuitEdge cedge, unsigned char radix){

	switch(radix){
		case 0: return cedge.e2;
		case 1: return cedge.e1;
		case 2: return cedge.c2;
	}
	return 0xFFFFFFFF;
}

/*
__global__ void sortCircuitGraphEdgeData3( unsigned int * cedgeOffset,
					   unsigned int * cedgeCount, 
					   unsigned int circuitCount, 
					   CircuitEdge * cedge,
					   unsigned short * mark,
					   unsigned int edgeCount,
					   unsigned char radix){


	unsigned int chunks=blockDim.x;
	unsigned int chunkSize=cedgeCount[bid]/chunks; //fix off by 1
	unsigned int offset=cedgeOffset[bid]+chunkSize*threadIdx.x;

	//now scan
	while(mark[offset]==0 && offset<cedgeCount[bid]) offset++;

	//__syncthreads();

	//everyone looking at its own chunk and we have to sort (mark[Offset] sized data)
	unsigned int count=mark[Offset];
	for(int i=0;i<count;i++){
		minIndex=offset+i;
		minValue=getValue(cedge[minIndex],radix);// cedge[minIndex].c2;
		for( j=offset+i+1;j<offset+count;j++){
			unsigned int nextValue=getValue(cedge[j],radix);
			if( minValue > nextValue){
				minIndex=j;
				minValue=nextValue;
			}
		}
		if(minIndex != offset+i){
			CircuitEdge temp=cedge[offset+i];
			cedge[offset+i]=cedge[minIndex];
			cedge[minIndex]=temp;
		}
	}
	mark[Offset]=0;
	offset+=count;
	//scan onemore time to count same
	
	


	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	unsigned int minIndex=0;
	unsigned int minValue=0;
	unsigned int i =0;
	unsigned int j=0;
	unsigned int count;
	unsigned int offset=0;

	if(tid<circuitCount){
		count=cedgeCount[tid];
		offset=cedgeOffset[tid];
		for (i=0;i<count;i++){
			minIndex=offset+i;
			minValue=getValue(cedge[minIndex],radix);// cedge[minIndex].c2;
			for( j=offset+i+1;j<offset+count;j++){
				unsigned int nextValue=getValue(cedge[j],radix);
				if( minValue > nextValue){
					minIndex=j;
					minValue=nextValue;
				}
			}
			if(minIndex != offset+i){
				CircuitEdge temp=cedge[offset+i];
				cedge[offset+i]=cedge[minIndex];
				cedge[minIndex]=temp;
			}
		}
	}
	
	
}
*/

__global__ void sortCircuitGraphEdgeData2( unsigned int * cedgeOffset,
					   unsigned int * cedgeCount, 
					   unsigned int circuitCount, 
					   CircuitEdge * cedge,unsigned char radix){

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	unsigned int minIndex=0;
	unsigned int minValue=0;
	unsigned int i =0;
	unsigned int j=0;
	unsigned int count;
	unsigned int offset=0;

	if(tid<circuitCount){
		count=cedgeCount[tid];
		offset=cedgeOffset[tid];
		for (i=0;i<count;i++){
			minIndex=offset+i;
			minValue=getValue(cedge[minIndex],radix);// cedge[minIndex].c2;
			for( j=offset+i+1;j<offset+count;j++){
				unsigned int nextValue=getValue(cedge[j],radix);
				if( minValue > nextValue){
					minIndex=j;
					minValue=nextValue;
				}
			}
			if(minIndex != offset+i){
				CircuitEdge temp=cedge[offset+i];
				cedge[offset+i]=cedge[minIndex];
				cedge[minIndex]=temp;
			}
		}
	}/*
		}
	*/
}

__global__ void sortCircuitGraphEdgeData( unsigned int * cedgeOffset,
					   unsigned int * cedgeCount, 
					   unsigned int circuitCount, 
					   CircuitEdge * cedge){

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	unsigned int minIndex=0;
	unsigned int minValue=0;
	unsigned int i =0;
	unsigned int j=0;
	unsigned int count;
	unsigned int offset=0;

	if(tid<circuitCount){
		count=cedgeCount[tid];
		offset=cedgeOffset[tid];
		for (i=0;i<count;i++){
			minIndex=offset+i;
			minValue=cedge[minIndex].c2;
			for( j=offset+i+1;j<offset+count;j++){
				if(minValue > cedge[j].c2){
					minIndex=j;
					minValue=cedge[j].c2;
				}else if( minValue == cedge[j].c2){
					if( cedge[minIndex].e1> cedge[j].e1){
						minIndex=j;
						minValue=cedge[j].c2;
					}else if(cedge[minIndex].e1 == cedge[j].e1){
						if(cedge[minIndex].e2 > cedge[j].e2) {
							minIndex=j;
							minValue=cedge[j].c2;
						}
					}
				}
			}
			if(minIndex != offset+i){
				CircuitEdge temp=cedge[offset+i];
				cedge[offset+i]=cedge[minIndex];
				cedge[minIndex]=temp;
			}
		}
	}/*
	if(tid<vCount && v[tid].ecount>0){
		index=v[tid].ep;
		maxIndex=index+v[tid].ecount-1;
		while (index<maxIndex   ){			
			unsigned int c1=map[D[e[index]]];
			unsigned int c2=map[D[e[index+1]]];
			if( c1 !=c2){
				unsigned int c=min(c1,c2);
				unsigned int t=max(c1,c2);
				unsigned int i=atomicDec(cedgeCount+c,ecount);
				i=i-1;
				cedge[cedgeOffset[c]+i].c1=c;
				cedge[cedgeOffset[c]+i].c2=t;
				cedge[cedgeOffset[c]+i].e1=e[index];
				cedge[cedgeOffset[c]+i].e2=e[index+1];
			}				
			index++;
		}
	}*/
}
__global__  void identifyContigStart( EulerEdge * e ,unsigned char * contigStart,unsigned int ecount){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;	
	if(tid<ecount){
		if(e[tid].s < ecount){
			contigStart[e[tid].s]=0;
			//atomicExch(contigStart+e[tid].s,0);
		}
	}
}

__global__ void  markSpanningEulerEdges(EulerEdge * ee, unsigned int * mark , unsigned int ecount,CircuitEdge * cg_edge,unsigned int cg_edgeCount,unsigned int * tree, unsigned int treeCount){

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;	
	if(tid < treeCount) {
		/*if(tree[tid]==1)*/{
			atomicExch(mark+min(cg_edge[tree[tid]].e1,cg_edge[tree[tid]].e2),1); // important: assumption if(mark[i]=1) means mark[i]and mark[i+1] are swipe
			//atomicExch(mark+cg_edge[tree[tid]].e2,1);
			
		}
	}
}

__global__ void executeSwipe(EulerVertex * ev,unsigned int * e, unsigned int vcount , EulerEdge * ee, unsigned int * mark,unsigned int ecount){

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;	
	unsigned int t;
	unsigned int index=0;
	unsigned int maxIndex;
	unsigned int s;
	if( tid< vcount){
		index=ev[tid].ep;
		maxIndex=index+ev[tid].ecount-1;
		while( index<maxIndex){

			if(mark[ee[e[index]].eid]==1){
				t=index;
				s=ee[e[index]].s;
				while(mark[ee[e[index]].eid]==1 && index < maxIndex){					
					ee[e[index]].s=ee[e[index+1]].s;
					index=index+1;
				}
				if(t!=index){
					ee[e[index]].s=s;
				}
			}
			index++;
		}

	}
}

 void executeSwipeHost(EulerVertex * ev,unsigned int * e, unsigned int vcount , EulerEdge * ee, unsigned int * mark,unsigned int ecount, unsigned int tid){

//	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;	
	unsigned int t;
	unsigned int index=0;
	unsigned int maxIndex;
	unsigned int s;
	if( tid< vcount){
		index=ev[tid].ep;
		maxIndex=index+ev[tid].ecount-1;
		while( index<maxIndex){

			if(mark[ee[e[index]].eid]==1){
				t=index;
				s=ee[e[index]].s;
				while(mark[ee[e[index]].eid]==1 && index < maxIndex){					
					ee[e[index]].s=ee[e[index+1]].s;
					index=index+1;
				}
				if(t!=index){
					ee[e[index]].s=s;
				}
			}
			index++;
		}

	}
}
void executeSwipeHostLaunch(EulerVertex * d_ev, unsigned int * d_e, unsigned int vcount, EulerEdge * d_ee, unsigned int * d_mark , unsigned int ecount){

		EulerVertex *  h_ev;
		unsigned int * h_e;
		EulerEdge * h_ee;
		unsigned int * h_mark;
		
		logMessage(LOG_LVL_DETAIL,"executeSwipeHostLaunch");
		h_ev=(EulerVertex *)malloc(vcount*sizeof(EulerVertex));	
		h_e=(unsigned int *) malloc(vcount * sizeof(unsigned int ));
		h_ee =(EulerEdge *) malloc(ecount * sizeof(EulerEdge));
		h_mark=(unsigned int *) malloc(ecount * sizeof(EulerEdge));
		
		checkCudaErrors( cudaMemcpy(h_ev,d_ev, vcount*sizeof(EulerVertex), cudaMemcpyDeviceToHost));
		checkCudaErrors( cudaMemcpy(h_e,d_e, vcount*sizeof(unsigned int), cudaMemcpyDeviceToHost));
		checkCudaErrors( cudaMemcpy(h_ee,d_ee, ecount*sizeof(EulerEdge), cudaMemcpyDeviceToHost));
		checkCudaErrors( cudaMemcpy(h_mark,d_mark, ecount*sizeof(unsigned int), cudaMemcpyDeviceToHost));
		
		for(unsigned tid =0;tid<vcount;tid++){
			executeSwipeHost(h_ev,h_e,vcount,h_ee,h_mark,ecount,tid);
		}
		
		free(h_ev);
		free(h_e);
		free(h_ee);
		free(h_mark);
		
}


extern "C"
void	markContigStart(EulerEdge * d_ee, unsigned char * d_contigStart, unsigned int ecount){
	
	dim3 grid;
	dim3 block;
	cudaMemset(d_contigStart,1,ecount);
	getOptimalLaunchConfiguration(ecount,&grid,&block);
	identifyContigStart<<<grid,block>>>(d_ee,d_contigStart,ecount);
	cudaThreadSynchronize();
	CheckCUDAError();

}
extern "C" 
void executeSwipeDevice(EulerVertex * d_ev,unsigned int * d_e, unsigned int vcount, EulerEdge * d_ee, unsigned int ecount, CircuitEdge * d_cg_edge,unsigned int cg_edgeCount , unsigned int * d_tree,unsigned int treeCount){
	dim3 grid ;
	dim3 block;

	unsigned int * d_mark;
	allocateMemory((void**) &d_mark, ecount* sizeof(unsigned int));

	cudaMemset(d_mark,1,ecount* sizeof(unsigned int));
	getOptimalLaunchConfiguration(treeCount,&grid,&block);
	logMessage(LOG_LVL_DETAIL,"kernel: markSpanningEulerEdges");
	markSpanningEulerEdges<<<grid,block>>>(d_ee, d_mark , ecount,d_cg_edge,cg_edgeCount,d_tree, treeCount);
	cudaThreadSynchronize();
	CheckCUDAError();

	//DEBUG_CALL(executeSwipeHostLaunch(d_ev,d_e,vcount,d_ee,d_mark,ecount));
	getOptimalLaunchConfiguration(vcount,&grid,&block);
	logMessage(LOG_LVL_DETAIL,"kernel: executeSwipe");
	executeSwipe<<<grid,block>>>(d_ev,d_e,vcount , d_ee, d_mark,ecount);
	cudaThreadSynchronize();
	CheckCUDAError();
	
	//printData(d_ev,vcount,d_ee,ecount);
	deallocateMemory(d_mark);

}

/**ok ! this is not something pleasent to the eyes :-\*/
inline bool edgeComp(CircuitEdge a, CircuitEdge b){
	if(a.c1<b.c1) { return true;}
	else if(a.c1==b.c1){
		if(a.c2<b.c2){	return true;	} 
		else if (a.c2==b.c2){
			if(a.e1<b.e1) {return true;}
			else if(a.e1==b.e1){
				if(a.e2<b.e2) {return true;}
				else return false;
			} else return false;
		} else return false;		
	}else return false;
}
extern "C"
void findEulerDevice(EulerVertex * d_ev,unsigned int * d_l, unsigned int * d_e, unsigned int vcount,EulerEdge * d_ee,unsigned int ecount,CircuitEdge ** d_cg_edge, unsigned int * cg_edgeCount,unsigned int * cg_vertexCount, unsigned int kmerLength){


	Vertex * d_v=NULL;
	unsigned int * d_D;
	unsigned int * d_C;
	unsigned int * d_cg_offset;
	unsigned int * d_cedgeCount;
	unsigned int * d_cv;
	unsigned int * d_cg_edge_start;


	
	dim3 grid;
	dim3 block;
	
	allocateMemory((void**) &d_v, ecount* sizeof(Vertex));
		
	//step 1:
	// assign sucessors
	getOptimalLaunchConfiguration(vcount,&grid,&block);
	logMessage(LOG_LVL_DETAIL,"kernel: assignSuccessor");
	assignSuccessor<<<grid,block>>>(d_ev,d_l,d_e,vcount,d_ee,ecount);
	cudaThreadSynchronize();
	CheckCUDAError();
	
	//validateSuccessors(d_ee,ecount);
	
	//printDebruijnGraph(d_ev,vcount,d_l,d_e,d_ee,ecount,kmerLength,0);

	//step 2 successor graph
	//constructSuccessorGraph P1
	getOptimalLaunchConfiguration(ecount,&grid,&block);
	logMessage(LOG_LVL_DETAIL,"kernel: constructSuccessorGraph P1");
	constructSuccessorGraphP1<<<grid,block>>>(d_ee,d_v,ecount);
	cudaThreadSynchronize();
	CheckCUDAError();
	// printSuccessorGraph( d_v , ecount);
	
	/* synchronize */
	logMessage(LOG_LVL_DETAIL,"kernel: constructSuccessorGraph P2");
	constructSuccessorGraphP2<<<grid,block>>>(d_ee,d_v,ecount);
	cudaThreadSynchronize();
	CheckCUDAError();

	 //printSuccessorGraph( d_v , ecount);

	//step 3findComponent
	allocateMemory((void**) &d_D, ecount * sizeof(unsigned int));
	findComponentDevice(d_v,&d_D,ecount);
	

	//step 4 circuit graph construction
	//step 4.a  vertex calculation
	allocateMemory((void**) &d_C, ecount * sizeof(unsigned int));
	getOptimalLaunchConfiguration(ecount,&grid,&block);
	logMessage(LOG_LVL_DETAIL,"kernel: calculateCircuitGraphVertexData");
	calculateCircuitGraphVertexData<<<grid,block>>>( d_D,d_C,ecount);
	cudaThreadSynchronize();
	CheckCUDAError();
	//printData(d_C,ecount);

	//step 4.b offset calculation .find prefix sum 
	CUDPPConfiguration config;
	config.op = CUDPP_ADD;
	config.datatype = CUDPP_UINT;
    	config.algorithm = CUDPP_SCAN;
    	config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;
    
    	CUDPPHandle scanplan = 0;
    	CUDPPResult result = cudppPlan(&scanplan, config,ecount, 1, 0);	
	
	// Run the scan
	allocateMemory((void**) &d_cg_offset, ecount * sizeof(unsigned int));
    	cudppScan(scanplan, d_cg_offset, d_C, ecount);
	cudppDestroyPlan(scanplan);
	
	//printData(d_cg_offset,ecount);

	//step 4.c create circuitGraph
	unsigned int buffer[2];
	readData(buffer,d_cg_offset+ecount-1,1,sizeof(unsigned int));
	readData(buffer+1,d_C+ecount-1,1,sizeof(unsigned int));
	unsigned int circuitVertexSize=buffer[0]+buffer[1];
	*cg_vertexCount=circuitVertexSize;
	logMessage(LOG_LVL_MSG,"#Circuit Graph Vertex : %d",circuitVertexSize);
	allocateMemory( (void**) &d_cv, circuitVertexSize * sizeof(unsigned int));
	getOptimalLaunchConfiguration(ecount,&grid,&block);
	logMessage(LOG_LVL_DETAIL,"kernel: constructCircuitGraphVertex");
	constructCircuitGraphVertex<<<grid,block>>>(d_C,d_cg_offset,ecount, d_cv, circuitVertexSize);
	cudaThreadSynchronize();
	CheckCUDAError();
//	printData(d_cv,circuitVertexSize);

	if(circuitVertexSize>1){
		//step 4.d calculate edge information 
		allocateMemory((void**) &d_cedgeCount, circuitVertexSize * sizeof(unsigned int ));
		getOptimalLaunchConfiguration(vcount,&grid,&block);
		calculateCircuitGraphEdgeData<<<grid,block>>>(d_ev,d_e,vcount , d_D,d_cg_offset, ecount, d_cedgeCount/*, circuitVertexSize*/);
		cudaThreadSynchronize();
		CheckCUDAError();
		
		//printData(d_cedgeCount,circuitVertexSize);

		//step 4.e calculate edge offsets
		config.op = CUDPP_ADD;
		config.datatype = CUDPP_UINT;
		config.algorithm = CUDPP_SCAN;
		config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;    
		scanplan = 0;
		result = cudppPlan(&scanplan, config,ecount, 1, 0);		
		// Run the scan
		allocateMemory((void**) &d_cg_edge_start, circuitVertexSize * sizeof(unsigned int));
		cudppScan(scanplan, d_cg_edge_start, d_cedgeCount, circuitVertexSize);
		cudppDestroyPlan(scanplan);	
		//printData(d_cg_edge_start,circuitVertexSize);

		//step 4.f construct edges
		readData(buffer,d_cg_edge_start+circuitVertexSize-1,1,sizeof(unsigned int));
		readData(buffer+1,d_cedgeCount+circuitVertexSize-1,1,sizeof(unsigned int));
		unsigned int circuitGraphEdgeCount=buffer[0]+buffer[1];
		*cg_edgeCount=circuitGraphEdgeCount;
		logMessage(LOG_LVL_MSG,"#Circuit Graph Edges : %d\n",circuitGraphEdgeCount);

		
		allocateMemory((void**) d_cg_edge, circuitGraphEdgeCount * sizeof(CircuitEdge));
		//unsigned int * h_cedgeCount=NULL;
	//	h_cedgeCount = (unsigned int *) malloc(circuitVertexSize*sizeof(unsigned int));
	//	checkCudaErrors( cudaMemcpy(h_cedgeCount, d_cedgeCount, circuitVertexSize*sizeof(unsigned int), cudaMemcpyDeviceToHost));
		getOptimalLaunchConfiguration(vcount,&grid,&block);
		logMessage(LOG_LVL_DETAIL,"kernel: assignCircuitGraphEdgeData");
		assignCircuitGraphEdgeData<<<grid,block>>>(d_ev,d_e, vcount,d_D,d_cg_offset,ecount, d_cg_edge_start,d_cedgeCount, circuitVertexSize, *d_cg_edge, circuitGraphEdgeCount);
		cudaThreadSynchronize();
		CheckCUDAError();
		
	//	checkCudaErrors( cudaMemcpy(d_cedgeCount,h_cedgeCount, circuitVertexSize*sizeof(unsigned int), cudaMemcpyHostToDevice));
	//	free(h_cedgeCount);

		/**try1***/
/*		getOptimalLaunchConfigCustomized(circuitVertexSize,&grid,&block,1);
		for(unsigned char radix=0;radix<3;radix++){
			sortCircuitGraphEdgeData2<<<grid,block>>>(d_cg_edge_start,d_cedgeCount, circuitVertexSize, *d_cg_edge,radix);
			cudaThreadSynchronize();
			CheckCUDAError();
		}
*/		 
		/**try 2***/ 
/*		getOptimalLaunchConfigCustomized(circuitVertexSize,&grid,&block,1);
		unsigned short * d_mark;
		unsigned short * d_t2;
		allocateMemory((void**),d_mark,circuitGraphEdgeCount*sizeof(unsigned short));
		allocateMemory((void**),d_t1,circuitGraphEdgeCount*sizeof(unsigned short));
		getOptimalLaunchConfiguration(CircuitVertexSize,&grid,&block);
		markSegments<<<grid,block>>>(d_mark,circuitGraphEdgeCount,d_cg_edge_start,d_cedgeCount,circuitVertexSize);
		unsigned int thread=1;
		for(unsigned char radix=2;radix>=0;radix--){
			getOptimalLaunchConfigCustomized(circuitVertexSize,&grid,&block,threads);
			sortCircuitGraphEdgeData3<<<grid,block>>>(d_cg_edge_start,d_cedgeCount, circuitVertexSize, *d_cg_edge,d_mark,circuitGraphEdgeCount,radix);
			cudaThreadSynchronize();
			CheckCUDAError();
			thread++;
		}
		
		deallocateMemory(d_mark);
		deallocateMemory(d_t2);*/


		//*** try 3 cpu sorting**/
		CircuitEdge * h_cg_edge=(CircuitEdge *)malloc(circuitGraphEdgeCount*sizeof(CircuitEdge));
		checkCudaErrors( cudaMemcpy(h_cg_edge,*d_cg_edge, circuitGraphEdgeCount*sizeof(CircuitEdge), cudaMemcpyDeviceToHost));
		std::sort(h_cg_edge,h_cg_edge+circuitGraphEdgeCount,edgeComp);
		checkCudaErrors( cudaMemcpy(*d_cg_edge,h_cg_edge, circuitGraphEdgeCount*sizeof(CircuitEdge), cudaMemcpyHostToDevice));		
		free(h_cg_edge);	
	//	printCircuitGraph(*d_cg_edge,circuitGraphEdgeCount);	
		deallocateMemory(d_cg_edge_start);
		deallocateMemory(d_cedgeCount);
	}
	else{
		*cg_edgeCount=0;
	}
	//printData(*d_cg_edge,circuitGraphEdgeCount);
	
	/*****/
	//step 6 swipe execution

	
	deallocateMemory(d_cv);
	deallocateMemory(d_cg_offset);
	deallocateMemory(d_C);
	deallocateMemory(d_D);
	deallocateMemory(d_v);

	

	/*****/
	//calcuate contig start
	/***/

	
}
/*
extern "C"
void findEulerGold(EulerVertex * h_ev,
			unsigned int * h_l, 
			unsigned int * h_e, 
			unsigned int vcount,
			EulerEdge * h_ee,
			unsigned int ecount,
			unsigned int kmerLength){




	//find start vertices.
	unsigned char vertexMap;

	vertexMap=(unsigned char *) malloc(sizeof(unsigned char)*vcount);
	
	for(unsigned int i=0; i<vcount; i++){
		if(h_ev[i].ecount>h_ev[i].lcount){
		 	vertexMap[i]=1; //SINK
		}
		else if(h_ev[i].lcount>h_ev[i].ecount){
			vertexMap[i]=2; //SOURCE
		}
		else vertexMap[i]=0;
		// USED=100;
	
	}
	

	free(vertexMap);
		
}

*/
