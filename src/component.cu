#include "Graph.h"
//#include <cutil_inline.h>
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_cuda.h"
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_timer.h"
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_functions.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "common.h"
#include "utils.h"
#include "stats.h"
/**  TODO
* Convert Vertex to EulerEdge
*
*/
__global__ void componentStepInit(Vertex * v, unsigned int * D,  unsigned int* Q, unsigned int length){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid <length){
		//v[tid].vid;
		D[tid]=tid;
		Q[tid]=0;
	}
}

__global__ void componentStepOne_ShortCuttingP1(Vertex * v, unsigned  int * prevD, unsigned  int * curD, unsigned int * Q, unsigned int length, int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid <length){
		curD[tid] =prevD[prevD[tid]];
	}
}
__global__ void componentStepOne_ShortCuttingP2(Vertex * v, unsigned  int * prevD, unsigned  int * curD, unsigned int * Q, unsigned int length, int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid <length){
		if(curD[tid]!=prevD[tid]){
			Q[curD[tid]]=s;
		}
	}
}
//for edge
__global__ void componentStepTwo(Vertex * v,  unsigned int * prevD,  unsigned int * curD, unsigned int * Q, unsigned int length, unsigned  int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	int a;
	int valIdx;
	int val;

	if( tid <length ){
		//it will done for each edge 1
		if(curD[tid] == prevD[tid] && v[tid].n1<length ){
			if(curD[v[tid].n1] < curD[tid]){
				a=curD[tid]; valIdx= v[tid].n1; val=curD[valIdx];
				__syncthreads();
				atomicMin(curD+a,val);
				atomicExch(Q+val,s);
			}
		}

		//it will done for each edge 2
		if(curD[tid] == prevD[tid] && v[tid].n2<length){
			if(curD[v[tid].n2] < curD[tid]){
				a=curD[tid]; valIdx= v[tid].n2; val=curD[valIdx];
				__syncthreads();
				atomicMin(curD+a,val);
				atomicExch(Q+val,s);
				
				
			}
		}

	}
}

__global__ void componentStepTwoP1(Vertex * v,  unsigned int * prevD,  unsigned int * curD, unsigned int * Q,unsigned int * t1,unsigned int *val1 ,unsigned int * t2,unsigned int * val2, unsigned int length, unsigned  int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	
	if( tid <length ){
		t1[tid]=length;t2[tid]=length;
		//it will done for each edge 1
		if(curD[tid] == prevD[tid] && v[tid].n1<length ){
			if(curD[v[tid].n1] < curD[tid]){
				t1[tid]=curD[tid]; 
				val1[tid]=curD[v[tid].n1];
				
			}
		}

		//it will done for each edge 2
		if(curD[tid] == prevD[tid] && v[tid].n2<length){
			if(curD[v[tid].n2] < curD[tid]){
				t2[tid]=curD[tid]; 
				val2[tid]=curD[v[tid].n2];				
			}
		}

	}
}
__global__ void componentStepTwoP2(Vertex * v,  unsigned int * prevD,  unsigned int * curD, unsigned int * Q,unsigned int * t1,unsigned int *val1 ,unsigned int * t2,unsigned int * val2, unsigned int length, unsigned  int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	
	int a;	
	int val;

	if( tid <length ){	
		//it will done for each edge 1
		if(t1[tid]<length){			
			a=t1[tid];
			val=val1[tid];				
			atomicMin(curD+a,val);
			atomicExch(Q+val,s);
			
		}

		//it will done for each edge 2
		if(t2[tid]<length){			
			a=t2[tid];
			val=val2[tid];				
			atomicMin(curD+a,val);
			atomicExch(Q+val,s);		
		}
	}
}
//for edge
__global__ void componentStepThree(Vertex * v, unsigned int * prevD,unsigned  int * curD,unsigned int * Q,unsigned int length,unsigned int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	int a;
	int valIdx;
	int val;
	if( tid< length) {
		//it will be done for each edge 1
		if(curD[tid]==curD[curD[tid]] && Q[curD[tid]] < s && v[tid].n1<length){
			if( curD[tid] != curD[v[tid].n1] ){
				a=curD[tid]; valIdx=v[tid].n1;val= curD[valIdx];
				__syncthreads();
				atomicMin(curD+a,val);
				//curD[curD[tid]]= curD[v[tid].n1];
			}
		}
		//it will be done for each edge 2
		if(curD[tid]==curD[curD[tid]] && Q[curD[tid]] < s && v[tid].n2<length){
			if( curD[tid] != curD[v[tid].n2] ){
				a=curD[tid]; valIdx=v[tid].n2;val= curD[valIdx];
				__syncthreads();
				atomicMin(curD+a,val);
				//curD[curD[tid]]= curD[v[tid].n2];
			}
		}
		
	}
}
__global__ void componentStepThreeP1(Vertex * v, unsigned int * prevD,unsigned  int * curD,unsigned int * Q,unsigned int * t1,unsigned int *val1 ,unsigned int * t2,unsigned int * val2,unsigned int length,unsigned int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid< length) {
		t1[tid]=length; t2[tid]=length;
		//it will be done for each edge 1
		if(curD[tid]==curD[curD[tid]] && Q[curD[tid]] < s && v[tid].n1<length){
			if( curD[tid] != curD[v[tid].n1] ){				
				t1[tid]=curD[tid];
				val1[tid]= curD[v[tid].n1];				
			}
		}
		//it will be done for each edge 2
		if(curD[tid]==curD[curD[tid]] && Q[curD[tid]] < s && v[tid].n2<length){
			if( curD[tid] != curD[v[tid].n2] ){
				t2[tid]=curD[tid];
				val2[tid]= curD[v[tid].n2];
			}
		}
		
	}
}

__global__ void componentStepThreeP2(Vertex * v, unsigned int * prevD,unsigned  int * curD,unsigned int * Q,unsigned int * t1,unsigned int *val1 ,unsigned int * t2,unsigned int * val2,unsigned int length,unsigned int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	int a;	
	int val;
	if( tid< length) {
		//it will be done for each edge 1
		if(t1[tid]<length){			
			a=t1[tid]; 
			val= val1[tid];			
			atomicMin(curD+a,val);
				
		}
		//it will be done for each edge 2
		if(t2[tid]<length){			
			a=t2[tid];
			val= val2[tid];			
			atomicMin(curD+a,val);
				
		}		
	}
}

__global__ void componentStepFour(Vertex * v, unsigned  int * curD,unsigned int length){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid < length){
		unsigned val=curD[curD[tid]];
		__syncthreads();
		curD[tid]= val;
	}
	//curD[tid]= curD[curD[tid]];
}
__global__ void componentStepFourP1(Vertex * v, unsigned  int * curD,unsigned int * val1,unsigned int length){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid < length){
		val1[tid]=curD[curD[tid]];
	}	
}
__global__ void componentStepFourP2(Vertex * v, unsigned  int * curD,unsigned int * val1,unsigned int length){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if( tid < length){
		curD[tid]= val1[tid];
	}
	
}
__global__ void componentStepFive(unsigned int * Q,unsigned int length,unsigned  int * sprimtemp,unsigned int s){
	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if(tid <length) {
		if(Q[tid]==s){		
			atomicExch(sprimtemp,1);
			//*sprime=*sprimtemp+1;
		}
	}
}

extern "C"
void findComponentDevice(Vertex *d_v,unsigned int ** d_D, unsigned int length){
	
	
	unsigned int * d_prevD;
	unsigned int * d_Q;	 
	unsigned int * d_t1;
	unsigned int * d_t2;
	unsigned int * d_val1;
	unsigned int * d_val2;
	unsigned int sp;
	unsigned int * sptemp;
	unsigned int * d_sptemp;
	unsigned int s;
	unsigned int * temp;
	StopWatchInterface *timer = NULL;

	dim3 grid;
	dim3 block;

//	cutilCheckError(cutCreateTimer(&timer));
//	cutilCheckError(cutStartTimer(timer));

	sdkCreateTimer(&timer);
	sdkResetTimer(&timer);
	sdkStartTimer(&timer);
	
	getOptimalLaunchConfiguration(length,&grid,&block);
	
	allocateMemory((void**) &d_Q, length* sizeof(int));
	allocateMemory((void**) &d_t1, length* sizeof(int));
	allocateMemory((void**) &d_t2, length* sizeof(int));
	allocateMemory((void**) &d_val1, length* sizeof(int));
	allocateMemory((void**) &d_val2, length* sizeof(int));
	allocateMemory((void**) &d_prevD, length* sizeof(int));
	cudaHostAlloc((void **)&sptemp, sizeof(int), cudaHostAllocMapped);	
	cudaHostGetDevicePointer((void **)&d_sptemp, (void *)sptemp, 0);	
	CheckCUDAError();
	//Initialize
	logMessage(LOG_LVL_DETAIL,"kernel: componentStepInit");
	componentStepInit<<<grid,block>>>(d_v,*d_D,d_Q,length);
	cudaThreadSynchronize();
	CheckCUDAError();

	s=1;
	sp=1;	
	while( s==sp)
	{		
		
		temp=*d_D;
		*d_D=d_prevD;
		d_prevD=temp;

		/**		componentStepOne_ShortCuttingP1		**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepOne_ShortCuttingP1");
		componentStepOne_ShortCuttingP1<<<grid,block>>>(d_v, d_prevD, *d_D, d_Q,length,s);
		cudaThreadSynchronize();		
		CheckCUDAError();

		/**		componentStepOne_ShortCuttingP2		**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepOne_ShortCuttingP2");
		componentStepOne_ShortCuttingP2<<<grid,block>>>(d_v, d_prevD, *d_D, d_Q,length,s);
		cudaThreadSynchronize();		
		CheckCUDAError();

		/**		componentStepTwoP1					**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepTwoP1");
		componentStepTwoP1<<<grid,block>>>(d_v, d_prevD, *d_D, d_Q,d_t1,d_val1,d_t2,d_val2,length,s);
		cudaThreadSynchronize(); 
		CheckCUDAError();

		/**		componentStepTwoP2					**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepTwoP2");
		componentStepTwoP2<<<grid,block>>>(d_v, d_prevD, *d_D, d_Q,d_t1,d_val1,d_t2,d_val2,length,s);
		cudaThreadSynchronize(); 
		CheckCUDAError();

		/**		componentStepThreeP1				**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepThreeP1");
		componentStepThreeP1<<<grid,block>>>(d_v, d_prevD, *d_D, d_Q,d_t1,d_val1,d_t2,d_val2,length,s);
		cudaThreadSynchronize(); 
		CheckCUDAError();

		/**		componentStepThreeP2			**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepThreeP2");
		componentStepThreeP2<<<grid,block>>>(d_v, d_prevD, *d_D, d_Q,d_t1,d_val1,d_t2,d_val2,length,s);
		cudaThreadSynchronize(); 
		CheckCUDAError();

		/**		componentStepFourP1				**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepFourP1");
		componentStepFourP1<<<grid,block>>>(d_v,  *d_D,d_val1, length);
		cudaThreadSynchronize(); 
		CheckCUDAError();
		
		/**		componentStepFourP2				**/
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepFourP2");
		componentStepFourP2<<<grid,block>>>(d_v,  *d_D,d_val1, length);
		cudaThreadSynchronize(); 
		CheckCUDAError();

		
		*sptemp=0;
		cudaThreadSynchronize(); 
		logMessage(LOG_LVL_DETAIL,"kernel: componentStepFive");
		componentStepFive<<<grid,block>>>(d_Q,length,d_sptemp,s);		
		cudaThreadSynchronize(); 
		CheckCUDAError();
		sp=sp+*sptemp;
		cudaThreadSynchronize(); 
		//printf("%d %d\n",s,*sptemp);
		s=s+1;
		//printData(*d_D,length);	
		
				
	}
//	printData(*d_D,length);
	deallocateMemory(d_t1);
	deallocateMemory(d_val1);
	deallocateMemory(d_t2);
	deallocateMemory(d_val2);
	deallocateMemory(d_Q);
	deallocateMemory(d_prevD);
	checkCudaErrors(cudaFreeHost(sptemp)); 	
	
//	cutilCheckError(cutStopTimer(timer));
//	setStatItem(TM_COMPONENT,cutGetTimerValue(timer));
//	cutilCheckError(cutDeleteTimer(timer));
	
	sdkStopTimer(&timer);
	float time = sdkGetTimerValue(&timer);
	setStatItem(TM_COMPONENT, time);
	sdkDeleteTimer(&timer);
}

extern "C"
void findComponent(Vertex *v,unsigned int * D, unsigned int length){
	
	unsigned int * d_curD;
	Vertex * d_v=NULL;

	checkCudaErrors( cudaMalloc( (void**) &d_v, length* sizeof(Vertex)) );	
	checkCudaErrors( cudaMalloc( (void**) &d_curD, length* sizeof(int) ));	
	
	checkCudaErrors( cudaMemcpy(d_v, v, length*(sizeof(Vertex)), cudaMemcpyHostToDevice) );
	findComponentDevice(d_v,&d_curD,length);

	checkCudaErrors( cudaMemcpy(D, d_curD, length*(sizeof(int)), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaFree(d_v) );	
	checkCudaErrors( cudaFree(d_curD) );	
}

