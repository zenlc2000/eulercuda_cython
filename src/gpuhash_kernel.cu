

#include "common.h"



__global__ void phase1(	KEY_PTR  keys,
			unsigned int * offset,
			unsigned int length,
			unsigned int* count,
			unsigned int bucketCount)
{

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if(tid<length)
	{
		KEY_T key=keys[tid];
		unsigned int bucket=hash_h(key,bucketCount);
		offset[tid]=atomicInc (count+bucket,MAX_INT);

	}
	__syncthreads();
}


__global__ void copyToBucket(	KEY_PTR keys,
				VALUE_PTR values,
				unsigned int * offset,
				unsigned int length,
				unsigned int* start,
				unsigned int bucketCount,
				KEY_PTR  bufferK,
				VALUE_PTR bufferV)
{

	unsigned tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;

	if(tid<length)
	{
		KEY_T key =keys[tid];
		unsigned int bucket=hash_h(key,bucketCount);
		VALUE_T value=values[tid];
		unsigned int index=start[bucket]+offset[tid];
		//index=(index * BUCKET_ITEM_SIZE);
		bufferK[index]=key;
		bufferV[index]=value;
		//*(BUFFER_ITEM_KEY_PTR(buffer,index))=key;
		//*(BUFFER_ITEM_VALUE_PTR(buffer,index))=value;
	}
}
__global__ void bucketSort(KEY_PTR   bufferK,VALUE_PTR bufferV, unsigned int * start,unsigned int * bucketSize,unsigned int bucketCount,KEY_PTR TK,VALUE_PTR TV){


		__shared__ KEY_T keys[MAX_BUCKET_ITEM];
		unsigned int keyCount[MAX_BUCKET_ITEM/32];
		//unsigned int keyCount=0;
		unsigned int blockOffset=start[blockIdx.x];
		unsigned int size=bucketSize[blockIdx.x];

		unsigned int chunks=size>>5;
		chunks= (chunks<<5==size)?chunks:chunks+1;
		for(unsigned int j=0;j<chunks;j++){
			if((j<<5)+threadIdx.x<size)
				keys[(j<<5)+threadIdx.x]=bufferK[blockOffset+(j<<5)+threadIdx.x];//
		}

		__syncthreads();
		for(unsigned int j=0;j<chunks;j++){
			if((j<<5)+threadIdx.x<size){
				keyCount[j]=0;
				for(int i=0; i<size; i++){
					//if( keys[(i<<5)+threadIdx.x]> keys[i] ) keyCount++;
					keyCount[j]=( keys[(j<<5)+threadIdx.x]> keys[i] )?keyCount[j]+1:keyCount[j];
				}
			}
		}
			__syncthreads();
		for(unsigned int j=0;j<chunks;j++){
			if((j<<5)+threadIdx.x<size){
				TK[GET_KEY_INDEX(blockIdx.x,keyCount[j])]=keys[(j<<5)+threadIdx.x];
				TV[GET_VALUE_INDEX(blockIdx.x,keyCount[j])]=bufferV[blockOffset+(j<<5)+threadIdx.x];
			}
		}
}
