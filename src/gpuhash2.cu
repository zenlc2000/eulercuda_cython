#define CUDPP_STATIC_LIB



//#include <cutil_inline.h>
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_cuda.h" // lib above replaced w/this one at CUDA 5.0
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_timer.h"
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_functions.h"
#include "cudpp.h"
#include "utils.h"
#include "gpuhash_device2.h"
#include <time.h>


#ifdef EULER_NDEBUG
#define DEBUG_GPUHASH2_CU(x)
#else
#define DEBUG_GPUHASH2_CU(x) x
#endif
#define DEBUG_CALL(x) DEBUG_GPUHASH2_CU(x)


__global__ void phase12(	KEY_PTR  keys,
			unsigned int * offset, 
			unsigned int length,
			unsigned int* count,
			unsigned int bucketCount){

	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if(tid<length){
		KEY_T key=keys[tid];
		unsigned int bucket=HASH_H(key,bucketCount);
		offset[tid]=atomicInc (count+bucket,MAX_INT);
		
	}
}
__global__ void copyToBucket2(	KEY_PTR keys,
				VALUE_PTR values,
				unsigned int * offset,
				unsigned int length,
				unsigned int* start,
				unsigned int bucketCount,
				KEY_PTR  bufferK,
				VALUE_PTR bufferV){

	unsigned tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;;

	if(tid<length){
		KEY_T key =keys[tid];
		unsigned int bucket=HASH_H(key,bucketCount);
		VALUE_T value=values[tid];
		unsigned int index=start[bucket]+offset[tid];
		//index=(index * BUCKET_ITEM_SIZE);
		bufferK[index]=key;
		bufferV[index]=value;
		//*(BUFFER_ITEM_KEY_PTR(buffer,index))=key;
		//*(BUFFER_ITEM_VALUE_PTR(buffer,index))=value;
	}
}

__global__ void phase22(KEY_PTR   bufferK,VALUE_PTR bufferV,
					unsigned int * start, unsigned int * count,
					unsigned int * bucketSeed,unsigned int bucketCount,
					KEY_PTR TK,VALUE_PTR TV,unsigned int * randomSeed,
					unsigned int seedCount){



		__shared__ unsigned int L[L2_ITEM_COUNT*3];
		__shared__ unsigned int pending;

		unsigned int blockCount=count[blockIdx.x];
		unsigned int key=bufferK[ ((start[blockIdx.x]+threadIdx.x)%blockCount)];
		unsigned int g[3];
		int tIdx=-1;
		bool unplaced=true;

		unsigned int iterations=0;
		unsigned int seedIdx=0;
		unsigned int seed=1;

		unsigned int id=threadIdx.x % blockCount;

		while ((unplaced) && seedIdx<seedCount){
			seed=*(randomSeed+blockIdx.x*seedCount+seedIdx);
			g[0]=HASH_G1(key,seed);
			g[1]=HASH_G2(key,seed);
			g[2]=HASH_G3(key,seed);
			tIdx=-1;
			do{
				atomicExch(&pending,0);
				__syncthreads();
				if(unplaced ) {
					tIdx= (tIdx+1)%3;
					atomicExch(L+tIdx*L2_ITEM_COUNT + g[tIdx],id);
				}
				__syncthreads();
				unplaced =(L[tIdx*L2_ITEM_COUNT + g[tIdx]] != id );
				if(unplaced )
					atomicExch(&pending,1);
				__syncthreads();
				iterations++;
			}while(pending && iterations<MAX_ITERATIONS);
			seedIdx++;
			__syncthreads();
		}
		__syncthreads();
		if(seedIdx>=seedCount){
					if (threadIdx.x==0) *(bucketSeed+blockIdx.x)=MAX_INT;
		}else {
			if (threadIdx.x==0) bucketSeed[blockIdx.x]=seed;
			g[0]=HASH_G1(key,seed);
			g[1]=HASH_G2(key,seed);
			g[2]=HASH_G3(key,seed);

			if(threadIdx.x<blockCount){
				TK[blockIdx.x * BLOCK_SIZE +tIdx*L2_ITEM_COUNT+g[tIdx]]=key;
				TV[blockIdx.x * BLOCK_SIZE +tIdx*L2_ITEM_COUNT+g[tIdx]]=bufferV[(start[blockIdx.x]+threadIdx.x)];
			}

		}

	}
unsigned int host_phase2(KEY_PTR   bufferK,VALUE_PTR bufferV,
		unsigned int * start, unsigned int * count,
		KEY_PTR TK,VALUE_PTR TV,
		unsigned int seedCount, unsigned int blockIdxx){

	unsigned int blockCount=count[blockIdxx];
	unsigned int blockStart=start[blockIdxx];
	unsigned int iterations=0;
	bool unplaced=true;
	unsigned int seed;
	unsigned int g[3];
	bool inserted=false;

	while(iterations<seedCount && unplaced){
		seed=rand();
		for(unsigned int j=0;j<BLOCK_SIZE;j++){	//set every entry of table to NULL
				TK[blockIdxx*BLOCK_SIZE+j]=MAX_INT;
				TV[blockIdxx*BLOCK_SIZE+j]=MAX_INT;
			}
		for(unsigned int i=0;i<blockCount;i++){
			KEY_T key=bufferK[blockStart+i];
			VALUE_T value=bufferV[blockStart+i];
			//insert
			inserted=false;
			unsigned int tries=0;
			unsigned int gidx=0;
			KEY_T tempKey;
			VALUE_T tempValue;
			while(!inserted && tries<MAX_ITERATIONS){
				g[0]=HASH_G1(key,seed);
				g[1]=HASH_G2(key,seed);
				g[2]=HASH_G3(key,seed);

				tempKey=TK[blockIdxx*BLOCK_SIZE+L2_ITEM_COUNT*gidx+g[gidx]];
				tempValue=TV[blockIdxx*BLOCK_SIZE+L2_ITEM_COUNT*gidx+g[gidx]];
				TK[blockIdxx*BLOCK_SIZE+L2_ITEM_COUNT*gidx+g[gidx]]=key;
				TV[blockIdxx*BLOCK_SIZE+L2_ITEM_COUNT*gidx+g[gidx]]=value;
				if(tempKey==MAX_INT && tempKey==MAX_INT){
						//empty

					inserted=true;
				}else {
					key=tempKey;
					value=tempValue;
					gidx=((gidx+1)%3);
				}
				tries ++;
			}
			if(!inserted) break;
		}
		if(inserted) unplaced=false;
		iterations++;
	}
	if(iterations>=seedCount)
		return MAX_INT;
	else return seed;

}
void phase2HostLaunch(KEY_PTR   bufferK,VALUE_PTR bufferV,
					unsigned int * start, unsigned int * count,
					unsigned int * bucketSeed,unsigned int bucketCount,
					KEY_PTR TK,VALUE_PTR TV,
					unsigned int seedCount,unsigned int length){


	KEY_PTR h_bufferK=(KEY_PTR)malloc(KEY_SIZE*length);
	VALUE_PTR h_bufferV=(VALUE_PTR)malloc(VALUE_SIZE*length);
	KEY_PTR h_TK=(KEY_PTR)malloc(KEY_SIZE*BLOCK_SIZE*bucketCount);
	VALUE_PTR h_TV=(VALUE_PTR)malloc(VALUE_SIZE*BLOCK_SIZE*bucketCount);

	unsigned int * h_start=(unsigned int*)malloc(sizeof(unsigned int)* bucketCount);
	unsigned int * h_count=(unsigned int*)malloc(sizeof(unsigned int)* bucketCount);
	unsigned int * h_bucketSeed=(unsigned int*)malloc(sizeof(unsigned int)* bucketCount);


	checkCudaErrors( cudaMemcpy(h_bufferK, bufferK, length * (KEY_SIZE), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_bufferV, bufferV, length * (VALUE_SIZE), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_start, start,sizeof(unsigned int)* bucketCount, cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_count, count, sizeof(unsigned int)* bucketCount, cudaMemcpyDeviceToHost));

	for(unsigned int blockIdx=0;blockIdx<bucketCount;blockIdx++){

		h_bucketSeed[blockIdx]=host_phase2(h_bufferK,h_bufferV,h_start,h_count,h_TK,h_TV,seedCount,blockIdx);
	}

	checkCudaErrors( cudaMemcpy( bucketSeed, h_bucketSeed, sizeof(unsigned int)* bucketCount,cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy( TK, h_TK, KEY_SIZE*BLOCK_SIZE*bucketCount,cudaMemcpyHostToDevice) );
	checkCudaErrors( cudaMemcpy( TV, h_TV, VALUE_SIZE*BLOCK_SIZE*bucketCount,cudaMemcpyHostToDevice) );


	free(h_bufferK);
	free(h_bufferV);
	free(h_TK);
	free(h_TV);
	free(h_start);
	free(h_count);
	free(h_bucketSeed);

}
void verifyBuffer2(KEY_PTR d_keys,VALUE_PTR d_values,unsigned int length, KEY_PTR d_bufferK, VALUE_PTR d_bufferV, unsigned int * d_bucketOffset,unsigned int * d_bucketSize, unsigned int bucketCount){

	KEY_PTR h_keys;
	VALUE_PTR h_values;
	KEY_PTR	h_bufferK;
	VALUE_PTR h_bufferV;
	unsigned int * h_bucketOffset;
	unsigned int * h_bucketSize;
	
	
	h_keys=(KEY_PTR) malloc(length * KEY_SIZE);
	h_values=(VALUE_PTR) malloc( length * VALUE_SIZE);
	h_bufferK= (KEY_PTR) malloc ( length *KEY_SIZE);
	h_bufferV= (VALUE_PTR) malloc( length * VALUE_SIZE);
	h_bucketOffset =(unsigned int *) malloc( bucketCount * sizeof(unsigned int));
	h_bucketSize = (unsigned int * ) malloc ( bucketCount * sizeof(unsigned int));	
	
	checkCudaErrors( cudaMemcpy(h_keys, d_keys, length * (KEY_SIZE), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_values, d_values, length * (VALUE_SIZE), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_bufferK, d_bufferK, length * (KEY_SIZE), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_bufferV, d_bufferV, length * (VALUE_SIZE), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_bucketOffset, d_bucketOffset, bucketCount* (sizeof(unsigned int)), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_bucketSize, d_bucketSize, bucketCount * (sizeof(unsigned int)), cudaMemcpyDeviceToHost));

	unsigned int found=0;
	unsigned int notfound=0;
	unsigned int correctValue=0;
	unsigned int incorrectValue=0;

	for (unsigned int i=0; i<length ; i++){
		unsigned int j=0;
		unsigned int bucket= HASH_H(h_keys[i],bucketCount);
		unsigned int offset=h_bucketOffset[bucket];
		unsigned int size = h_bucketSize[bucket];
	
		while( j<size && h_bufferK[offset+j]!=h_keys[i]) j++;
		if( j<size) {
			found ++;
			if( h_values[i]== h_bufferV[offset+j]) {
				correctValue++;
			}else {
				incorrectValue++;
			}
		}
		else {
			notfound++;
			incorrectValue++;
		}
		
	}
	printf("found:[%u], notfound:[%u] , correct:[%u] ,incorrect:[%u]\n",found,notfound,correctValue,incorrectValue);
	free(h_keys);
	free(h_values);
	free(h_bufferK);
	free(h_bufferV);
	free(h_bucketOffset);
	free(h_bucketSize);
}
void verifyHashTable2(KEY_PTR  d_keys, VALUE_PTR d_values, unsigned int length, KEY_PTR  d_TK,VALUE_PTR d_TV ,unsigned int tableLength, unsigned int * d_bucketSize, unsigned int bucketCount){

	KEY_PTR		 h_keys;
	VALUE_PTR	 h_values;
	KEY_PTR		 h_TK;
	VALUE_PTR	 h_TV;
	unsigned int * 	 h_bucketSize;
	
	unsigned int *	 bCount;
	unsigned int b;
	
	h_keys=(KEY_PTR) malloc(length * (KEY_SIZE));
	h_values=(VALUE_PTR) malloc(length * (VALUE_SIZE));
	h_TK= (KEY_PTR ) malloc( BUCKET_KEY_SIZE* bucketCount);
	h_TV= (VALUE_PTR ) malloc( BUCKET_VALUE_SIZE* bucketCount);
	h_bucketSize= (unsigned int *) malloc ( bucketCount * sizeof(unsigned int));
	bCount= (unsigned int * ) malloc( bucketCount *sizeof(unsigned int));
	
	checkCudaErrors( cudaMemcpy(h_keys, d_keys, length * (KEY_SIZE), cudaMemcpyDeviceToHost));
	checkCudaErrors( cudaMemcpy(h_values, d_values, length * (VALUE_SIZE), cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy(h_TK, d_TK,  BUCKET_KEY_SIZE*bucketCount, cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy(h_TV, d_TV,  BUCKET_VALUE_SIZE*bucketCount, cudaMemcpyDeviceToHost) );
	checkCudaErrors( cudaMemcpy(h_bucketSize, d_bucketSize, bucketCount * sizeof(unsigned int), cudaMemcpyDeviceToHost) );
	
	CheckCUDAError();
	/*for(int j=0; j<bucketCount; j++){
		printf(" [%u]:%u ,",h_bucketSeed[j],j);
		}
		printf("\n");*/

	unsigned int correct=0;
	unsigned int incorrect=0;
	
	memset(bCount,0,bucketCount * sizeof(unsigned int) );
	for (int i=0;i <length; i++){
		if( host_getHashValue2(h_keys[i],h_TK,h_TV,h_bucketSize,bucketCount,&b) == h_values[i]) {
			correct++;
			bCount[b]++;
		}
		else incorrect++;
	}
	printf("total verified %u, incorrect %u\n",correct,incorrect);
	
	//for(int  i =0; i<bucketCount; i++){ printf(" Count Bucket-%d=%u\n",i,bCount[i]);}
	free(h_keys);
	free(h_values);
	free(h_TK);
	free(h_TV);
	free(h_bucketSize);
	free(bCount);
	
}
extern "C"
void createHashTable2(KEY_PTR d_keys,VALUE_PTR d_values, unsigned int length, KEY_PTR *  d_TK,VALUE_PTR * d_TV,unsigned int * tableLength, unsigned int ** d_bucketSeed,unsigned int * bucketCount){


	unsigned int * d_count;
	unsigned int * d_offset;	
	unsigned int * d_start;
	
	KEY_PTR d_bufferK;
	VALUE_PTR d_bufferV;
	unsigned int * d_randomSeed;
	unsigned int * h_randomSeed;

//	unsigned int timer = 0;
//	cutilCheckError(cutCreateTimer(&timer));
	StopWatchInterface *timer = NULL;
	sdkCreateTimer(&timer);
	sdkResetTimer(&timer);
	sdkStartTimer(&timer);
	srand ( time(NULL) );
	*bucketCount=(length /409)+1; //ceil
	unsigned int dataSize=length*sizeof(unsigned int);
	unsigned int bucketDataSize=*bucketCount*sizeof(unsigned int);

	checkCudaErrors( cudaMalloc( (void**) &d_offset, dataSize));
	//allocate count 
	checkCudaErrors( cudaMalloc( (void**) &d_count, bucketDataSize));
	
	
	//initialize offset to zero
	checkCudaErrors( cudaMemset(d_offset,0,dataSize));
	//initialize count to zero
	checkCudaErrors( cudaMemset(d_count,0,bucketDataSize));
	
	
	/**********Initiating Phase 1*********/
//	cutilCheckError(cutStartTimer(timer));
	sdkStartTimer(&timer);
	//launch phase 1 , bucket allocation
	phase12<<<length/512+1,512>>>(d_keys,d_offset,length,d_count,*bucketCount);
	CheckCUDAError();
//	cutilCheckError(cutStopTimer(timer));
	sdkStopTimer(&timer);

	/************  Calculating Start of each bucket (prefix sum of Count) **********/
	//allocate and initiazlie start 
	checkCudaErrors( cudaMalloc( (void**) &d_start, bucketDataSize));
	checkCudaErrors( cudaMemset(d_start,0,bucketDataSize));

	//find prefix sum 
	CUDPPConfiguration config;
    config.op = CUDPP_ADD;
	config.datatype = CUDPP_UINT;
    config.algorithm = CUDPP_SCAN;
    config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;

    CUDPPHandle scanplan = 0;
    CUDPPResult result = cudppPlan(&scanplan, config, *bucketCount, 1, 0);
//	cutilCheckError(cutStartTimer(timer));
	sdkStartTimer(&timer);
	// Run the scan
    cudppScan(scanplan, d_start, d_count, *bucketCount);
    CheckCUDAError();
//	cutilCheckError(cutStopTimer(timer));
	sdkStopTimer(&timer);
	cudppDestroyPlan(scanplan);
	

	/************* Copying to buffer **************/

	
	//allocate buffer
	checkCudaErrors( cudaMalloc( (void**) &d_bufferK, length*KEY_SIZE));
	checkCudaErrors( cudaMalloc( (void**) &d_bufferV, length*VALUE_SIZE));
//	cutilCheckError(cutStartTimer(timer));
	sdkStartTimer(&timer);
	//copy to buckets
	copyToBucket2<<<length/512+1,512>>>(d_keys,d_values,d_offset,length,d_start,*bucketCount,d_bufferK,d_bufferV);
	CheckCUDAError();
//	cutilCheckError(cutStopTimer(timer));
	sdkStopTimer(&timer);
	
	//free up some resources
	checkCudaErrors(cudaFree(d_offset));

	/***************     Cuckoo Hashing        ******************/
	checkCudaErrors( cudaMalloc( (void**) d_bucketSeed, bucketDataSize));
	checkCudaErrors( cudaMalloc( (void**) d_TK, (*bucketCount)*BUCKET_KEY_SIZE));
	checkCudaErrors( cudaMalloc( (void**) d_TV, (*bucketCount)*BUCKET_VALUE_SIZE));

	unsigned int randomSeedSize=*bucketCount*MAX_SEED_COUNT* sizeof(int);
	h_randomSeed=( unsigned int *)  malloc(randomSeedSize);
	checkCudaErrors( cudaMalloc( (void**) &d_randomSeed,randomSeedSize));
	for(unsigned int i=0;i<*bucketCount;i++){
		for(int j=0;j<MAX_SEED_COUNT;j++){
			*(h_randomSeed+i*MAX_SEED_COUNT+j)=rand();
		}
	}
	checkCudaErrors( cudaMemcpy( d_randomSeed, h_randomSeed, randomSeedSize,cudaMemcpyHostToDevice) );
	free(h_randomSeed);


	//phase2<<<*bucketCount,512>>>(d_buffer,d_start,d_count,*d_bucketSeed,d_bucketState,*bucketCount,*d_T,d_randomSeed,MAX_SEED_COUNT);
	//phase22<<<*bucketCount,512>>>(d_bufferK,d_bufferV,d_start,d_count,*d_bucketSeed,*bucketCount,*d_TK,*d_TV,d_randomSeed,MAX_SEED_COUNT);
	phase2HostLaunch(d_bufferK,d_bufferV,
			d_start, d_count,
			*d_bucketSeed,*bucketCount,
			*d_TK,*d_TV,
			MAX_SEED_COUNT,length);
	CheckCUDAError();
//	cutilCheckError(cutStopTimer(timer));
	sdkStopTimer(&timer);
	cudaError_t err=cudaGetLastError();
	if(cudaSuccess != err ){
		printf("%s\n",cudaGetErrorString(err));
	}
	checkCudaErrors(cudaFree(d_start));
		checkCudaErrors(cudaFree(d_count));

	checkCudaErrors(cudaFree(d_randomSeed));
	checkCudaErrors(cudaFree(d_bufferK));
	checkCudaErrors(cudaFree(d_bufferV));

	
	
	*tableLength=*bucketCount*BLOCK_SIZE;

	DEBUG_CALL(verifyHashTable2(d_keys,d_values,length,*d_TK,*d_TV,*tableLength,*d_bucketSeed,*bucketCount));
//	cutilCheckError(cutDeleteTimer(timer));
	sdkDeleteTimer(&timer);
}
