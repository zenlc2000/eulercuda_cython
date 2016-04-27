#define CUDPP_STATIC_LIB



#include "cutil_inline.h"
#include "cudpp.h"
#include "utils.h"
#include "gpuhash_device.h"
#include "gpuhash_kernel.cu"
#include <time.h>


#ifdef EULER_NDEBUG
#define DEBUG_GPUHASH_CU(x)
#else
#define DEBUG_GPUHASH_CU(x) x
#endif
#define DEBUG_CALL(x) DEBUG_GPUHASH_CU(x)


void bucketSortHost(KEY_PTR   bufferK,VALUE_PTR bufferV, unsigned int * start,unsigned int * bucketSize,unsigned int bucketCount,KEY_PTR TK,VALUE_PTR TV,unsigned int bid,unsigned int tid){

		KEY_T keys[512];
		unsigned int keyCount=0;
		unsigned int blockOffset=start[bid];
		unsigned int size=bucketSize[bid];
		if (tid< size) 
			keys[tid]=bufferK[blockOffset+tid];// 
//		__syncthreads();
		if(tid<size){
			for(int i=0; i<size; i++){
				if( keys[tid]> keys[i] ) keyCount++;
			}
		}
		//__syncthreads();
		
		if(tid<size) {
			TK[GET_KEY_INDEX(bid,keyCount)]=keys[tid];
			TV[GET_VALUE_INDEX(bid,keyCount)]=bufferV[blockOffset+tid];
		}
}


void verifyBuffer(KEY_PTR d_keys,VALUE_PTR d_values,unsigned int length, KEY_PTR d_bufferK, VALUE_PTR d_bufferV, unsigned int * d_bucketOffset,unsigned int * d_bucketSize, unsigned int bucketCount){

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
	
	cutilSafeCall( cudaMemcpy(h_keys, d_keys, length * (KEY_SIZE), cudaMemcpyDeviceToHost));
	cutilSafeCall( cudaMemcpy(h_values, d_values, length * (VALUE_SIZE), cudaMemcpyDeviceToHost));
	cutilSafeCall( cudaMemcpy(h_bufferK, d_bufferK, length * (KEY_SIZE), cudaMemcpyDeviceToHost));
	cutilSafeCall( cudaMemcpy(h_bufferV, d_bufferV, length * (VALUE_SIZE), cudaMemcpyDeviceToHost));
	cutilSafeCall( cudaMemcpy(h_bucketOffset, d_bucketOffset, bucketCount* (sizeof(unsigned int)), cudaMemcpyDeviceToHost));
	cutilSafeCall( cudaMemcpy(h_bucketSize, d_bucketSize, bucketCount * (sizeof(unsigned int)), cudaMemcpyDeviceToHost));

	unsigned int found=0;
	unsigned int notfound=0;
	unsigned int correctValue=0;
	unsigned int incorrectValue=0;

	for (unsigned int i=0; i<length ; i++){
		unsigned int j=0;
		unsigned int bucket= host_hash_h(h_keys[i],bucketCount);
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
}/*
unsigned int  host_hash_h(KEY_T key, unsigned int bucketCount){
	return ((C0+C1*key)% LARGE_PRIME )% bucketCount;
}*/
unsigned int host_hash_g1(KEY_T key,unsigned int seed){
	return ((C10^seed+(C11^seed)*key)% LARGE_PRIME )%L2_SIZE;
}
unsigned int host_hash_g2(KEY_T key,unsigned int seed){
	return ((C20^seed+(C21^seed)*key)% LARGE_PRIME )%L2_SIZE;
}
unsigned int host_hash_g3(KEY_T  key,unsigned int seed){
	return ((C30^seed+(C31^seed)*key)% LARGE_PRIME )%L2_SIZE;
}
/*

VALUE_T host_getHashValue(KEY_T key, KEY_PTR TK,VALUE_PTR TV,unsigned int * bucketSize, unsigned int bucketCount,unsigned int * bucket){
	*bucket=host_hash_h(key,bucketCount);

	unsigned int l=0;
	unsigned int r=bucketSize[*bucket];
	unsigned int mid;
	while(l<r){
		mid =l+((r-l)/2);
		//if( GET_HASH_KEY(T,(*bucket),mid) <key) {
		if( TK[GET_KEY_INDEX(*bucket,mid)] <key) {
			l=mid+1;
		}else {
			r=mid;
		}
	}
	//if(l < bucketSize[*bucket] && (GET_HASH_KEY(T,(*bucket),l))==key){ 
	if(l < bucketSize[*bucket] && TK[GET_KEY_INDEX(*bucket,l)]==key){ 
		//return GET_HASH_VALUE(T,(*bucket),l);
		return TV[GET_VALUE_INDEX(*bucket,l)];
	}else {
		printf("value not found\nprinting bucket data");
		for(unsigned int i =0;i<bucketSize[*bucket];i++){
			printf("[%u]:{%lu}=>%u\t\t",i,(unsigned long )(TK[GET_KEY_INDEX((*bucket),i)]),(unsigned int )(TV[GET_VALUE_INDEX((*bucket),i)]));
		}
		printf("\n");
		return MAX_INT;
	}
}
*/
void verifyHashTable(KEY_PTR  d_keys, VALUE_PTR d_values, unsigned int length, KEY_PTR  d_TK,VALUE_PTR d_TV ,unsigned int tableLength, unsigned int * d_bucketSize, unsigned int bucketCount){

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
	
	cutilSafeCall( cudaMemcpy(h_keys, d_keys, length * (KEY_SIZE), cudaMemcpyDeviceToHost));
	cutilSafeCall( cudaMemcpy(h_values, d_values, length * (VALUE_SIZE), cudaMemcpyDeviceToHost) );
	cutilSafeCall( cudaMemcpy(h_TK, d_TK,  BUCKET_KEY_SIZE*bucketCount, cudaMemcpyDeviceToHost) );
	cutilSafeCall( cudaMemcpy(h_TV, d_TV,  BUCKET_VALUE_SIZE*bucketCount, cudaMemcpyDeviceToHost) );
	cutilSafeCall( cudaMemcpy(h_bucketSize, d_bucketSize, bucketCount * sizeof(unsigned int), cudaMemcpyDeviceToHost) );
	
	CheckCUDAError();
	/*for(int j=0; j<bucketCount; j++){
		printf(" [%u]:%u ,",h_bucketSeed[j],j);
		}
		printf("\n");*/
	for(int j=0; j<bucketCount; j++){
		if(h_bucketSize[j] >=500) printf("possible invalid bucket [%u] size[%u]\n",j,h_bucketSize[j]);
		/*printf("BUCKET -> [%u]\n",j);
		for (int l=0; l<3; l++){
			printf("TABLE : T%d\n",l);
			for (int k=0; k< L2_SIZE; k++){
				printf("[%u]={%u}, ",h_T[j*BLOCK_SIZE + L2_SIZE*2*l +k*2],h_T[j*BLOCK_SIZE + L2_SIZE*2*l +k*2+1]);
			}
			printf("\n");
		}*/
	}
	unsigned int correct=0;
	unsigned int incorrect=0;
	
	memset(bCount,0,bucketCount * sizeof(unsigned int) );
	for (int i=0;i <length; i++){
		if( host_getHashValue(h_keys[i],h_TK,h_TV,h_bucketSize,bucketCount,&b) == h_values[i]) { 
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

/**	Create Hash table with linear probing sorting the buckets and performing Binary Search on Lookup****/
extern "C"
void createHashTable(KEY_PTR d_keys,VALUE_PTR d_values, unsigned int length, KEY_PTR *  d_TK,VALUE_PTR * d_TV,unsigned int * tableLength, unsigned int ** d_bucketSize,unsigned int * bucketCount){
		
	unsigned int * d_offset;	
	unsigned int * d_start;
	
	KEY_PTR d_bufferK;
	VALUE_PTR d_bufferV;
	

	unsigned int timer = 0;
	cutilCheckError(cutCreateTimer(&timer));
	srand ( time(NULL) );
	*bucketCount=(length /409)+1; //ceil
	unsigned int dataSize=length*sizeof(unsigned int);
	unsigned int bucketDataSize=*bucketCount*sizeof(unsigned int);

	cutilSafeCall( cudaMalloc( (void**) &d_offset, dataSize));
	//allocate count 
	cutilSafeCall( cudaMalloc( (void**) d_bucketSize, bucketDataSize));
	
	
	//initialize offset to zero
	cutilSafeCall( cudaMemset(d_offset,0,dataSize));
	//initialize count to zero
	cutilSafeCall( cudaMemset(*d_bucketSize,0,bucketDataSize));
	
	
	/**********Initiating Phase 1*********/
	cutilCheckError(cutStartTimer(timer));
	//launch phase 1 , bucket allocation
	phase1<<<length/512+1,512>>>(d_keys,d_offset,length,*d_bucketSize,*bucketCount);
	CheckCUDAError();
	cutilCheckError(cutStopTimer(timer));
	/************  Calculating Start of each bucket (prefix sum of Count) **********/
	//allocate and initiazlie start 
	cutilSafeCall( cudaMalloc( (void**) &d_start, bucketDataSize));
	cutilSafeCall( cudaMemset(d_start,0,bucketDataSize));

	//find prefix sum 
	CUDPPConfiguration config;
    	config.op = CUDPP_ADD;
	config.datatype = CUDPP_UINT;
    	config.algorithm = CUDPP_SCAN;
    	config.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;
  
	CUDPPHandle scanplan = 0;
    	CUDPPResult result = cudppPlan(&scanplan, config, *bucketCount, 1, 0);
	cutilCheckError(cutStartTimer(timer));
	// Run the scan
    	cudppScan(scanplan, d_start, *d_bucketSize, *bucketCount);
    	CheckCUDAError();
	cutilCheckError(cutStopTimer(timer));
	cudppDestroyPlan(scanplan);
	

	/************* Copying to buffer **************/

	
	//allocate buffer
	cutilSafeCall( cudaMalloc( (void**) &d_bufferK, length*KEY_SIZE));
	cutilSafeCall( cudaMalloc( (void**) &d_bufferV, length*VALUE_SIZE));
	cutilCheckError(cutStartTimer(timer));
	//copy to buckets
	copyToBucket<<<length/512+1,512>>>(d_keys,d_values,d_offset,length,d_start,*bucketCount,d_bufferK,d_bufferV);
	CheckCUDAError();
	cutilCheckError(cutStopTimer(timer));
	
	//verifyBuffer(d_keys,d_values,length,d_bufferK,d_bufferV,d_start,*d_bucketSize,*bucketCount);

	
	//free up some resources
	cutilSafeCall(cudaFree(d_offset));
	
	/***************     Cuckoo Hashing        ******************/
	
	cutilSafeCall( cudaMalloc( (void**) d_TK, (*bucketCount)*BUCKET_KEY_SIZE));
	cutilSafeCall( cudaMalloc( (void**) d_TV, (*bucketCount)*BUCKET_VALUE_SIZE));
	bucketSort<<<*bucketCount,32>>>(d_bufferK,d_bufferV,d_start,*d_bucketSize,*bucketCount,*d_TK,*d_TV);
	CheckCUDAError();
	cutilCheckError(cutStopTimer(timer));
	cudaError_t err=cudaGetLastError();
	if(cudaSuccess != err ){
		printf("%s\n",cudaGetErrorString(err));
	}
	
//	verifyBuffer(d_keys,d_values,length,d_bufferK,d_bufferV,d_start,*d_bucketSize,*bucketCount);
	cutilSafeCall(cudaFree(d_start));
	cutilSafeCall(cudaFree(d_bufferK));
	cutilSafeCall(cudaFree(d_bufferV));
	
	*tableLength=*bucketCount*MAX_BUCKET_ITEM ;

	DEBUG_CALL(verifyHashTable(d_keys,d_values,length,*d_TK,*d_TV,*tableLength,*d_bucketSize,*bucketCount));
	cutilCheckError(cutDeleteTimer(timer));
	
}
