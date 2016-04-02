#ifndef GPUHASH_DEVICE_H
#define GPUHASH_DEVICE_H

#define GPUHASH_IMPL //for gpuhash.h
#include "gpuhash.h"
/*

*/
__forceinline__ __device__ unsigned int hash_h(KEY_T  key, unsigned int bucketCount){
	return ((C0+C1*key)% LARGE_PRIME )% bucketCount;
}

__forceinline__ __device__ unsigned int hash_g1(KEY_T key,unsigned int seed){
	return ((C10^seed+(C11^seed)*key)% LARGE_PRIME )%L2_SIZE;
}
__forceinline__ __device__ unsigned int hash_g2(KEY_T key,unsigned int seed){
	return ((C20^seed+(C21^seed)*key)% LARGE_PRIME )%L2_SIZE;
}
__forceinline__ __device__ unsigned int hash_g3(KEY_T key,unsigned int seed){
	return ((C30^seed+(C31^seed)*key)% LARGE_PRIME )%L2_SIZE;
}



__forceinline__ __device__ VALUE_T getHashValue(KEY_T key,KEY_PTR TK,VALUE_PTR TV,unsigned int *bucketSize, unsigned int bucketCount){

	unsigned int bucket=hash_h(key,bucketCount);
	unsigned int l=0;
	unsigned int r=bucketSize[bucket];
	unsigned int mid;
	while(l<r){
		mid =l+((r-l)/2);
		//if( (GET_HASH_KEY(T,bucket,mid)) < key) {
		if( TK[GET_KEY_INDEX(bucket,mid)] <  key) {
			l=mid+1;
		}else {
			r=mid;
		}
	}
	//if(l < bucketSize[bucket] && GET_HASH_KEY(T,bucket,l)==key){ 
	if(l < bucketSize[bucket] && TK[GET_KEY_INDEX(bucket,l)]==key){ 
	//	return GET_HASH_VALUE(T,bucket,l);
		return TV[GET_VALUE_INDEX(bucket,l)];
	}else {
		return MAX_INT;
	}
}
#endif //GPUHASH_DEVICE
