#ifndef GPUHASH_DEVICE2_H
#define GPUHASH_DEVICE2_H

#define GPUHASH2_IMPL //for gpuhash.h
#include "gpuhash2.h"
/*


__device__ unsigned int hash_h(KEY_T  key, unsigned int bucketCount){
	return ((C0+C1*key)% LARGE_PRIME )% bucketCount;
}

__device__ unsigned int hash_g1(KEY_T key,unsigned int seed){
	return ((C10^seed+(C11^seed)*key)% LARGE_PRIME )%L2_SIZE;
}
__device__ unsigned int hash_g2(KEY_T key,unsigned int seed){
	return ((C20^seed+(C21^seed)*key)% LARGE_PRIME )%L2_SIZE;
}
__device__ unsigned int hash_g3(KEY_T key,unsigned int seed){
	return ((C30^seed+(C31^seed)*key)% LARGE_PRIME )%L2_SIZE;
}

*/

__device__ VALUE_T getHashValue2(KEY_T key,KEY_PTR TK,VALUE_PTR TV,unsigned int *bucketSeed, unsigned int bucketCount){
	unsigned int bucket=HASH_H(key,bucketCount);
		unsigned int seed=bucketSeed[bucket];


		if ( TK[bucket * BLOCK_SIZE + T1_OFFSET + HASH_G1(key,seed) ]==key)
			return TV[bucket * BLOCK_SIZE + T1_OFFSET + HASH_G1(key,seed) ];
		if ( TK[bucket * BLOCK_SIZE + T2_OFFSET + HASH_G2(key,seed) ]==key)
			return TV[bucket * BLOCK_SIZE + T2_OFFSET + HASH_G2(key,seed) ];
		if ( TK[bucket * BLOCK_SIZE + T3_OFFSET + HASH_G3(key,seed) ]==key)
					return TV[bucket * BLOCK_SIZE + T3_OFFSET + HASH_G3(key,seed) ];

		return MAX_INT;

}



#endif //GPUHASH_DEVICE
