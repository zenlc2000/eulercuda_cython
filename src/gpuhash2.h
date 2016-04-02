#ifndef GPUHASH2_H
#define GPUHASH2_H
#include "common.h"



#define LARGE_PRIME 1900813

#define KEY_SIZE (sizeof(KEY_T))
#define VALUE_SIZE (sizeof(VALUE_T))

#define L2_ITEM_COUNT 192
#define L2_KEY_SIZE (KEY_SIZE*L2_ITEM_COUNT)
#define L2_VALUE_SIZE (VALUE_SIZE* L2_ITEM_COUNT)
#define MAX_ITERATIONS 100
#define MAX_INT 0xffffffff
#define MAX_SEED_COUNT 25
#define C0  0x01010101
#define C1	0x12345678
#define C10 0xABCDEFAB
#define C11 0xCDEFABCD
#define C20 0xEFABCDEF
#define C21 0xBAFEDCBA
#define C30 0xFEDCBAFE
#define C31 0xDCBAFEDC

#define T1_OFFSET 0
#define T2_OFFSET L2_ITEM_COUNT
#define T3_OFFSET (L2_ITEM_COUNT << 1)

#define BLOCK_SIZE ((L2_ITEM_COUNT<<1) +(L2_ITEM_COUNT ))


#define MAX_BUCKET_ITEM (BLOCK_SIZE)
//#define KEY_OFFSET (0)
//#define VALUE_OFFSET (0)
//#define ITEM_SIZE (KEY_SIZE+VALUE_SIZE)
#define BUCKET_KEY_SIZE (KEY_SIZE * MAX_BUCKET_ITEM)
#define BUCKET_VALUE_SIZE (VALUE_SIZE * MAX_BUCKET_ITEM)
/*
#define GET_HASH_ITEM_OFFSET(ptrTable,blockIdx,itemIdx) ((ptrTable)+(blockIdx)*BUCKET_SIZE + (itemIdx)*BUCKET_ITEM_SIZE )
#define HASH_VALUE_PTR(ptrTable,blockIdx,itemIdx) ((VALUE_PTR ) (GET_HASH_ITEM_OFFSET(ptrTable,blockIdx,itemIdx)+VALUE_OFFSET))
#define GET_HASH_VALUE(ptrTable,blockIdx,itemIdx) (*(HASH_VALUE_PTR(ptrTable,blockIdx,itemIdx)))
//#define SET_HASH_VALUE(ptrTable,blockIdx,itemIdx,value)
#define HASH_KEY_PTR(ptrTable, blockIdx,itemIdx) ((KEY_PTR)(GET_HASH_ITEM_OFFSET(ptrTable,blockIdx,itemIdx)+KEY_OFFSET))
#define GET_HASH_KEY(ptrTable,blockIdx,itemIdx) (*(HASH_KEY_PTR(ptrTable,blockIdx,itemIdx)))
*/
#define GET_KEY_INDEX(blockIdx,itemIdx) ((blockIdx)*MAX_BUCKET_ITEM+(itemIdx))
#define GET_VALUE_INDEX(blockIdx,itemIdx) ((blockIdx)*MAX_BUCKET_ITEM+(itemIdx))

#define HASH_H(key,bucketCount) (((C0+C1*key)% LARGE_PRIME )% bucketCount)

#define HASH_G1(key,seed)  (((C10^seed+(C11^seed)*key)% LARGE_PRIME )%L2_ITEM_COUNT)
#define HASH_G2(key,seed)  (((C20^seed+(C21^seed)*key)% LARGE_PRIME )%L2_ITEM_COUNT)
#define HASH_G3(key,seed)  (((C30^seed+(C31^seed)*key)% LARGE_PRIME )%L2_ITEM_COUNT)


#ifndef GPUHASH2_IMPL

extern "C"
void createHashTable2(KEY_PTR d_keys,VALUE_PTR d_values, unsigned int length, KEY_PTR * d_TK,VALUE_PTR * d_TV,unsigned int * tableLength, unsigned int ** d_seed,unsigned int * bucketCount);

#endif//GPUHASH_IMPL




inline VALUE_T host_getHashValue2(KEY_T key, KEY_PTR TK,VALUE_PTR TV,unsigned int * bucketSeed, unsigned int bucketCount,unsigned int * bucket){
	*bucket=HASH_H(key,bucketCount);
	unsigned int seed=bucketSeed[*bucket];
	if ( TK[*bucket * BLOCK_SIZE + T1_OFFSET + HASH_G1(key,seed) ]==key)
		return TV[*bucket * BLOCK_SIZE + T1_OFFSET + HASH_G1(key,seed) ];
	if ( TK[*bucket * BLOCK_SIZE + T2_OFFSET + HASH_G2(key,seed) ]==key)
		return TV[*bucket * BLOCK_SIZE + T2_OFFSET + HASH_G2(key,seed) ];
	if ( TK[*bucket * BLOCK_SIZE + T3_OFFSET + HASH_G3(key,seed) ]==key)
		return TV[*bucket * BLOCK_SIZE + T3_OFFSET + HASH_G3(key,seed) ];

	return MAX_INT;
}

#endif //GPUHASH_H
