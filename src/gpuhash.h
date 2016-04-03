#ifndef GPUHASH_H
#define GPUHASH_H
#include "common.h"



#define LARGE_PRIME 1900813
#define L2_SIZE 192
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
#define T2_OFFSET L2_SIZE 
#define T3_OFFSET (L2_SIZE << 1)

#define BLOCK_SIZE ((L2_SIZE<<1) +(L2_SIZE <<2))
typedef struct HASH_TABLE{
	KEY_PTR keys;
	VALUE_PTR values;
	unsigned int bucketCount;
	unsigned int * bucketLength;
} HASH_TABLE_T;


#define MAX_BUCKET_ITEM (520)
#define KEY_SIZE (sizeof(KEY_T))
#define VALUE_SIZE (sizeof(VALUE_T))
#define KEY_OFFSET (0)
#define VALUE_OFFSET (0)
#define ITEM_SIZE (KEY_SIZE+VALUE_SIZE)
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


#ifndef GPUHASH_IMPL

//extern "C"
void createHashTable(KEY_PTR d_keys,VALUE_PTR d_values, unsigned int length, KEY_PTR * d_TK,VALUE_PTR * d_TV,unsigned int * tableLength, unsigned int ** d_bucketSize,unsigned int * bucketCount);

#endif//GPUHASH_IMPL


inline unsigned int  host_hash_h(KEY_T key, unsigned int bucketCount){
	return ((C0+C1*key)% LARGE_PRIME )% bucketCount;
}


inline VALUE_T host_getHashValue(KEY_T key, KEY_PTR TK,VALUE_PTR TV,unsigned int * bucketSize, unsigned int bucketCount,unsigned int * bucket){
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
		printf("value for key:{%u}not found in bucket %u \nprinting bucket data",key,*bucket);
		for(unsigned int i =0;i<MAX_BUCKET_ITEM;i++){
			printf("[%u]:{%lu}=>%u\t\t",i,(unsigned long )(TK[GET_KEY_INDEX((*bucket),i)]),(unsigned int )(TV[GET_VALUE_INDEX((*bucket),i)]));
		}
		printf("\n");
		return MAX_INT;
	}
}

#endif //GPUHASH_H
