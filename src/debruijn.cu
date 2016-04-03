#define CUDPP_STATIC_LIB
// #include "../nvidia_sdk/C/common/inc/cutil_inline.h"
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_cuda.h" // lib above replaced w/this one at CUDA 5.0
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include <stdio.h>      /* defines printf for tests */
#include <time.h>       /* defines time_t for timings in the test */
#include <math.h>

#ifdef linux
# include <endian.h>    /* attempt to define endianness */
#endif

//#include "debruijn.h"
#include "kmer.h"
#include "graph.h"
#include "cudpp.h"
#include "utils.h"
#include "common.h"
#include "gpuhash.h"
#include "gpuhash_device.h"

//#include "utils.cpp"

#if (defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN) && \
	__BYTE_ORDER == __LITTLE_ENDIAN) || \
	(defined(i386) || defined(__i386__) || defined(__i486__) || \
	defined(__i586__) || defined(__i686__) || defined(vax) || defined(MIPSEL))
# define HASH_LITTLE_ENDIAN 1
# define HASH_BIG_ENDIAN 0
#elif (defined(__BYTE_ORDER) && defined(__BIG_ENDIAN) && \
	__BYTE_ORDER == __BIG_ENDIAN) || \
	(defined(sparc) || defined(POWERPC) || defined(mc68000) || defined(sel))
# define HASH_LITTLE_ENDIAN 0
# define HASH_BIG_ENDIAN 1
#else
# define HASH_LITTLE_ENDIAN 0
# define HASH_BIG_ENDIAN 0
#endif

/*
 unsigned int _host_hash_h(KEY_T key, unsigned int bucketCount){
 return ((C0+C1*key)% LARGE_PRIME )% bucketCount;
 }


 VALUE_T getHashValue2Host(KEY_T key, TABLE_PTR T,unsigned int *bucketSize, unsigned int bucketCount){

 unsigned int bucket=_host_hash_h(key,bucketCount);
 unsigned int l=0;
 unsigned int r=bucketSize[bucket];
 unsigned int offset=bucket * BUCKET_SIZE;
 unsigned int mid=(l+r)>>1;
 while(l<r){
 mid =l+((r-l)/2);
 if( T[offset+(mid<<1)] <key) {
 l=mid+1;
 }else {
 r=mid;
 }
 }
 if(l < bucketSize[bucket] && T[offset+(l<<1)]==key){
 return T[offset+(l<<1)+1];
 }else {
 return MAX_INT;
 }
 }
 */

/***
 * Inline Printing Routine for l and e structures
 */
inline void printData(unsigned int * d_lstart, unsigned int * d_lcount,
		unsigned int * d_estart, unsigned int * d_ecount, unsigned int length) {

	unsigned int * h_lstart;
	unsigned int * h_lcount;
	unsigned int * h_estart;
	unsigned int * h_ecount;

	h_lstart = (unsigned int *) malloc(sizeof(unsigned int) * length);
	h_lcount = (unsigned int *) malloc(sizeof(unsigned int) * length);
	h_estart = (unsigned int *) malloc(sizeof(unsigned int) * length);
	h_ecount = (unsigned int *) malloc(sizeof(unsigned int) * length);

	checkCudaErrors(
			cudaMemcpy(h_lstart, d_lstart, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_lcount, d_lcount, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_estart, d_estart, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_ecount, d_ecount, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));

	for (unsigned int i = 0; i < length; i++) {
		printf("[%4u]\t es:%4u\t ec:%4u\t ls:%4u\t lc:%4u\n", i, h_estart[i],
				h_ecount[i], h_lstart[i], h_lcount[i]);
	}

	free(h_lstart);
	free(h_lcount);
	free(h_estart);
	free(h_ecount);

}
//__global__ void debruijnCount(unsigned int k,unsigned long kmerCount, unsigned char * idata,unsigned int * icount,unsigned int * vcount, unsigned int * lcount,unsigned int * ecount,unsigned int  validBitMask){

/*
 *  This kernel works on each l-mer ,counting edges of the graph.
 */
__global__ void debruijnCount(KEY_PTR lmerKeys, /* lmer keys	*/
VALUE_PTR lmerValues, /* lmer frequency */
unsigned int lmerCount, /* total lmers */
KEY_PTR TK, /* Keys' pointer for Hash table*/
VALUE_PTR TV, /* Value pointer for Hash table*/
unsigned int * bucketSeed, /* bucketSize: size of each bucket (it should be renamed to bucketSize)*/
unsigned int bucketCount, /* total buckets */
unsigned int * lcount, /* leaving edge count array : OUT */
unsigned int * ecount, /* entering edge count array: OUT */
KEY_T validBitMask /* bit mask for K length encoded bits*/
) {

	unsigned int tid = (blockDim.x * blockDim.y * gridDim.x * blockIdx.y)
			+ (blockDim.x * blockDim.y * blockIdx.x)
			+ (blockDim.x * threadIdx.y) + threadIdx.x;
	if (tid < lmerCount) {
		KEY_T lmer = lmerKeys[tid];
		VALUE_T lmerValue = lmerValues[tid];
		KEY_T prefix = (lmer & (validBitMask << 2)) >> 2;
		KEY_T suffix = (lmer & validBitMask);

		KEY_T lomask = 3;
		VALUE_T prefixIndex = getHashValue(prefix, TK, TV, bucketSeed,
				bucketCount);
		VALUE_T suffixIndex = getHashValue(suffix, TK, TV, bucketSeed,
				bucketCount);
		KEY_T transitionTo = (lmer & lomask);
		KEY_T transitionFrom = ((lmer >> __popcll(validBitMask)) & lomask);
		//atomicAdd(lcount+(prefixIndex<<2 )+transition,lmerValue);
		//atomicAdd(ecount+(suffixIndex<<2)+transition,lmerValue);
		lcount[(prefixIndex << 2) + transitionTo] = lmerValue;
		ecount[(suffixIndex << 2) + transitionFrom] = lmerValue;
	}
}

/**
 * This is cpu version for same kernel. for Debugging purpose only
 */
void debruijnCountHost(KEY_PTR lmerKeys, VALUE_PTR lmerValues,
		unsigned int lmerCount, KEY_PTR TK, VALUE_PTR TV,
		unsigned int * bucketSeed, unsigned int bucketCount,
		unsigned int * lcount, unsigned int * ecount, KEY_T validBitMask,
		unsigned int bitCount, unsigned int tid) {

//	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if (tid < lmerCount) {
		KEY_T lmer = lmerKeys[tid];
		VALUE_T lmerValue = lmerValues[tid];
		KEY_T prefix = (lmer & (validBitMask << 2)) >> 2;
		KEY_T suffix = (lmer & validBitMask);

		KEY_T lomask = 3;
		unsigned int b;
		VALUE_T prefixIndex = host_getHashValue(prefix, TK, TV, bucketSeed,
				bucketCount, &b);
		VALUE_T suffixIndex = host_getHashValue(suffix, TK, TV, bucketSeed,
				bucketCount, &b);
		KEY_T transitionTo = (lmer & lomask);
		KEY_T transitionFrom = ((lmer >> bitCount) & lomask);
		//atomicAdd(lcount+(prefixIndex<<2 )+transition,lmerValue);
		//atomicAdd(ecount+(suffixIndex<<2)+transition,lmerValue);
		if (lcount[(prefixIndex << 2) + transitionTo] > 0) {
			lcount[(prefixIndex << 2) + transitionTo] = lmerValue;
		} else {
			lcount[(prefixIndex << 2) + transitionTo] = lmerValue;
		}
		if (ecount[(suffixIndex << 2) + transitionFrom] > 0) {
			ecount[(suffixIndex << 2) + transitionFrom] = lmerValue;
		} else {
			ecount[(suffixIndex << 2) + transitionFrom] = lmerValue;
		}

	}
}

/*
 * stub for debruijnCountHost for debugging purpose
 */
void verifyDebruijnCountHost(KEY_PTR d_lmerKeys, VALUE_PTR d_lmerValues,
		unsigned int lmerCount, KEY_PTR d_TK, VALUE_PTR d_TV,
		unsigned int * d_bucketSeed, unsigned int bucketCount,
		unsigned int * d_lcount, unsigned int * d_ecount, KEY_T validBitMask,
		unsigned int kmerCount) {

	KEY_PTR h_lmerKeys;
	VALUE_PTR h_lmerValues;
	KEY_PTR h_TK;
	VALUE_PTR h_TV;
	unsigned int * h_bucketSeed;
	unsigned int * h_lcount;
	unsigned int * h_ecount;
	unsigned int * hq_lcount;
	unsigned int * hq_ecount;

	h_lmerKeys = (KEY_PTR) malloc(lmerCount * sizeof(KEY_T));
	h_lmerValues = (VALUE_PTR) malloc(lmerCount * sizeof(VALUE_T));
	h_TK = (KEY_PTR) malloc(bucketCount * BUCKET_KEY_SIZE);
	h_TV = (VALUE_PTR) malloc(bucketCount * BUCKET_VALUE_SIZE);
	h_bucketSeed = (unsigned int *) malloc(bucketCount * sizeof(unsigned int));
	h_lcount = (unsigned int *) malloc(4 * kmerCount * sizeof(unsigned int));
	h_ecount = (unsigned int *) malloc(4 * kmerCount * sizeof(unsigned int));
	hq_lcount = (unsigned int *) malloc(4 * kmerCount * sizeof(unsigned int));
	hq_ecount = (unsigned int *) malloc(4 * kmerCount * sizeof(unsigned int));

	checkCudaErrors(
			cudaMemcpy(h_lmerKeys, d_lmerKeys, lmerCount * KEY_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_lmerValues, d_lmerValues, lmerCount * VALUE_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_TK, d_TK, bucketCount * BUCKET_KEY_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_TV, d_TV, bucketCount * BUCKET_VALUE_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_bucketSeed, d_bucketSeed,
					bucketCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));

	checkCudaErrors(
			cudaMemcpy(hq_lcount, d_lcount,
					4 * kmerCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(hq_ecount, d_ecount,
					4 * kmerCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));

	memset(h_lcount, 0, 4 * kmerCount * sizeof(unsigned int));
	memset(h_ecount, 0, 4 * kmerCount * sizeof(unsigned int));
	unsigned int bitCount = 0;
	KEY_T bit = 1;
	while (bit != 0) {
		if (bit & validBitMask)
			bitCount++;
		bit = bit << 1;
	}
	unsigned int edgesCount = 0;
	for (unsigned int i = 0; i < lmerCount; i++) {
		debruijnCountHost(h_lmerKeys, h_lmerValues, lmerCount, h_TK, h_TV,
				h_bucketSeed, bucketCount, h_lcount, h_ecount, validBitMask,
				bitCount, i);
		edgesCount += h_lmerValues[i];
	}
	unsigned int esum = 0;
	unsigned int qesum = 0;
	unsigned int lsum = 0;
	unsigned int qlsum = 0;
	unsigned int ei = 0;
	unsigned int li = 0;
	unsigned int qei = 0;
	unsigned int qli = 0;
	for (int j = 0; j < 4 * kmerCount; j++) {
		esum += h_ecount[j];
		lsum += h_lcount[j];
		qesum += hq_ecount[j];
		qlsum += hq_lcount[j];
		if (esum > edgesCount && ei < 1)
			ei = j;
		if (lsum > edgesCount && li < 1)
			li = j;
		if (qesum > edgesCount && qei < 1)
			qei = j;
		if (qlsum > edgesCount && qli < 1)
			qli = j;
	}
	printf(
			"lmerCount: %u, esum: %u, lsum: %u, ei: %u li:%u \n qesum:%u , qlsum:%u \n",
			edgesCount, esum, lsum, ei, li, qesum, qlsum);

	unsigned int enc = 0;
	unsigned int lnc = 0;
	for (unsigned int k = 0; k < 4 * kmerCount; k++) {
		if (h_lcount[k] != hq_lcount[k])
			lnc++;
		if (h_ecount[k] != hq_ecount[k])
			enc++;
	}
	printf("enc: %u,  lnc:%u \n", enc, lnc);
	free(h_lmerValues);
	free(h_TK);
	free(h_TV);
	free(h_bucketSeed);
	free(h_lcount);
	free(h_ecount);
	free(hq_lcount);
	free(hq_ecount);

}

/**
 * CPU prefix scan
 */
void prefixScan(unsigned int * h_out, unsigned int * h_in, unsigned int length,
		bool inclusive) {

	memset(h_out, 0, length * sizeof(unsigned int));
	/*calculate gold*/
	if (inclusive) {
		h_out[0] = h_in[0];
	} else {
		h_out[0] = 0;
	}
	for (unsigned int i = 1; i < length; i++) {
		h_out[i] = h_out[i - 1] + h_in[i - (inclusive ? 0 : 1)];
	}

}
/*
 * prefix sum validator
 **/
void validatePrefixScan(unsigned int * d_output, unsigned int * d_input,
		unsigned int length, bool inclusive) {

	unsigned int * h_input;
	unsigned int * h_output;
	unsigned int * hq_output;

	h_input = (unsigned int *) malloc(length * sizeof(unsigned int));
	h_output = (unsigned int *) malloc(length * sizeof(unsigned int));
	hq_output = (unsigned int *) malloc(length * sizeof(unsigned int));

	checkCudaErrors(
			cudaMemcpy(h_input, d_input, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(hq_output, d_output, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));

	/*
	 memset(h_output,0,length*sizeof(unsigned int));
	 if(inclusive)
	 {	h_output[0]=h_input[0];}
	 else
	 {	h_output[0]=0;	}
	 for(unsigned int i=1;i<length;i++){
	 h_output[i]=h_output[i-1]+h_input[i-(inclusive?0:1)];
	 }
	 */
	prefixScan(h_output, h_input, length, inclusive);
	/*compare*/
	for (unsigned int j = 0; j < length; j++) {
		if (h_output[j] != hq_output[j]) {
			printf("differnce at index:%u is gold:%u, cudpp\n", j, h_output[j], hq_output[j]);
		}
	}

	free(h_input);
	free(h_output);
	free(hq_output);
}
/*
 *  This kernel works on a k-mer (l-1mer) which are vertices of the graph.
 */
__global__ void setupVertices(KEY_PTR kmerKeys, unsigned int kmerCount,
		KEY_PTR TK, VALUE_PTR TV, unsigned int * bucketSeed,
		unsigned int bucketCount, EulerVertex * ev, unsigned int * lcount,
		unsigned int * loffset, unsigned int * ecount, unsigned int * eoffset) {
	unsigned int tid = (blockDim.x * blockDim.y * gridDim.x * blockIdx.y)
			+ (blockDim.x * blockDim.y * blockIdx.x)
			+ (blockDim.x * threadIdx.y) + threadIdx.x;
	if (tid < kmerCount) {
		KEY_T key = kmerKeys[tid];
		VALUE_T index = getHashValue(key, TK, TV, bucketSeed, bucketCount);
		;
		ev[index].vid = key;
		ev[index].lp = loffset[(index << 2)];
		ev[index].lcount = lcount[(index << 2)] + lcount[(index << 2) + 1]
				+ lcount[(index << 2) + 2] + lcount[(index << 2) + 3];
		ev[index].ep = eoffset[(index << 2)];
		ev[index].ecount = ecount[(index << 2)] + ecount[(index << 2) + 1]
				+ ecount[(index << 2) + 2] + ecount[(index << 2) + 3];
	}
}
void setupVerticesHost(KEY_PTR kmerKeys, unsigned int kmerCount, KEY_PTR TK,
		VALUE_PTR TV, unsigned int * bucketSeed, unsigned int bucketCount,
		EulerVertex * ev, unsigned int * lcount, unsigned int * loffset,
		unsigned int * ecount, unsigned int * eoffset, unsigned int tid) {
//	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if (tid < kmerCount) {
		KEY_T key = kmerKeys[tid];
		unsigned int bucket;
		VALUE_T index = host_getHashValue(key, TK, TV, bucketSeed, bucketCount,
				&bucket);
		;
		ev[index].vid = key;
		ev[index].lp = loffset[(index << 2)];
		ev[index].lcount = lcount[(index << 2)] + lcount[(index << 2) + 1]
				+ lcount[(index << 2) + 2] + lcount[(index << 2) + 3];
		ev[index].ep = eoffset[(index << 2)];
		ev[index].ecount = ecount[(index << 2)] + ecount[(index << 2) + 1]
				+ ecount[(index << 2) + 2] + ecount[(index << 2) + 3];
	}
}

/* 
 *  This kernel works on an l-mer, which represents an edge
 *  in the debruijn Graph.
 */
__global__ void setupEdges( KEY_PTR  lmerKeys,  VALUE_PTR  lmerValues,
		 unsigned int *  lmerOffsets, const unsigned int lmerCount,
		 KEY_PTR  TK, VALUE_PTR  TV, unsigned int *  bucketSeed,
		const unsigned int bucketCount, unsigned int *  l,
		 unsigned int *  e, EulerEdge *  ee,
		 unsigned int *  loffsets, unsigned int *  eoffsets,
		const KEY_T validBitMask) {

	unsigned int tid = (blockDim.x * blockDim.y * gridDim.x * blockIdx.y)
			+ (blockDim.x * blockDim.y * blockIdx.x)
			+ (blockDim.x * threadIdx.y) + threadIdx.x;
	if (tid < lmerCount) {
		KEY_T lmer = lmerKeys[tid];
		VALUE_T lmerValue = lmerValues[tid];
		KEY_T prefix = (lmer & (validBitMask << 2)) >> 2;
		KEY_T suffix = (lmer & validBitMask);
		KEY_T lomask = 3;
		//prefix and suffix index must be less than kmer count
		VALUE_T prefixIndex = getHashValue(prefix, TK, TV, bucketSeed,
				bucketCount);
		VALUE_T suffixIndex = getHashValue(suffix, TK, TV, bucketSeed,
				bucketCount);
		KEY_T transitionTo = (lmer & lomask);
		KEY_T transitionFrom = ((lmer >> __popcll(validBitMask)) & lomask);
		unsigned int loffset = loffsets[(prefixIndex << 2) + transitionTo];
		unsigned int eoffset = eoffsets[(suffixIndex << 2) + transitionFrom];

		unsigned int lmerOffset = lmerOffsets[tid];
		for (unsigned int i = 0; i < lmerValue; i++) {

			ee[lmerOffset].eid =lmerOffset;
			ee[lmerOffset].v1 = prefixIndex;
			ee[lmerOffset].v2 = suffixIndex;
			// lmerOffset;
			ee[lmerOffset].s = lmerValues[lmerCount - 1]
					+ lmerOffsets[lmerCount - 1];

			l[loffset] = lmerOffset;
			e[eoffset] = lmerOffset;
			loffset++;
			eoffset++;
			lmerOffset++;
		}
	}
}
void setupEdgesHost(KEY_PTR const lmerKeys, VALUE_PTR const lmerValues,
		unsigned int * const lmerOffsets, const unsigned int lmerCount,
		KEY_PTR const TK, VALUE_PTR const TV, unsigned int * const bucketSeed,
		const unsigned int bucketCount, unsigned int * const l,
		unsigned int * const e, EulerEdge * const ee,
		unsigned int * const loffsets, unsigned int * const eoffsets,
		const KEY_T validBitMask, const unsigned int tid) {

//	unsigned int tid=(blockDim.x*blockDim.y * gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	if (tid < lmerCount) {
		KEY_T lmer = lmerKeys[tid];
		VALUE_T lmerValue = lmerValues[tid];
		KEY_T prefix = (lmer & (validBitMask << 2)) >> 2;
		KEY_T suffix = (lmer & validBitMask);
		KEY_T lomask = 3;
		unsigned int bucket;
		//prefix and suffix index must be less than kmer count
		VALUE_T prefixIndex = host_getHashValue(prefix, TK, TV, bucketSeed,
				bucketCount, &bucket);
		VALUE_T suffixIndex = host_getHashValue(suffix, TK, TV, bucketSeed,
				bucketCount, &bucket);
		KEY_T transitionTo = (lmer & lomask);
		KEY_T transitionFrom = ((lmer >> 16) & lomask);
		unsigned int loffset = loffsets[(prefixIndex << 2) + transitionTo];
		unsigned int eoffset = eoffsets[(suffixIndex << 2) + transitionFrom];

		unsigned int lmerOffset = lmerOffsets[tid];
		for (int i = 0; i < lmerValue; i++) {
			ee[lmerOffset].eid = lmerOffset;
			ee[lmerOffset].v1 = prefixIndex;
			ee[lmerOffset].v2 = suffixIndex;
			ee[lmerOffset].s = lmerValues[lmerCount - 1]
					+ lmerOffsets[lmerCount - 1];

			l[loffset] = lmerOffset;
			e[eoffset] = lmerOffset;
			loffset++;
			eoffset++;
			lmerOffset++;
		}
	}
}

void verifySetupEdges(KEY_PTR d_lmerKeys, VALUE_PTR d_lmerValues,
		unsigned int * d_lmerOffsets, const unsigned int lmerCount,
		KEY_PTR d_TK, VALUE_PTR d_TV, unsigned int * d_bucketSeed,
		const unsigned int bucketCount, unsigned int * d_l, unsigned int * d_e,
		EulerEdge * d_ee, unsigned int * d_lcount, unsigned int * d_loffsets,
		unsigned int * d_ecount, unsigned int * d_eoffsets,
		unsigned int kmerCount, unsigned int ecount, const KEY_T validBitMask) {

	KEY_PTR h_lmerKeys;
	VALUE_PTR h_lmerValues;
	unsigned int * h_lmerOffsets;
	KEY_PTR h_TK;
	VALUE_PTR h_TV;
	unsigned int * h_bucketSeed;
	unsigned int * h_l;
	unsigned int * h_e;
	EulerEdge * h_ee;
	unsigned int * h_loffsets;
	unsigned int * h_lcount;
	unsigned int * h_eoffsets;
	unsigned int * h_ecount;

	h_lmerKeys = (KEY_PTR) malloc(lmerCount * KEY_SIZE);
	h_lmerValues = (VALUE_PTR) malloc(lmerCount * VALUE_SIZE);
	h_lmerOffsets = (unsigned int *) malloc(lmerCount * sizeof(unsigned int));
	h_TK = (KEY_PTR) malloc(bucketCount * BUCKET_KEY_SIZE);
	h_TV = (VALUE_PTR) malloc(bucketCount * BUCKET_VALUE_SIZE);
	h_bucketSeed = (unsigned int *) malloc(bucketCount * sizeof(unsigned int));
	h_bucketSeed = (unsigned int *) malloc(bucketCount * sizeof(unsigned int));
	h_l = (unsigned int *) malloc(ecount * sizeof(unsigned int));
	h_e = (unsigned int *) malloc(ecount * sizeof(unsigned int));
	h_ee = (EulerEdge *) malloc(ecount * sizeof(EulerEdge));
	h_loffsets = (unsigned int *) malloc(kmerCount * 4 * sizeof(unsigned int));
	h_eoffsets = (unsigned int *) malloc(kmerCount * 4 * sizeof(unsigned int));
	h_lcount = (unsigned int *) malloc(kmerCount * 4 * sizeof(unsigned int));
	h_ecount = (unsigned int *) malloc(kmerCount * 4 * sizeof(unsigned int));

	checkCudaErrors(
			cudaMemcpy(h_lmerKeys, d_lmerKeys, lmerCount * KEY_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_lmerValues, d_lmerValues, lmerCount * VALUE_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_lmerOffsets, d_lmerOffsets,
					lmerCount * sizeof(unsigned int), cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_TK, d_TK, bucketCount * BUCKET_KEY_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_TV, d_TV, bucketCount * BUCKET_VALUE_SIZE,
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_bucketSeed, d_bucketSeed,
					bucketCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_ee, d_ee, ecount * sizeof(EulerEdge),
					cudaMemcpyDeviceToHost));

	checkCudaErrors(
			cudaMemcpy(h_loffsets, d_loffsets,
					kmerCount * 4 * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_eoffsets, d_eoffsets,
					kmerCount * 4 * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_lcount, d_lcount, kmerCount * 4 * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_ecount, d_ecount, kmerCount * 4 * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_l, d_l, ecount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_e, d_e, ecount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));

	printf(".....diff....\n");
	for (unsigned int j = 0; j < (4 * kmerCount) - 1; j++) {
		if (h_lcount[j] != h_loffsets[j + 1] - h_loffsets[j]) {
			printf(" lcount mismatch j:[%u] lcount:[%u] diff:[%u]\n", j,
					h_lcount[j], h_loffsets[j + 1] - h_loffsets[j]);
		}
		if (h_ecount[j] != h_eoffsets[j + 1] - h_eoffsets[j]) {
			printf(" ecount mismatch j:[%u] ecount:[%u] diff:[%u]\n", j,
					h_ecount[j], h_eoffsets[j + 1] - h_eoffsets[j]);
		}
	}

	/*	for(unsigned int k=0;k<4*kmerCount;k++){
	 printf("[%u]:  loffset[%u] ,lcount[%u] ,eoffset[%u], ecount[%u]\n",k,h_loffsets[k],h_lcount[k],h_eoffsets[k],h_ecount[k]);
	 }
	 */
	for (int i = 0; i < lmerCount; i++) {
		setupEdgesHost(h_lmerKeys, h_lmerValues, h_lmerOffsets, lmerCount, h_TK,
				h_TV, h_bucketSeed, bucketCount, h_l, h_e, h_ee, h_loffsets,
				h_eoffsets, validBitMask, i);
	}

	free(h_lmerKeys);
	free(h_lmerValues);
	free(h_lmerOffsets);
	free(h_TK);
	free(h_TV);
	free(h_bucketSeed);
	free(h_l);
	free(h_e);
	free(h_ee);
	free(h_loffsets);
	free(h_eoffsets);
	free(h_lcount);
	free(h_ecount);

}

void verifyleOffsets(unsigned int * d_lOffsets, unsigned int * d_lcount,
		unsigned int * d_eOffsets, unsigned int * d_ecount, unsigned int length,
		unsigned int ecount) {

	unsigned int * h_lOffsets;
	unsigned int * h_eOffsets;
	unsigned int * h_lcount;
	unsigned int * h_ecount;

	h_lOffsets = (unsigned int*) malloc(length * sizeof(unsigned int));
	h_eOffsets = (unsigned int *) malloc(length * sizeof(unsigned int));
	h_lcount = (unsigned int*) malloc(length * sizeof(unsigned int));
	h_ecount = (unsigned int *) malloc(length * sizeof(unsigned int));

	checkCudaErrors(
			cudaMemcpy(h_lOffsets, d_lOffsets, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_eOffsets, d_eOffsets, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_lcount, d_lcount, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_ecount, d_ecount, length * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	unsigned int esum = 0;
	unsigned int lsum = 0;
	for (unsigned int t = 0; t < length; t++) {
		esum += h_ecount[t];
		lsum += h_lcount[t];
	}
	printf("esum : %u , lsum : %u \n", esum, lsum);
	unsigned int incorrectTotal = 0;
	for (unsigned int i = 0; i < length; i++) {
		if (h_lOffsets[i] > ecount || h_lOffsets[i] + h_lcount[i] > ecount) {
			incorrectTotal++;
			printf("incorrect l @ %u,  value %u\n",i,h_lOffsets[i]);
		}
		if (h_eOffsets[i] > ecount || h_eOffsets[i] + h_ecount[i] > ecount) {
			incorrectTotal++;
			printf("incorrect e @ %u,  value %u\n",i,h_eOffsets[i]);
		}
	}
	free(h_lOffsets);
	free(h_eOffsets);
	free(h_lcount);
	free(h_ecount);
}
extern "C" void constructDebruijnGraphGold(unsigned int * ecount,
		KEY_PTR h_lmerKeys, //in lmer keys
		VALUE_PTR h_lmerValues, //in lmer values
		unsigned int lmerCount, //in total lmers
		KEY_PTR h_kmerKeys, //in
		unsigned long kmerCount, //in  total kmers
		unsigned int l, //in k
		KEY_PTR h_TK, VALUE_PTR h_TV, unsigned int * h_bucketSeed,
		unsigned int bucketCount, EulerVertex ** h_ev, //out
		unsigned int ** h_l, //out
		unsigned int ** h_e, //out
		EulerEdge ** h_ee //out
		) { //out

	dim3 grid;
	dim3 block;

	unsigned int * h_lcount;
	unsigned int * h_lstart;
	unsigned int * h_ecount;
	unsigned int * h_estart;
	unsigned int * h_lmerOffsets;

	unsigned int memsize;
	KEY_T validBitMask = 0;
	//unsigned int timerGPU = 0;	
	unsigned int k = l - 1;

	//cutilCheckError(cutCreateTimer(&timerGPU));

	memsize = (kmerCount) * sizeof(unsigned int) * 4; // 4-tuple for each kmer

	h_lcount = (unsigned int *) malloc(memsize);
	h_lstart = (unsigned int *) malloc(memsize);
	h_estart = (unsigned int *) malloc(memsize);
	h_ecount = (unsigned int *) malloc(memsize);
	h_lmerOffsets = (unsigned int*) malloc(lmerCount * VALUE_SIZE);

	for (unsigned int i = 0; i < k * 2; i++) {
		validBitMask = (validBitMask << 1) | 1;
	}

	unsigned int bitCount = 0;
	KEY_T bit = 1;
	while (bit != 0) {
		if (bit & validBitMask)
			bitCount++;
		bit = bit << 1;
	}
	memset(h_lcount, 0, sizeof(unsigned int) * 4 * kmerCount);
	memset(h_ecount, 0, sizeof(unsigned int) * 4 * kmerCount);

	//verifyDebruijnCountHost(d_lmerKeys,d_lmerValues,lmerCount,d_TK,d_TV,d_bucketSeed,bucketCount,d_lcount,d_ecount,validBitMask,kmerCount);
	for (unsigned int tid = 0; tid < lmerCount; tid++) {
		debruijnCountHost(h_lmerKeys, h_lmerValues, lmerCount, h_TK, h_TV,
				h_bucketSeed, bucketCount, h_lcount, h_ecount, validBitMask,
				bitCount, tid);
	}

	/* we need to perform pre-fix scan on , lcount, ecount, lmerValues,
	 * lcount and ecount has equal number of elements ,4*kmercount
	 * lmer has lmerCount elements, choose whichever is larger
	 */

	memset(h_lstart, 0, sizeof(unsigned int) * 4 * kmerCount);
	memset(h_estart, 0, sizeof(unsigned int) * 4 * kmerCount);
	memset(h_lmerOffsets, 0, sizeof(unsigned int) * lmerCount);

	prefixScan(h_lstart, h_lcount, 4 * kmerCount, false);
	prefixScan(h_estart, h_ecount, 4 * kmerCount, false);

	prefixScan(h_lmerOffsets, h_lmerValues, lmerCount, false);

	/*
	 unsigned int buffer[2];
	 readData(buffer,d_lmerOffsets+lmerCount-1,1,sizeof(unsigned int));
	 readData(buffer+1,d_lmerValues+lmerCount-1,1,sizeof(unsigned int));
	 *ecount=buffer[0]+buffer[1];
	 */
	*ecount = h_lmerOffsets[lmerCount - 1] + h_lmerValues[lmerCount - 1];

	*h_ev = (EulerVertex *) malloc(sizeof(EulerVertex) * (kmerCount));
	*h_l = (unsigned int *) malloc(sizeof(unsigned int) * (*ecount));
	*h_e = (unsigned int *) malloc(sizeof(unsigned int) * (*ecount));
	*h_ee = (EulerEdge *) malloc(sizeof(EulerEdge) * (*ecount));
	memset(*h_e, 0, sizeof(unsigned int) * (*ecount));
	memset(*h_l, 0, sizeof(unsigned int) * (*ecount));

//	getOptimalLaunchConfiguration(kmerCount,&grid,&block);
	for (unsigned int tid = 0; tid < kmerCount; tid++) {
		setupVerticesHost(h_kmerKeys, kmerCount, h_TK, h_TV, h_bucketSeed,
				bucketCount, *h_ev, h_lcount, h_lstart, h_ecount, h_estart,
				tid);
	}

	//getOptimalLaunchConfiguration(lmerCount,&grid,&block);
	for (unsigned int tid = 0; tid < lmerCount; tid++) {
		setupEdgesHost(h_lmerKeys, h_lmerValues, h_lmerOffsets, lmerCount, h_TK,
				h_TV, h_bucketSeed, bucketCount, *h_l, *h_e, *h_ee, h_lstart,
				h_estart, validBitMask, tid);
	}

	free(h_lmerOffsets);
	free(h_lcount);
	free(h_lstart);
	free(h_estart);
	free(h_ecount);

}
//extern "C" 
void constructDebruijnGraphDevice(unsigned int * ecount,
		KEY_PTR d_lmerKeys, //in lmer keys
		VALUE_PTR d_lmerValues, //in lmer values
		unsigned int lmerCount, //in total lmers
		KEY_PTR d_kmerKeys, //in
		unsigned long kmerCount, //in  total kmers
		unsigned int l, //in k
		KEY_PTR d_TK, 
		VALUE_PTR d_TV, 
		unsigned int * d_bucketSeed,
		unsigned int bucketCount, 
		EulerVertex ** d_ev, //out
		unsigned int ** d_l, //out
		unsigned int ** d_e, //out
		EulerEdge ** d_ee //out
		) { //out

	dim3 grid;
	dim3 block;

	unsigned int * d_lcount;
	unsigned int * d_lstart;
	unsigned int * d_ecount;
	unsigned int * d_estart;
	unsigned int * d_lmerOffsets;

	unsigned int mem_size;
	KEY_T validBitMask = 0;
	//unsigned int timerGPU = 0;	
	unsigned int k = l - 1;

	//cutilCheckError(cutCreateTimer(&timerGPU));

	mem_size = (kmerCount) * sizeof(unsigned int) * 4; // 4-tuple for each kmer

	allocateMemory((void**) &d_lcount, mem_size);
	allocateMemory((void**) &d_lstart, mem_size);
	allocateMemory((void**) &d_estart, mem_size);
	allocateMemory((void**) &d_ecount, mem_size);
	allocateMemory((void**) &d_lmerOffsets, lmerCount * VALUE_SIZE);

	for (unsigned int i = 0; i < k * 2; i++) {
		validBitMask = (validBitMask << 1) | 1;
	}

	logMessage(LOG_LVL_DETAIL,"deb bit mask %lu\n",validBitMask);
	logMessage(LOG_LVL_DETAIL, "kernel: debruijnCount");
	getOptimalLaunchConfiguration(lmerCount, &grid, &block);
	//verifyDebruijnCountHost(d_lmerKeys,d_lmerValues,lmerCount,d_TK,d_TV,d_bucketSeed,bucketCount,d_lcount,d_ecount,validBitMask,kmerCount);
	debruijnCount<<<grid,block>>>(d_lmerKeys,d_lmerValues,lmerCount,d_TK,d_TV,d_bucketSeed,bucketCount,d_lcount,d_ecount,validBitMask);
	CheckCUDAError();
	//verifyDebruijnCountHost(d_lmerKeys,d_lmerValues,lmerCount,d_TK,d_TV,d_bucketSeed,bucketCount,d_lcount,d_ecount,validBitMask,kmerCount);

	/* we need to perform pre-fix scan on , lcount, ecount, lmerValues,
	 * lcount and ecount has equal number of elements ,4*kmercount
	 * lmer has lmerCount elements, choose whichever is larger
	 */

//	unsigned int maxElements=(lmerCount>4*kmerCount)?lmerCount:4*kmerCount;
	CUDPPConfiguration configKmer;
	configKmer.op = CUDPP_ADD;
	configKmer.datatype = CUDPP_UINT;
	configKmer.algorithm = CUDPP_SCAN;
	configKmer.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;

	CUDPPHandle scanplanKmer = 0;
	cudppPlan(&scanplanKmer, configKmer, 4 * kmerCount, 1, 0);
	CheckCUDAError();

	cudppScan(scanplanKmer, d_lstart, d_lcount, 4 * kmerCount);
	cudppScan(scanplanKmer, d_estart, d_ecount, 4 * kmerCount);
	cudppDestroyPlan(scanplanKmer);

	CUDPPConfiguration configLmer;
	configLmer.op = CUDPP_ADD;
	configLmer.datatype = CUDPP_UINT;
	configLmer.algorithm = CUDPP_SCAN;
	configLmer.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_EXCLUSIVE;

	CUDPPHandle scanplanLmer = 0;
	CUDPPResult result = cudppPlan(&scanplanLmer, configLmer, lmerCount, 1, 0);
	CheckCUDAError();

	cudppScan(scanplanLmer, d_lmerOffsets, d_lmerValues, lmerCount);
	cudppDestroyPlan(scanplanLmer);

	//validatePrefixScan(d_lstart,d_lcount,4*kmerCount,false);
	//validatePrefixScan(d_estart,d_ecount,4*kmerCount,false);
	//validatePrefixScan(d_lmerOffsets,d_lmerValues,lmerCount,false);

	unsigned int buffer[2];
	readData(buffer, d_lmerOffsets + lmerCount - 1, 1, sizeof(unsigned int));
	readData(buffer + 1, d_lmerValues + lmerCount - 1, 1, sizeof(unsigned int));
	*ecount = buffer[0] + buffer[1];

	logMessage(LOG_LVL_MSG, "debruijn vertex count:%d \ndebruijn edge count:%d",
			kmerCount, *ecount);

	allocateMemory((void**) d_ev, sizeof(EulerVertex) * (kmerCount));
	allocateMemory((void**) d_l, sizeof(unsigned int) * (*ecount));
	allocateMemory((void**) d_e, sizeof(unsigned int) * (*ecount));
	allocateMemory((void**) d_ee, sizeof(EulerEdge) * (*ecount));
	CheckCUDAError();
	cudaMemset(*d_e, 0, sizeof(unsigned int) * (*ecount));
	cudaMemset(*d_l, 0, sizeof(unsigned int) * (*ecount));
	CheckCUDAError();

	logMessage(LOG_LVL_DETAIL, "kernel: setupVertices");
	getOptimalLaunchConfiguration(kmerCount, &grid, &block);
	//cutilCheckError(cutStartTimer(timerGPU));
	setupVertices<<<grid,block>>>(d_kmerKeys,kmerCount,d_TK,d_TV,d_bucketSeed,bucketCount,*d_ev,d_lcount,d_lstart,d_ecount,d_estart);
	CheckCUDAError();

	///*DEBUG*/verifyleOffsets(d_lstart,d_lcount,d_estart,d_ecount,4*kmerCount,*ecount);

	getOptimalLaunchConfiguration(lmerCount, &grid, &block);
	//verifySetupEdges(d_lmerKeys,d_lmerValues,d_lmerOffsets,lmerCount, d_TK,d_TV,d_bucketSeed,bucketCount,*d_l,*d_e,*d_ee,d_lcount,d_lstart,d_ecount,d_estart,kmerCount,*ecount,validBitMask);
	logMessage(LOG_LVL_DETAIL,"kernel: setupEdges");
	setupEdges<<<grid,block>>>(d_lmerKeys,d_lmerValues,d_lmerOffsets,lmerCount, d_TK,d_TV,d_bucketSeed,bucketCount,*d_l,*d_e,*d_ee,d_lstart,d_estart,validBitMask);

	CheckCUDAError();

	//cutilCheckError(cutStopTimer(timerGPU));
	//logMessage(LOG_LVL_MSG,"CPU Time : %f",cutGetTimerValue(timerGPU));

	//constructDebruijnGold( d_idata, d_icount, kmerCount,kmerLength,totalVertices,validBitMask);
	//printDebruijnGraph(*d_ev, kmerCount, *d_l, *d_e, *d_ee, *ecount, k, 0); // may not need it

	//printDebruijnGraph(*d_ev,kmerCount,*d_l,*d_e,*d_ee,*ecount,k,1);
	//printData(*d_ev,*vcount);
	//printData(*d_ee,*ecount);

	//cutilCheckError(cutDeleteTimer(timerGPU));

	deallocateMemory(d_lmerOffsets);
	deallocateMemory(d_lcount);
	deallocateMemory(d_lstart);
	deallocateMemory(d_estart);
	deallocateMemory(d_ecount);

}

