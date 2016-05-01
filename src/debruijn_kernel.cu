



/*
 *  This kernel works on each l-mer ,counting edges of the graph.
 */
__global__ void debruijnCount(
    KEY_PTR lmerKeys, /* lmer keys	*/
    VALUE_PTR lmerValues, /* lmer frequency */
    unsigned int lmerCount, /* total lmers */
    KEY_PTR TK, /* Keys' pointer for Hash table*/
    VALUE_PTR TV, /* Value pointer for Hash table*/
    unsigned int * bucketSeed, /* bucketSize: size of each bucket (it should be renamed to bucketSize)*/
    unsigned int bucketCount, /* total buckets */
    unsigned int * lcount, /* leaving edge count array : OUT */
    unsigned int * ecount, /* entering edge count array: OUT */
    KEY_T validBitMask /* bit mask for K length encoded bits*/
    )
{

	unsigned int tid = (blockDim.x * blockDim.y * gridDim.x * blockIdx.y)
			+ (blockDim.x * blockDim.y * blockIdx.x)
			+ (blockDim.x * threadIdx.y) + threadIdx.x;
	if (tid < lmerCount)
	{
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

/*
 *  This kernel works on a k-mer (l-1mer) which are vertices of the graph.
 */

__global__ void setupVertices(KEY_PTR kmerKeys, unsigned int kmerCount,
		KEY_PTR TK, VALUE_PTR TV, unsigned int * bucketSeed,
		unsigned int bucketCount, EulerVertex * ev, unsigned int * lcount,
		unsigned int * loffset, unsigned int * ecount, unsigned int * eoffset)
{
	unsigned int tid = (blockDim.x * blockDim.y * gridDim.x * blockIdx.y)
			+ (blockDim.x * blockDim.y * blockIdx.x)
			+ (blockDim.x * threadIdx.y) + threadIdx.x;
	if (tid < kmerCount)
	{
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
		const KEY_T validBitMask)
{

	unsigned int tid = (blockDim.x * blockDim.y * gridDim.x * blockIdx.y)
			+ (blockDim.x * blockDim.y * blockIdx.x)
			+ (blockDim.x * threadIdx.y) + threadIdx.x;
	if (tid < lmerCount)
	{
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
		for (unsigned int i = 0; i < lmerValue; i++)
		{

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