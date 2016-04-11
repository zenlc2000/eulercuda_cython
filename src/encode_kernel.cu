__device__ __constant__ KEY_T lmerMask[] ={
    0x0000000000000003, 0x000000000000000F, 0x000000000000003F, 0x00000000000000FF, // 0   1   2   3
    0x00000000000003FF, 0x0000000000000FFF, 0x0000000000003FFF, 0x000000000000FFFF, // 4   5   6   7
    0x000000000003FFFF, 0x00000000000FFFFF, 0x00000000003FFFFF, 0x0000000000FFFFFF, // 8   9   10  11
    0x0000000003FFFFFF, 0x000000000FFFFFFF, 0x000000003FFFFFFF, 0x00000000FFFFFFFF, // 12  13  14  15
    0x00000003FFFFFFFF, 0x0000000FFFFFFFFF, 0x0000003FFFFFFFFF, 0x000000FFFFFFFFFF, // 16  17  18  19
    0x000003FFFFFFFFFF, 0x00000FFFFFFFFFFF, 0x00003FFFFFFFFFFF, 0x0000FFFFFFFFFFFF, // 20  21  22  23
    0x0003FFFFFFFFFFFF, 0x000FFFFFFFFFFFFF, 0x003FFFFFFFFFFFFF, 0x00FFFFFFFFFFFFFF, // 24  25  26  27
    0x03FFFFFFFFFFFFFF, 0x0FFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF // 28  29  30  31
};

__device__ __constant__ unsigned char shifter[4] [4]=
{
		{0,0,0,0},
		{1,4,16,64},
		{2,8,32,128},
		{3,12,48,192},
};

__device__ __constant__ char  codeF[]={0,0,0,1,3,0,0,2};
__device__ __constant__ char  codeR[]={0,3,0,2,0,0,0,1};
__global__ void encodeLmerDevice(char  * buffer,
			//	const unsigned int buffSize,
			//	const unsigned int readLength,
				KEY_PTR lmers,
				const unsigned int lmerLength
				)
{
//    printf("in GPU");
	extern __shared__ char dnaRead[]; // MB: changed from 'read' to solve compile error
	const unsigned int tid=threadIdx.x;
	const unsigned int rOffset=(blockDim.x*blockDim.y*gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x) + (blockDim.x*threadIdx.y);
	KEY_T lmer=0;

	dnaRead[tid]=buffer[rOffset+tid];
	__syncthreads();

	for (unsigned int i = 0; i < 8; i++)    //calculate lmer
	{
	    lmer= (lmer<< 8) |	((KEY_T)(shifter[codeF[dnaRead[threadIdx.x+i*4]& 0x07]][3] |
							shifter[codeF[dnaRead[threadIdx.x+i*4+1]& 0x07]][2] |
							shifter[codeF[dnaRead[threadIdx.x+i*4+2]& 0x07]][1] |
							codeF[dnaRead[threadIdx.x+i*4+3] & 0x07]) ) ;
	}
	lmer = (lmer >> ((32 - lmerLength) << 1)) & lmerMask[lmerLength-1];
//    printf("%llu", lmer);
	lmers[rOffset+tid]=lmer;


}

__global__ void computeKmerDevice( 	KEY_PTR lmers,
				KEY_PTR pkmers,
				KEY_PTR skmers,
				KEY_T validBitMask
			){

	const unsigned int tid=(blockDim.x*blockDim.y*gridDim.x*blockIdx.y) +(blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
	KEY_T lmer;
	//fetch lmer
	lmer=lmers[tid];
	//find prefix
	pkmers[tid]=LMER_PREFIX(lmer,validBitMask);
	//find suffix
	skmers[tid] = LMER_SUFFIX(lmer,validBitMask);
}

__global__ void encodeLmerComplementDevice(	char  * buffer,
				const unsigned int buffSize,
				const unsigned int readLength,
				KEY_PTR lmers,
				const unsigned int lmerLength
				)
{

	extern __shared__ char dnaRead[];//have to fix it
	const unsigned int tid=threadIdx.x;
	const unsigned int rOffset=(blockDim.x*blockDim.y*gridDim.x*blockIdx.y) +(blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y);
	KEY_T lmer=0;
	KEY_T temp=0;

	dnaRead[tid]=buffer[rOffset+tid];
	__syncthreads();
	dnaRead[tid]=codeR[dnaRead[tid] & 0x07];
	__syncthreads();
	for(unsigned int i =0; i< lmerLength; i++)
	{
		temp=((KEY_T)dnaRead[(tid+i)%blockDim.x]);
		lmer = (temp<<(i<<1)) | lmer;
	}
	lmers[rOffset+tid]=lmer;

}