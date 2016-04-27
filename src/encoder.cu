#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define CUDPP_STATIC_LIB


// #include <cutil_inline.h>
#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_cuda.h" 
#include <stdio.h>      /* defines printf for tests */
#include <time.h>       /* defines time_t for timings in the test */
#include <math.h>

#ifdef linux
# include <endian.h>    /* attempt to define endianness */
#endif


#include "Kmer.h"
#include "Graph.h"
#include "cudpp.h"
#include "utils.h"
#include "common.h"
#include "encode_kernel.cu"

// idea ! have all the cuda code as macro to be used in cpu code as well cuda code

/*
one read per block
|R| threads
copy each Ri to shared mem
|lmers|=|R|
l-1 dummy enteris or l-1 thread stall/branch

max overhead 256 bytes per read.

1) use dummy enteris
2) stall extra threads.
*/
//A=65=41=0100-0001
//C=67=43=0100-0011
//T=84=54=0101-0100
//G=71=47=0100-1111
//0A0CT00G
//00013002
//0T0GA00C
//03020001
//__device__ __constant__ KEY_T lmerMask[] ={
//    0x0000000000000003, 0x000000000000000F, 0x000000000000003F, 0x00000000000000FF, // 0   1   2   3
//    0x00000000000003FF, 0x0000000000000FFF, 0x0000000000003FFF, 0x000000000000FFFF, // 4   5   6   7
//    0x000000000003FFFF, 0x00000000000FFFFF, 0x00000000003FFFFF, 0x0000000000FFFFFF, // 8   9   10  11
//    0x0000000003FFFFFF, 0x000000000FFFFFFF, 0x000000003FFFFFFF, 0x00000000FFFFFFFF, // 12  13  14  15
//    0x00000003FFFFFFFF, 0x0000000FFFFFFFFF, 0x0000003FFFFFFFFF, 0x000000FFFFFFFFFF, // 16  17  18  19
//    0x000003FFFFFFFFFF, 0x00000FFFFFFFFFFF, 0x00003FFFFFFFFFFF, 0x0000FFFFFFFFFFFF, // 20  21  22  23
//    0x0003FFFFFFFFFFFF, 0x000FFFFFFFFFFFFF, 0x003FFFFFFFFFFFFF, 0x00FFFFFFFFFFFFFF, // 24  25  26  27
//    0x03FFFFFFFFFFFFFF, 0x0FFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF // 28  29  30  31
//};
//
//__device__ __constant__ unsigned char shifter[4] [4]=
//{
//		{0,0,0,0},
//		{1,4,16,64},
//		{2,8,32,128},
//		{3,12,48,192},
//};
//
//__device__ __constant__ char  codeF[]={0,0,0,1,3,0,0,2};
//__device__ __constant__ char  codeR[]={0,3,0,2,0,0,0,1};
//__global__ void encodeLmerDevice(char  * buffer,
//			//	const unsigned int buffSize,
//			//	const unsigned int readLength,
//				KEY_PTR lmers,
//				const unsigned int lmerLength
//				)
//{
//    printf("%s", buffer);
//	extern __shared__ char dnaRead[]; // MB: changed from 'read' to solve compile error
//	const unsigned int tid=threadIdx.x;
//	const unsigned int rOffset=(blockDim.x*blockDim.y*gridDim.x*blockIdx.y) + (blockDim.x*blockDim.y*blockIdx.x) + (blockDim.x*threadIdx.y);
//	KEY_T lmer=0;
//
//	dnaRead[tid]=buffer[rOffset+tid];
//	__syncthreads();
//
//	for (unsigned int i = 0; i < 8; i++)    //calculate lmer
//	{
//	    lmer= (lmer<< 8) |	((KEY_T)(shifter[codeF[dnaRead[threadIdx.x+i*4]& 0x07]][3] |
//							shifter[codeF[dnaRead[threadIdx.x+i*4+1]& 0x07]][2] |
//							shifter[codeF[dnaRead[threadIdx.x+i*4+2]& 0x07]][1] |
//							codeF[dnaRead[threadIdx.x+i*4+3] & 0x07]) ) ;
//	}
//	lmer = (lmer >> ((32 - lmerLength) << 1)) & lmerMask[lmerLength-1];
//
//	lmers[rOffset+tid]=lmer;
//
//
//}
//__global__ void encodeLmerComplementDevice(	char  * buffer,
//				const unsigned int buffSize,
//				const unsigned int readLength,
//				KEY_PTR lmers,
//				const unsigned int lmerLength
//				){
//
//	extern __shared__ char dnaRead[];//have to fix it
//	const unsigned int tid=threadIdx.x;
//	const unsigned int rOffset=(blockDim.x*blockDim.y*gridDim.x*blockIdx.y) +(blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y);
//	KEY_T lmer=0;
//	KEY_T temp=0;
//
//	dnaRead[tid]=buffer[rOffset+tid];
//	__syncthreads();
//	dnaRead[tid]=codeR[dnaRead[tid] & 0x07];
//	__syncthreads();
//	for(unsigned int i =0; i< lmerLength; i++){
//		temp=((KEY_T)dnaRead[(tid+i)%blockDim.x]);
//		lmer = (temp<<(i<<1)) | lmer;
//	}
//	lmers[rOffset+tid]=lmer;
//
//}


//__global__ void computeKmerDevice( 	KEY_PTR lmers,
//				KEY_PTR pkmers,
//				KEY_PTR skmers,
//				KEY_T validBitMask
//			)
//{
//
//	const unsigned int tid=(blockDim.x*blockDim.y*gridDim.x*blockIdx.y) +(blockDim.x*blockDim.y*blockIdx.x)+(blockDim.x*threadIdx.y)+threadIdx.x;
//	KEY_T lmer;
//	//fetch lmer
//	lmer=lmers[tid];
//	//find prefix
//	pkmers[tid]=LMER_PREFIX(lmer,validBitMask);
//	//find suffix
//	skmers[tid] = LMER_SUFFIX(lmer,validBitMask);
//}

extern "C"
void encodeLmer(
		char * d_buffer,
		const unsigned int bufferSize,
		const unsigned int readLength,
		KEY_PTR d_lmers,
		const unsigned int lmerLength,
		const unsigned int entriesCount
		)
{
	dim3 grid, block;
	char * d_reads = NULL;

	cudaError_t err = cudaMalloc((void**) &d_reads, bufferSize);
    err = cudaMemcpy(d_reads, d_buffer, bufferSize,cudaMemcpyHostToDevice);
    if (err != 0)
        printf("cudaMemcpy error encodeLmer: %d\n",err);
	getOptimalLaunchConfigCustomized(entriesCount,&grid,&block,readLength);

	encodeLmerDevice<<<grid,block,readLength+31>>>(d_reads,d_lmers,lmerLength);

	CheckCUDAError();
		

}
extern "C"
void encodeLmerComplement(
		char * d_buffer,
		const unsigned int bufferSize,
		const unsigned int readLength,
		KEY_PTR d_lmers,
		const unsigned int lmerLength,
		const unsigned int entriesCount
		)
{
	dim3 grid, block;
	char * d_reads = NULL;

	cudaError_t err = cudaMalloc((void**) &d_reads, bufferSize);
    err = cudaMemcpy(d_reads, d_buffer, bufferSize,cudaMemcpyHostToDevice);
    if (err != 0)
        printf("cudaMemcpy error encodeLmerComplement: %d\n",err);
	getOptimalLaunchConfigCustomized(entriesCount,&grid,&block,readLength);		
	encodeLmerComplementDevice<<<grid,block,readLength>>>(d_reads,bufferSize,readLength,d_lmers,lmerLength);
	CheckCUDAError();
}

extern "C"
void computeKmer(	KEY_PTR d_lmers,
			KEY_PTR d_pkmers,
			KEY_PTR d_skmers,
			KEY_PTR h_lmers,
			KEY_PTR h_pkmers,
			KEY_PTR h_skmers,
			KEY_T validBitMask,
			const unsigned int readLength,
			const unsigned int entriesCount
		)
{
	dim3 grid, block;
	unsigned int ebSize = entriesCount * sizeof(KEY_T);
	getOptimalLaunchConfigCustomized(entriesCount,&grid,&block,readLength);
	computeKmerDevice<<<grid,block>>>(d_lmers,d_pkmers,d_skmers,validBitMask);

    cudaError_t err1 = cudaMemcpy(h_lmers, d_lmers, ebSize, cudaMemcpyDeviceToHost);
    cudaError_t err2 = cudaMemcpy(h_pkmers, d_pkmers, ebSize,cudaMemcpyDeviceToHost);
    cudaError_t err3 = cudaMemcpy(h_skmers, d_skmers, ebSize,cudaMemcpyDeviceToHost);

    if ((err1 != 0) || (err2 != 0) || (err3 != 0))
        printf("err1 = %d err2 = %d err3 = %d\n",err1, err2, err3);

	CheckCUDAError();
			
}

//extern "C"
//void computeKmerComplement(	KEY_PTR d_lmers,
//			KEY_PTR d_pkmers,
//			KEY_PTR d_skmers,
//			KEY_PTR h_lmers,
//			KEY_PTR h_pkmers,
//			KEY_PTR h_skmers,
//			KEY_T validBitMask,
//			const unsigned int readLength,
//			const unsigned int entriesCount
//		)
//{
//	dim3 grid, block;
//	unsigned int ebSize = entriesCount * sizeof(KEY_T);
//	getOptimalLaunchConfigCustomized(entriesCount,&grid,&block,readLength);
//	computeKmerDevice<<<grid,block>>>(d_lmers,d_pkmers,d_skmers,validBitMask);
//
//    cudaError_t err1 = cudaMemcpy(h_lmers, d_lmers, ebSize, cudaMemcpyDeviceToHost);
//    cudaError_t err2 = cudaMemcpy(h_pkmers, d_pkmers, ebSize,cudaMemcpyDeviceToHost);
//    cudaError_t err3 = cudaMemcpy(h_skmers, d_skmers, ebSize,cudaMemcpyDeviceToHost);
//
//    if ((err1 != 0) || (err2 != 0) || (err3 != 0))
//        printf("err1 = %d err2 = %d err3 = %d\n",err1, err2, err3);
//
//	CheckCUDAError();
//
//}
