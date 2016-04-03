#include "/Volumes/Macintosh HD/Developer/NVIDIA/CUDA-7.5/samples/common/inc/helper_cuda.h" // lib above replaced w/this one at CUDA 5.0
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdarg.h>
#include "Graph.h"
#include "common.h"
 #include <execinfo.h>
#define DEFAULT_BLOCK_SIZE 512

#ifdef EULER_NDEBUG
#define DEBUG_UTIL_CPP(x)
#else
#define DEBUG_UTIL_CPP(x) x
#endif
#define DEBUG_CALL(x) DEBUG_UTIL_CPP(x)

#define LOG_ENABLED 
typedef struct MemNode {
	void * d_ptr;
	MemNode * next;
} MemNode;

/********** Globals *****/
MemNode * head;
int LOG_LEVEL = 0;
int blockSize = DEFAULT_BLOCK_SIZE;
FILE * logFile;

/* function for loggin purpose**/


extern "C" void logMessage(int logLevel, const char * format, ...) {

#if defined LOG_ENABLED
		if( logLevel<=LOG_LEVEL) {
			va_list ap;
			int r;
			va_start (ap, format);
			fprintf(stderr,"<!-- ");
			//printf ("P[%d]:",rank);
		r = vfprintf (stderr,format, ap);
		va_end (ap);
		fprintf(stderr," --> \n");
		//	fprintf(stderr,"\n");

	}
#endif
	}
extern "C" void logMessageNL(int logLevel, const char * format, ...) { //this is actually logMessage No LineFeed

#if defined LOG_ENABLED
		if( logLevel<=LOG_LEVEL) {
			va_list ap;
			int r;
			va_start (ap, format);
			fprintf(stderr,"<!-- ");
			//printf ("P[%d]:",rank);
		r = vfprintf (stderr,format, ap);
		fprintf(stderr," -->");
		va_end (ap);

	}
#endif
	}

void logCaller(){
	 void *buffer[6];
	 char **strings;

	 unsigned int nptrs = backtrace(buffer, 5);
	strings = backtrace_symbols(buffer, nptrs);
	if (strings == NULL) {
		perror("backtrace_symbols");

	}else {
		logMessage(LOG_LVL_DEBUG, "Caller:%s  , Callee:%s", strings[2],strings[1]);
		 free(strings);
	}


}
extern "C" void CheckCUDAError() {
	cudaThreadSynchronize();
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err) {
		logMessage(LOG_LVL_ERROR, "%s\n", cudaGetErrorString(err));
	}else{
		logMessage(LOG_LVL_DETAIL, "CUDA call succeeded");
	}
}

extern "C" void initMemList() {
	head = (MemNode *) malloc(sizeof(MemNode));
	head->d_ptr = NULL;
	head->next = NULL;
}

extern "C" void addMemNode(void * d_ptr) {
	MemNode * tmp = (MemNode*) malloc(sizeof(MemNode));
	tmp->d_ptr = d_ptr;
	tmp->next = head->next;
	head->next = tmp;
}

extern "C" void freeMemList() {
	MemNode * tmp;
	while (head->next != NULL) {
		tmp = head->next;
		head->next = head->next->next;
		cudaFree(tmp->d_ptr);
		free(tmp);
		CheckCUDAError();
	}
}
extern "C" void cleanupMemList() {
	freeMemList();
	free(head);
}
extern "C" void printData(unsigned int * d_buffer, int length) {

	unsigned int *h_buffer = (unsigned int *) malloc(length * sizeof(int));
	checkCudaErrors(
			cudaMemcpy(h_buffer, d_buffer, length * sizeof(int),
					cudaMemcpyDeviceToHost));
	for (int i = 0; i < length; i++) {
		logMessageNL(LOG_LVL_MSG, "@[%d] %u \n", i, h_buffer[i]);
	}
	logMessage(LOG_LVL_DETAIL, "");
	free(h_buffer);

}
/*
 extern "C"
 void printData(void * d_buffer,int length,int width){

 unsigned int *h_buffer= (unsigned int *)malloc(length* sizeof(int)*width);
 checkCudaErrors( cudaMemcpy( h_buffer, d_buffer, length* sizeof(int)*width,cudaMemcpyDeviceToHost) );
 for( int i=0;i<length;i++){
 for(int j=0;j< width;j++) {
 printf("[%d]:%u ",j,h_buffer[i*length+j]);
 }
 printf("\n");
 }
 printf("\n");
 free(h_buffer);

 }*/
/*
 extern "C"
 void getOptimalLaunchConfiguration(unsigned int threadCount,unsigned int * gridx,unsigned int * gridy,unsigned int * threads){

 *threads=32;
 *gridy=threadCount/(*threads);
 if(threadCount%(*threads) >0) (*gridy)++;
 (*gridx)=(*gridy)/65535+1;
 (*gridy)=(*gridy)%65535;
 }
 */

extern "C" void setBlockSize(int newSize) {
	blockSize = newSize;
}

extern "C" void getOptimalLaunchConfigCustomized(unsigned int threadCount,
		dim3 * grid, dim3 * block, unsigned int threadPerBlock) {

	*block = make_uint3(threadPerBlock, 1, 1);
	*grid = make_uint3(1, 1, 1);

	/*grid->y=threadCount/(block->x);
	 if(threadCount%(block->x) >0) grid->y++;
	 grid->x=grid->y / 65535 +1;
	 grid->y=grid->y % 65535;
	 grid->z=1;*/
	if (threadCount > block->x) {
		grid->y = threadCount / (block->x);
		if (threadCount % (block->x) > 0)
			grid->y++;
		grid->x = grid->y / 65535 + 1;
		grid->y = (grid->y > 65535 ) ? 65535 : grid->y;
		grid->z = 1;
	}
}
extern "C" void getOptimalLaunchConfiguration(unsigned int threadCount,
		dim3 * grid, dim3 * block) {
	getOptimalLaunchConfigCustomized(threadCount, grid, block, blockSize);
}

extern "C" void readData(void * h_out, void * d_in, int length, int width) {
	checkCudaErrors(
			cudaMemcpy(h_out, d_in, length * width, cudaMemcpyDeviceToHost));
}

extern "C" void allocateMemory(void ** d_buffer, unsigned int memSize) {
	size_t  free, total;
	DEBUG_CALL(logCaller());
	cuMemGetInfo(&free, &total);
	logMessage(LOG_LVL_DEBUG, "\t\t\tMemory Requested %d bytes", memSize);
	logMessage(
			LOG_LVL_DEBUG,
			"\t\t\tMemory Status Before Alloc :: total:[%u]\t used:[%u]\t free:[%u]",
			total, total - free, free);
	if (total - free != 0 && free < memSize) {
		freeMemList();
	}
	checkCudaErrors( cudaMalloc((void**) d_buffer, memSize));
	;
	CheckCUDAError();
	logMessage(LOG_LVL_DEBUG, "\t\t\tmemory address %u", *d_buffer);
	cudaMemset(*d_buffer, 0, memSize);
	cuMemGetInfo(&free, &total);
	logMessage(
			LOG_LVL_DEBUG,
			"\t\t\tMemory Status After Alloc :: total:[%u]\t used:[%u]\t free:[%u]",
			total, total - free, free);
}

extern "C" void deallocateMemory(void * d_buffer) {

	size_t freeBefore, total, freeAfter;
	DEBUG_CALL(logCaller());
	logMessage(LOG_LVL_DEBUG, "\t\t\tReleasing Memory");

	cuMemGetInfo(&freeBefore, &total);
	logMessage(
			LOG_LVL_DEBUG,
			"\t\t\tMemory Status Before Releasing :: total:[%u]\t used:[%u]\t free:[%u]",
			total, total - freeBefore, freeBefore);
	//cudaFree(d_buffer) ;//
	checkCudaErrors(cudaFree(d_buffer));
	//addMemNode(d_buffer);
	logMessage(LOG_LVL_DEBUG, "\t\t\t adding memory to junk list %u\n",
			d_buffer);
	CheckCUDAError();
	cuMemGetInfo(&freeAfter, &total);
	logMessage(
			LOG_LVL_DEBUG,
			"\t\t\tMemory Status After Releasing :: total:[%u]\t used:[%u]\t free:[%u]\t",
			total, total - freeAfter, freeAfter);
	logMessage(LOG_LVL_DEBUG, "\t\t\tMemory Freed %d bytes\n",
			freeAfter - freeBefore);

}
char _translate(int i) {
	if (i == 0)
		return 'A';
	if (i == 1)
		return 'C';
	if (i == 2)
		return 'G';
	if (i == 3)
		return 'T';
	return '.';
}
void _getString(char * kmer, int length, unsigned int value) {

	unsigned int currentValue = value;
	for (int i = 1; i <= length; i++) {
		kmer[length - i] = _translate(currentValue % 4);
		currentValue = currentValue / 4;
	}
}

extern "C" void printDebruijnGraphLongFmt(EulerVertex * d_ev, int vertexCount,
		unsigned int * d_l, unsigned int * d_e, EulerEdge * d_ee, int edgeCount,
		unsigned int kmerLength) {
	EulerVertex *h_ev = (EulerVertex *) malloc(
			vertexCount * sizeof(EulerVertex));
	EulerEdge *h_ee = (EulerEdge *) malloc(edgeCount * sizeof(EulerEdge));
	unsigned int * h_l = (unsigned int *) malloc(
			edgeCount * sizeof(unsigned int));
	unsigned int * h_e = (unsigned int *) malloc(
			edgeCount * sizeof(unsigned int));
	checkCudaErrors(
			cudaMemcpy(h_ev, d_ev, vertexCount * sizeof(EulerVertex),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_ee, d_ee, edgeCount * sizeof(EulerEdge),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_e, d_e, edgeCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_l, d_l, edgeCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));

	for (int i = 0; i < vertexCount; i++) {
		logMessage(LOG_LVL_DEBUG,
				"$[%d]:vid[%lu], ep[%u], eend[%u],  lp[%u], lend[%u]", i,
				(unsigned long) h_ev[i].vid, h_ev[i].ep, h_ev[i].ecount,
				h_ev[i].lp, h_ev[i].lcount);
		logMessageNL(LOG_LVL_DEBUG, "$e: [");
		for (unsigned int j = 0; j < h_ev[i].ecount; j++) {
			logMessageNL(LOG_LVL_DEBUG, " %d", h_e[h_ev[i].ep + j]);
		}
		logMessage(LOG_LVL_DEBUG, "]");
		logMessageNL(LOG_LVL_DEBUG, "$l: [");
		for (unsigned int j = 0; j < h_ev[i].lcount; j++) {
			logMessageNL(LOG_LVL_DEBUG, " %d", h_l[h_ev[i].lp + j]);
		}
		logMessage(LOG_LVL_DEBUG, "]");

	}

	logMessageNL(LOG_LVL_DEBUG, "$e: [");
	for (int i = 0; i < edgeCount; i++) {
		logMessageNL(LOG_LVL_DEBUG, " %u", h_e[i]);
	}
	logMessage(LOG_LVL_DEBUG, "]");
	logMessageNL(LOG_LVL_DEBUG, "$l: [");
	for (int i = 0; i < edgeCount; i++) {
		logMessageNL(LOG_LVL_DEBUG, " %u", h_l[i]);
	}
	logMessage(LOG_LVL_DEBUG, "]");

	logMessage(LOG_LVL_DEBUG, "$edges...");
	for (int i = 0; i < edgeCount; i++) {
		logMessage(LOG_LVL_DEBUG, "$[%d]: eid[%u], v1[%u], v2[%u], s[%u]\n", i,
				h_ee[i].eid, h_ee[i].v1, h_ee[i].v2, h_ee[i].s);
	}

	free(h_l);
	free(h_e);
	free(h_ev);
	free(h_ee);
}
extern "C" void printDebruijnGraphVizFmt(EulerVertex * d_ev, int vertexCount,
		unsigned int * d_l, unsigned int * d_e, EulerEdge * d_ee, int edgeCount,
		unsigned int kmerLength) {

	EulerVertex *h_ev = (EulerVertex *) malloc(
			vertexCount * sizeof(EulerVertex));
	EulerEdge *h_ee = (EulerEdge *) malloc(edgeCount * sizeof(EulerEdge));
	unsigned int * h_l = (unsigned int *) malloc(
			edgeCount * sizeof(unsigned int));
	unsigned int * h_e = (unsigned int *) malloc(
			edgeCount * sizeof(unsigned int));
	checkCudaErrors(
			cudaMemcpy(h_ev, d_ev, vertexCount * sizeof(EulerVertex),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_ee, d_ee, edgeCount * sizeof(EulerEdge),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_e, d_e, edgeCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));
	checkCudaErrors(
			cudaMemcpy(h_l, d_l, edgeCount * sizeof(unsigned int),
					cudaMemcpyDeviceToHost));

	char * v1 = (char *) malloc(kmerLength * sizeof(char));
	v1[kmerLength - 1] = '\0';

	char * v2 = (char *) malloc(kmerLength * sizeof(char));
	v2[kmerLength - 1] = '\0';

	logMessage(LOG_LVL_DEBUG, "$digraph G{");
	for (int i = 0; i < vertexCount; i++) {
		//_getString(v2,kmerLength-1,h_ev[i].vid);
		_getString(v1, kmerLength - 1, h_ev[i].vid);
		for (unsigned int j = 0; j < h_ev[i].ecount; j++) {
			//_getString(v1,kmerLength-1,h_ev[h_ee[h_e[h_ev[i].ep+j]].v1].vid);
			_getString(v2, kmerLength - 1,
					h_ev[h_ee[h_l[h_ev[i].lp + j]].v2].vid);
			//printf("$\t%s -> %s [ label= %s]\n",v1,v2,v2+kmerLength-2);
			logMessageNL(LOG_LVL_DEBUG, "$\t%s -> %s [ label= %s]\n", v1, v2,
					v2 + kmerLength - 2);
		}
	}
	logMessage(LOG_LVL_DEBUG, "$}");

	free(h_l);
	free(h_e);
	free(h_ev);
	free(h_ee);
	free(v1);
	free(v2);
}
extern "C" void printDebruijnGraph(EulerVertex * d_ev, unsigned int vertexCount,
		unsigned int * d_l, unsigned int * d_e, EulerEdge * d_ee,
		unsigned int edgeCount, unsigned int kmerLength, int format) {
	if (LOG_LEVEL >= LOG_LVL_DEBUG) {
		switch (format) {
		case 0:
			printDebruijnGraphLongFmt(d_ev, vertexCount, d_l, d_e, d_ee,
					edgeCount, kmerLength);
			break;
		case 1:
			printDebruijnGraphVizFmt(d_ev, vertexCount, d_l, d_e, d_ee,
					edgeCount, kmerLength);
			break;
		}
	}
}
extern "C" void printDeviceInfo(int argc, char** argv) {

	int devID;
	cudaDeviceProp deviceProp;

	// get number of SMs on this GPU
	checkCudaErrors(cudaGetDevice(&devID));
	checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

	if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
		// printf("[%s] - there is no device supporting CUDA.\n", sSDKsample);
		cudaThreadExit();
//		cutilExit(argc, argv);
	} else {
		printf("#> Device %d: \"%s\"\n", devID, deviceProp.name);
		printf("#> SM Capability %d.%d detected:\n", deviceProp.major,
				deviceProp.minor);
	}
	printf("#Device %d: \"%s\"\n", 0, deviceProp.name);
	printf("#  CUDA Capability Major revision number:         %d\n",
			deviceProp.major);
	printf("#  CUDA Capability Minor revision number:         %d\n",
			deviceProp.minor);
	printf("#  Total amount of global memory:                 %lu bytes\n",
			(unsigned long) deviceProp.totalGlobalMem);
	// #if CUDART_VERSION >= 2000
	printf("#  Number of multiprocessors:                     %d\n",
			deviceProp.multiProcessorCount);
	printf("#  Number of cores:                               %d\n",
			8 * deviceProp.multiProcessorCount);
	// #endif
	printf("#  Total amount of constant memory:               %lu bytes\n",
			(unsigned long) deviceProp.totalConstMem);
	printf("#  Total amount of shared memory per block:       %lu bytes\n",
			(unsigned long) deviceProp.sharedMemPerBlock);
	printf("#  Total number of registers available per block: %d\n",
			deviceProp.regsPerBlock);
	printf("#  Warp size:                                     %d\n",
			deviceProp.warpSize);
	printf("#  Maximum number of threads per block:           %d\n",
			deviceProp.maxThreadsPerBlock);
	printf("#  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
			deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
			deviceProp.maxThreadsDim[2]);
	printf("#  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
			deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
			deviceProp.maxGridSize[2]);
	printf("#  Maximum memory pitch:                          %lu bytes\n",
			(unsigned long) deviceProp.memPitch);
	printf("#  Texture alignment:                             %lu bytes\n",
			(unsigned long) deviceProp.textureAlignment);
	printf("#  Clock rate:                                    %.2f GHz\n",
			deviceProp.clockRate * 1e-6f);
	// #if CUDART_VERSION >= 2000
	printf("#  Concurrent copy and execution:                 %s\n",
			deviceProp.deviceOverlap ? "Yes" : "No");
	//  #endif
	//  #if CUDART_VERSION >= 2020
	printf("#  Run time limit on kernels:                     %s\n",
			deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
	printf("#  Integrated:                                    %s\n",
			deviceProp.integrated ? "Yes" : "No");
	printf("#  Support host page-locked memory mapping:       %s\n",
			deviceProp.canMapHostMemory ? "Yes" : "No");
	printf(
			"#  Compute mode:                                  %s\n",
			deviceProp.computeMode == cudaComputeModeDefault ?
					"Default (multiple host threads can use this device simultaneously)" :
			deviceProp.computeMode == cudaComputeModeExclusive ?
					"Exclusive (only one host thread at a time can use this device)" :
			deviceProp.computeMode == cudaComputeModeProhibited ?
					"Prohibited (no host thread can use this device)" :
					"Unknown");
	// #endif
}

extern "C"
void initDevice() {
	checkCudaErrors(cudaThreadExit());

	checkCudaErrors(cudaSetDevice(gpuGetMaxGflopsDeviceId() ));
	checkCudaErrors(cudaSetDeviceFlags(cudaDeviceMapHost));
	CheckCUDAError();
}

extern "C"
void resetDevice(){

}
