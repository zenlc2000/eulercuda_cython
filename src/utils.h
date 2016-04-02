#ifndef UTILS_H
#define UTILS_H

#include "Graph.h"

#define DEFAULT_BLOCK_SIZE 512

extern "C"
void printData(unsigned int * d_buffer,int length);

extern int LOG_LEVEL;

extern "C"
void logMessage(int logLevel,const char *  format,...);
extern "C"
void logMessageNL(int logLevel,const char *  format,...);

extern "C" 
void cleanupMemList();
extern "C"
void freeMemList();
extern "C"
void addMemNode(void * d_ptr);
extern "C" 
void initMemList();
/*
extern "C" 
void getOptimalLaunchConfiguration(unsigned int threadCount,unsigned int * gridx,unsigned int * gridy,unsigned int * threads);
*/

extern "C"
void setBlockSize(int newSize);

extern "C"
void getOptimalLaunchConfigCustomized(unsigned int threadCount,dim3 * grid,dim3 * block,unsigned int threadPerBlock);

extern "C" 
void getOptimalLaunchConfiguration(unsigned int threadCount,dim3 * grid,dim3 * block);

extern "C"
void CheckCUDAError();

extern "C" 
void readData(void * h_out, void * d_in, int length,int width);

extern "C"
void printDeviceInfo(int argc, char** argv);


extern "C"
void deallocateMemory(void * d_buffer);

extern "C"
void allocateMemory(void ** d_buffer,unsigned int memSize);

extern "C"
void printDebruijnGraphLongFmt(EulerVertex * d_ev, int vertexCount, unsigned int * d_l, unsigned int * d_e, EulerEdge * d_ee,int edgeCount,unsigned int kmerLength);

extern "C"
void printDebruijnGraphVizFmt(EulerVertex * d_ev, int vertexCount, unsigned int * d_l, unsigned int * d_e, EulerEdge * d_ee,int edgeCount,unsigned int kmerLength);

extern "C"
void printDebruijnGraph(EulerVertex * d_ev,unsigned int vertexCount, unsigned int * d_l, unsigned int * d_e, EulerEdge * d_ee,unsigned int edgeCount,unsigned int kmerLength,int format);

extern "C"
void initDevice();

extern "C"
void resetDevice();
#endif //UTILS_H //UTILS_H
