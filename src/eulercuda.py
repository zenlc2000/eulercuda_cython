import sys
import argparse
import pyencoder
import pygpuhash
#import pydebruijn
import numpy as np
from numba import cuda
import cython

def parse_fastq(filename):
    """
    Read fastq formatted <filename> and return a dictionary of
    read_name : read
    """
    file = open(filename)
    result = {}
    current_name = None
    for i, line in enumerate(file):
        if i % 4 == 0:
            current_name = line.rstrip('\n')
        if i % 4 == 1:
            result[current_name] = line.rstrip('\n')
        if i % 4 == 3:
            print('quality is ' + line.rstrip('\n'))
    return result


def read_fastq(filename):
    """
    Read fastq formatted <filename> and return a list of reads
    """
    with open(filename, "r") as infile:
        result = []
        for i, line in enumerate(infile):
            if i % 4 == 1:
                result.append(line.rstrip('\n'))
        return result


def doErrorCorrection(readBuffer, readCount, ec_tuple_size, max_ec_pos):
    return readCount


c_lmerMask = [
    0x0000000000000003, 0x000000000000000F, 0x000000000000003F, 0x00000000000000FF,  # 0   1   2   3
    0x00000000000003FF, 0x0000000000000FFF, 0x0000000000003FFF, 0x000000000000FFFF,  # 4   5   6   7
    0x000000000003FFFF, 0x00000000000FFFFF, 0x00000000003FFFFF, 0x0000000000FFFFFF,  # 8   9   10  11
    0x0000000003FFFFFF, 0x000000000FFFFFFF, 0x000000003FFFFFFF, 0x00000000FFFFFFFF,  # 12  13  14  15
    0x00000003FFFFFFFF, 0x0000000FFFFFFFFF, 0x0000003FFFFFFFFF, 0x000000FFFFFFFFFF,  # 16  17  18  19
    0x000003FFFFFFFFFF, 0x00000FFFFFFFFFFF, 0x00003FFFFFFFFFFF, 0x0000FFFFFFFFFFFF,  # 20  21  22  23
    0x0003FFFFFFFFFFFF, 0x000FFFFFFFFFFFFF, 0x003FFFFFFFFFFFFF, 0x00FFFFFFFFFFFFFF,  # 24  25  26  27
    0x03FFFFFFFFFFFFFF, 0x0FFFFFFFFFFFFFFF, 0x3FFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF  # 28  29  30  31
]


c_codeF = [0, 0, 0, 1, 3, 0, 0, 2]
c_codeR = [0, 3, 0, 2, 0, 0, 0, 1]

c_shifter = [
    [0, 0, 0, 0],
    [1, 4, 16, 64],
    [2, 8, 32, 128],
    [3, 12, 48, 192],
]

lmerMask = np.array(c_lmerMask, dtype=np.uint64)
codeF = np.array(c_codeF)
codeR = np.array(c_codeR)
shifter = np.array(c_shifter)


# @cuda.jit #(argtypes=[i1[:],i1,i1[:],i1[:]])  # NUMBA_ENABLE_CUDASIM =0 = GPU, 1 = CPU/debug
# def py_encode_lmer(buffer,lmerLength, dlmers,read, lmerMask, codeF, shifter):
#
#     tid = cuda.threadIdx.x
#     rOffset = (cuda.blockDim.x * cuda.blockDim.y * cuda.blockIdx.y) + (cuda.blockDim.x * cuda.blockDim.y * cuda.blockIdx.x) + (cuda.blockDim.x * cuda.threadIdx.y);
#     lmer = 0
# #    lmers = []
#
#     read[tid] = buffer[rOffset + tid]
#     cuda.syncthreads()
#
#     for i in range(7):
#         lmer = (lmer >> 8) | (shifter[codeF[read[cuda.threadIdx.x+i*4]& 0x07]][3] |
#                 shifter[codeF[read[cuda.threadIdx.x+i*4+1]& 0x07]][2] |
#                 shifter[codeF[read[cuda.threadIdx.x+i*4+2]& 0x07]][1] |
#                 codeF[read[cuda.threadIdx.x+i*4+3] & 0x07])
#
#     lmer = (lmer << ((32 - lmerLength) >> 1))
#     u_lmer = np.uint64(lmer)
#     mask = lmerMask[lmerLength - 1]
#     lmer = u_lmer & mask
#
#     dlmers[rOffset + tid] = lmer

def readLmersKmersCuda(readBuffer, readLength, readCount, lmerLength, lmerKeys, lmerValues, lmerCount, kmerKeys,
                       kmerValues, kmerCount):
    """
    char * d_reads = NULL;
    KEY_PTR h_lmersF = NULL;
    KEY_PTR h_lmersR = NULL;
    KEY_PTR d_lmers = NULL;
    KEY_PTR h_pkmersF = NULL;
    KEY_PTR h_pkmersR = NULL;
    KEY_PTR h_skmersF = NULL;
    KEY_PTR h_skmersR = NULL;
    KEY_PTR d_pkmers = NULL;
    KEY_PTR d_skmers = NULL;
    unsigned int readProcessed=0;
    unsigned int kmerGPUEncTimer = 0;
    unsigned int kmerExtractTimer=0;

    typedef dense_hash_map<KEY_T, VALUE_T> map;
    map kmerMap(readLength*readCount);
    map lmerMap(readLength*readCount);
    """

   # bufferSize = buffer.size
    kmerMap = {}
    lmerMap = {}
    d_lmers = np.zeros(len(readBuffer), dtype = 'uint64')
    d_pkmers = np.empty_like(d_lmers)
    d_skmers = np.empty_like(d_lmers)
    h_lmersF = np.empty_like(d_lmers)
    h_pkmersF = np.empty_like(d_lmers)
    h_skmersF = np.empty_like(d_lmers)
    h_lmersR = np.empty_like(d_lmers)
    h_pkmersR = np.empty_like(d_lmers)
    h_skmersR = np.empty_like(d_lmers)

    CUDA_NUM_READS = 1024 * 32
    if readCount < CUDA_NUM_READS:
        readToProcess = readCount
    else:
        readToProcess = CUDA_NUM_READS
    kmerBitMask = 0


    bufferSize = sys.getsizeof(np.uint8) * readLength * readToProcess

    entriesCount = readLength * readCount

    for _ in range(0, (lmerLength - 1) * 2):
        kmerBitMask = (kmerBitMask << 1) | 1
       # print("kmerBitMask = " + str(kmerBitMask))
    readProcessed = 0
    # Originally a loop slicing readBuffer into chunks then process each chunk
   # d_reads = readBuffer[:]  # copy list by value
    d_lmerMask = cuda.to_device(lmerMask)
    d_codeF = cuda.to_device(codeF)
    d_codeR = cuda.to_device(codeR)
    d_shifter = cuda.to_device(shifter)
    validLmerCount = readLength - lmerLength + 1
    while readProcessed < readCount:
        buffer = np.fromstring(''.join(readBuffer), dtype='uint8')
        # d_buffer = cuda.to_device(buffer)
        # read = np.empty_like(buffer)
        # d_read = cuda.to_device(read)
        # dlmers = cuda.to_device(d_lmers)
        pyencoder.encode_lmer(buffer, bufferSize, readLength, d_lmers, lmerLength, entriesCount)

        # extract kmer
        pyencoder.compute_kmer(d_lmers, d_pkmers, d_skmers, h_lmersF, h_pkmersF, h_skmersF, kmerBitMask, readLength, entriesCount)



        # reverse
        pyencoder.encode_lmer_complement(buffer, bufferSize, readLength, d_lmers, lmerLength, entriesCount)

        # extract kmer
        pyencoder.compute_kmer(d_lmers, d_pkmers, d_skmers, h_lmersR, h_pkmersR, h_skmersR, kmerBitMask, readLength, entriesCount)



        lmerEmpty, kmerEmpty = 0, 0

        # Here he fills the kmerMap and lmerMap with a nested for loop
        # for j in range(readToProcess):
        #     for i in range(validLmerCount):
        for index in range(readToProcess):
                # index = j * readLength + i
            kmerMap[h_pkmersF[index]] = 1
            kmerMap[h_skmersF[index]] = 1
            kmerMap[h_pkmersR[index]] = 1
            kmerMap[h_skmersR[index]] = 1

            if h_lmersF[index] == 0:
                lmerEmpty += 1
            else:
                if lmerMap.get(h_lmersF[index]) == None:
                    lmerMap[h_lmersF[index]] = 1
                else:
                    lmerMap[h_lmersF[index]] += 1
            if h_lmersR[index] == 0:
                lmerEmpty += 1
            else:
                if lmerMap.get(h_lmersR[index]) == None:
                    lmerMap[h_lmersR[index]] = 1
                else:
                    lmerMap[h_lmersR[index]] += 1
        readProcessed += readToProcess
        readToProcess = readCount - readToProcess
        if readCount < CUDA_NUM_READS:
            readToProcess = readCount
        else:
            readToProcess = CUDA_NUM_READS
        # End of chunking loop

    kmerCount = len(kmerMap) + kmerEmpty
    # TODO: Log message with kmer count
    print('kmer count = ' + kmerCount)

    # kmerKeys = []
    # kmerValues = []
    #
    for index, k, v in enumerate(kmerMap.items()):
        kmerKeys[index] = k
        kmerValues[index] = index

    # original code has if below, but I don't see that it will ever be executed
    #     if (kmerEmpty > 0){
    #     ( * kmerKeys)[index]=0;
    #     ( * kmerValues)[index]=index;
    #     }

    lmerCount = len(lmerMap) + lmerEmpty
    # lmerKeys = []
    # lmerValues = []
    #
    lmerKeys = lmerMap.keys()
    lmerValues = lmerMap.values()

    if lmerEmpty > 0:
        lmerKeys[len(lmerMap)] = 0
        lmerValues[len(lmerMap)] = lmerEmpty

    return lmerCount, kmerCount
def constructDebruijnGraph(readBuffer, readCount, readLength, lmerLength, evList, eeList, levEdgeList, entEdgeList,
                           edgeCountList, vertexCountList):
    """
    ///variables

    KEY_PTR			h_lmerKeys =NULL;
    VALUE_PTR 		h_lmerValues= NULL;
    KEY_PTR			d_lmerKeys =NULL;
    VALUE_PTR 		d_lmerValues= NULL;
    unsigned int 	lmerCount=0;
    KEY_PTR 		h_kmerKeys=NULL;
    VALUE_PTR 		h_kmerValues=NULL;
    KEY_PTR 		d_kmerKeys=NULL;
    VALUE_PTR 		d_kmerValues=NULL;
    unsigned int 	kmerCount=0;
    KEY_PTR			d_TK=NULL;
    VALUE_PTR		d_TV=NULL;
    unsigned int 	tableLength=0;
    unsigned int	bucketCount=0;
    unsigned int *	d_bucketSize=NULL;

    unsigned int coverage =20;

    EulerVertex * d_ev=NULL;
    EulerEdge 	* d_ee=NULL;
    unsigned int * d_levEdge=NULL;
    unsigned int * d_entEdge=NULL;

    """
    h_lmerKeys = []
    h_lmerValues = []

    lmerCount = 0
    h_kmerKeys = []
    h_kmerValues = []
    d_kmerKeys = []
    d_kmerValues = []
    kmerCount = 0
    d_TK = []
    d_TV = []
    tableLength = 0
    bucketCount = 0
    d_bucketSize = []

    coverage = 20
    d_ev = []
    d_ee = []
    d_levEdge = []
    d_entEdge = []

    py_buffer = '\n'.join(readBuffer)

    # May need to return unpacked tuple of integer variables
    lmerCount, kmerCount = readLmersKmersCuda(readBuffer, readLength, readCount, lmerLength, h_lmerKeys, h_lmerValues, lmerCount, h_kmerKeys,
                       h_kmerValues, kmerCount)
    # initDevice()
    # setStatItem(NM_LMER_COUNT, lmerCount);
    # setStatItem(NM_KMER_COUNT, kmerCount);

    # lots of memory movement to Device
    d_lmerKeys = np.empty(lmerCount * cython.sizeof(cython.unsignedlonglong))
    d_lmerValues = []
    # from gpuhash & gpuhash2

    # Need to return these back out...
    # d_TK,  d_TV, tableLength,  d_bucketSize,  bucketCount)
    pygpuhash.create_hash_table(d_kmerKeys, d_kmerValues, kmerCount)

#    constructDebruijnGraphDevice(d_lmerKeys, d_lmerValues, lmerCount,
#        d_kmerKeys, kmerCount, l, d_TK, d_TV, d_bucketSize, bucketCount,
#        & d_ev, & d_levEdge, & d_entEdge, & d_ee, edgeCount);



    def findEulerTour(evList, eeList, levEdgeList, entEdgeList, edgeCountList, vertexCountLsit, lmerLength, outfile):
    pass


def assemble2(infile, outfile, lmerLength, errorCorrection, max_ec_pos, ec_tuple_size):
    """
    Do the assemble
    """
    # TODO: figure out logging
    # TODO: Unit testing

    # for performance reasons, may want to make these Numpy arrays

    # char * 		readBuffer=NULL;
    # EulerVertex * 	ev=NULL;
    # EulerEdge 	* ee=NULL;#
    # unsigned int * 	levEdge=NULL;
    # unsigned int * 	entEdge=NULL;
    # unsigned int  	edgeCount=0;
    # unsigned int 	vertexCount=0;
    # unsigned int 	readCount=0;

    readBuffer = read_fastq(infile)

    readCount = len(readBuffer)
    readLength = len(readBuffer[0])
    evList = []
    eeList = []
    levEdgeList = []
    entEdgeList = []
    edgeCountList = []
    vertexCountList = []

    if readCount > 0:
        if errorCorrection:
            readCount = doErrorCorrection(readBuffer, readCount, ec_tuple_size, max_ec_pos)
        constructDebruijnGraph(readBuffer, readCount, readLength,
                               lmerLength, evList, eeList, levEdgeList, entEdgeList, edgeCountList, vertexCountList)
        findEulerTour(evList, eeList, levEdgeList, entEdgeList, edgeCountList, vertexCountList, lmerLength, outfile)


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', action='store', dest='input_filename',
                        help='Input Fie Name')
    parser.add_argument('-o', action='store', dest='output_filename',
                        help='Output File Name')
    parser.add_argument('-k', action='store', dest='k', type=int,
                        help='kmer size')
    parser.add_argument('-d', action='store_true', default=False,
                        help='Use DDFS')
    results = parser.parse_args()
    # Need to process commandline args. Probably just copy=paste from disco3
    if results.input_filename == '':
        fname = '../data/read_1.fq'
    else:
        fname = results.input_filename
    readBuffer = read_fastq(fname)
    # assemble2(inputFileName, outputFileName, readLength, assemble,
    # lmerLength, coverage,errorCorrection, max_ec_pos,ec_tuple_size);

    # void assemble2(	const char * filename, 	//input filename
    # 		const char * output, 	//output filename
    # 		unsigned int readLength,	//readLength
    # 		bool assemble,
    # 		unsigned  int l,		//lmer length
    # 		unsigned int coverage,	//coverage M
    # 		bool errorCorrection,
    # 		unsigned int max_ec_pos,	//ec positions
    # 		unsigned int ec_tuple_size	//ec tuple size
    # 		){

    assemble2(results.input_filename, results.output_filename, 17, False, 20, 20)
