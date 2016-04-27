from libc.stdlib cimport malloc, free
cimport cython
from cpython.string cimport PyString_AsString

cimport pyencoder
from pyencoder cimport encodeLmer, computeKmer,encodeLmerComplement

import numpy as np
cimport numpy as np

def encode_lmer(char[::1] py_buffer,
                const unsigned int buffSize,
                const unsigned int readLength,
                KEY_T[::1] py_lmers,
                const unsigned int lmerLength,
                const unsigned int entriesCount):

    return encodeLmer(&py_buffer[0], buffSize, readLength,&py_lmers[0],lmerLength,entriesCount)

def compute_kmer(KEY_T[::1] py_lmers,
                KEY_T[::1] py_pkmers,
                KEY_T[::1] py_skmers,
                KEY_T[::1] h_lmers,
			    KEY_T[::1] h_pkmers,
			    KEY_T[::1] h_skmers,
                unsigned int validBitMask,
                const unsigned int readLength,
                const unsigned int entriesCount):

    return computeKmer(&py_lmers[0], &py_pkmers[0], &py_skmers[0], &h_lmers[0], &h_pkmers[0], &h_skmers[0], validBitMask,readLength, entriesCount)


def encode_lmer_complement(char[::1] py_buffer,
                const unsigned int buffSize,
                const unsigned int readLength,
                KEY_T[::1] py_lmers,
                const unsigned int lmerLength,
                const unsigned int entriesCount):


    return encodeLmerComplement(&py_buffer[0], buffSize, readLength,&py_lmers[0],lmerLength,entriesCount)


