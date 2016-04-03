from libc.stdlib cimport malloc, free

cimport pyencoder


def encodeLmer(char[:] d_buffer, 
                const unsigned int buffSize,
                const unsigned int readLength,
                int[:] d_lmers,
                const unsigned int lmerLength,
                const unsigned int entriesCount):
        return pyencoder.encodeLmer(d_buffer, buffSize, readLength,d_lmers,lmerLength,entriesCount)

def encodeLmerComplement(char[:] d_buffer, 
                const unsigned int buffSize,
                const unsigned int readLength,
                int[:] d_lmers,
                const unsigned int lmerLength,
                const unsigned int entriesCount):
        return pyencoder.encodeLmerComplement(d_buffer, buffSize, readLength,d_lmers,lmerLength,entriesCount)


def computeKmer(int[:] d_lmers,
                int[:] d_pkmers,
                int[:] d_skmers,
                int validBitMask,
                const unsigned int readLength,
                const unsigned int entriesCount):
        return pyencoder.computeKmer(d_lmers, d_pkmers, d_skmers, validBitMask,readLength, entriesCount)
