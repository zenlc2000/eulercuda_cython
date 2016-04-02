from libc.stdlib cimport malloc, free

cimport encoder


def encodeLmer(char[:] d_buffer, 
        const unsigned int buffSize,
        const unsigned int readLength, 
        int[:] d_lmers, 
        const unsigned int lmerLength,
        const unsigned int entriesCount):
        return encoder.encodeLmer(<char *> &d_buffer[0], buffSize, readLength,<KEY_PTR> &d_lmers[0],lmerLength,entriesCount)

def encodeLmerComplement(char[:] d_buffer, 
        const unsigned int buffSize,
        const unsigned int readLength, 
        int[:] d_lmers, 
        const unsigned int lmerLength,
        const unsigned int entriesCount):
        return encoder.encodeLmerComplement(<char *> &d_buffer[0], buffSize, readLength,<KEY_PTR> &d_lmers[0],lmerLength,entriesCount)


def computeKmer(int[:] d_lmers,
        int[:] d_pkmers,
        int[:] d_skmers,
        int validBitMask,
        const unsigned int readLength,
        const unsigned int entriesCount):
        return encoder.computeKmer(<KEY_PTR> &d_lmers[0], <KEY_PTR> &d_pkmers[0], <KEY_PTR> &d_skmers[0], validBitMask,readLength, entriesCount)
