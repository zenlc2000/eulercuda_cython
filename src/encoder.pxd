#ifndef ENCODER_H
#define ENCODER_H

cdef extern from "common.h":
	#define KEY_TYPE_VALUE_TYPE
	ctypedef int KEY_T 
	ctypedef int * KEY_PTR 
	ctypedef int VALUE_T
	ctypedef int * VALUE_PTR 


cdef void encodeLmer(char  * d_buffer, 
		const unsigned int buffSize,
		const unsigned int readLength, 
		KEY_PTR d_lmers, 
		const unsigned int lmerLength,
		const unsigned int entriesCount
		)

cdef void encodeLmerComplement(	char  * d_buffer, 
				const unsigned int buffSize,
				const unsigned int readLength, 
				KEY_PTR d_lmers, 
				const unsigned int lmerLength,
				const unsigned int entriesCount
				)


cdef void computeKmer(	KEY_PTR d_lmers,
			KEY_PTR d_pkmers,
			KEY_PTR d_skmers,
			KEY_T validBitMask,
			const unsigned int readLength,
			const unsigned int entriesCount
			)



#endif
