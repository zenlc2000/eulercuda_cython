#ifndef ENCODER_H
#define ENCODER_H

extern "C"
void encodeLmer(char  * d_buffer, 
		const unsigned int buffSize,
		const unsigned int readLength, 
		KEY_PTR d_lmers, 
		const unsigned int lmerLength,
		const unsigned int entriesCount
		);
extern "C"
 void encodeLmerComplement(	char  * d_buffer, 
				const unsigned int buffSize,
				const unsigned int readLength, 
				KEY_PTR d_lmers, 
				const unsigned int lmerLength,
				const unsigned int entriesCount
				);

extern "C"
void computeKmer(	KEY_PTR d_lmers,
			KEY_PTR d_pkmers,
			KEY_PTR d_skmers,
			KEY_T validBitMask,
			const unsigned int readLength,
			const unsigned int entriesCount
			);

//compress reads

#endif
