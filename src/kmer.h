#ifndef KMER_H
#define KMER_H

#define THRD 128
#define GRID 1
//#define KMER_LENGTH 10
#define GRID_INIT GRID 
#define GRID_END GRID
#define MIN_THRD THRD-1
#define MAX_THRD THRD-1

typedef struct KmerData{
	//unsigned char  kmer[KMER_LENGTH];		//this can be compressed to use 2 bits only, at the moment using fixed length allocation	
	unsigned int	count;
	unsigned char neighbours; // HO incoming  LO outgoing TGCA : TGCA  XXXXYYYY
} KmerData;


#endif