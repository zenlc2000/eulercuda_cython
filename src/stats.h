#ifndef STATS_H
#define STATS_H

#include <stdio.h>
/**This is not a thread safe library ;use with caution-- faraz*/

typedef unsigned int SIID_T;	//stat item ID _type

#define FMT_XML		0;	//
#define FMT_CSV		1;	//Not Implemented
#define FMT_TAB		2;	//Not Implemented
#define FMT_SMP		3;	//NOT Implemented


#define TM_DEBRUIJN_CONSTRUCTION	0
#define TM_EULER_TOUR_P1		1
#define TM_EULER_TOUR_P2		2
#define TM_EULER_TOUR			3
#define TM_COMPONENT			4
#define TM_SPANNING_TREE		5
#define TM_SWIPE_EXECUTION		6
#define TM_CONTIG_GENERATION		7
#define TM_TOTAL_TIME			8
#define NM_LMER_LENGTH			9
#define NM_BLOCK_SIZE			10
#define NM_LMER_COUNT			11
#define NM_KMER_COUNT			12
#define NM_DEBRUIJN_VTX			13
#define NM_DEBRUIJN_EDG			14
#define NM_CIRCUIT_VTX			15
#define NM_CIRCUIT_EDG			16
#define NM_TOTAL_MEMORY_ALLOCATED	17
#define NM_PEAK_MEMORY_ALLOCATED	18
#define TM_KMER_EXTRACTION_TIME		19
#define TM_HASHTABLE_CONSTRUCTION	20
#define TM_KMER_GPUENC_TIME			21
#define TM_PROCESS_READ_TIME		22
#define NM_READ_COUNT				23
#define TM_ERROR_CORRECTION_TIME	24
#define NM_CORRECTED_READ_COUNT		25
#define NM_MODIFIED_READ_COUNT		26
#define TM_MUTATION_GPU_TIME		27
#define TM_ACCUMULATE_GPU_TIME		28
#define TM_SELECT_MUTATION_GPU_TIME	29
#define TM_SELECT_POSITION_GPU_TIME	30


#define _MAX_STAT_ITEM			31

#define _ATT_VAL			0
#define _ATT_NAME			1
#define _ATT_FMT			2
#define _MAX_ATT			3

//typedef tuple<double ,std::string ,std::string> StatItem;

//#ifndef STATS_IMPL 
#define _EXTERN_ extern "C"
//#else 
//#define _EXTERN_
//#endif

_EXTERN_ 
void  setStatItem(SIID_T statId, double value); //only using predefine Ids.
/*
 for future ,add custom stat item,
e.g
SIID_T createStatItem(string msg, string fmt);

**/
_EXTERN_
double getStatItem(SIID_T statId);

_EXTERN_
void writeStatItem(FILE * file, SIID_T statId,unsigned int fmt);

_EXTERN_
void writeStat(FILE * file,unsigned int fmt);

#endif //STATS_H
