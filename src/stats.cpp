/* maintain different statistics */
#define STATS_IMPL


#include "stats.h"
#undef STATS_IMPL

typedef struct StatItem{
	double	 value;
	const char * name;
	unsigned int fmtIdx;
} StatItem;

const char *fmts[]={ "<item name='%s' value='%.0f'/>\n","<item name='%s' value='%f'/>\n"};
#define _FMT_INDEX_NM			0
#define _FMT_INDEX_TM			1

StatItem _valmap[_MAX_STAT_ITEM]={
	{-1,"Debruijn Construction Time",_FMT_INDEX_TM},	//0
	{-1,"Euler Tour Phase 1 Time",_FMT_INDEX_TM},		//1
	{-1,"Euler Tour Phase 2 Time",_FMT_INDEX_TM},		//2
	{-1,"Euler Tour Total Time",_FMT_INDEX_TM},		//3
	{-1,"Component Label Time",_FMT_INDEX_TM},		//4
	{-1,"Spanning Tree Const. Time",_FMT_INDEX_TM},		//5
	{-1,"Swipe Execution Time",_FMT_INDEX_TM},		//6
	{-1,"Contig Generation Time",_FMT_INDEX_TM},		//7
	{-1,"Total Time",_FMT_INDEX_TM},			//8
	{-1,"l-mer Length",_FMT_INDEX_NM},			//9
	{-1,"Block Size",_FMT_INDEX_NM},			//10
	{-1,"l-mer Count",_FMT_INDEX_NM},			//11
	{-1,"k-mer Count",_FMT_INDEX_NM},			//12
	{-1,"Debruijn Vertices",_FMT_INDEX_NM},			//13
	{-1,"Debruijn Edges",_FMT_INDEX_NM},			//14
	{-1,"Circuit Vertices",_FMT_INDEX_NM},			//15
	{-1,"Circuit Edges",_FMT_INDEX_NM},			//16
	{-1,"Total Memory Allocated",_FMT_INDEX_NM},		//17
	{-1,"Peak Memory Allocated",_FMT_INDEX_NM},		//18
	{-1,"k-mer Extraction Time",_FMT_INDEX_TM},		//19
	{-1,"Hash table Construction Time",_FMT_INDEX_TM},	//20
	{-1,"k-mer GPU Encoding Time",_FMT_INDEX_TM},		//21
	{-1,"Input Read Processing Time",_FMT_INDEX_TM},		//22
	{-1,"Original Read Count",_FMT_INDEX_NM},		//23
	{-1,"ERROR CORRECTION TIME",_FMT_INDEX_TM},		//24
	{-1,"Corrected Read Count",_FMT_INDEX_NM},		//25
	{-1,"Modified Read Count",_FMT_INDEX_NM},		//26
	{-1,"Mutation GPU Time",_FMT_INDEX_TM},		//27
	{-1,"Mutation Score Accumulation Time",_FMT_INDEX_TM},		//28
	{-1,"Select Mutation Time",_FMT_INDEX_TM},		//29
	{-1,"Select Position Time",_FMT_INDEX_TM}		//30
	};

extern "C"
void setStatItem(SIID_T statId,double value){
	
	if(statId<_MAX_STAT_ITEM){
		_valmap[statId].value=value;
	}
}

extern "C"
double getStatItem(SIID_T statId){

	return _valmap[statId].value;
}

extern "C"
void writeStatItem(FILE * file ,SIID_T statId, unsigned int fmt) {
	//right now using only xml formatting
	
	fprintf(file,fmts[_valmap[statId].fmtIdx],_valmap[statId].name,_valmap[statId].value);
}

extern "C"
void writeStat(FILE * file, unsigned int fmt){
	for (int i =0; i< _MAX_STAT_ITEM;i++){
		if(_valmap[i].value>0){
			writeStatItem(file,i,fmt);
		}
	}
}
