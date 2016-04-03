#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include "graph.h"
#include "gpuhash.h"

//extern "C"
void constructDebruijnGraphDevice(	unsigned int * ecount,
					KEY_PTR d_lmerKeys,		//in lmer keys
					VALUE_PTR d_lmerValues,	//in lmer values
					unsigned int lmerCount,			//in total lmers
					KEY_PTR d_kmerKeys,		//in
					unsigned long kmerCount,		//in  total kmers 
					unsigned int l,					//in k
					KEY_PTR d_TK,
					VALUE_PTR d_TV,
					unsigned int * d_bucketSeed,
					unsigned int bucketCount,	
					EulerVertex ** d_ev,			//out
					unsigned int ** d_l,			//out
					unsigned int ** d_e,			//out									
					EulerEdge ** d_ee			//out
					);


//extern "C"
void constructDebruijnGraphGold(	KEY_PTR h_lmerKeys,		//in lmer keys
					VALUE_PTR h_lmerValues,	//in lmer values
					unsigned int lmerCount,			//in total lmers
					KEY_PTR h_kmerKeys,		//in
					unsigned long kmerCount,		//in  total kmers 
					unsigned int l,					//in k
					KEY_PTR h_TK,VALUE_PTR h_TV,
					unsigned int * h_bucketSeed,
					unsigned int bucketCount,	
					EulerVertex ** h_ev,			//out
					unsigned int ** h_l,			//out
					unsigned int ** h_e,			//out									
					EulerEdge ** h_ee,			//out
					unsigned int * ecount);		//out


#endif
