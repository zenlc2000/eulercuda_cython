#ifndef COMMON_H
#define COMMON_H //common.h

#ifndef LOG_LEVELS
#define LOG_LEVELS

#define LOG_LVL_OFF		0
#define LOG_LVL_ERROR		2
#define LOG_LVL_WARN		4
#define LOG_LVL_MSG		5
#define LOG_LVL_INFO		6
#define	LOG_LVL_DETAIL		7
#define LOG_LVL_DEBUG		8
#define LOG_LVL_ALL		9

#endif

/*
#ifndef STATS_COUNTER
#define STATS_COUNTER

#define	DEBRUIJN_CONSTRUCTION
#define EULER_TOUR
#define COMPONENT
#define SPANNING_TREE
#define SWIPE_EXECUTION
#define TOTAL_TIME
#define DEBRUIJN_VERTICES
#define DEBRUIJN_EDGES
#define CIRCUIT_VERTICES
#define CIRCUITS
#define CIRCUIT_EDGES
#define TOTAL_MEMORY_ALLOCATED
#define PEAK_MEMORY_ALLOCATED

#endif

*/


#ifndef KEY_TYPE_VALUE_TYPE
#define KEY_TYPE_VALUE_TYPE
typedef unsigned  long long KEY_T ;
typedef KEY_T * KEY_PTR ;
typedef unsigned int VALUE_T;
typedef VALUE_T * VALUE_PTR ;

#define LMER_PREFIX(lmer,bitMask) ((lmer & (bitMask<<2))>>2)
#define LMER_SUFFIX(lmer,bitMask) ((lmer & bitMask))
#define CUDA_NUM_READS 1024* 32
#endif 


#endif
