#ifndef GRAPH_H
#define GRAPH_H
#include "common.h"

typedef struct Vertex{
	unsigned int vid;
	unsigned int n1;
	unsigned int n2;
	
} Vertex;

typedef struct EulerVertex{
	KEY_T	vid;
	unsigned int  ep;
	unsigned int  ecount;
	unsigned int  lp;
	unsigned int  lcount;
	
}EulerVertex;
typedef struct EulerEdge{
	KEY_T eid;
	unsigned int v1;
	unsigned int v2;
	unsigned int s;
	unsigned int pad;
	//unsigned int readId;
}EulerEdge;
/*
typedef struct CircuitVertex{
	unsigned int cid;
	unsigned int estart;
	unsigned int eend;
}CircuitVertex;*/
typedef struct CircuitEdge{
	unsigned int ceid;
	unsigned e1;
	unsigned e2;
	unsigned c1;
	unsigned c2;
}CircuitEdge;
#endif
