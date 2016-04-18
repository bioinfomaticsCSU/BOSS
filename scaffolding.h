#ifndef SCAFFOLDING_H_INCLUDED 
#define SCAFFOLDING_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "scaffoldgraph.h"

using namespace std;

typedef struct BFSGraphNode{
    long int contigIndex;
    BFSGraphNode * next;
    BFSGraphNode(){
        contigIndex = -1;
        next = NULL;
    }
}BFSGraphNode;

typedef struct BFSGraph{
    BFSGraphNode * node;
    BFSGraph(){
        node = NULL;
    }
}BFSGraph;

typedef struct BFSNode{
    long int orientation;
    long int contigIndex;
    long int gapDistance;
    BFSNode(){
        orientation = -1;
        contigIndex = -1;
        gapDistance = 0;
    }
}BFSNode;

typedef struct ShortContig{
    long int contigIndex;
    long int gapDistance;
    int orientation;
    ShortContig * next;
    ShortContig(){
        contigIndex = -1;
        next = NULL;
        gapDistance = 0;
        orientation = -2;
    }
}ShortContig;

typedef struct MergeIndex{
    long int contigIndex;
    long int rightPosIndex;
    long int leftPosIndex;
    MergeIndex * next;
    MergeIndex(){
        contigIndex = -1;
        rightPosIndex = -1;
        leftPosIndex = -1;
        next = NULL;
    }
}MergeIndex;

typedef struct MergeIndexHead{
    long int count;
    long int rightStart;
    long int leftStart;
    int orientation;
    MergeIndex * head;
    MergeIndexHead(){
        count = 0;
        rightStart = -1;
        leftStart = -1;
        orientation = -1;
        head = NULL;
    }
}MergeIndexHead;

typedef struct StringArray{
    char * str;
    StringArray * next;
    StringArray(){
        str = NULL;
        next = NULL;
    }
}StringArray;

typedef struct InsertSize{
    long int insertsize;
    InsertSize * next;
    InsertSize(){
        insertsize = 0;
        next = NULL;
    }
}InsertSize;

typedef struct SubGraphNode{
    long int index;
    SubGraphNode * next;
    SubGraphNode(){
        index = -1;
        next = NULL;
    }
}SubGraphNode;

typedef struct SubGraph{
    SubGraphNode * startNode;
    long int left;
    long int right;
    long int weight;
    double fitNumber;
    long int mateIndex;
    SubGraph * next;
    SubGraph(){
        left = 0;
        right = 0;
        weight = 0;
        mateIndex = -1;
        fitNumber = 0;
        startNode = NULL;
        next = NULL;
    }
}SubGraph;

ScaffoldSet * InitScaffoldSet(ContigSet * contigSet, long int contigCount);


ShortContig * ReverseShortContig(ShortContig * shortContig);
MergeIndexHead * DetermineMergeBetweenTwoScaffold1(ShortContig * right, ShortContig * left1, int index);

int SearchScaffoldEdge(long int index, ContigSequence * contigSequence);
int InsertShortContigBetweenScaffoldSet(ScaffoldGraph * scaffoldGraph, ScaffoldSet * scaffoldSet, long int contigCount, long int contigCutOff);
int InsertShortContigInScaffoldSet(ScaffoldGraph * scaffoldGraph, ScaffoldSet * scaffoldSet, long int contigCount, long int contigCutOff);
ScaffoldSet * GetScaffoldSet(ScaffoldGraph * scaffoldGraph, long int contigCount, long int contigCutOff);
void OutPutScaffoldTag(ScaffoldSet * scaffoldSet);
void OutPutScaffoldTag(ScaffoldSet * scaffoldSet, char * result);
void OutPutScaffoldSet(ScaffoldSet * scaffoldSet, ContigSet * contigSet, long int contigCount, char * result); 
void OutPutScaffoldSet(ScaffoldGraph * scaffoldGraph, ContigSet * contigSet, double cutoff, long int contigCount, char * suffix);
void OutPutScaffoldSetAllNumber(ScaffoldGraph * scaffoldGraph, ContigSet * contigSet, int cutoff, long int contigCount);
void OutPutScaffoldSetAll(ScaffoldGraph * scaffoldGraph, ContigSet * contigSet, double min, long int contigCount, char * outputFile);
long int GetScaffoldGraphEdgeNumber(ScaffoldGraph * scaffoldGraph, long int contigCount);

ScaffoldSet * OptimizeScaffoldSet(ContigSet * contigSet, ScaffoldSet * scaffoldSet, ScaffoldGraph * scaffoldGraph, long int & contigCount, long int realContigCount, long int * contigLength, long int insertsize, long int std, long int lambda);

BFSNode * InsertShortContigOfTail(ScaffoldGraph * scaffoldGraph, ContigSequence * previousBfsContigSequence, ContigSequence * bfsContigSequence, ScaffoldSet * curretScaffold, bool * visited, long int contigCutOff);
BFSNode * InsertShortContigOfHead(ScaffoldGraph * scaffoldGraph, ContigSequence * previousBfsContigSequence, ContigSequence * bfsContigSequence, ScaffoldSet * currentScaffold, bool * visited, long int contigCutOff);
ContigSequence * FindAndMergeTwoScaffold(BFSNode * bfsNode, ScaffoldSet * currentScaffold, ScaffoldSet * scaffoldSet, bool token);
void BFSScaffold(ScaffoldSet * scaffoldSet, ScaffoldSet * currentScaffold, ScaffoldGraph * scaffoldGraph, long int contigCount, long int contigCutOff);
int BFSScaffolding(ScaffoldGraph * scaffoldGraph, ScaffoldSet * scaffoldSet, long int contigCount, long int contigCutOff);

int AddShortContigToScaffoldSet(ScaffoldSet * scaffoldSet, long int contigCount, bool * index);

long int DetermineOrientationOfContigs(ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation, bool ** fixIndex, double minScore, bool last);
long int * IterativeDetermineOrderOfContigs(ContigSet * contigSet, ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation, long int * tempContigOrder, long int * contigPosition, bool ** fixIndex, long int insertsize, long int std, long int lambda, double minScore);

double * weightOrder(ScaffoldGraph * scaffoldGraph, long int contigCount);

#endif
