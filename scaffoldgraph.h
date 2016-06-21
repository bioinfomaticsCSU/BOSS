#ifndef BUILDSCAFFOLDGRAPH_H_INCLUDED 
#define BUILDSCAFFOLDGRAPH_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

extern int lambda;
extern double scoreThreshold;
extern int minMappedPairedReadNumber;

typedef struct InputArg{
    char * bamFileName1;
    char * bamFileName2;
    long int readLength;
    long int insertsize;
    long int std;
    double minEdgeWeight;
    double minEdgeLinkNumber;
    double minRepetitiveCov;
    bool pairedRead;
    bool edgeMeanMethod;
    InputArg(){
        bamFileName1 = NULL;
        bamFileName2 = NULL;
        readLength = 0;
        insertsize = 0;
        std = 0;
        minEdgeWeight = 0;
        minEdgeLinkNumber = 0;
        minRepetitiveCov = 0;
        pairedRead = false;
        edgeMeanMethod = false;
    }
}InputArg;


typedef struct ContigSet{
    char * contig;
    ContigSet * next;
    ContigSet(){
        contig = NULL;
        next = NULL;
    }
}ContigSet;

typedef struct ReadMapPosition{
    long int * index;
    long int * leftIndex;
    long int * rightIndex;
    long int * noMapLeftIndex;
    long int * noMapRightIndex;
    long int * repeatIndex;
    long int * leftReadCoverage;
    long int * rightReadCoverage;
    double avg;
    ReadMapPosition(){
        index = NULL;
        leftIndex = NULL;
        rightIndex = NULL;
        noMapLeftIndex = NULL;
        noMapRightIndex = NULL;
        repeatIndex = NULL;
        leftReadCoverage = NULL;
        rightReadCoverage = NULL;
        avg = 0;
    }
}ReadMapPosition;


typedef struct MapPosition{
    long int pos;
    long int matePos;
    bool orientation;
    bool mateOrientation;
    MapPosition * next;
    MapPosition(){
        pos = -1;
        matePos = -1;
        orientation = false;
        mateOrientation = false;
        next = NULL;
    }
}MapPosition;

typedef struct ScaffoldEdge{
    long int contigIndex;
    double fitDistribution;
    double fitNumber;
    double readDistribution;
    long int gapDistance;
    MapPosition * mapPosition;
    ScaffoldEdge * next;
    long int edgeWeight;
    ScaffoldEdge(){
        contigIndex = -1;
        fitDistribution = 0;
        fitNumber = 0;
        readDistribution = 0;
        mapPosition = NULL;
        edgeWeight = 0;
        next = NULL;
    }
}ScaffoldEdge;


typedef struct ScaffoldGraph{
    long int length;
    double leftReadCoverage;
    double rightReadCoverage;
    ScaffoldEdge * outLink;
    ScaffoldEdge * inLink;
    ScaffoldGraph(){
        length = 0;
        leftReadCoverage = 0;
        rightReadCoverage = 0;
        outLink = NULL;
        inLink = NULL;
    }
}ScaffoldGraph;

typedef struct ScaffoldGraphHead{
    
    ScaffoldGraph * scaffoldGraph;
    ScaffoldGraphHead(){
        scaffoldGraph = NULL;
    }
        
}ScaffoldGraphHead;


typedef struct ContigSequence{
    long int index;
    int orientation;
    long int gapDistance;
    ContigSequence * next;
    ContigSequence(){
        index = -1;
        orientation = -1;
        gapDistance = 0;
        next = NULL;
    }
}ContigSequence;

typedef struct ScaffoldSet{
    char * scaffold;
    long int length;
    ContigSequence * contigSequence;
    ScaffoldSet * next;
    ScaffoldSet(){
        scaffold = NULL;
        length = 0;
        contigSequence = NULL;
        next = NULL;
    }
}ScaffoldSet;

typedef struct ScaffoldToContig{
    
    long int contigIndex;
    long int start;
    bool orientation;
    long int length;
    ScaffoldToContig(){
        contigIndex = -1;
        start = -1;
        orientation = 0;
        length = 0;
    }
    
}ScaffoldToContig;

typedef struct PairedReadMappedData{
    long int leftPosition;
    bool orientationLeft;
    long int rightPosition;
    bool orientationRight;
    long int leftReference;
    long int rightReference;
    PairedReadMappedData(){
        leftPosition = -1;
        orientationLeft = 0;
        rightPosition = -1;
        orientationRight = 0;
        leftReference = -1;
        rightReference = -1;
    }
}PairedReadMappedData;

double AddPositionIndex(ScaffoldGraph * scaffoldGraph, ReadMapPosition * readMapPosition, long int scaffoldCount);

long int * GetScaffoldSetLength(ScaffoldSet * temp, long int scaffoldCount);
void GetScaffoldSetLength(ScaffoldSet * scaffoldSet, long int & scaffoldCount, long int * contigLength);

ScaffoldToContig * GetScaffoldToContig(ScaffoldSet * scaffoldSet, long int * contigLength, long int contigCount, long int & scaffoldCount);
long int GetScaffoldofContigPosition(ScaffoldToContig * scaffoldToContig, long int contigIndex, long int * contigLength, long int position, long int readLength);
bool GetScaffoldofContigOrientation(ScaffoldToContig * scaffoldToContig, long int contigIndex, bool orientation);

ContigSet * GetContigSet(char * str, long int & contigNum);
int BuildScaffoldGraphFromTwoBam(ScaffoldSet * scaffoldSet, long int * scaffoldLength, ContigSet * contigSet, long int * contigLength, ScaffoldToContig * scaffoldToContig, long int scaffoldCount, ScaffoldGraphHead * scaffoldGraphHead, InputArg * inputArg, char * result);
bool ReverseComplement(char * temp1, char * temp2);
void AppendRight( char * temp, char * temp1, char * temp2,long int kmerLength);
void CopyContig(char *temp, char * temp1);
double f(double x);
double F(double a,double b,double ep=1e-6);
double MapFitDistribution(long int * index, long int * noMapIndex, MapPosition * tempMapPosition, long int * distance, double correctP, long int length, long int length1, long int readLength, long int gap, long int insertsize, long int std, double & cgap, double & fitN, double & fitD, double & mapP, int token, double avgP);
double MapFitDistribution1(long int * index, long int * noMapIndex, MapPosition * tempMapPosition, long int * distance, double correctP, long int length, long int length1, long int readLength, long int gap, long int insertsize, long int std, double & cgap, double & fitN, double & fitD, double & mapP, int token, double avgP);
int AddScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int refId, bool refStrand, long int pos, long int mateRefId, bool mateRefStrand, long int matePos);
int DeleteScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, bool orientation, long int cutOff,long int maxIndex);
int DeleteScaffoldEdgeShort(ScaffoldGraph * scaffoldGraph, long int index, bool orientation, long int cutOff, long int maxIndex);
int DeleteScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, bool orientation);
int DeleteSpecialScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, long int index1);
double ComputeReadDistribution(MapPosition * read, long int * readMapPosition, long int * distance, long int contigLength, long int contigLength1, long int readLength, long int gap, long int insertsize, long int std1, double p);
int ComputeTruePairedMapping(BamAlignment * alignment, long int len, long int len1, long int readLength, long int insertsize, long int std);
int ComputeTruePairedMapping(PairedReadMappedData * alignment, long int readLength, long int insertsize, long int std);
MapPosition * OptimizeMapPosition(MapPosition * mapPosition, long int & edgeWeight);
long int GetGapDistance(long int * distance, long int count, long int insertsize);
int OptimizeScaffoldGraph(ScaffoldGraph * scaffoldGraph, long int contigCount, long int readLength, long int insertsize, long int std, long int contigCutOff, bool edgeWeightMethod);
void * EqualScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int contigCount);

ScaffoldEdge * SearchScaffoldGraphEdge(ScaffoldEdge * scaffoldEdge, long int contigIndex);
long int DeleteScaffoldGraphEdge(ScaffoldGraph * scaffoldGraph, long int index, ScaffoldEdge * scaffoldEdge);

int KMPIndexOfContig(char * contig, char * pattern, long int * next);
int MergeSubContig(ContigSet * contigSet, long int contigCount);
void KMPGetNext(char * pattern, long int * next);
double ComputeVar(long int * distance, long int count);

int GetOutAndInLinkNumber(ScaffoldGraph * scaffoldGraph, long int contigCount, long int min, long int * outLinkNumber, long int * inLinkNumber);
int PrintScaffoldGraph(ScaffoldGraph * scaffoldGraph, long int contigCount, char * result);

int MergeTwoScaffoldSet(ScaffoldSet * previousScaffoldSet, long int & contigCount, ScaffoldSet * scaffoldSet);
ContigSequence * ReverseContigSequence(ContigSequence * contigSequence);

long int FindScaffoldGraphIndexOfScaffoldSet(ScaffoldSet * scaffoldSet, long int index, bool & orientation);
long int GetInLinkNumberOfScaffoldNode(ScaffoldGraph * scaffoldGraph, long int index);
long int GetOutLinkNumberOfScaffoldNode(ScaffoldGraph * scaffoldGraph, long int index);
#endif
