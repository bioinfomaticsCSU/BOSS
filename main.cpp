#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "scaffoldgraph.h"
#include "scaffolding.h"
#include "lp/lp_lib.h"

using namespace std;

int main(int argc, char *argv[]){

    time_t timep;
    time(&timep);
    ofstream ocout;
    ocout.open("LNSG_time_memory.txt");
    
    for(int i = 0; i<argc; i++){
        ocout<<argv[i]<<" ";
    }
    ocout<<endl;
    
    ocout<<"LNSG start time:"<<asctime(gmtime(&timep))<<endl;
      
    if((argc-3)%8!=0){
        cout<<"Please Input Correct Arguments!"<<endl;
        exit(0);
    }
    int libraryCount = (argc-3)/8;

    InputArg * inputArg = new InputArg[libraryCount];
    for(int i = 0; i < libraryCount; i++){
        inputArg[i].bamFileName = argv[2+8*i];
        inputArg[i].readLength = atoi(argv[2+8*i+1]);
        inputArg[i].insertsize = atoi(argv[2+8*i+2]);
        inputArg[i].std = atoi(argv[2+8*i+3]);
        inputArg[i].minEdgeWeight = atof(argv[2+8*i+4]);
        inputArg[i].minEdgeLinkNumber = atoi(argv[2+8*i+5]);
        inputArg[i].minRepetitiveCov = atof(argv[2+8*i+6]);
        inputArg[i].pairedRead = atoi(argv[2+8*i+7]);
    }
    
    char * contigSetFileName = argv[1];
    long int contigCount = 0;
    ContigSet * contigSet = GetContigSet(contigSetFileName, contigCount); 
    
    long int * scaffoldLength = NULL;
    long int * contigLength = new long int[contigCount];
    for(long int i = 0; i<contigCount; i++){
        contigLength[i] = strlen(contigSet[i].contig);
    }
    ScaffoldSet * scaffoldSet = InitScaffoldSet(contigSet, contigCount);
    ScaffoldToContig * scaffoldToContig = NULL;
    long int scaffoldCount = contigCount;
    ScaffoldGraphHead * scaffoldGraphHead = new ScaffoldGraphHead[libraryCount];
    
    for(int i = 0; i < libraryCount; i++){
        
        scaffoldToContig = GetScaffoldToContig(scaffoldSet, contigLength, contigCount, scaffoldCount); 
        
        scaffoldLength = GetScaffoldSetLength(scaffoldSet, scaffoldCount);
        
        BuildScaffoldGraphFromTwoBam(scaffoldSet, scaffoldLength, contigSet, contigLength, scaffoldToContig, scaffoldCount, scaffoldGraphHead + i, inputArg + i, argv[argc-1]);
        
        OptimizeScaffoldGraph((scaffoldGraphHead+i)->scaffoldGraph, scaffoldCount, inputArg[i].readLength, inputArg[i].insertsize, inputArg[i].std, inputArg[i].insertsize+lambda*inputArg[i].std, 0.4);
        
        scaffoldSet = OptimizeScaffoldSet(scaffoldSet, (scaffoldGraphHead+i)->scaffoldGraph, scaffoldCount, contigCount, contigLength, inputArg[i].insertsize, inputArg[i].std, lambda);
    }
    OutPutScaffoldSet(scaffoldSet, contigSet, contigCount, argv[argc-1]);
    
    time(&timep);
    ocout<<"LNSG end time:"<<asctime(gmtime(&timep))<<endl; 
    
    
    char * file = new char[100];
    long int pid = getpid();
    sprintf(file,"/proc/%ld/status",pid);
    ifstream icin;
    icin.open(file);
    char * line = new char[1000];
    while(icin.getline(line,1000)){
        ocout<<line<<endl;
    }
    

}

