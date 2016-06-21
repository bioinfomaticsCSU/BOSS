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
    
    if(argc == 1){
        cout<<"Command line:"<<endl;
        cout<<"./boss <contigs.fa> <bamfile_left.bam> <bamfile_right.bam> <read_length>"<<endl; 
        cout<<"<insert_size> <std> <min_weight> <min_number> <is_paired_end>"<<endl;
        cout<<"<edge_weight_method> <scaffold_file_name>"<<endl;
        cout<<endl;
        cout<<"Parameters:"<<endl;
	    cout<<"<contigs.fa>:"<<endl;
        cout<<'\t'<<"The file includes contigs produced by one assembler."<<endl;
	    cout<<"<bamfile_left.bam>:"<<endl;
        cout<<'\t'<<"Before scaffolding, users should using one mapping"<<endl;
        cout<<'\t'<<"tool(bowtie2, bwa or bowtie) to map left read library"<<endl; 
        cout<<'\t'<<"to contigs, and this will produce the file bamfile_left.bam."<<endl;
        cout<<"<bamfile_right.bam>:"<<endl;
        cout<<'\t'<<"Before scaffolding, users should using one mapping"<<endl;
        cout<<'\t'<<"tool(bowtie2, bwa or bowtie) to map right read library"<<endl; 
        cout<<'\t'<<"to contigs, and this will produce the file bamfile_right.bam."<<endl;
	    cout<<"<read_length>:"<<endl;
        cout<<'\t'<<"The length of read."<<endl;
	    cout<<"<insert_size>:"<<endl;
        cout<<'\t'<<"The insert size of read library."<<endl;
	    cout<<"<std_percentage>:"<<endl;
        cout<<'\t'<<"The percentage of standard deviation to insert size,"<<endl;
        cout<<'\t'<<"std = insert_size*std_percentage. In default,"<<endl;
        cout<<'\t'<<"std_percentage = 0.07."<<endl;
	    cout<<"<min_weight>:"<<endl;
        cout<<'\t'<<"One cutoff for removing suprious edgs in"<<endl;
        cout<<'\t'<<"the scaffold graph. In default, min_weight = 0.2."<<endl;
	    cout<<"<min_number>:"<<endl;
        cout<<'\t'<<"The minimum number of links between contigs. In default,"<<endl;
        cout<<'\t'<<"min_number = 2."<<endl;
	    cout<<"<is_paired_end>:"<<endl;
        cout<<'\t'<<"It is equal to 0 or 1, 1 represents that the read library is"<<endl;
        cout<<'\t'<<"paired-end, 0 represents that the read library is mate-paired."<<endl;
	    cout<<"<edge_weight_method>:"<<endl; 
        cout<<'\t'<<"It is equal to 0 or 1, 0 represents that the edge weight calculated"<<endl;
        cout<<'\t'<<"by arithmetic mean, 1 represents that the edge weight calculated"<<endl;
        cout<<'\t'<<"by geometric mean. In default, edge_weight_method = 0."<<endl;
	    cout<<"<scaffold_file_name>:"<<endl; 
        cout<<'\t'<<"Output file name, the file includes scaffolds produced by BOSS."<<endl;
        cout<<endl;
        cout<<"Example:"<<endl;
        cout<<'\t'<<"./boss contigs.fa left.bam right.bam 76 650 0.07 0.2 2 0 0 ecoli"<<endl;
        cout<<'\t'<<"This command will produce scaffolding file, ecoli_ScaffoldSet.fa."<<endl;
        exit(0);
    }
      
    time_t timep;
    time(&timep);
    ofstream ocout;
    char * bossInforFileName = new char[20];
    strcpy(bossInforFileName, "BOSS_infor.txt");
    ocout.open(bossInforFileName);
    
    for(int i = 0; i<argc; i++){
        ocout<<argv[i]<<" ";
    }
    ocout<<endl;
    
    ocout<<"BOSS start time:"<<asctime(gmtime(&timep))<<endl;
      
    if((argc-3)%9!=0){
        cout<<"Please Input Correct Arguments!"<<endl;
        exit(0);
    }
    int libraryCount = (argc-3)/9;

    InputArg * inputArg = new InputArg[libraryCount];
    for(int i = 0; i < libraryCount; i++){
        inputArg[i].bamFileName1 = argv[2+9*i];
        inputArg[i].bamFileName2 = argv[2+9*i+1];
        inputArg[i].readLength = atoi(argv[2+9*i+2]);
        inputArg[i].insertsize = atoi(argv[2+9*i+3]);
        inputArg[i].std = atof(argv[2+9*i+4])*inputArg[i].insertsize;
        inputArg[i].minEdgeWeight = atof(argv[2+9*i+5]);
        inputArg[i].minEdgeLinkNumber = atoi(argv[2+9*i+6]);
        inputArg[i].pairedRead = atoi(argv[2+9*i+7]);
        inputArg[i].edgeMeanMethod = atoi(argv[2+9*i+8]);
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
    
    long int ss = 0;
    
    for(int i = 0; i < libraryCount; i++){
       
        scaffoldToContig = GetScaffoldToContig(scaffoldSet, contigLength, contigCount, scaffoldCount); 
        
        scaffoldLength = GetScaffoldSetLength(scaffoldSet, scaffoldCount);
        
        BuildScaffoldGraphFromTwoBam(scaffoldSet, scaffoldLength, contigSet, contigLength, scaffoldToContig, scaffoldCount, scaffoldGraphHead + i, inputArg + i, argv[argc-1]);
        
        OptimizeScaffoldGraph((scaffoldGraphHead+i)->scaffoldGraph, scaffoldCount, inputArg[i].readLength, inputArg[i].insertsize, inputArg[i].std, inputArg[i].insertsize+lambda*inputArg[i].std, inputArg[i].edgeMeanMethod);
        
        scaffoldSet = OptimizeScaffoldSet(contigSet, scaffoldSet, (scaffoldGraphHead+i)->scaffoldGraph, scaffoldCount, contigCount, contigLength, inputArg[i].insertsize, inputArg[i].std, lambda);
        
    }
    
    OutPutScaffoldSet(scaffoldSet, contigSet, contigCount, argv[argc-1]);
    time_t timep1;
    time(&timep1);
    ocout<<"BOSS end time:"<<asctime(gmtime(&timep1))<<endl; 
    ocout<<"AllTime:"<<difftime(timep1, timep)<<endl;
    
}

