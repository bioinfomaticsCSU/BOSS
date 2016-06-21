#ifndef BUILDSCAFFOLDGRAPH_CPP_INCLUDED 
#define BUILDSCAFFOLDGRAPH_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <exception>

#include "scaffoldgraph.h"
#include "scaffolding.h"
#include "lp/lp_lib.h"


using namespace std;


int lambda = 3;
double scoreThreshold = 0.3;
int minMappedPairedReadNumber = 4;
double allAverageReadCoverage = 0;
double repeativeCutOff = 3;

long int * GetScaffoldSetLength(ScaffoldSet * temp, long int scaffoldCount){

    long int * scaffoldLength = new long int[scaffoldCount];
    
    for(long int j = 0; j<scaffoldCount; j++){
        if(temp->contigSequence == NULL){
            temp = temp->next;
            j--;
            continue;
        }
        scaffoldLength[j] = temp->length;
        temp = temp->next;
    }
    return scaffoldLength;

}

void GetScaffoldSetLength(ScaffoldSet * scaffoldSet, long int & scaffoldCount, long int * contigLength){

    ContigSequence * temp = NULL;
    long int i = 0;
    scaffoldCount = 0;
    while(scaffoldSet!=NULL){
        temp = scaffoldSet->contigSequence;
        scaffoldSet->length = 0;
        if(temp!=NULL){
            scaffoldCount++;
        }
        while(temp!=NULL){
            scaffoldSet->length = scaffoldSet->length + contigLength[temp->index] + temp->gapDistance;
            temp = temp->next;
        }
        scaffoldSet = scaffoldSet->next;
        i++;
    }

}

ScaffoldToContig * GetScaffoldToContig(ScaffoldSet * scaffoldSet, long int * contigLength, long int contigCount, long int & scaffoldCount){
    
    long int start = 0;
    long int len = 0;
    scaffoldCount = -1;
    ScaffoldToContig * scaffoldToContig = new ScaffoldToContig[contigCount];
    
    ContigSequence * contigSequence = NULL;
    while(scaffoldSet != NULL){
        contigSequence = scaffoldSet->contigSequence;
        start = 0;
        len = 0;
        if(contigSequence!=NULL){
            scaffoldCount++;
        }
        while(contigSequence!=NULL){
            scaffoldToContig[contigSequence->index].contigIndex = scaffoldCount;
            scaffoldToContig[contigSequence->index].start = start + len;
            len = contigLength[contigSequence->index] + contigSequence->gapDistance;
            start = scaffoldToContig[contigSequence->index].start;
            scaffoldToContig[contigSequence->index].orientation = contigSequence->orientation;
            contigSequence = contigSequence->next;
        }
        scaffoldSet = scaffoldSet->next;
    }
    
    scaffoldCount++;
    
    return scaffoldToContig;
    
}

long int GetScaffoldofContigPosition(ScaffoldToContig * scaffoldToContig, long int contigIndex, long int * contigLength, long int position, long int readLength){
    if(scaffoldToContig[contigIndex].orientation == 0){
        return position + scaffoldToContig[contigIndex].start;
    }
    if(scaffoldToContig[contigIndex].orientation == 1){
        return contigLength[contigIndex] - position - readLength + scaffoldToContig[contigIndex].start;
    }
}

bool GetScaffoldofContigOrientation(ScaffoldToContig * scaffoldToContig, long int contigIndex, bool orientation){
    if(scaffoldToContig[contigIndex].orientation == 0){
        return orientation;
    }else{
        return !orientation;
    }
}

bool ReverseComplement(char * temp1, char * temp2){
    long int len = strlen(temp1);
    long int i = 0;
    long int j = 0;
    for(i=0;i<len;i++){
        if(temp1[i]=='A'){
            temp2[len-1-i]='T';
        }else if(temp1[i]=='T'){
            temp2[len-1-i]='A';
        }else if(temp1[i]=='G'){
            temp2[len-1-i]='C';
        }else if(temp1[i]=='C'){
            temp2[len-1-i]='G';
        }else if(temp1[i]=='a'){
            temp2[len-1-i]='t';
        }else if(temp1[i]=='t'){
            temp2[len-1-i]='a';
        }else if(temp1[i]=='g'){
            temp2[len-1-i]='c';
        }else if(temp1[i]=='c'){
            temp2[len-1-i]='g';
        }else{
            temp2[len-1-i]='N';
        }
    }
    temp2[len]='\0';
    return true;
}

void AppendRight( char * temp, char * temp1, char * temp2,long int kmerLength){
     long int i = 0;
     long int j = 0;
     long int len1 = strlen(temp1);
     long int len2 = strlen(temp2);
     for(i =0;i<len1+len2-kmerLength+1;i++){
         if(i<len1){
             temp[i] = temp1[i];
         }else{
             temp[i] = temp2[kmerLength-1+j];
             j++;
         }
     }
     temp[i] = '\0';
}

void CopyContig(char *temp, char * temp1){
     long int i = 0;
     long int j = 0;
     long int len1 = strlen(temp1);
     for(i =0;i<len1;i++){
         temp[i] = temp1[i];
     }
     temp[i] = '\0';
}

double f(double x){
    return exp(-x*x/2);
}

double F(double a,double b,double ep){
    
    double h,s1=0,s2=(b-a)*(f(a)+f(b))/2;
    int n,k;
    for(int n=1;fabs(s1-s2)>ep;n*=2)
    {
        h=(b-a)/n;
        s1 = s2;
        s2 = 0;
        for(int k=0;k<n;++k)
        {
            s2 += h*f(a+(k+0.5)*h);
        }
        s2 = (s1+s2)/2;        
    }
    return s2*sqrt(1/(8*atan(1.0)));
}

void KMPGetNext(char * pattern, long int * next){
    next[0] = -1;
    long int k = -1;
    long int j = 0;
    while(pattern[j]!= '\0'){
        if(k!=-1  && pattern[k]!= pattern[j]){
            k=next[k];
        }        
        ++j;
        ++k;
        if(pattern[k]==pattern[j]){
            next[j]=next[k];
        }else{
            next[j]=k;
        }
    }
}

int KMPIndexOfContig(char * contig, char * pattern, long int * next){
    long int i = 0;
    long int j = 0;
    long int index = 0;
    long int len1 = strlen(contig);
    long int len2 = strlen(pattern);
    if(len1<len2){
        return -1;
    }
    KMPGetNext(pattern, next);
    while(i<len1&&j<len2){
        if(contig[i]==pattern[j]){
            i++;
            j++;
        }else{
            index += j-next[j];
            if(next[j]!=-1){
               j=next[j];
            }else{
                j=0;
                ++i;
            }
        }
    }
    if(pattern[j]=='\0'){
       return index;
    }else{
       return -1;  
    }
} 

int MergeSubContig(ContigSet * contigSet, long int contigCount){
    long int i = 0;
    long int j = 0;
    int count = 0;
    
    long int * contigLength = new long int[contigCount];
    for(i = 0; i<contigCount; i++){
        contigLength[i] = strlen(contigSet[i].contig);
    }
    
    i = 0;
    while(i<contigCount){
        if(contigSet[i].contig==NULL){
            i++;
            continue;
        }
        j = 0;
        while(j<contigCount){
            if(contigSet[j].contig==NULL || i==j || contigLength[i]<contigLength[j]){
                j++;
                continue;
            }
            long int len = contigLength[j];
            long int * next = new long int[len+1];
            long int p = KMPIndexOfContig(contigSet[i].contig,contigSet[j].contig,next);
            if(p!=-1){
                delete []contigSet[j].contig;
                contigSet[j].contig = NULL;
                delete []next;
                next = NULL;
                j++;
                continue;
            }                      
      
            char * tempReverseContig = new char[len+1];
            ReverseComplement(contigSet[j].contig,tempReverseContig);
            p = KMPIndexOfContig(contigSet[i].contig,tempReverseContig,next);
            if(p!=-1){
                delete []contigSet[j].contig;
                contigSet[j].contig = NULL;
                delete []next;
                next = NULL;
                delete []tempReverseContig;
                tempReverseContig = NULL;
                j++;
                continue;
            }
                        
            delete []tempReverseContig;
            delete []next;
            next = NULL; 
            j++;          
        }
        i++;

    }
    return 1;
}

ContigSet * GetContigSet(char * str, long int & contigNum){
    long int i = 0;
    long int j = 0;
    long int maxContigLength = 3000000;
    
    ifstream icin;
    icin.open(str);    
    if(icin==NULL){
        cout<<"Open File Wrong!"<<endl;
        exit(0);
    }    
    char * temp = new char[maxContigLength];  
    while(icin.getline(temp,maxContigLength)){
        if(temp[0]=='>'){
            contigNum++;
        }
    }  
    icin.close();
    
    ContigSet * contigSet = new ContigSet[contigNum];    
    ifstream icin1;
    icin1.open(str);
    i = -1;
    while(icin1.getline(temp,maxContigLength)){
        if(temp[0]=='>'){
            i++;
            continue;
        }
        long int len = strlen(temp);
        if(contigSet[i].contig==NULL){
            contigSet[i].contig = new char[len+1];
            CopyContig(contigSet[i].contig,temp);
            continue;
        }
        len = strlen(contigSet[i].contig) + len;
        
        char * temp1 = new char[len + 1];
        AppendRight(temp1,contigSet[i].contig,temp,1);
        delete []contigSet[i].contig;
        contigSet[i].contig = temp1;
        temp1 = NULL;  
    }
    
    delete [] temp;
    temp = NULL;
    return contigSet;
    
}

double ComputeVar(long int * distance, long int count){
    long int i = 0;
    long int all = 0;
    double avg = 0;
    for(i = 0; i<count; i++){
        all = distance[i] + all;
    }
    avg = (double)all/(double)count;
    
    double var = 0;
    
    for(i = 0; i<count; i++){
        var = var + pow((distance[i]-avg),2);
    }
    var = var/count;
    return sqrt(var)/avg;
    
}

double AddPositionIndex(ScaffoldGraph * scaffoldGraph, ReadMapPosition * readMapPosition, long int scaffoldCount){

    long int i = 0;
    long int j = 0;
    
    for(i=0;i<scaffoldCount;i++){
        
        ScaffoldEdge * scaffoldOutEdge = scaffoldGraph[i].outLink;
        while(scaffoldOutEdge!=NULL){
            if(scaffoldOutEdge->edgeWeight<minMappedPairedReadNumber){
                scaffoldOutEdge = scaffoldOutEdge->next;
                continue;
            }
            MapPosition * tempMapPosition = scaffoldOutEdge->mapPosition;
            while(tempMapPosition != NULL){
                if(tempMapPosition->orientation == 1){
                    readMapPosition[i].rightIndex[tempMapPosition->pos]++;
                }else{
                    readMapPosition[i].leftIndex[tempMapPosition->pos]++;
                }
                
                tempMapPosition = tempMapPosition->next;
            } 
            scaffoldOutEdge = scaffoldOutEdge->next;
        }
     
    }
    
    
    
}


double MapFitDistribution(long int * index, long int * noMapIndex, MapPosition * tempMapPosition, long int * distance, double correctP, long int length, long int length1, long int readLength, long int gap, long int insertsize, long int std, double & cgap, double & fitN, double & fitD, double & mapP, int token, double avgP){

    long int i = 0;
    long int j = 0;
    
    double p = 0;
             
    double max = 4;
    double min = -4;
    
    long int continousMaxGap = 0;
    long int tempMaxGap = 0;
     
    mapP = 0;
    
    i = 0;
    long int edgeCount = 0;
    
    
    long int * realCount = new long int[insertsize + lambda*std+1];
    double * expectCount = new double[insertsize + lambda*std+1];
    double * expectPro = new double[insertsize + lambda*std+1];
    double * realPro = new double[insertsize + lambda*std+1];
    
    for(i = 0; i<insertsize +lambda*std; i++){
        realCount[i] = 0;
        expectCount[i] = 0;
        expectPro[i] = 0;
        realPro[i] = 0;
    }
    
    
    i = 0;
    MapPosition * temp11 = tempMapPosition;
    while(tempMapPosition!=NULL){
        if(gap+distance[i]<=insertsize+lambda*std && gap+distance[i]>=insertsize-lambda*std){
            realCount[length - tempMapPosition->pos - 1]++;
            edgeCount++;
        }
        i++;
        tempMapPosition = tempMapPosition->next;
    }
    
    if(edgeCount == 0){
        return 0;
    }
    
    long int * distanceT = new long int[edgeCount+1];
    distanceT[edgeCount] = 0;
    i = 0;
    j = 0; 
    long int avgDistanceT = 0;
    while(temp11!=NULL){
        if(gap+distance[i]<=insertsize+lambda*std && gap+distance[i]>=insertsize-lambda*std){
            distanceT[j] = distance[i];
            avgDistanceT = distanceT[j] + avgDistanceT;
            j++;
        }
        i++;
        temp11 = temp11->next;
    }
    
    avgDistanceT = avgDistanceT/j;
    
    long int posMap = 0;
    long int posMapCount = 0;
    
    
    double bP = 1;
    double tempBP = 0;
    double eBP = 0;
    double stdBP = 0;
    double varBP = 0;
    double allVarBP = 0;
    
    double aa = 0;
    double bb = 0;
    
    long int allSingleMapped = 0;
    double p1 = 0;
    double minVarBP = 0;
    
    long int allNoMapped = 0;
       
    i = readLength -1;
    long int AllEdgeCount = edgeCount;
    AllEdgeCount = 0;
    edgeCount = 0;
    for(i = i; i<insertsize+lambda*std-readLength -gap && i<=length -1-readLength; i++){
        
        edgeCount = edgeCount + realCount[i];
        long int tempNumber = 0;
        if(token == 0){
            tempNumber = index[length-1-i] + realCount[i];
            
        }else{
            tempNumber = avgP*realCount[i];
            
        }
        if(tempNumber == 0){
            tempMaxGap++;
            continue;
        }else{
            
        }
        AllEdgeCount++;
        if(realCount[i]!=0 && realCount[i]+noMapIndex[length-1-i]!=index[length-1-i]){
            aa++;
        }else if(realCount[i]!=0){
            bb++;
        }
          
        max = 4;
        min = -4;
        double end = (double)(i+gap+length1-insertsize)/(double)std;
        double start = (double)(i + gap + readLength -insertsize)/(double)std;
        if(max<end){
            end = max;
        }
        if(min>start){
            start = min;
        }
        double ndP = 0;
        if(start<end){
            ndP = (double)F(start,end);
        }
        
        
        expectCount[i] = (double)(1-correctP)*(index[length-1-i])*ndP;
        
        p1 = p1 + (double)(realCount[i])*ndP;
        allNoMapped = allNoMapped + noMapIndex[length-1-i];
        
        p = p + expectCount[i];
        
        allVarBP = allVarBP + expectCount[i]*(1-(double)(1-correctP)*ndP);
   
        
    }
    

    long int realGap = labs(p - edgeCount);
    
    mapP = tempMaxGap;
    
    if(p>edgeCount){
        p1 = p1/(double)edgeCount;
        p = edgeCount/p;
    }else{
        //p1 = p1/(double)edgeCount;
        if(p == 0 && edgeCount == 0){
            p = 0;
        }else{
            p = p/edgeCount;
        }
    }
    
    fitD = 1;
    
    fitN = p;
    
    

    delete [] realCount;
    delete [] expectCount;
    delete [] expectPro;
    delete [] realPro;
    delete [] distanceT;
   
}

double MapFitDistribution1(long int * index, long int * noMapIndex, MapPosition * tempMapPosition, long int * distance, double correctP, long int length, long int length1, long int readLength, long int gap, long int insertsize, long int std, double & cgap, double & fitN, double & fitD, double & mapP, int token, double avgP){
    long int i = 0;
    long int j = 0;
    
    double p = 0;
       
    
       
    double max = 4;
    double min = -4;
    
    long int continousMaxGap = 0;
    long int tempMaxGap = 0;
    
    mapP = 0;
    
    
    long int * realCount = new long int[insertsize + lambda*std];
    double * expectCount = new double[insertsize + lambda*std];
    double * expectPro = new double[insertsize + lambda*std];
    double * realPro = new double[insertsize + lambda*std];
    
    for(i = 0; i<insertsize +lambda*std; i++){
        realCount[i] = 0;
        expectCount[i] = 0;
        expectPro[i] = 0;
        realPro[i] = 0;
        
    }
    
    long int edgeCount = 0;
    i = 0;
    MapPosition * temp11 = tempMapPosition;
    while(tempMapPosition!=NULL){
        if(gap+distance[i]<=insertsize+lambda*std && gap+distance[i]>=insertsize-lambda*std){
            realCount[tempMapPosition->pos]++;
            edgeCount++;
        }
        i++;
        tempMapPosition = tempMapPosition->next;
    }
    
    if(edgeCount == 0){
        return 0;
    }
    
    long int * distanceT = new long int[edgeCount+1];
    distanceT[edgeCount] = 0;
    i = 0;
    j = 0; 
    long int avgDistanceT = 0;
    
    while(temp11!=NULL){
        if(gap+distance[i]<=insertsize+lambda*std && gap+distance[i]>=insertsize-lambda*std){
            distanceT[j] = distance[i];
            avgDistanceT = avgDistanceT + distanceT[j];
            j++;
        }
        i++;
        temp11 = temp11->next;
    }
   
    avgDistanceT = avgDistanceT/j;
    
    i = 0;
    long int posMap = 0;
    long int posMapCount = 0;
    
    double tempBP = 0;
    double eBP = 0;
    double stdBP = 0;
    double varBP = 0;
    double allVarBP = 0;
    
    double aa = 0;
    double bb = 0;
    
    long int allSingleMapped = 0;
    double p1 = 0;
    
    double minVarBP = 0;
    long int allNoMapped = 0;
    
    i = 0;
    long int AllEdgeCount = edgeCount;
    AllEdgeCount = 0;
    edgeCount = 0;
    for(i = i; i<insertsize+lambda*std-readLength -gap && i<=length -1 -readLength; i++){
        
        edgeCount = edgeCount + realCount[i];
        long int tempNumber = 0;
        if(token == 0){
            tempNumber = index[i] +avgP*realCount[i];
            
        }else{
            tempNumber = avgP*realCount[i];
            
        }
        
        if(tempNumber== 0){
            tempMaxGap++;
            continue;
        }else{
            
        }
        AllEdgeCount++;
        if(realCount[i]!=0 && realCount[i]+noMapIndex[i]!=index[i]){
            aa++;
        }else if(realCount[i]!=0){
            bb++;
        }
          
        max = 4;
        min = -4;
        double end = (double)(i+gap+readLength + length1-insertsize)/(double)std;
        double start = (double)(i + gap + readLength -insertsize)/(double)std;
        if(max<end){
            end = max;
        }
        if(min>start){
            start = min;
        }
        double ndP = 0;
        if(start<end){
            ndP = (double)F(start,end);
        }
        
        expectCount[i] = (double)(1-correctP)*(index[i])*ndP;
       
        p1 = p1 + (double)(realCount[i])*ndP;
        allNoMapped = allNoMapped + noMapIndex[i];
       
        p = p + expectCount[i];
        
        allVarBP = allVarBP + expectCount[i]*(1-(double)(1-correctP)*ndP);
        
        
    }
    
    mapP = tempMaxGap;

    if(p>edgeCount){
        p1 = p1/(double)edgeCount;
        p = edgeCount/p;
    }else{
        //p1 = p1/(double)edgeCount;
        if(p == 0 && edgeCount == 0){
            p = 0;
        }else{
            p = p/edgeCount;
        }
    }
    
    fitD = 1; 
    fitN = p;
    

    delete [] realCount;
    delete [] expectCount;
    delete [] expectPro;
    delete [] realPro;
    delete [] distanceT;
}



int AddScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int refId, bool refStrand, long int pos, long int mateRefId, bool mateRefStrand, long int matePos){
    
    ScaffoldEdge * tempScaffoldEdge = scaffoldGraph[refId].outLink;
    
    while(tempScaffoldEdge!=NULL){
        if(tempScaffoldEdge->contigIndex==mateRefId){
            break;
        }
        tempScaffoldEdge = tempScaffoldEdge->next;
    }
    
    if(tempScaffoldEdge!=NULL){
        MapPosition * tempMapPosition = tempScaffoldEdge->mapPosition;
        while(tempMapPosition!=NULL){
            if(tempMapPosition->pos==pos && tempMapPosition->matePos==matePos){
                break;
            }
            tempMapPosition = tempMapPosition->next;
        }
        if(tempMapPosition==NULL){
            tempMapPosition = new MapPosition;
            tempMapPosition->pos = pos;
            tempMapPosition->matePos= matePos;
            tempMapPosition->orientation = refStrand;
            tempMapPosition->mateOrientation = mateRefStrand;
            tempMapPosition->next = tempScaffoldEdge->mapPosition->next;
            tempScaffoldEdge->mapPosition->next = tempMapPosition;  
            tempScaffoldEdge->edgeWeight++;         
        }
    }else{
        
        tempScaffoldEdge = new ScaffoldEdge;
        tempScaffoldEdge->contigIndex = mateRefId;
        tempScaffoldEdge->mapPosition = new MapPosition;
        tempScaffoldEdge->mapPosition->pos = pos;
        tempScaffoldEdge->mapPosition->matePos= matePos;
        tempScaffoldEdge->mapPosition->orientation = refStrand;
        tempScaffoldEdge->mapPosition->mateOrientation = mateRefStrand;
        
        if(scaffoldGraph[refId].outLink!=NULL){
            tempScaffoldEdge->next = scaffoldGraph[refId].outLink->next;
            scaffoldGraph[refId].outLink->next = tempScaffoldEdge;
        }else{
            scaffoldGraph[refId].outLink = tempScaffoldEdge;
        }

        tempScaffoldEdge->edgeWeight++;
        
    }

}

void DeleteMapPositionOfEdge(ScaffoldEdge * edge){
    
    MapPosition * tempMapPosition = edge->mapPosition;
    edge->mapPosition = NULL;
    MapPosition * pre = NULL;
    while(tempMapPosition != NULL){
        pre = tempMapPosition;
        tempMapPosition = tempMapPosition->next;
        pre->next = NULL;
        delete pre;
    }
    
}

int DeleteScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, bool orientation, long int cutOff, long int maxIndex){
    
    ScaffoldEdge * temp = scaffoldGraph[index].outLink;
    ScaffoldEdge * tempMate = NULL;
    ScaffoldEdge * previous = NULL;
    ScaffoldEdge * matePrevious = NULL;
    long int mateIndex = -1;
    
    while(temp!=NULL){        
        if(temp->mapPosition->orientation == orientation && scaffoldGraph[temp->contigIndex].length > cutOff && temp->contigIndex != maxIndex){
            
            mateIndex = temp->contigIndex;
            tempMate = scaffoldGraph[mateIndex].outLink;
            matePrevious = NULL;
            while(tempMate != NULL){
                if(tempMate->contigIndex == index){
                    if(matePrevious == NULL){
                        scaffoldGraph[mateIndex].outLink = tempMate->next;
                    }else{
                        matePrevious->next = tempMate->next;
                    }
                    DeleteMapPositionOfEdge(tempMate);
                    tempMate->next = NULL;
                    delete tempMate;
                    tempMate = NULL;
                    break;
                }
                matePrevious = tempMate;
                tempMate = tempMate->next;
            }
            
                        
            if(previous==NULL){
                scaffoldGraph[index].outLink = temp->next;
                DeleteMapPositionOfEdge(temp);
                temp->next = NULL;
                delete temp;
                temp = scaffoldGraph[index].outLink;
                continue;                  
            }
            previous->next = temp->next;
            temp->next = NULL;
            DeleteMapPositionOfEdge(temp);
            delete temp;
            temp = previous->next;
            continue;
        }
        previous = temp;
        temp = temp->next;
    } 
    
}

int DeleteScaffoldEdgeShort(ScaffoldGraph * scaffoldGraph, long int index, bool orientation, long int cutOff, long int maxIndex){
    
    ScaffoldEdge * temp = scaffoldGraph[index].outLink;
    ScaffoldEdge * tempMate = NULL;
    ScaffoldEdge * previous = NULL;
    ScaffoldEdge * matePrevious = NULL;
    long int mateIndex = -1;
    
    while(temp!=NULL){        
        if(temp->mapPosition->orientation == orientation && scaffoldGraph[temp->contigIndex].length < cutOff && temp->contigIndex != maxIndex){
          
            mateIndex = temp->contigIndex;
            tempMate = scaffoldGraph[mateIndex].outLink;
            matePrevious = NULL;
            while(tempMate != NULL){
                if(tempMate->contigIndex == index){
                    if(matePrevious == NULL){
                        scaffoldGraph[mateIndex].outLink = tempMate->next;
                    }else{
                        matePrevious->next = tempMate->next;
                    }
                    DeleteMapPositionOfEdge(tempMate);
                    tempMate->next = NULL;
                    delete tempMate;
                    tempMate = NULL;
                    break;
                }
                matePrevious = tempMate;
                tempMate = tempMate->next;
            }
            
                        
            if(previous==NULL){
                scaffoldGraph[index].outLink = temp->next;
                DeleteMapPositionOfEdge(temp);
                temp->next = NULL;
                delete temp;
                temp = scaffoldGraph[index].outLink;
                continue;                  
            }
            previous->next = temp->next;
            DeleteMapPositionOfEdge(temp);
            temp->next = NULL;
            delete temp;
            temp = previous->next;
            continue;
        }
        previous = temp;
        temp = temp->next;
    } 
    
}

int DeleteScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, bool orientation){
    
    ScaffoldEdge * temp = scaffoldGraph[index].outLink;
    ScaffoldEdge * tempMate = NULL;
    ScaffoldEdge * previous = NULL;
    ScaffoldEdge * matePrevious = NULL;
    long int mateIndex = -1;
    
    while(temp!=NULL){        
        if(temp->mapPosition->orientation == orientation){
            
            mateIndex = temp->contigIndex;
            tempMate = scaffoldGraph[mateIndex].outLink;
            matePrevious = NULL;
            while(tempMate != NULL){
                if(tempMate->contigIndex == index){
                    if(matePrevious == NULL){
                        scaffoldGraph[mateIndex].outLink = tempMate->next;
                    }else{
                        matePrevious->next = tempMate->next;
                    }
                    DeleteMapPositionOfEdge(tempMate);
                    tempMate->next = NULL;
                    delete tempMate;
                    tempMate = NULL;
                    break;
                }
                matePrevious = tempMate;
                tempMate = tempMate->next;
            }
            
                        
            if(previous==NULL){
                scaffoldGraph[index].outLink = temp->next;
                DeleteMapPositionOfEdge(temp);
                temp->next = NULL;
                delete temp;
                temp = scaffoldGraph[index].outLink;
                continue;                  
            }
            previous->next = temp->next;
            DeleteMapPositionOfEdge(temp);
            temp->next = NULL;
            delete temp;
            temp = previous->next;
            continue;
        }
        previous = temp;
        temp = temp->next;
    } 
    
}

int DeleteSpecialScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int index, long int index1){
    
    ScaffoldEdge * temp = scaffoldGraph[index].outLink;
    ScaffoldEdge * previous = NULL;
    while(temp!=NULL){        
        if(temp->contigIndex == index1){
            if(previous==NULL){
                scaffoldGraph[index].outLink = temp->next;
                DeleteMapPositionOfEdge(temp);
                temp->next = NULL;
                delete temp;
                temp = scaffoldGraph[index].outLink;
                continue;                  
            }
            previous->next = temp->next;
            DeleteMapPositionOfEdge(temp);
            temp->next = NULL;
            delete temp;
            temp = previous->next;
            continue;
        }
        previous = temp;
        temp = temp->next;
    } 
    
    temp = scaffoldGraph[index].inLink;
    previous = NULL;
    
    while(temp!=NULL){        
        if(temp->contigIndex == index1){
            if(previous==NULL){
                scaffoldGraph[index].inLink = temp->next;
                DeleteMapPositionOfEdge(temp);
                temp->next = NULL;
                delete temp;
                temp = scaffoldGraph[index].inLink;
                continue;                  
            }
            previous->next = temp->next;
            DeleteMapPositionOfEdge(temp);
            temp->next = NULL;
            delete temp;
            temp = previous->next;
            continue;
        }
        previous = temp;
        temp = temp->next;
    } 
    
}


double ComputeReadDistribution(MapPosition * read, long int * readMapPosition, long int * distance, long int contigLength, long int contigLength1, long int readLength, long int gap, long int insertsize, long int std1, double p){
    
    long int i = 0;
    long int j = 0;
    
    int * index = new int[contigLength];
    
    for(i=0;i<contigLength;i++){
        index[i] = 0;
    }
    
    long int count = 0;
    long int start = insertsize - lambda*std1 - gap;
    long int end = insertsize + lambda*std1 - gap;
    if(start>contigLength1){
        start = start - contigLength1;
    }else{
        start = 0;
    }
    if(end>contigLength){
        end = contigLength;
    }
    count = end - start -readLength;
    
    double expectation = p*count;
    double std = sqrt(p*(1-p)*count);
    
    long int readNumber = 0;
    MapPosition * tempRead = read;
    
    j = 0;
    i = 0;
    while(tempRead!=NULL){
        
        if(index[tempRead->pos] == 0){
            index[tempRead->pos]++;
            readNumber++;
            i++;
        }            
          
        j++;       
        tempRead = tempRead->next;
    }
    
    for(i = 0 ; i< end; i++){
        if(read->orientation == 1){
            if(readMapPosition[contigLength - 1 - i] >0){
                readNumber++;
            }
        }else{
            if(readMapPosition[i] >0){
                readNumber++;
            }
        }
    }
    
    long int variation = readNumber-expectation;
    if(variation<0){
        variation = -variation;
    }
    
    double score = (double)variation/(double)(lambda*std);    
    return score;

}

int ComputeTruePairedMapping(BamAlignment * alignment, long int len, long int len1, long int readLength, long int insertsize, long int std){
    
    if(alignment->IsReverseStrand()==0&&alignment->IsMateReverseStrand()==1){
        if(alignment->Position - alignment->MatePosition + readLength > insertsize + lambda*std || alignment->Position - alignment->MatePosition + readLength < insertsize - lambda*std){
            return 0;
        }
        return 1;
    }
            
    if(alignment->IsReverseStrand()==1&&alignment->IsMateReverseStrand()==0){
        if(alignment->MatePosition - alignment->Position + readLength > insertsize + lambda*std || alignment->MatePosition - alignment->Position + readLength < insertsize - lambda*std){
            return 0;
        }
        return 1;
    }
    return 0;
}

int ComputeTruePairedMapping(PairedReadMappedData * alignment, long int readLength, long int insertsize, long int std){
    
    if(alignment->orientationRight==0&&alignment->orientationLeft==1){
        if(alignment->rightPosition - alignment->leftPosition + readLength > insertsize + lambda*std || alignment->rightPosition - alignment->leftPosition + readLength < insertsize - lambda*std){
            return 0;
        }
        return 1;
    }
            
    if(alignment->orientationRight==1&&alignment->orientationLeft==0){
        if(alignment->leftPosition - alignment->rightPosition + readLength > insertsize + lambda*std || alignment->leftPosition - alignment->rightPosition + readLength < insertsize - lambda*std){
            return 0;
        }
        return 1;
    }
    
    return 0;
}

MapPosition * OptimizeMapPosition(MapPosition * mapPosition, long int & edgeWeight, char * fileName){
    long int i = 0;
    long int j = 0;
    long int oriLeft0 = 0;
    long int oriLeft1 = 0;
    long int oriRight0 = 0;
    long int oriRight1 = 0;
    
    MapPosition * tempMapPosition = mapPosition;            
    while(tempMapPosition!=NULL){                  
        if(tempMapPosition->orientation==0 && tempMapPosition->mateOrientation==1){
            oriLeft1++;
        }
        if(tempMapPosition->orientation==0 && tempMapPosition->mateOrientation==0){
            oriLeft0++;
        }
        if(tempMapPosition->orientation==1 && tempMapPosition->mateOrientation==1){
            oriRight1++;
        }
        if(tempMapPosition->orientation==1 && tempMapPosition->mateOrientation==0){
            oriRight0++;
        }
        tempMapPosition = tempMapPosition->next;
    }
    
    
    long int maxLeft = 0;
    long int maxRight = 1;
    long int max = oriLeft1;
    if(max<oriLeft0){
        max = oriLeft0;
        maxLeft = 0;
        maxRight = 0;
    }
    if(max<oriRight1){
        max = oriRight1;
        maxLeft = 1;
        maxRight = 1;
    }
    if(max<oriRight0){
        max = oriRight0;
        maxLeft = 1;
        maxRight = 0;
    }

    MapPosition * previous = NULL;
    MapPosition * first = mapPosition;
    tempMapPosition = mapPosition;
    long int ss = 0;
    while(tempMapPosition!=NULL){                  
        if(tempMapPosition->orientation!=maxLeft || tempMapPosition->mateOrientation!=maxRight){
            if(previous==NULL){
                tempMapPosition = tempMapPosition->next;
                first->next = NULL;
                first = tempMapPosition;
            }else{
                previous->next = tempMapPosition->next;
                tempMapPosition->next = NULL;
                tempMapPosition = previous->next;
            }
            
            edgeWeight--;           
        }else{
            previous = tempMapPosition;
            tempMapPosition = tempMapPosition->next;
        }      
        ss++;
        
    }
    
    return first;
    
}

long int GetGapDistance(long int * distance, long int count, long int insertsize){
    
    long int i = 0;
    long int j = 0;
    long int all = 0;
    
    for(i=0;i<count;i++){
        all = all + distance[i];
    }
    all = labs(insertsize-all/count);
    return all;
}


void * EqualScaffoldEdge(ScaffoldGraph * scaffoldGraph, long int contigCount, bool edgeWeightMethod){
    long int i = 0;
    long int j = 0;
    long int index = -1;
    
    bool * visited = new bool[contigCount];
    for(i=0;i<contigCount;i++){
        visited[i] = false;
    }
    
    
    ScaffoldEdge * temp = NULL;
    ScaffoldEdge * temp1 = NULL;
    for(i=0;i<contigCount;i++){
        temp = scaffoldGraph[i].outLink;
        while(temp!=NULL){
            index = temp->contigIndex;
            temp1 = scaffoldGraph[index].outLink;
            if(visited[index]==true){
                temp = temp->next;
                continue;
            }
            while(temp1!=NULL){
                if(temp1->contigIndex==i){
                    if(edgeWeightMethod == false){
                        temp1->fitNumber = (temp1->fitNumber + temp->fitNumber)/2;
                        temp->fitNumber = temp1->fitNumber;
                    }else{
                        temp1->fitNumber = sqrt(temp1->fitNumber*temp->fitNumber);
                        temp->fitNumber = temp1->fitNumber;
                    }
                    break;
                }
                temp1 = temp1->next;
            }
            temp = temp->next;
        }
        visited[i] = true;
    }
    
    delete [] visited;
    visited = NULL;
}

int GetOutAndInLinkNumber(ScaffoldGraph * scaffoldGraph, long int contigCount, long int min, long int * outLinkNumber, long int * inLinkNumber){
    ScaffoldEdge * temp = NULL;
    long int i = 0;
    for(i = 0; i<contigCount; i++){
        temp = scaffoldGraph[i].outLink;
        while(temp!=NULL){
            if(temp->mapPosition->orientation == 1 && temp->edgeWeight >= min){
                outLinkNumber[i]++;
            }else if(temp->mapPosition->orientation == 0 && temp->edgeWeight >= min){
                inLinkNumber[i]++;
            }
            temp = temp->next;
        }
    }
    return 1;
    
}

int PrintScaffoldGraph(ScaffoldGraph * scaffoldGraph, long int contigCount, char * result){
    long int i = 0;
    long int j = 0;
    ofstream ocout;
    ocout.open(result);
    ScaffoldEdge * tempEdge = NULL;
    for(i = 0; i<contigCount; i++){
        tempEdge = scaffoldGraph[i].outLink;
        ocout<<i<<":";
        while(tempEdge!=NULL){
            ocout<<tempEdge->contigIndex<<","<<tempEdge->fitNumber<<","<<tempEdge->edgeWeight<<","<<tempEdge->gapDistance<<"--";
            tempEdge = tempEdge->next;
        }
        tempEdge = scaffoldGraph[i].inLink;
        ocout<<"--in--";
        while(tempEdge!=NULL){
            ocout<<tempEdge->contigIndex<<","<<tempEdge->fitNumber<<","<<tempEdge->edgeWeight<<","<<tempEdge->gapDistance<<"--";
            tempEdge = tempEdge->next;
        }
        ocout<<endl;
    }
}

ContigSequence * ReverseContigSequence(ContigSequence * contigSequence){

    long int i = 0;
    ContigSequence * temp = NULL;
    ContigSequence * previous = NULL;
    long int gapDistance = 0;
    long int gapDistance1 = 0;
    while(contigSequence!=NULL){
        
        if(contigSequence->orientation==0){
            contigSequence->orientation = 1;
        }else{
            contigSequence->orientation = 0;
        }
        
        temp = contigSequence->next;
        contigSequence->next = previous;
        gapDistance1 = contigSequence->gapDistance;
        contigSequence->gapDistance = gapDistance;
        gapDistance = gapDistance1;
        previous = contigSequence;       
        contigSequence = temp;
                        
    }

    return previous;

}


int MergeTwoScaffoldSet(ScaffoldSet * previousScaffoldSet, long int & contigCount, ScaffoldSet * scaffoldSet){
    long int i = 0;
    long int j = 0;
    
    bool * visited = new bool[contigCount];
    for(i = 0; i<contigCount; i++){
        visited[i] = false;
    }
    ScaffoldSet * last = NULL;
    j = 0;
    while(scaffoldSet != NULL){

        ContigSequence * contigSequence = scaffoldSet->contigSequence;
        
        contigSequence = scaffoldSet->contigSequence;
        
        ContigSequence * previous = NULL;
        
        while(contigSequence!=NULL){
            i = 0;
            ScaffoldSet * temp = previousScaffoldSet;
            visited[contigSequence->index] = true;
            while(temp!=NULL){
                if(i == contigSequence->index){
                    break;
                }
                i++;
                temp = temp->next;
            }
            
            if(contigSequence->orientation == 1){
                temp->contigSequence = ReverseContigSequence(temp->contigSequence);    
            }
            
            ContigSequence * tempContigSequence = temp->contigSequence;
            while(tempContigSequence->next!=NULL){
                tempContigSequence = tempContigSequence->next;
            }
            
            tempContigSequence->gapDistance = contigSequence->gapDistance;
            if(previous!=NULL){
                tempContigSequence->next = contigSequence->next;
                previous->next = temp->contigSequence;
            }else{
                tempContigSequence->next = contigSequence->next;
                scaffoldSet->contigSequence = temp->contigSequence;
            }
            previous = tempContigSequence;
            contigSequence = tempContigSequence->next;
        }
        j++;
        last = scaffoldSet;
        scaffoldSet = scaffoldSet->next;
    }
    i = 0;
    while(previousScaffoldSet!=NULL){
        if(visited[i] == false && previousScaffoldSet->contigSequence!=NULL){
            last->next = previousScaffoldSet;
            last = last->next;
            previousScaffoldSet = last->next;
            last->next = NULL;
            j++;
            i++;
            continue;
        }
        i++;
        previousScaffoldSet = previousScaffoldSet->next;
    }
    contigCount = j;
    
    delete [] visited;
    visited = NULL;
}

long int GetOutLinkNumberOfScaffoldNode(ScaffoldGraph * scaffoldGraph, long int index){
    long int count = 0;
    ScaffoldEdge * edge = scaffoldGraph[index].outLink;
    while(edge != NULL){
        count++;
        edge = edge->next;
    }
    return count;
}

long int GetInLinkNumberOfScaffoldNode(ScaffoldGraph * scaffoldGraph, long int index){
    long int count = 0;
    ScaffoldEdge * edge = scaffoldGraph[index].inLink;
    while(edge != NULL){
        count++;
        edge = edge->next;
    }
    return count;
}

void SetScaffoldToContigIndex(ScaffoldToContig * scaffoldToContig, long int index, long int newIndex, long int start, int orientation, long int contigCount){
    long int i = 0;
    for(i = 0; i<contigCount;i++){
        if(scaffoldToContig[i].contigIndex == index){
            scaffoldToContig[i].contigIndex = newIndex;
            scaffoldToContig[i].start = start + scaffoldToContig[i].start;
            scaffoldToContig[i].orientation = orientation * scaffoldToContig[i].orientation;
        }
    }
}

long int FindScaffoldGraphIndexOfScaffoldSet(ScaffoldSet * scaffoldSet, long int index, bool & orientation){
    long int i = 0;
    while(scaffoldSet != NULL){
        ContigSequence * temp = scaffoldSet->contigSequence;
        while(temp!=NULL){
            if(temp->index == index){
                orientation = temp->orientation;
                return i;
            }
            temp = temp->next;
        }
        scaffoldSet = scaffoldSet->next;
        i++;
    }
    return -1;
}


int OptimizeScaffoldGraph(ScaffoldGraph * scaffoldGraph, long int contigCount, long int readLength, long int insertsize, long int std, long int contigCutOff, bool edgeWeightMethod){
    
    long int i = 0;
    long int j = 0;
    
    long int * leftNumber = new long int[contigCount];
    long int * rightNumber = new long int[contigCount];
    long int * leftIndex = new long int[contigCount];
    long int * rightIndex = new long int[contigCount];
    bool * leftOrientation = new bool[contigCount];
    bool * rightOrientation = new bool[contigCount];
    
    
    for(i=0;i<contigCount;i++){       
        leftNumber[i]=0;
        rightNumber[i]=0;
        leftIndex[i]=-1;
        rightIndex[i]=-1;
        leftOrientation[i] = false;
        rightOrientation[i] = false;       
    }
      
    EqualScaffoldEdge(scaffoldGraph,contigCount,edgeWeightMethod);
       
    for(i=0;i<contigCount;i++){
        ScaffoldEdge * temp = scaffoldGraph[i].outLink;
        ScaffoldEdge * previous = NULL;
        
        long int allRightEdgeWeight = 0;
        long int allLeftEdgeWeight = 0;
        
        bool lengthToken = false;
        if(scaffoldGraph[i].length<contigCutOff){
            lengthToken = true;
        }
        
        while(temp!=NULL){
            if(temp->edgeWeight>=minMappedPairedReadNumber && temp->mapPosition->orientation == 0){
                allLeftEdgeWeight = allLeftEdgeWeight + temp->edgeWeight;      
            }
            if(temp->edgeWeight>=minMappedPairedReadNumber && temp->mapPosition->orientation == 1){
                allRightEdgeWeight = allRightEdgeWeight + temp->edgeWeight;
            }
            temp = temp->next;
        }
        
        temp = scaffoldGraph[i].outLink;
        while(temp!=NULL){
            
            
            bool numIndex = false;
            
            if(temp->mapPosition->orientation == 0){
                temp->fitDistribution = (double)(temp->edgeWeight)/(double)(allLeftEdgeWeight);
            }else{
                temp->fitDistribution = (double)(temp->edgeWeight)/(double)(allRightEdgeWeight);
            }
            
            
            if(temp->fitDistribution == 1 && temp->fitNumber > 0.1){
                numIndex = true;
                if(temp->fitNumber < scoreThreshold){
                    temp->fitNumber = scoreThreshold + 0.00001; 
                }
            }
            
            if((temp->fitNumber < scoreThreshold && numIndex == false) || temp->edgeWeight < minMappedPairedReadNumber || temp->fitNumber < 0){
                
                long int index = temp->contigIndex;
                
                ScaffoldEdge * temp1 = scaffoldGraph[index].outLink;
                ScaffoldEdge * previous1 = NULL;
                while(temp1!=NULL){
                    if(temp1->contigIndex==i){
                        if(previous1==NULL){
                            scaffoldGraph[index].outLink = temp1->next;
                            temp1->next = NULL;
                            temp1 = scaffoldGraph[index].outLink;
                            break;
                        }
                        previous1->next = temp1->next;
                        temp1->next = NULL;
                        break;
                    }
                    previous1 = temp1;
                    temp1 = temp1->next;
                }
                if(previous==NULL){
                    scaffoldGraph[i].outLink = temp->next;
                    temp->next = NULL;
                    temp = scaffoldGraph[i].outLink;
                    continue;                  
                }
                previous->next = temp->next;
                temp->next = NULL;
                temp = previous->next;
                continue;
            }
            
            previous = temp;
            temp = temp->next;           
        }               
    }

    
    
    double * edgeWeightLeft = new double[contigCount]; //the largest score of all edges;
    double * edgeWeightRight = new double[contigCount];
    double * edgeWeightLeft1 = new double[contigCount];// the second largest score of all edges;
    double * edgeWeightRight1 = new double[contigCount];
    
    for(i=0;i<contigCount;i++){
        
        ScaffoldEdge * temp = scaffoldGraph[i].outLink;
        
        leftNumber[i]=0;
        rightNumber[i]=0;
        leftIndex[i]=-1;
        rightIndex[i]=-1;
        leftOrientation[i] = false;
        rightOrientation[i] = false; 
        
        edgeWeightLeft[i] = 0; //the largest score of all edges;
        edgeWeightRight[i] = 0;
        edgeWeightLeft1[i] = 0;// the second largest score of all edges;
        edgeWeightRight1[i] = 0;
        
        //Find the largest score
        while(temp!=NULL){           
            
            if(scaffoldGraph[temp->contigIndex].length<contigCutOff){
                temp = temp->next;
                continue;
            } 
                 
            if(temp->mapPosition->orientation==0){
                leftNumber[i]++;
                if(edgeWeightLeft[i] < temp->fitNumber){
                    leftOrientation[i] = temp->mapPosition->mateOrientation;
                    leftIndex[i] = temp->contigIndex;
                    edgeWeightLeft[i] = temp->fitNumber;
                }
            }else{
                rightNumber[i]++;
                if(edgeWeightRight[i] < temp->fitNumber){
                    rightOrientation[i] = temp->mapPosition->mateOrientation;
                    rightIndex[i] = temp->contigIndex;
                    edgeWeightRight[i] = temp->fitNumber;
                }            
            }                  
            temp = temp->next;
        }
        
        //Find the second largest score
        temp = scaffoldGraph[i].outLink;
        while(temp!=NULL){      
            if(temp->mapPosition->orientation==0){
                if(edgeWeightLeft1[i] < temp->fitNumber && temp->contigIndex != leftIndex[i]){
                    edgeWeightLeft1[i] = temp->fitNumber;
                }
            }else{
                if(edgeWeightRight1[i] < temp->fitNumber && temp->contigIndex != rightIndex[i]){
                    edgeWeightRight1[i] = temp->fitNumber;
                }            
            }                  
            temp = temp->next;
        }
        
    } 
    
    
    for(i=0;i<contigCount;i++){
        
        ScaffoldEdge * temp = scaffoldGraph[i].outLink;
                
        if(temp == NULL){
            continue;
        }
        
        if(rightNumber[i]>1 && scaffoldGraph[i].length >= contigCutOff){
            
            DeleteScaffoldEdge(scaffoldGraph,i,1,contigCutOff,rightIndex[i]);
        }

        if(leftNumber[i]>1 && scaffoldGraph[i].length >= contigCutOff){
            
            DeleteScaffoldEdge(scaffoldGraph,i,0,contigCutOff,leftIndex[i]); 
        }

        
    }   
    
    for(i=0;i<contigCount;i++){
        
        ScaffoldEdge * temp = scaffoldGraph[i].outLink;
        scaffoldGraph[i].inLink = NULL;
        ScaffoldEdge * previous = NULL;
                
        j = 0;
        
        if(temp==NULL){
            continue;
        }
        while(temp!=NULL){
                          
            if(temp->mapPosition->orientation == 0){
                if(previous==NULL){
                    scaffoldGraph[i].outLink = temp->next;
                    temp->next = scaffoldGraph[i].inLink;
                    scaffoldGraph[i].inLink = temp;
                    temp = scaffoldGraph[i].outLink;                   
                }else{
                    previous->next = temp->next;                
                    temp->next = scaffoldGraph[i].inLink;
                    scaffoldGraph[i].inLink = temp;
                    temp = previous->next;
                }                   
                continue;
            }
            
            previous = temp;
            temp = temp->next;
        }
        
    } 
    
    
    delete [] leftNumber;
    delete [] rightNumber;
    delete [] leftIndex;
    delete [] rightIndex;
    delete [] leftOrientation;
    delete [] rightOrientation;
    
    delete [] edgeWeightLeft;
    delete [] edgeWeightRight;
    delete [] edgeWeightLeft1;
    delete [] edgeWeightRight1;
    
    
    
}


ScaffoldEdge * SearchScaffoldGraphEdge(ScaffoldEdge * scaffoldEdge, long int contigIndex){
    while(scaffoldEdge!=NULL){
        if(scaffoldEdge->contigIndex == contigIndex ){
            return scaffoldEdge;
        }
        scaffoldEdge = scaffoldEdge->next;
    }
    return NULL;
}

long int DeleteScaffoldGraphEdge(ScaffoldGraph * scaffoldGraph, long int index, ScaffoldEdge * scaffoldEdge){
    
    long int mateIndex = scaffoldEdge->contigIndex;
    bool orientation = scaffoldEdge->mapPosition->orientation;
    bool mateOrientation = scaffoldEdge->mapPosition->mateOrientation;
    
    ScaffoldEdge * temp = NULL;
    ScaffoldEdge * temp1 = NULL;
    if(orientation == 1 && mateOrientation == 0){
        temp = scaffoldGraph[index].outLink;
        temp1 = scaffoldGraph[mateIndex].inLink;
    }
    if(orientation == 1 && mateOrientation == 1){
        temp = scaffoldGraph[index].outLink;
        temp1 = scaffoldGraph[mateIndex].outLink;
    }
    if(orientation == 0 && mateOrientation == 1){
        temp = scaffoldGraph[index].inLink;
        temp1 = scaffoldGraph[mateIndex].outLink;
    }
    if(orientation == 0 && mateOrientation == 0){
        temp = scaffoldGraph[index].inLink;
        temp1 = scaffoldGraph[mateIndex].inLink;
    }
    
    ScaffoldEdge * previous = NULL;
    while(temp!=NULL){
        if(temp->contigIndex == mateIndex){
            
            if(previous == NULL){
                if(orientation == 0){
                    scaffoldGraph[index].inLink = temp->next;
                }else{
                    scaffoldGraph[index].outLink = temp->next;
                }
            }else{
                previous->next = temp->next;
            }
            temp->next = NULL;
            break;
        }
        previous = temp;
        temp = temp->next;
    }
    previous = NULL;
    while(temp1!=NULL){
        if(temp1->contigIndex == index){
            
            if(previous == NULL){
                if(mateOrientation == 0){
                    scaffoldGraph[mateIndex].inLink = temp1->next;
                }else{
                    scaffoldGraph[mateIndex].outLink = temp1->next;
                }
            }else{
                previous->next = temp1->next;
            }
            temp1->next = NULL;
            break;
        }
        previous = temp1;
        temp1 = temp1->next;
    }
    
}


int BuildScaffoldGraphFromTwoBam(ScaffoldSet * scaffoldSet, long int * scaffoldLength, ContigSet * contigSet, long int * contigLength, ScaffoldToContig * scaffoldToContig, long int scaffoldCount, ScaffoldGraphHead * scaffoldGraphHead, InputArg * inputArg, char * result){
    
    const std::string bamFileName1 = inputArg->bamFileName1;
    const std::string bamFileName2 = inputArg->bamFileName2;
    long int insertsize = inputArg->insertsize;
    double std = inputArg->std; 
    long int contigLengthCutOff = insertsize + lambda*std;
    long int readLength = inputArg->readLength;
    
    scoreThreshold = inputArg->minEdgeWeight;
    minMappedPairedReadNumber = inputArg->minEdgeLinkNumber;
    repeativeCutOff = inputArg->minRepetitiveCov;
    
    double isPairedRead = inputArg->pairedRead;
    
    long int allContigLength = 0;
    long int noMappedReadCount = 0;
    long int readCount = 0;
    double noMappedReadRate = 0;
    long int mappedReadCount = 0;
        
    int pairedMappingIndex = 0;
    
    long int RefLength = 0;
    long int MateRefLength = 0;
    long int Position = 0;
    long int MatePosition = 0;
    bool Orientation = 0;
    bool MateOrientation = 0;
    long int RefID = -1;
    long int MateRefID = -1;
    
    long int allMappedReadNumber = 0;
    long int allReadNumber = 0;
    
    BamReader bamReader;
    BamReader bamReaderLeft;
    BamReader bamReaderRight;
    string bamFileNameRight = bamFileName1;
    string bamFileNameLeft = bamFileName2; 
    
    bamReader.Open(bamFileNameLeft);
    bamReaderLeft.Open(bamFileNameLeft);
    bamReaderRight.Open(bamFileNameRight);
    BamAlignment alignment;
    BamAlignment alignmentLeft;
    BamAlignment alignmentRight;
    
    long int contigCount = bamReaderLeft.GetReferenceCount();
    
    ScaffoldGraph * scaffoldGraph = new ScaffoldGraph[scaffoldCount];
    scaffoldGraphHead->scaffoldGraph = scaffoldGraph;
    
    const RefVector & refVector = bamReaderLeft.GetReferenceData();    
    
    long int count = 0;
    long int i = 0;
    long int j = 0;
    ReadMapPosition * readMapPosition = new ReadMapPosition[scaffoldCount];
    
    bool * printContigIndex = new bool[scaffoldCount];
    
    
    long int * peMapped = new long int[scaffoldCount];
    long int * noPeMapped = new long int[scaffoldCount];
    long int * noDistancePeMapped = new long int[scaffoldCount];
    double * mappedRate = new double[scaffoldCount];
    
    for(i=0;i<scaffoldCount;i++){
        RefLength = scaffoldLength[i]; 
        readMapPosition[i].index = new long int[RefLength];
        readMapPosition[i].leftIndex = new long int[RefLength];
        readMapPosition[i].rightIndex = new long int[RefLength];
        readMapPosition[i].noMapRightIndex = new long int[RefLength];
        readMapPosition[i].noMapLeftIndex = new long int[RefLength];
        readMapPosition[i].leftReadCoverage = new long int[RefLength];
        readMapPosition[i].rightReadCoverage = new long int[RefLength];
        for(j=0;j<RefLength;j++){
            readMapPosition[i].index[j] = 0;
            readMapPosition[i].leftIndex[j] = 0;
            readMapPosition[i].rightIndex[j] = 0;
            readMapPosition[i].noMapRightIndex[j] = 0;
            readMapPosition[i].noMapLeftIndex[j] = 0;
            readMapPosition[i].leftReadCoverage[j] = 0;
            readMapPosition[i].rightReadCoverage[j] = 0;
        }
        
        scaffoldGraph[i].length = RefLength;
        allContigLength = allContigLength + RefLength;
        printContigIndex[i] = false;
        peMapped[i] = 0;
        noPeMapped[i] = 0;
        noDistancePeMapped[i] = 0;
        mappedRate[i] = 0;
    }
    
    
    
    std::vector<int> clipSizes;
    std::vector<int> readPositions;
    std::vector<int> genomePositions;
    bool usePadded = false;
    
    
    i = 0;
    j = 0;
    
    
    long int possionLeftAllMappedCount = 0;
    long int possionRightAllMappedCount = 0;
    long int possionAllMappedCount = 0;
    
    while(bamReader.GetNextAlignmentCore(alignment)){
        
        if((alignment.AlignmentFlag & 0x900) != 0){
            continue;
        } 
        allReadNumber++;
    }
    bamReader.Close();
    
    
    long int aa = 0;
    long int bb = 0;
    long int cc = 0;
    long int dd = 0;
    PairedReadMappedData * allPairedReadMappedData = new PairedReadMappedData[allReadNumber];
    long int * leftReadMapNumber = new long int[allReadNumber];
    long int * rightReadMapNumber = new long int[allReadNumber];
    
    string leftReadName;
    string rightReadName;
    bool repetiveToken = false;
    
    bamReader.Open(bamFileNameLeft);
    
    i = -1;
    
    i = 0;
    while(bamReaderLeft.GetNextAlignmentCore(alignmentLeft) && bamReaderRight.GetNextAlignmentCore(alignmentRight)){
        
        while((alignmentLeft.AlignmentFlag & 0x900) != 0){
            bamReaderLeft.GetNextAlignmentCore(alignmentLeft);
            continue;
        } 
        while((alignmentRight.AlignmentFlag & 0x900) != 0){
            bamReaderRight.GetNextAlignmentCore(alignmentRight);
            continue;
        }
        
        
        if(alignmentLeft.IsMapped()){                                    
            allPairedReadMappedData[i].leftPosition = GetScaffoldofContigPosition(scaffoldToContig, alignmentLeft.RefID, contigLength, alignmentLeft.Position, readLength);
            allPairedReadMappedData[i].orientationLeft = GetScaffoldofContigOrientation(scaffoldToContig, alignmentLeft.RefID, alignmentLeft.IsReverseStrand());
            if(isPairedRead == true){
                allPairedReadMappedData[i].orientationLeft = !allPairedReadMappedData[i].orientationLeft;
            }
            allPairedReadMappedData[i].leftReference = scaffoldToContig[alignmentLeft.RefID].contigIndex;
            
            aa++;
        }
        if(alignmentRight.IsMapped()){
            allPairedReadMappedData[i].rightPosition = GetScaffoldofContigPosition(scaffoldToContig, alignmentRight.RefID, contigLength, alignmentRight.Position, readLength);
            allPairedReadMappedData[i].orientationRight = GetScaffoldofContigOrientation(scaffoldToContig, alignmentRight.RefID, alignmentRight.IsReverseStrand());
            if(isPairedRead == true){
                allPairedReadMappedData[i].orientationRight = !allPairedReadMappedData[i].orientationRight;
            }
            allPairedReadMappedData[i].rightReference = scaffoldToContig[alignmentRight.RefID].contigIndex;
            bb++;
        }
        
        if(allPairedReadMappedData[i].leftPosition>=0 && allPairedReadMappedData[i].rightPosition>=0){
            cc++;
            if(allPairedReadMappedData[i].leftReference != allPairedReadMappedData[i].rightReference){
                dd++;
            }
        }
        i++; 
    }
    
    i = 0;
    j = 0;
    long int falsePaired = 0;
    while(i<allReadNumber){
        long int token = 1;
        if(allPairedReadMappedData[i].leftPosition>=0 && allPairedReadMappedData[i].rightPosition>=0 && allPairedReadMappedData[i].leftReference == allPairedReadMappedData[i].rightReference){
            token = ComputeTruePairedMapping(allPairedReadMappedData + i, readLength, insertsize, std);
        } 
        if(token == 0){
            allPairedReadMappedData[i].leftPosition = -1;
            allPairedReadMappedData[i].rightPosition = -1;
            falsePaired++;
        }
        
        if(allPairedReadMappedData[i].leftPosition>=0){
            if(allPairedReadMappedData[i].leftPosition<insertsize+lambda*std){
                scaffoldGraph[allPairedReadMappedData[i].leftReference].leftReadCoverage++;
            }
            if(allPairedReadMappedData[i].leftPosition>scaffoldGraph[allPairedReadMappedData[i].leftReference].length - insertsize - lambda*std){
                scaffoldGraph[allPairedReadMappedData[i].leftReference].rightReadCoverage++;
            }
            for(long int pp = 0; pp<readLength && allPairedReadMappedData[i].leftPosition+pp<scaffoldGraph[allPairedReadMappedData[i].leftReference].length; pp++){
                readMapPosition[allPairedReadMappedData[i].leftReference].leftReadCoverage[allPairedReadMappedData[i].leftPosition+pp]++;
            }
            
        }
        if(allPairedReadMappedData[i].rightPosition>=0){
            if(allPairedReadMappedData[i].rightPosition<insertsize+lambda*std){
                scaffoldGraph[allPairedReadMappedData[i].rightReference].leftReadCoverage++;
            }
            if(allPairedReadMappedData[i].rightPosition>scaffoldGraph[allPairedReadMappedData[i].rightReference].length - insertsize - lambda*std){
                scaffoldGraph[allPairedReadMappedData[i].rightReference].rightReadCoverage++;
            }
            for(long int pp = 0; pp<readLength && allPairedReadMappedData[i].rightPosition+pp<scaffoldGraph[allPairedReadMappedData[i].rightReference].length; pp++){
                readMapPosition[allPairedReadMappedData[i].rightReference].rightReadCoverage[allPairedReadMappedData[i].rightPosition+pp]++;
            }
        }
        i++;
    }
    
    double pp = 0;
    for(i=0;i<scaffoldCount;i++){
        long int tempLength = scaffoldGraph[i].length - readLength;
        if(tempLength > insertsize+lambda*std - readLength ){
            tempLength = insertsize+lambda*std - readLength;
        }
        scaffoldGraph[i].leftReadCoverage = scaffoldGraph[i].leftReadCoverage/(double)(tempLength);
        scaffoldGraph[i].rightReadCoverage = scaffoldGraph[i].rightReadCoverage/(double)(tempLength);
        pp = pp + scaffoldGraph[i].leftReadCoverage + scaffoldGraph[i].rightReadCoverage;
    }
    allAverageReadCoverage = (double)(aa+bb-2*falsePaired)*readLength/(double)allContigLength;
    double stdOfReadCoverage = 0;
    //cout<<"tt"<<endl;
    for(i=0;i<allReadNumber;i++){
        if(allPairedReadMappedData[i].leftPosition>=0){
            stdOfReadCoverage = stdOfReadCoverage + pow((double)(readMapPosition[allPairedReadMappedData[i].leftReference].leftReadCoverage[allPairedReadMappedData[i].leftPosition] + readMapPosition[allPairedReadMappedData[i].leftReference].rightReadCoverage[allPairedReadMappedData[i].leftPosition])-allAverageReadCoverage,2);
        }
        if(allPairedReadMappedData[i].rightPosition>=0){
            stdOfReadCoverage = stdOfReadCoverage + pow((double)(readMapPosition[allPairedReadMappedData[i].rightReference].rightReadCoverage[allPairedReadMappedData[i].rightPosition] + readMapPosition[allPairedReadMappedData[i].rightReference].rightReadCoverage[allPairedReadMappedData[i].rightPosition])-allAverageReadCoverage,2);
        }
    }
    
    stdOfReadCoverage = sqrt(stdOfReadCoverage/allReadNumber);
    //cout<<"stdOfReadCoverage:"<<stdOfReadCoverage<<"--allAverageReadCoverage:"<<allAverageReadCoverage<<endl;
    
    long int repeatRemoveNumber = 0;
    for(i=0;i<allReadNumber;i++){
        if(allPairedReadMappedData[i].leftPosition>=0){
            double tempStd = ((double)(readMapPosition[allPairedReadMappedData[i].leftReference].leftReadCoverage[allPairedReadMappedData[i].leftPosition] + readMapPosition[allPairedReadMappedData[i].leftReference].rightReadCoverage[allPairedReadMappedData[i].leftPosition]) - allAverageReadCoverage)/stdOfReadCoverage;
            tempStd = fabs(tempStd);
            if(tempStd > 2){
                allPairedReadMappedData[i].leftPosition = -1;
                allPairedReadMappedData[i].rightPosition = -1;
                repeatRemoveNumber++;
                continue;
            }
        }
        if(allPairedReadMappedData[i].rightPosition>=0){
            double tempStd = ((double)(readMapPosition[allPairedReadMappedData[i].rightReference].rightReadCoverage[allPairedReadMappedData[i].rightPosition] + readMapPosition[allPairedReadMappedData[i].rightReference].rightReadCoverage[allPairedReadMappedData[i].rightPosition]) - allAverageReadCoverage)/stdOfReadCoverage;
            tempStd = fabs(tempStd);
            if(tempStd > 2){
                allPairedReadMappedData[i].leftPosition = -1;
                allPairedReadMappedData[i].rightPosition = -1;
                repeatRemoveNumber++;
                continue;
            }
        }
    }  

    i = 0;
    
    while(i<allReadNumber){
        if(allPairedReadMappedData[i].leftPosition>=0){           
            allMappedReadNumber++;
            RefID = allPairedReadMappedData[i].leftReference;
            Position = allPairedReadMappedData[i].leftPosition;
            Orientation = allPairedReadMappedData[i].orientationLeft;
            readMapPosition[RefID].index[Position]++; 
            possionAllMappedCount++;        
            if(Orientation==0){
                if(allPairedReadMappedData[i].rightPosition<0){
                    readMapPosition[RefID].leftIndex[Position]++;
                    readMapPosition[RefID].noMapLeftIndex[Position]++;
                }else{
                    MateRefID = allPairedReadMappedData[i].rightReference;
                    MateOrientation =allPairedReadMappedData[i].orientationRight;
                    if(RefID == MateRefID && MateOrientation==1){
                        readMapPosition[RefID].leftIndex[Position]++;
                    }
                    
                }
                
                possionLeftAllMappedCount++;
            }else{
                
                if(allPairedReadMappedData[i].rightPosition<0){
                    readMapPosition[RefID].rightIndex[Position]++;
                    readMapPosition[RefID].noMapRightIndex[Position]++;
                }else{
                    MateRefID = allPairedReadMappedData[i].rightReference;
                    MateOrientation = allPairedReadMappedData[i].orientationRight;
                    if(RefID == MateRefID && MateOrientation==0){
                        readMapPosition[RefID].rightIndex[Position]++;
                    }
                }
                
                possionRightAllMappedCount++;
            }   
        }
        
        if(allPairedReadMappedData[i].rightPosition>=0){                   
            allMappedReadNumber++;
            RefID = allPairedReadMappedData[i].rightReference;
            Position = allPairedReadMappedData[i].rightPosition;
            Orientation =allPairedReadMappedData[i].orientationRight;
            readMapPosition[RefID].index[Position]++; 
            
            possionAllMappedCount++;  
                    
            if(Orientation==0){
                if(allPairedReadMappedData[i].leftPosition<0){
                    readMapPosition[RefID].leftIndex[Position]++;
                    readMapPosition[RefID].noMapLeftIndex[Position]++;
                }else{
                    MateRefID = allPairedReadMappedData[i].leftReference;
                    MateOrientation = allPairedReadMappedData[i].orientationLeft;
                    if(RefID == MateRefID && MateOrientation==1){
                        readMapPosition[RefID].leftIndex[Position]++;
                    }
                }
                
                possionLeftAllMappedCount++;
            }else{
                
                if(allPairedReadMappedData[i].leftPosition<0){
                    readMapPosition[RefID].rightIndex[Position]++;
                    readMapPosition[RefID].noMapRightIndex[Position]++;
                }else{
                    MateRefID = allPairedReadMappedData[i].leftReference;
                    MateOrientation = allPairedReadMappedData[i].orientationLeft;
                    if(RefID == MateRefID && MateOrientation==0){
                        readMapPosition[RefID].rightIndex[Position]++;
                    }
                }
                
                possionRightAllMappedCount++;
            }   
        }    
        

        if(allPairedReadMappedData[i].rightPosition>=0 && allPairedReadMappedData[i].leftPosition>=0 && allPairedReadMappedData[i].rightReference != allPairedReadMappedData[i].leftReference){   
            Position = allPairedReadMappedData[i].rightPosition;
            MatePosition = allPairedReadMappedData[i].leftPosition;
            MateOrientation = allPairedReadMappedData[i].orientationLeft;
            Orientation = allPairedReadMappedData[i].orientationRight;
            MateRefID = allPairedReadMappedData[i].leftReference;
            RefID = allPairedReadMappedData[i].rightReference;
            
            MateRefLength = scaffoldLength[MateRefID];
            RefLength = scaffoldLength[RefID]; 
            
            bool tag = false;
            
            if(Orientation == 0 && MateOrientation == 1 && tag != true){
                if(MateRefLength - MatePosition + Position + readLength <= insertsize + lambda*std){
                    readMapPosition[RefID].leftIndex[Position]++;
                    readMapPosition[MateRefID].rightIndex[MatePosition]++;
                }else{
                    tag = true;
                }
            }
            
            if(Orientation == 1 && MateOrientation == 0 && tag != true){
                if(RefLength - Position + MatePosition + readLength <= insertsize + lambda*std){
                    readMapPosition[RefID].leftIndex[Position]++;
                    readMapPosition[MateRefID].rightIndex[MatePosition]++;
                }else{
                    tag = true;
                }       
            }
            
            if(Orientation == 0 && MateOrientation == 0 && tag != true){
                if(MatePosition + Position +2*readLength <= insertsize + lambda*std){
                    readMapPosition[RefID].leftIndex[Position]++;
                    readMapPosition[MateRefID].rightIndex[MatePosition]++;
                }else{
                    tag = true;
                }           
            }
            
            if(Orientation == 1 && MateOrientation == 1 && tag != true){
                if(MateRefLength + RefLength - MatePosition - Position <= insertsize + lambda*std){
                    readMapPosition[RefID].leftIndex[Position]++;
                    readMapPosition[MateRefID].rightIndex[MatePosition]++;
                }else{
                    tag = true;
                }
            }
            
            if(tag == false){
                AddScaffoldEdge(scaffoldGraph, RefID, Orientation, Position, MateRefID, MateOrientation, MatePosition);
                AddScaffoldEdge(scaffoldGraph, MateRefID, MateOrientation, MatePosition, RefID, Orientation, Position);
                count++;
            }     
        }
              
        i++;
        
    }

    long int posMatch = 0;
    long int posMatchCount = 0;
    
    long int posMatch1 = 0;
    long int posMatchCount1 = 0;
    
    double rightMapP = 0;//(double)posMatchCount/(double)posMatch;
    double leftMapP = 0;//(double)posMatchCount1/(double)posMatch1;


    noMappedReadRate = 1-(double)(aa+bb-2*falsePaired)/(double)(2*allReadNumber);
    
    //cout<<1-noMappedReadRate<<"--"<<(double)(aa+bb)/(double)(2*allReadNumber)<<endl;
    //exit(0);
    
    for(i=0;i<scaffoldCount;i++){
        
        if(scaffoldGraph[i].outLink!=NULL){
            ScaffoldEdge * temp = scaffoldGraph[i].outLink;   
                     
            while(temp!=NULL){
                
                if(temp->edgeWeight<minMappedPairedReadNumber){
                    temp->fitNumber = 0;
                    temp->fitDistribution = 0;
                    temp = temp->next;
                    continue;
                }
                
                long int * tempDistance = new long int[temp->edgeWeight + 1];
                long int * tempMateDistance = new long int[temp->edgeWeight + 1];
                long int * tempAllDistance = new long int[temp->edgeWeight + 1];
                int p = 0;
                for(p = 0; p<temp->edgeWeight+1; p++){
                    tempDistance[p] = 0;
                    tempMateDistance[p] = 0;
                    tempAllDistance[p] = 0;
                }
                
                MapPosition * tempMapPosition = temp->mapPosition;
                p = 0;

                temp->mapPosition = OptimizeMapPosition(temp->mapPosition, temp->edgeWeight, NULL);
                
                long int maxPosition = 0;
                long int minPosition = scaffoldLength[i];
                long int maxMatePosition = 0;
                long int minMatePosition = scaffoldLength[temp->contigIndex];
                

                tempMapPosition = temp->mapPosition;
                while(tempMapPosition!=NULL){
                    MateRefID = temp->contigIndex;
                    MateRefLength = scaffoldLength[MateRefID];
                    RefLength = scaffoldLength[i]; 
                    
                    if(tempMapPosition->orientation==0 && tempMapPosition->mateOrientation==1){
                        tempDistance[p] = tempMapPosition->pos + readLength;
                        tempMateDistance[p] = MateRefLength - tempMapPosition->matePos;
                        tempAllDistance[p] = tempDistance[p] + tempMateDistance[p];
                        p++;
                    }
                    if(tempMapPosition->orientation==1 && tempMapPosition->mateOrientation==0){
                        tempDistance[p] = RefLength - tempMapPosition->pos;
                        tempMateDistance[p] = tempMapPosition->matePos + readLength;
                        tempAllDistance[p] = tempDistance[p] + tempMateDistance[p];
                        p++;
                    }
                    if(tempMapPosition->orientation==0 && tempMapPosition->mateOrientation==0){
                        tempDistance[p] = tempMapPosition->pos + readLength;
                        tempMateDistance[p] = tempMapPosition->matePos + readLength;
                        tempAllDistance[p] = tempDistance[p] + tempMateDistance[p];
                        p++;
                    }
                    if(tempMapPosition->orientation==1 && tempMapPosition->mateOrientation==1){
                        tempDistance[p] = RefLength - tempMapPosition->pos;
                        tempMateDistance[p] = MateRefLength - tempMapPosition->matePos;
                        tempAllDistance[p] = tempDistance[p] + tempMateDistance[p];
                        p++;
                    } 
                    
                    tempMapPosition = tempMapPosition->next;
                                      
                }
                
                
                tempMapPosition = temp->mapPosition;
                long int tempGap = GetGapDistance(tempAllDistance, p, insertsize);
                double tempFD = 0;
                
                
                temp->fitDistribution = sqrt(temp->fitDistribution);
                long int gap = tempGap;
                if(gap <=0 ){
                    gap = 10;
                } 
                temp->gapDistance = gap;
                
                long int tempEdgeWeight = temp->edgeWeight;
                long int temp88 = tempEdgeWeight;
                
                double continousGap = 0;
                
                for(p = 0; p<temp->edgeWeight; p++){
                    if(gap+tempAllDistance[p]>insertsize+lambda*std || gap+tempAllDistance[p]<insertsize-lambda*std){
                        tempEdgeWeight--;
                    }
                }
                
                
                double mapP = 0;
                MateRefID = temp->contigIndex;
                MateRefLength = scaffoldLength[MateRefID];
                RefLength = scaffoldLength[i]; 
                double avgP = 0;

                if(temp->mapPosition->orientation==1){
                    
                    avgP = scaffoldGraph[i].rightReadCoverage/allAverageReadCoverage;
                    avgP = 1;
                    
                    MapFitDistribution(readMapPosition[i].rightIndex,readMapPosition[i].noMapRightIndex,tempMapPosition,tempAllDistance,noMappedReadRate,RefLength,MateRefLength,readLength,gap,insertsize, std, continousGap,temp->fitNumber,temp->fitDistribution,rightMapP,0,avgP);
                    
                }else{
                    
                    avgP = scaffoldGraph[i].leftReadCoverage/allAverageReadCoverage;
                    avgP = 1;
                    MapFitDistribution1(readMapPosition[i].leftIndex,readMapPosition[i].noMapLeftIndex,tempMapPosition,tempAllDistance,noMappedReadRate,RefLength,MateRefLength,readLength,gap,insertsize, std, continousGap, temp->fitNumber,temp->fitDistribution,leftMapP,0,avgP);  
                    
                }              
                
                
                temp = temp->next;

                
                delete [] tempDistance;
                delete [] tempMateDistance;
                delete [] tempAllDistance;
                
            }
            
        }
    }
    
    
    delete [] peMapped;
    delete [] noPeMapped;
    delete [] noDistancePeMapped;
    delete [] mappedRate;
    
}





#endif
