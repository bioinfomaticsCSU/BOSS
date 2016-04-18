#ifndef SCAFFOLDING_CPP_INCLUDED 
#define SCAFFOLDING_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <typeinfo>

#include "scaffolding.h"
#include "scaffoldgraph.h"
#include "lp/lp_lib.h"

using namespace std;
long int fileIndex = 0;

long int allContigLength = 0;

ScaffoldSet * InitScaffoldSet(ContigSet * contigSet, long int contigCount){
    ScaffoldSet * first = NULL;
    ScaffoldSet * previous = NULL;
    for(long int i = 0; i<contigCount; i++){
        ScaffoldSet * temp = new ScaffoldSet;
        ContigSequence * temp1 = new ContigSequence;
        temp1->index = i;
        temp1->gapDistance = 0;
        temp1->orientation = 0;
        temp->contigSequence = temp1;
        temp->length = strlen(contigSet[i].contig);
        if(previous!=NULL){
            previous->next = temp;
            previous = temp;
        }else{
            first = temp;
            previous = temp;
        }
    }
    return first;
}



ShortContig * ReverseShortContig(ShortContig * shortContig){
    
    long int i = 0;
    long int j = 0;
    
    ShortContig * temp = NULL;
    ShortContig * previous = NULL;
    long int gapDistance = 0;
    long int gapDistance1 = 0;
    while(shortContig!=NULL){
        
        if(shortContig->orientation==0){
            shortContig->orientation = 1;
        }else{
            shortContig->orientation = 0;
        }
        
        temp = shortContig->next;
        shortContig->next = previous;
        gapDistance1 = shortContig->gapDistance;
        shortContig->gapDistance = gapDistance;
        gapDistance = gapDistance1;
        previous = shortContig;       
        shortContig = temp;
                        
    }
    
    return previous;
}



MergeIndexHead * DetermineMergeBetweenTwoScaffold1(ShortContig * right, ShortContig * left1, int index){
    
    long int i = 0;
    long int j = 0;
    
    MergeIndexHead * head = new MergeIndexHead;
    MergeIndex * tempFirst = NULL;
    
    ShortContig * temp = NULL;
    ShortContig * temp1 = NULL;
    
    temp = right;
    temp1 = left1;
    
    long int leftNumber = 0;
    while(temp1!=NULL){
        leftNumber++;
        temp1 = temp1->next;
    }
    temp1 = left1;
    long int * maxNumber = new long int[leftNumber];
    long int * rightStart = new long int[leftNumber];
    for(i=0;i<leftNumber;i++){
        maxNumber[i] = 0;
        rightStart[i] = -1;
    }
    i = 0;
    
    ShortContig * tempReverse = new ShortContig;
    
    if(index == 1){
       
        while(temp!=NULL){
            ShortContig * tempShortContig = new ShortContig;
            tempShortContig->contigIndex = temp->contigIndex;
            
            tempShortContig->orientation = temp->orientation;
            tempShortContig->next = tempReverse->next;
            tempReverse->next = tempShortContig;
            temp = temp->next;
        }
        temp = tempReverse->next;
    } 
    
    if(index == 2){
        while(temp1!=NULL){
            ShortContig * tempShortContig = new ShortContig;
            tempShortContig->contigIndex = temp1->contigIndex;
            tempShortContig->orientation = temp1->orientation;
            tempShortContig->next = tempReverse->next;
            tempReverse->next = tempShortContig;
            temp1 = temp1->next;
        }
        temp1 = tempReverse->next;
    } 
    
    
    ShortContig * right1 = right;
    while(temp1!=NULL){
        i = 0;
        if(index == 0 || index == 2){
            temp = right;
        }else{
            temp = tempReverse->next;
        }
        while(temp!=NULL){
            if(temp1->contigIndex!=temp->contigIndex){
                i++;
                temp=temp->next;
            }else{      
                               
                ShortContig * temp2 = temp1;
                rightStart[j] = i;
                while(temp2!=NULL && temp!=NULL){
                    if(temp2->contigIndex!=temp->contigIndex){
                        break;
                    }else{
                        maxNumber[j]++;
                    }
                    temp2 = temp2->next;
                    temp = temp->next;                    
                }
                break;
            }
        }

        temp1 = temp1->next;  
        j++;  
         
    }
    
    long int max = 0;
    for(i = 0; i<leftNumber; i++){
        if(maxNumber[i]>max){
            j = i;
            max = maxNumber[i];
        }
    }
    
    
    head->count = max;
    head->rightStart = rightStart[j];
    head->leftStart = j;
    
    return head;   
    
}

BFSNode * InsertShortContigOfTail(ScaffoldGraph * scaffoldGraph, ContigSequence * previousBfsContigSequence, ContigSequence * bfsContigSequence, ScaffoldSet * currentScaffold, bool * visited, long int contigCutOff){
    
    long int i = 0;
    long int j = 0;
    long int longContigCount = 0;
    long int indexValue = -3;
    long int orientation = -1;
    long int preGapDistance = 0;
    long int gapDistance = 0;
    
    long int longGapDistance = 0;
    
    if(bfsContigSequence == NULL){
        return NULL;
    }
    
    ContigSequence * contigSequence = currentScaffold->contigSequence;
    while(contigSequence->next!=NULL){
        
        contigSequence = contigSequence->next;
    }
    
    ScaffoldEdge * bfsScaffoldEdge = NULL;   
    if(bfsContigSequence->orientation==1){
        bfsScaffoldEdge = scaffoldGraph[bfsContigSequence->index].inLink;
    }else{
        bfsScaffoldEdge = scaffoldGraph[bfsContigSequence->index].outLink;
    }
    
    ShortContig * shortContigHead = new ShortContig;   
    shortContigHead->next = NULL; 
    while(bfsScaffoldEdge!=NULL){
        
        if(visited[bfsScaffoldEdge->contigIndex] == false){
            
            ShortContig * temp = new ShortContig;
            temp->contigIndex = bfsScaffoldEdge->contigIndex;
            temp->gapDistance = bfsScaffoldEdge->gapDistance;
                
            if(bfsContigSequence->orientation==0){
                if(bfsScaffoldEdge->mapPosition->mateOrientation==1){
                    temp->orientation = 1;
                }else{
                    temp->orientation = 0;
                }
            }else{
                if(bfsScaffoldEdge->mapPosition->mateOrientation==0){
                    temp->orientation = 0;
                }else{
                    temp->orientation = 1;
                }
            }
            if(scaffoldGraph[bfsScaffoldEdge->contigIndex].length>=contigCutOff){
                longContigCount++;
                indexValue = bfsScaffoldEdge->contigIndex;
                orientation = temp->orientation;    
                longGapDistance =  temp->gapDistance;          
            }
                
            temp->next = shortContigHead->next;
            shortContigHead->next = temp;
            shortContigHead->gapDistance++;
            visited[bfsScaffoldEdge->contigIndex] = true;
            
                    
        }
        
        bfsScaffoldEdge = bfsScaffoldEdge->next;
    }

    if(longContigCount>1){
        return NULL;
    }
    
    if(longContigCount == 1){
        ScaffoldEdge * insertShortNode = NULL;
        if(orientation == 0){
            insertShortNode = scaffoldGraph[indexValue].inLink;
        }else{
            insertShortNode = scaffoldGraph[indexValue].outLink;
        }   
        while(insertShortNode != NULL){
            if(visited[insertShortNode->contigIndex]==false && scaffoldGraph[insertShortNode->contigIndex].length < contigCutOff){
                ShortContig * temp = new ShortContig;
                temp->contigIndex = insertShortNode->contigIndex;
                temp->gapDistance = longGapDistance - insertShortNode->gapDistance;
                
                if(orientation == 0){
                    if(insertShortNode->mapPosition->mateOrientation == 0){
                        temp->orientation = 1;
                    }else{
                        temp->orientation = 0;
                    }
                }else{
                    if(insertShortNode->mapPosition->mateOrientation == 0){
                        temp->orientation = 0;
                    }else{
                        temp->orientation = 1;
                    }
                }    
                
                visited[insertShortNode->contigIndex]=true;
                temp->next = shortContigHead->next;
                shortContigHead->next = temp;
                shortContigHead->gapDistance++;
            }
            insertShortNode = insertShortNode->next;
        }   
    }  
    
    ShortContig * temp = shortContigHead->next;
    ShortContig * temp1 = NULL;
    
    while(temp!=NULL){
        temp1 = temp->next;
        
        while(temp1!=NULL){
            if(temp->gapDistance > temp1->gapDistance){
                        
                long int tempGap = temp->gapDistance;
                temp->gapDistance = temp1->gapDistance;
                temp1->gapDistance = tempGap;
                        
                tempGap = temp->contigIndex;
                temp->contigIndex = temp1->contigIndex;
                temp1->contigIndex = tempGap;
                    
                tempGap = temp->orientation;
                temp->orientation = temp1->orientation;
                temp1->orientation = tempGap;
                        
            }
            temp1=temp1->next;
        }
        temp = temp->next;
    }
    
    temp = shortContigHead->next;
    temp1 = shortContigHead;
    if(longContigCount == 1){
        while(temp!=NULL){
            if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                preGapDistance = temp1->gapDistance;
                gapDistance = temp->gapDistance;
                temp1->next = NULL;
            }
            temp1 = temp;
            temp = temp->next;
        }
        
    }
    temp = shortContigHead->next;
    ContigSequence * firstContigSequence = NULL;
    
    while(temp!=NULL){ 
        ContigSequence * insertShortContigSequence = new ContigSequence;
        insertShortContigSequence->index = temp->contigIndex;
        insertShortContigSequence->orientation = temp->orientation;
        if(firstContigSequence!=NULL){
            contigSequence->gapDistance = abs(temp->gapDistance - firstContigSequence->gapDistance - scaffoldGraph[firstContigSequence->index].length);
        }else{
            if(temp->gapDistance < 0){
                previousBfsContigSequence->gapDistance = labs(previousBfsContigSequence->gapDistance + temp->gapDistance);
                previousBfsContigSequence->next = insertShortContigSequence;
                insertShortContigSequence->gapDistance = -temp->gapDistance;
                insertShortContigSequence->next = contigSequence;
                previousBfsContigSequence = insertShortContigSequence;
                temp=temp->next;
                continue;
            }
            contigSequence->gapDistance = temp->gapDistance;
        }
        
        contigSequence->next = insertShortContigSequence;
        contigSequence = insertShortContigSequence;
        firstContigSequence = insertShortContigSequence;
        temp = temp->next;
    }
    
    contigSequence = currentScaffold->contigSequence;
    while(contigSequence->next!=NULL){
        contigSequence = contigSequence->next;
    }
    
    
    if(longContigCount == 1){
        BFSNode * bfsNode = new BFSNode;
        bfsNode->orientation = orientation;
        bfsNode->contigIndex = indexValue;
        if(firstContigSequence != NULL){
            bfsNode->gapDistance = abs(gapDistance - preGapDistance - scaffoldGraph[contigSequence->index].length);
        }else{
            bfsNode->gapDistance = gapDistance;
        }
        
        return bfsNode;
    }
    
    return NULL;
    
    
}

BFSNode * InsertShortContigOfHead(ScaffoldGraph * scaffoldGraph, ContigSequence * previousBfsContigSequence, ContigSequence * bfsContigSequence, ScaffoldSet * currentScaffold, bool * visited, long int contigCutOff){
    
    long int i = 0;
    long int j = 0;
    long int longContigCount = 0;
    long int indexValue = -3;
    long int orientation = -1;
    long int preGapDistance = 0;
    long int gapDistance = 0;
    
    if(bfsContigSequence == NULL){
        return NULL;
    }
    
    ScaffoldEdge * bfsScaffoldEdge = NULL;   
    if(bfsContigSequence->orientation==0){
        bfsScaffoldEdge = scaffoldGraph[bfsContigSequence->index].inLink;
    }else{
        bfsScaffoldEdge = scaffoldGraph[bfsContigSequence->index].outLink;
    }
    
    ShortContig * shortContigHead = new ShortContig;
    
    long int longGapDistance = 0;
            
    while(bfsScaffoldEdge!=NULL){
            
        if(visited[bfsScaffoldEdge->contigIndex] == false){
            
            ShortContig * temp = new ShortContig;
            temp->contigIndex = bfsScaffoldEdge->contigIndex;
            temp->gapDistance = bfsScaffoldEdge->gapDistance;
                
            if(bfsContigSequence->orientation==0){
                if(bfsScaffoldEdge->mapPosition->mateOrientation==0){
                    temp->orientation = 1;
                }else{
                    temp->orientation = 0;
                }
            }else{
                if(bfsScaffoldEdge->mapPosition->mateOrientation==1){
                    temp->orientation = 0;
                }else{
                    temp->orientation = 1;
                }
            }
            
            if(scaffoldGraph[bfsScaffoldEdge->contigIndex].length>=contigCutOff){
                longContigCount++;
                indexValue = bfsScaffoldEdge->contigIndex;
                orientation = temp->orientation;
                longGapDistance = temp->gapDistance;
            }
                
            temp->next = shortContigHead->next;
            shortContigHead->next = temp;
            shortContigHead->gapDistance++;
            visited[bfsScaffoldEdge->contigIndex] = true;
                    
        }
                
        bfsScaffoldEdge = bfsScaffoldEdge->next;
    }
    if(longContigCount>1){
        return NULL;
    }
    
    if(longContigCount == 1){
        ScaffoldEdge * insertShortNode = NULL;
        if(orientation == 1){
            insertShortNode = scaffoldGraph[indexValue].inLink;
        }else{
            insertShortNode = scaffoldGraph[indexValue].outLink;
        }   
        while(insertShortNode != NULL){
            if(visited[insertShortNode->contigIndex]==false && scaffoldGraph[insertShortNode->contigIndex].length < contigCutOff){
                ShortContig * temp = new ShortContig;
                temp->contigIndex = insertShortNode->contigIndex;
                
                temp->gapDistance = longGapDistance - insertShortNode->gapDistance;
               
                if(orientation == 0){
                    if(insertShortNode->mapPosition->mateOrientation == 0){
                        temp->orientation = 0;
                       
                    }else{
                        temp->orientation = 1;
                        
                    }
                }else{
                    if(insertShortNode->mapPosition->mateOrientation == 0){
                        temp->orientation = 1;
                        
                    }else{
                        temp->orientation = 0;
                        
                    }
                }
                
                visited[insertShortNode->contigIndex]=true;
                temp->next = shortContigHead->next;
                shortContigHead->next = temp;
                shortContigHead->gapDistance++;
            }
            insertShortNode = insertShortNode->next;
        }   
    }
    
       
    ShortContig * temp = shortContigHead->next;
    ShortContig * temp1 = NULL;
    while(temp!=NULL){
        temp1 = temp->next;
        while(temp1!=NULL){
            if(temp->gapDistance > temp1->gapDistance){
                        
                long int tempGap = temp->gapDistance;
                temp->gapDistance = temp1->gapDistance;
                temp1->gapDistance = tempGap;
                        
                tempGap = temp->contigIndex;
                temp->contigIndex = temp1->contigIndex;
                temp1->contigIndex = tempGap;
                    
                tempGap = temp->orientation;
                temp->orientation = temp1->orientation;
                temp1->orientation = tempGap;
                        
            }
            temp1=temp1->next;
        }
        temp = temp->next;
    }
    
    temp = shortContigHead->next;
    temp1 = shortContigHead;
    if(longContigCount == 1){
        while(temp!=NULL){
            if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                preGapDistance = temp1->gapDistance;
                gapDistance = temp->gapDistance;
                temp1->next = NULL;
            }
            temp1 = temp;
            temp = temp->next;
        }
    }
        
    temp = shortContigHead->next;
    ContigSequence * firstContigSequence = NULL;
    ContigSequence * tempContigSequence = currentScaffold->contigSequence;
    while(temp!=NULL){ 
        ContigSequence * insertShortContigSequence = new ContigSequence;
        insertShortContigSequence->index = temp->contigIndex;
        insertShortContigSequence->orientation = temp->orientation;
        if(firstContigSequence!=NULL){
            insertShortContigSequence->gapDistance = abs(temp->gapDistance - firstContigSequence->gapDistance - scaffoldGraph[firstContigSequence->index].length);
        }else{
            if(temp->gapDistance < 0){
                insertShortContigSequence->gapDistance = labs(currentScaffold->contigSequence->gapDistance + temp->gapDistance);
                insertShortContigSequence->next = currentScaffold->contigSequence->next;
                currentScaffold->contigSequence->gapDistance = -temp->gapDistance;
                currentScaffold->contigSequence->next = insertShortContigSequence;
                temp=temp->next;
                continue;
            }
            insertShortContigSequence->gapDistance = temp->gapDistance;
        }
        currentScaffold->contigSequence = insertShortContigSequence;
        insertShortContigSequence->next = tempContigSequence;
        tempContigSequence = currentScaffold->contigSequence;
        firstContigSequence = currentScaffold->contigSequence;
        temp = temp->next;
    }
    
    if(longContigCount == 1){
        BFSNode * bfsNode = new BFSNode;
        bfsNode->orientation = orientation;
        bfsNode->contigIndex = indexValue;
        if(firstContigSequence != NULL){
            bfsNode->gapDistance = abs(gapDistance - preGapDistance - scaffoldGraph[currentScaffold->contigSequence->index].length);
        }else{
            bfsNode->gapDistance = gapDistance;
        }
        
        return bfsNode;
    }
    
    return NULL;
    
}

ContigSequence * FindAndMergeTwoScaffold(BFSNode * bfsNode, ScaffoldSet * currentScaffold, ScaffoldSet * scaffoldSet, bool token){

    long int i = 0;
    long int j = 0;
    
    long int first = -1;
    long int last = -1;
    
    long int startIndex = currentScaffold->contigSequence->index;
    
    ContigSequence * contigSequence = currentScaffold->contigSequence;
    while(contigSequence->next !=NULL){
        contigSequence = contigSequence->next;
    }
    
    long int endIndex = contigSequence->index;
    
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    ContigSequence * tempContigSequence = NULL;
    while(tempScaffoldSet != NULL){
        tempContigSequence = tempScaffoldSet->contigSequence;
        if(tempContigSequence == NULL){
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        
        
        if(tempContigSequence->index==bfsNode->contigIndex && tempContigSequence->index != startIndex && tempContigSequence->index != endIndex && tempContigSequence->next == NULL){
            first = 1;
            last = 1;
            
            break;
        }
        
        
        if(tempContigSequence->index==bfsNode->contigIndex && tempContigSequence->index != startIndex && tempContigSequence->index != endIndex){
            first = 1;
            break;
        }
        while(tempContigSequence->next!=NULL){
            tempContigSequence = tempContigSequence->next;
        }
        if(tempContigSequence->index==bfsNode->contigIndex && tempContigSequence->index != startIndex && tempContigSequence->index != endIndex){
            last = 1;
            break;
        }
        tempScaffoldSet = tempScaffoldSet->next;
    }
    
    if(tempScaffoldSet == NULL){
        return NULL;
    }
    
    if(last == 1 && first == last){
        if(token == 1 && bfsNode->orientation == 0){
            last = 0;
        }
        if(token == 1 && bfsNode->orientation == 1){
            first = 0;
        }
        if(token == 0 && bfsNode->orientation == 0){
            first = 0;
        }
        if(token == 0 && bfsNode->orientation == 1){
            last = 0;
        }
        
    }

    if(token == 1){
        
        if(first == 1 && bfsNode->orientation == tempContigSequence->orientation){
            contigSequence->gapDistance = bfsNode->gapDistance;
            contigSequence->next = tempContigSequence;
            tempScaffoldSet->contigSequence = NULL;
        }
        if(last == 1 && bfsNode->orientation != tempContigSequence->orientation){
            tempScaffoldSet->contigSequence = ReverseContigSequence(tempScaffoldSet->contigSequence);
            contigSequence->gapDistance = bfsNode->gapDistance;
            contigSequence->next = tempScaffoldSet->contigSequence;
            tempScaffoldSet->contigSequence = NULL;
        }
        while(contigSequence->next!=NULL){
            contigSequence = contigSequence->next;
        }
        return contigSequence;
         
    }
    if(token == 0){
        if(last == 1){
            tempContigSequence->gapDistance = bfsNode->gapDistance;
            tempContigSequence->next = currentScaffold->contigSequence;
            currentScaffold->contigSequence = tempScaffoldSet->contigSequence;
            tempScaffoldSet->contigSequence = NULL;
            
        }
        
        if(first == 1){
            tempScaffoldSet->contigSequence = ReverseContigSequence(tempScaffoldSet->contigSequence);
            tempContigSequence = tempScaffoldSet->contigSequence;
            while(tempContigSequence->next!=NULL){
                tempContigSequence = tempContigSequence->next;
            }
            tempContigSequence->gapDistance = bfsNode->gapDistance;
            tempContigSequence->next = currentScaffold->contigSequence;
            currentScaffold->contigSequence = tempScaffoldSet->contigSequence;
            tempScaffoldSet->contigSequence = NULL;
            
        }
        
        return currentScaffold->contigSequence;
        
    }

}

void BFSScaffold(ScaffoldSet * scaffoldSet, ScaffoldSet * currentScaffold, ScaffoldGraph * scaffoldGraph, long int contigCount, long int contigCutOff, bool * visited){
    
    long int i = 0;
    long int j = 0;
    
    long int firstContigIndex = 0;
    long int firstOrientation = 0;
    long int lastContigIndex = 0;
    long int lastOrientation = 0;

    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    ContigSequence * tempContigSequence = currentScaffold->contigSequence;
    
    firstContigIndex = tempContigSequence->index;
    firstOrientation = tempContigSequence->orientation;
    while(tempContigSequence->next!=NULL){
        tempContigSequence = tempContigSequence->next;
    }
    lastContigIndex = tempContigSequence->index;
    lastOrientation = tempContigSequence->orientation;
    
    ContigSequence * bfsContigSequence = tempContigSequence;
    ContigSequence * previousBfsContigSequence = NULL;
    long int indexValue = -1;
    
    while(bfsContigSequence != NULL){        
        BFSNode * bfsNode = InsertShortContigOfTail(scaffoldGraph, previousBfsContigSequence, bfsContigSequence, currentScaffold, visited, contigCutOff);
       
        if(bfsNode != NULL){
            bfsContigSequence = FindAndMergeTwoScaffold(bfsNode, currentScaffold, scaffoldSet, 1);
            
            continue;
        }
        previousBfsContigSequence = bfsContigSequence;
        bfsContigSequence = bfsContigSequence->next;
        
    }
    
    bfsContigSequence = currentScaffold->contigSequence;
    previousBfsContigSequence = NULL;
    while(bfsContigSequence != NULL){        
        BFSNode * bfsNode = InsertShortContigOfHead(scaffoldGraph, previousBfsContigSequence, bfsContigSequence, currentScaffold, visited, contigCutOff);
        
        if(bfsNode != NULL){
            
            bfsContigSequence = FindAndMergeTwoScaffold(bfsNode, currentScaffold, scaffoldSet, 0);
            tempContigSequence = bfsContigSequence;
            
            continue;
        }
        long int index0 = bfsContigSequence->index;
        tempContigSequence = currentScaffold->contigSequence;
        while(tempContigSequence->next!=NULL){
            if(tempContigSequence->next->index == bfsContigSequence->index){
                previousBfsContigSequence = bfsContigSequence;
                bfsContigSequence = tempContigSequence;
                break;
            }
            tempContigSequence = tempContigSequence->next;
        }
        if(index0 == bfsContigSequence->index){
            break;
        }
        
    }
    
    
}

int BFSScaffolding(ScaffoldGraph * scaffoldGraph, ScaffoldSet * scaffoldSet, long int contigCount, long int contigCutOff){
    
    long int i = 0;
    long int j = 0;
    
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    ContigSequence * tempContigSequence = NULL;
    
    bool * visited = new bool[contigCount];
    for(i = 0; i<contigCount; i++){
        visited[i] = false;
    }
    
    tempScaffoldSet = scaffoldSet;
    while(tempScaffoldSet!=NULL){
        tempContigSequence = tempScaffoldSet->contigSequence;
        while(tempContigSequence != NULL){
            if(scaffoldGraph[tempContigSequence->index].length < contigCutOff){
                visited[tempContigSequence->index] = true;
            }
            tempContigSequence = tempContigSequence->next;
        }
        tempScaffoldSet = tempScaffoldSet->next;
    }
   
    i = 0;
    tempScaffoldSet = scaffoldSet;
    while(tempScaffoldSet!=NULL){
        if(tempScaffoldSet->contigSequence==NULL){
            i++;
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        
        BFSScaffold(scaffoldSet, tempScaffoldSet, scaffoldGraph, contigCount, contigCutOff, visited);
        tempScaffoldSet = tempScaffoldSet->next;
        i++;
    }
    
    
    for(i = 0; i<contigCount; i++){
        visited[i] = false;
    }
    
    tempScaffoldSet = scaffoldSet;
    ScaffoldSet * pre = NULL;
    while(tempScaffoldSet!=NULL){
        tempContigSequence = tempScaffoldSet->contigSequence;
        while(tempContigSequence != NULL){
            visited[tempContigSequence->index] = true;
            tempContigSequence = tempContigSequence->next;
        }
        pre = tempScaffoldSet;
        tempScaffoldSet = tempScaffoldSet->next;
    }
    
    bool * visited1 = new bool[contigCount];
    for(i = 0; i<contigCount; i++){
        visited1[i] = false;
    }
    AddShortContigToScaffoldSet(scaffoldSet, contigCount, visited1);
    
    tempScaffoldSet = pre->next;
    i = 0;
    while(tempScaffoldSet!=NULL){
        if(tempScaffoldSet->contigSequence==NULL){
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        if(visited[tempScaffoldSet->contigSequence->index] == true){
            tempScaffoldSet->contigSequence = NULL;
            i++;
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        
        visited[tempScaffoldSet->contigSequence->index] = true;
        BFSScaffold(scaffoldSet, tempScaffoldSet, scaffoldGraph, contigCount, contigCutOff, visited);
        
        tempScaffoldSet = tempScaffoldSet->next;
        i++;
    }
    
}


int MergeScaffoldSet(ScaffoldGraph * scaffoldGraph, ScaffoldSet * scaffoldSet, long int contigCutOff){
    long int i = 0;
    long int j = 0;
    long int p = 0;
    
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    ScaffoldSet * tempScaffoldSet1 = scaffoldSet;
    ContigSequence * tempContigSequence = NULL;
    ContigSequence * firstContigSequence = NULL;
    ScaffoldEdge * tempScaffoldEdge = NULL;
    
    long int scaffoldCount = 0;
    
    while(tempScaffoldSet1!=NULL){
        scaffoldCount++;
        tempScaffoldSet1 = tempScaffoldSet1->next;
    }
    tempScaffoldSet1 = tempScaffoldSet;
    
    ShortContig * outShortContigHead = new ShortContig[scaffoldCount];
    ShortContig * inShortContigHead = new ShortContig[scaffoldCount];
    
    while(tempScaffoldSet1!=NULL){
        
        tempContigSequence = tempScaffoldSet1->contigSequence;
        j = 0;
        ShortContig * firstShortContig = NULL;
        ShortContig * lastShortContig = NULL;
        while(tempContigSequence!=NULL){
            
            if(j==0&&scaffoldGraph[tempContigSequence->index].length<contigCutOff){
            
                ShortContig * tempShortContig = new ShortContig;
                tempShortContig->contigIndex = tempContigSequence->index;
                tempShortContig->gapDistance = tempContigSequence->gapDistance;
                tempShortContig->orientation = tempContigSequence->orientation;
                
                if(lastShortContig == NULL){
                    inShortContigHead[i].next = tempShortContig;
                    lastShortContig = tempShortContig;
                }else{
                    lastShortContig->next = tempShortContig;   
                    lastShortContig = tempShortContig;
                }
                inShortContigHead[i].gapDistance++;
                         
            }else{
                j++;
            }
            if(j>0&&scaffoldGraph[tempContigSequence->index].length<contigCutOff){
                
                ContigSequence * temp22 = tempContigSequence;
                while(temp22!=NULL){
                    if(scaffoldGraph[temp22->index].length>=contigCutOff){
                        break;
                    }
                    temp22 = temp22->next;
                }
                if(temp22!=NULL){
                    tempContigSequence = tempContigSequence->next;
                    continue;
                }
                
                ShortContig * tempShortContig = new ShortContig;
                tempShortContig->contigIndex = tempContigSequence->index;
                tempShortContig->gapDistance = tempContigSequence->gapDistance;
                tempShortContig->orientation = tempContigSequence->orientation;
                if(firstShortContig == NULL){
                    outShortContigHead[i].next = tempShortContig;
                    firstShortContig = tempShortContig;
                }else{
                    firstShortContig->next = tempShortContig;   
                    firstShortContig = tempShortContig;
                } 
                outShortContigHead[i].gapDistance++;           
            }
            
            
            tempContigSequence = tempContigSequence->next;
        }        
        
        i++;
        tempScaffoldSet1 = tempScaffoldSet1->next;
    }
    
    tempScaffoldSet1 = tempScaffoldSet;
    i = 0;
    j = 0;
    while(tempScaffoldSet!=NULL){
        j = 0;
        
        tempScaffoldSet1 = scaffoldSet;
        
        if(tempScaffoldSet->contigSequence==NULL){
            i++;
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        
        while(tempScaffoldSet1!=NULL){
            
            if(i == j || tempScaffoldSet1->contigSequence==NULL){
                tempScaffoldSet1 = tempScaffoldSet1->next;
                j++;
                continue;
            }
            
            long int token = 0;
            
            long int t = 0;
            long int t1 = 0;
            long int t2 = 0;
            long int t3 = 0;
            
            MergeIndexHead * tempHead = DetermineMergeBetweenTwoScaffold1(outShortContigHead[i].next, inShortContigHead[j].next, 0);
            MergeIndexHead * tempHead1 = DetermineMergeBetweenTwoScaffold1(inShortContigHead[i].next, inShortContigHead[j].next, 1);
            MergeIndexHead * tempHead2 = DetermineMergeBetweenTwoScaffold1(outShortContigHead[i].next, outShortContigHead[j].next, 2);
            MergeIndexHead * tempHead3 = DetermineMergeBetweenTwoScaffold1(outShortContigHead[j].next, inShortContigHead[i].next, 0);
            t = tempHead->count;
            t1 = tempHead1->count;
            t2 = tempHead2->count;
            t3 = tempHead3->count;

            if((t2>0||t1>0)&&(t>0||t3>0)){
                tempScaffoldSet1 = tempScaffoldSet1->next;
                j++;
                continue;
            }
            
            long int max = 0;
            long int max1 = 0;
            
            if(t2>t1){
                max = t2;
            }else{
                max = t1;
            }
                        
            if(t>t3){
                max1 = t;
            }else{
                max1 = t3;
            }

            t = tempHead->count;
            
            if(t>0 && t == max1){
                
                
                p = 1;
                while(p<tempHead->leftStart + t){
                    tempScaffoldSet1->contigSequence = tempScaffoldSet1->contigSequence->next;
                    p++;
                }
                tempContigSequence = tempScaffoldSet->contigSequence;
                ContigSequence * pp = tempContigSequence;
                long int rightNumber = 0;
                while(pp!=NULL){
                    pp = pp->next;
                    rightNumber++;
                }
                p = 1;
                while(p < rightNumber - outShortContigHead[i].gapDistance + tempHead->rightStart + t){
                    tempContigSequence = tempContigSequence->next;
                    p++;
                }
                
                tempContigSequence->gapDistance = tempScaffoldSet1->contigSequence->gapDistance;
                tempContigSequence->next = tempScaffoldSet1->contigSequence->next;
                
                if(tempHead->orientation == 1){
                    
                    while(tempScaffoldSet1->contigSequence!=NULL){
                        if(tempScaffoldSet1->contigSequence->orientation == 0){
                            tempScaffoldSet1->contigSequence->orientation = 1;
                        }else{
                            tempScaffoldSet1->contigSequence->orientation = 0;
                        }                        
                        tempScaffoldSet1->contigSequence = tempScaffoldSet1->contigSequence->next;
                    }                    
                }
                token = 1;
                tempScaffoldSet1->contigSequence = NULL;
                outShortContigHead[i].next = outShortContigHead[j].next;
                outShortContigHead[i].gapDistance = outShortContigHead[j].gapDistance;
            }
            
            if(token == 1){
                tempScaffoldSet1 = tempScaffoldSet1->next;
                j++;
                continue;
            }
            
            
            
            t = tempHead1->count;
            
            p = 0;
            
            if(t>0 && t == max1){
                
                tempScaffoldSet->contigSequence = ReverseContigSequence(tempScaffoldSet->contigSequence);
                
                ContigSequence * temp22 = tempScaffoldSet1->contigSequence;
                while(temp22!=NULL){
                    
                    temp22 = temp22->next;
                }
                
                p = 1;
                while(p<tempHead1->leftStart + t){
                    tempScaffoldSet1->contigSequence = tempScaffoldSet1->contigSequence->next;
                    p++;
                }
                
                tempContigSequence = tempScaffoldSet->contigSequence;
                ContigSequence * pp = tempContigSequence;
                long int rightNumber = 0;
                while(pp!=NULL){
                    pp = pp->next;
                    rightNumber++;
                }
                p = 1;
                
                while(p < rightNumber - inShortContigHead[i].gapDistance + tempHead1->rightStart + t){
                    tempContigSequence = tempContigSequence->next;
                    p++;
                }
                
                tempContigSequence->gapDistance = tempScaffoldSet1->contigSequence->gapDistance;
                tempContigSequence->next = tempScaffoldSet1->contigSequence->next;
                
                
                if(tempHead->orientation == 1){
                    
                    while(tempScaffoldSet1->contigSequence!=NULL){
                        if(tempScaffoldSet1->contigSequence->orientation == 0){
                            tempScaffoldSet1->contigSequence->orientation = 1;
                        }else{
                            tempScaffoldSet1->contigSequence->orientation = 0;
                        }                        
                        tempScaffoldSet1->contigSequence = tempScaffoldSet1->contigSequence->next;
                    }                    
                }
                token = 1;
                tempScaffoldSet1->contigSequence = NULL;
                
                inShortContigHead[i].next = ReverseShortContig(outShortContigHead[i].next);
                inShortContigHead[i].gapDistance = outShortContigHead[i].gapDistance;
                outShortContigHead[i].next = NULL;
                
                outShortContigHead[i].next = outShortContigHead[j].next;
                outShortContigHead[i].gapDistance = outShortContigHead[j].gapDistance;
                outShortContigHead[j].next = NULL;
                
                 
            }
            
            if(token == 1){
                tempScaffoldSet1 = tempScaffoldSet1->next;
                j++;
                continue;
            }
            
            
            
            t = tempHead2->count;
            p = 0;
            
            if(t>0 && t == max1){
                
                tempScaffoldSet1->contigSequence = ReverseContigSequence(tempScaffoldSet1->contigSequence);
                
                p = 1;
                while(p<tempHead2->leftStart + t){
                    tempScaffoldSet1->contigSequence = tempScaffoldSet1->contigSequence->next;
                    p++;
                }
                tempContigSequence = tempScaffoldSet->contigSequence;
                ContigSequence * pp = tempContigSequence;
                long int rightNumber = 0;
                while(pp!=NULL){
                    pp = pp->next;
                    rightNumber++;
                }
                p = 1;
                while(p < rightNumber - outShortContigHead[i].gapDistance + tempHead2->rightStart + t){
                    tempContigSequence = tempContigSequence->next;
                    p++;
                }
                
                tempContigSequence->gapDistance = tempScaffoldSet1->contigSequence->gapDistance;
                tempContigSequence->next = tempScaffoldSet1->contigSequence->next;
                
                if(tempHead->orientation == 1){
                    
                    while(tempScaffoldSet1->contigSequence!=NULL){
                        if(tempScaffoldSet1->contigSequence->orientation == 0){
                            tempScaffoldSet1->contigSequence->orientation = 1;
                        }else{
                            tempScaffoldSet1->contigSequence->orientation = 0;
                        }                        
                        tempScaffoldSet1->contigSequence = tempScaffoldSet1->contigSequence->next;
                    }                    
                }
                token = 1;
                tempScaffoldSet1->contigSequence = NULL;
                outShortContigHead[i].next = ReverseShortContig(inShortContigHead[j].next);
                outShortContigHead[i].gapDistance = inShortContigHead[j].gapDistance;
                inShortContigHead[j].next = NULL;
            }
            
            tempScaffoldSet1 = tempScaffoldSet1->next;
            j++;
            
            
        }
        i++;
        tempScaffoldSet = tempScaffoldSet->next;
    }
    
    
}

int SearchScaffoldEdge(long int index, ContigSequence * contigSequence){
    
    long int i = 0;
    
    while(contigSequence!=NULL){
        if(contigSequence->index==index){
            return 0;
        }
        contigSequence = contigSequence->next;
    }
    
    return 1;
    
}

int InsertShortContigBetweenScaffoldSet(ScaffoldGraph * scaffoldGraph, ScaffoldSet * scaffoldSet, long int contigCount, long int contigCutOff){
    
    long int i = 0;
    long int j = 0;
    
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    ContigSequence * tempContigSequence = NULL;
    ContigSequence * firstContigSequence = NULL;
    ScaffoldEdge * tempScaffoldEdge = NULL;
    
    while(tempScaffoldSet != NULL){
        
        tempContigSequence = tempScaffoldSet->contigSequence;
        
        if(tempContigSequence==NULL){
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        
        if(tempContigSequence->orientation==0){
            tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].inLink;
        }else{
            tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].outLink;
        }
        
        ShortContig * shortContigHead = new ShortContig;
            
        while(tempScaffoldEdge!=NULL){
            
            int cc = SearchScaffoldEdge(tempScaffoldEdge->contigIndex, tempScaffoldSet->contigSequence);
            
            if(scaffoldGraph[tempScaffoldEdge->contigIndex].length<contigCutOff && cc!=0){
                
                ShortContig * temp = new ShortContig;
                temp->contigIndex = tempScaffoldEdge->contigIndex;
                temp->gapDistance = tempScaffoldEdge->gapDistance;
                
                if(tempContigSequence->orientation==0){
                    if(tempScaffoldEdge->mapPosition->mateOrientation==0){
                        temp->orientation = 1;
                    }else{
                        temp->orientation = 0;
                    }
                }else{
                    if(tempScaffoldEdge->mapPosition->mateOrientation==1){
                        temp->orientation = 0;
                    }else{
                        temp->orientation = 1;
                    }
                }
                
                temp->next = shortContigHead->next;
                shortContigHead->next = temp;
                shortContigHead->gapDistance++;
                    
            }
                
            tempScaffoldEdge = tempScaffoldEdge->next;
        }
        
        ShortContig * temp = shortContigHead->next;
        ShortContig * temp1 = NULL;
        while(temp!=NULL){
            temp1 = temp->next;
            while(temp1!=NULL){
                if(temp->gapDistance > temp1->gapDistance){
                        
                    long int tempGap = temp->gapDistance;
                    temp->gapDistance = temp1->gapDistance;
                    temp1->gapDistance = tempGap;
                        
                    tempGap = temp->contigIndex;
                    temp->contigIndex = temp1->contigIndex;
                    temp1->contigIndex = tempGap;
                    
                    tempGap = temp->orientation;
                    temp->orientation = temp1->orientation;
                    temp1->orientation = tempGap;
                        
                }
                temp1=temp1->next;
            }
            temp = temp->next;
        }
        
        temp = shortContigHead->next;
        firstContigSequence = NULL;
        while(temp!=NULL){
                
            ContigSequence * insertShortContigSequence = new ContigSequence;
            insertShortContigSequence->index = temp->contigIndex;
            insertShortContigSequence->orientation = temp->orientation;
            if(firstContigSequence!=NULL){
                insertShortContigSequence->gapDistance = abs(temp->gapDistance - firstContigSequence->gapDistance - scaffoldGraph[firstContigSequence->index].length);
            }else{
                insertShortContigSequence->gapDistance = temp->gapDistance;
            }
            tempScaffoldSet->contigSequence = insertShortContigSequence;
            insertShortContigSequence->next = tempContigSequence;
            tempContigSequence = tempScaffoldSet->contigSequence;
            firstContigSequence = tempScaffoldSet->contigSequence;
            temp = temp->next;
            
        }
        
        
        while(tempContigSequence!=NULL && tempContigSequence->next!=NULL){
            tempContigSequence = tempContigSequence->next;
        }
        
        if(tempContigSequence==NULL){
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        
        
        if(tempContigSequence->orientation==1){
            tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].inLink;
        }else{
            tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].outLink;
        }
        
        shortContigHead = new ShortContig;
            
        while(tempScaffoldEdge!=NULL){
            
            int cc = SearchScaffoldEdge(tempScaffoldEdge->contigIndex, tempScaffoldSet->contigSequence);
            if(scaffoldGraph[tempScaffoldEdge->contigIndex].length<contigCutOff && cc!=0){
                ShortContig * temp = new ShortContig;
                temp->contigIndex = tempScaffoldEdge->contigIndex;
                temp->gapDistance = tempScaffoldEdge->gapDistance;
                
                if(tempContigSequence->orientation==0){
                    if(tempScaffoldEdge->mapPosition->mateOrientation==1){
                        temp->orientation = 1;
                    }else{
                        temp->orientation = 0;
                    }
                }else{
                    if(tempScaffoldEdge->mapPosition->mateOrientation==0){
                        temp->orientation = 0;
                    }else{
                        temp->orientation = 1;
                    }
                }
                
                temp->next = shortContigHead->next;
                shortContigHead->next = temp;
                shortContigHead->gapDistance++;
                    
            }
                
            tempScaffoldEdge = tempScaffoldEdge->next;
        }
        
        temp = shortContigHead->next;
        temp1 = NULL;
        while(temp!=NULL){
            temp1 = temp->next;
            while(temp1!=NULL){
                if(temp->gapDistance > temp1->gapDistance){
                        
                    long int tempGap = temp->gapDistance;
                    temp->gapDistance = temp1->gapDistance;
                    temp1->gapDistance = tempGap;
                        
                    tempGap = temp->contigIndex;
                    temp->contigIndex = temp1->contigIndex;
                    temp1->contigIndex = tempGap;
                    
                    tempGap = temp->orientation;
                    temp->orientation = temp1->orientation;
                    temp1->orientation = tempGap;
                        
                }
                temp1=temp1->next;
            }
            temp = temp->next;
        }
        
        temp = shortContigHead->next;
        firstContigSequence = NULL;
        while(temp!=NULL){
        
            ContigSequence * insertShortContigSequence = new ContigSequence;
            insertShortContigSequence->index = temp->contigIndex;
            insertShortContigSequence->orientation = temp->orientation;
            if(firstContigSequence!=NULL){
                tempContigSequence->gapDistance = abs(temp->gapDistance - firstContigSequence->gapDistance - scaffoldGraph[firstContigSequence->index].length);
            }else{
                tempContigSequence->gapDistance = temp->gapDistance;
            }
            tempContigSequence->next = insertShortContigSequence;
            tempContigSequence = insertShortContigSequence;
            firstContigSequence = insertShortContigSequence;
            temp = temp->next;
            
        }        
        
        tempScaffoldSet = tempScaffoldSet->next;       
    }

} 

int InsertShortContigInScaffoldSet(ScaffoldGraph * scaffoldGraph, ScaffoldSet * scaffoldSet, long int contigCount, long int contigCutOff){
    
    long int i = 0;
    long int j = 0;
    
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    ContigSequence * tempContigSequence = NULL;
    ContigSequence * firstContigSequence = NULL;
    ScaffoldEdge * tempScaffoldEdge = NULL;
    
    bool * visited = new bool[contigCount];
    for(i=0;i<contigCount;i++){
        visited[i] = false;
    }
    
    
    while(tempScaffoldSet != NULL){
        
        tempContigSequence = tempScaffoldSet->contigSequence;
        
        while(tempContigSequence!=NULL && tempContigSequence->next!=NULL){
            
            firstContigSequence = tempContigSequence;
            
            if(tempContigSequence->orientation==0){
                tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].outLink;
            }else{
                tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].inLink;
            }
            
            ShortContig * shortContigHead = new ShortContig;
            
            while(tempScaffoldEdge!=NULL){
                
                if(scaffoldGraph[tempScaffoldEdge->contigIndex].length<contigCutOff && visited[tempScaffoldEdge->contigIndex] !=true){
                   
                    ShortContig * temp = new ShortContig;
                    temp->contigIndex = tempScaffoldEdge->contigIndex;
                    temp->gapDistance = tempScaffoldEdge->gapDistance;
                    
                    if(tempContigSequence->orientation==0){
                        if(tempScaffoldEdge->mapPosition->mateOrientation==1){
                            temp->orientation = 1;
                        }else{
                            temp->orientation = 0;
                        }
                    }else{
                        if(tempScaffoldEdge->mapPosition->mateOrientation==0){
                             temp->orientation = 0;
                        }else{
                             temp->orientation = 1;
                        }
                    }
                    
                    temp->next = shortContigHead->next;
                    shortContigHead->next = temp;
                    shortContigHead->gapDistance++;
                    
                }
                
                tempScaffoldEdge = tempScaffoldEdge->next;
            }
            
            long int longGapDistance = tempContigSequence->gapDistance;
            tempContigSequence = tempContigSequence->next;
            if(tempContigSequence->orientation==1){
                tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].outLink;
            }else{
                tempScaffoldEdge = scaffoldGraph[tempContigSequence->index].inLink;
            }
            
            ShortContig * shortContigHead1 = new ShortContig;
            
            while(tempScaffoldEdge!=NULL){
                
                if(scaffoldGraph[tempScaffoldEdge->contigIndex].length<contigCutOff && visited[tempScaffoldEdge->contigIndex] !=true){
                    
                    ShortContig * temp = new ShortContig;
                    temp->contigIndex = tempScaffoldEdge->contigIndex;
                    temp->gapDistance = tempScaffoldEdge->gapDistance;
                    temp->orientation = tempScaffoldEdge->mapPosition->orientation;
                    temp->next = shortContigHead1->next;
                    shortContigHead1->next = temp;
                    shortContigHead1->gapDistance++;    
                }
                
                tempScaffoldEdge = tempScaffoldEdge->next;
            }
            
            ShortContig * temp1 = shortContigHead->next;
            ShortContig * temp = shortContigHead1->next;
            
            while(temp!=NULL){
                long int tempContigIndex = temp->contigIndex;
                temp1 = shortContigHead1->next;
                while(temp1!=NULL){
                    if(temp1->contigIndex==tempContigIndex){
                        break;
                    }
                    temp1= temp1->next;
                }
                
                if(temp1==NULL){
                    ShortContig * temp2 = new ShortContig;
                    temp2->contigIndex = temp->contigIndex;
                    temp2->gapDistance = longGapDistance - temp->gapDistance;
                    temp2->orientation = temp->orientation;
                    temp2->next = shortContigHead->next;
                    shortContigHead->next = temp2;
                }
                
                temp = temp->next;
            }
            if(temp!=NULL){
                continue;
            }
            temp = shortContigHead->next;
            while(temp!=NULL){
                temp1 = temp->next;
                while(temp1!=NULL){
                    if(temp->gapDistance > temp1->gapDistance){
                        
                        long int tempGap = temp->gapDistance;
                        temp->gapDistance = temp1->gapDistance;
                        temp1->gapDistance = tempGap;
                        
                        tempGap = temp->contigIndex;
                        temp->contigIndex = temp1->contigIndex;
                        temp1->contigIndex = tempGap;
                        
                        tempGap = temp->orientation;
                        temp->orientation = temp1->orientation;
                        temp1->orientation = tempGap;
                        
                    }
                    temp1=temp1->next;
                }
                temp = temp->next;
            }
            
            temp = shortContigHead->next;
            while(temp!=NULL){
                
                visited[temp->contigIndex] = true;
                ContigSequence * insertShortContigSequence = new ContigSequence;
                insertShortContigSequence->index = temp->contigIndex;
                insertShortContigSequence->orientation = temp->orientation;
                insertShortContigSequence->gapDistance = abs(firstContigSequence->gapDistance - temp->gapDistance - scaffoldGraph[temp->contigIndex].length);
                firstContigSequence->next = insertShortContigSequence;
                insertShortContigSequence->next = tempContigSequence;
                firstContigSequence->gapDistance = temp->gapDistance;
                firstContigSequence = insertShortContigSequence;
                temp = temp->next;
            
            }
            
        }
        
        
        
        tempScaffoldSet = tempScaffoldSet->next;
        
    }
    

}

long int DetermineOrientationOfContigsHeuristic(ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation, bool ** fixIndex, double minScore){

}




long int DetermineOrientationOfContigs(ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation, bool ** fixIndex, double minScore, bool last){
    
    long int i = 0;
    long int j = 0;
    long int p = 1;
    long int c = 100;
    
    bool * contigIndex = new bool[contigCount];
    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
    }
    long int edgeNumber = 0;
    long int constraintNumber = 0;
    ScaffoldEdge * tempEdge = NULL;
    
    for(i=0;i<contigCount;i++){
        tempEdge = scaffoldGraph[i].outLink;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
        tempEdge = scaffoldGraph[i].inLink;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
        contigIndex[i] = false;
    }
    
    edgeNumber = edgeNumber/2;
    constraintNumber = 0;
    
    lprec *lp;
    int Ncol, *colno = NULL, ret = 0;
    REAL *row = NULL;
    
    Ncol = edgeNumber + contigCount;
    lp = make_lp(0, Ncol);
    if(lp == NULL){
        printf("couldn't construct a new linear programming model");
        exit(0);
    }
    double * weight = new double[edgeNumber]; 
    long int * edgeLeftNode = new long int[edgeNumber];
    long int * edgeRightNode = new long int[edgeNumber];
    
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL)){
        printf("couldn't new colno and row");
        exit(0);
    }
    
    set_add_rowmode(lp, TRUE);
    
    
    for(i=0;i<contigCount;i++){
       
        for(int q=0;q<2;q++){
            if(q==0){
                tempEdge = scaffoldGraph[i].outLink;
            }else{
                tempEdge = scaffoldGraph[i].inLink;
            }
            while(tempEdge!=NULL){
                
                bool token = false;
                if(tempEdge->fitNumber>=minScore && fixIndex[i][tempEdge->contigIndex]==true){
                    token = true;
                }
                if(tempEdge->fitNumber<minScore){
                    tempEdge = tempEdge->next;
                    continue;
                }
                
                if(index[i][tempEdge->contigIndex] == false && index[tempEdge->contigIndex][i] == false){
                    if(tempEdge->mapPosition->orientation == tempEdge->mapPosition->mateOrientation){
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = 1;
                        colno[j] = tempEdge->contigIndex+1; 
                        row[j++] = 1;
                        if(token == false){
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c+1)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }else{
                            if(!add_constraintex(lp, j, row, colno, LE, 1)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                            if(!add_constraintex(lp, j, row, colno, GE, 1)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }

                        j = 0;
                        colno[j] = i+1;
                        row[j++] = -1;
                        colno[j] = tempEdge->contigIndex+1; 
                        row[j++] = -1;
                        if(token == false){
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c-1)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }
                        
                        constraintNumber = constraintNumber+2;
                        
                    }else{
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = 1;
                        colno[j] = tempEdge->contigIndex+1; 
                        row[j++] = -1;
                        if(token == false){
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }else{
                            if(!add_constraintex(lp, j, row, colno, LE, 0)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                            if(!add_constraintex(lp, j, row, colno, GE, 0)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }
                        
                    
                        j = 0;
                        colno[j] = tempEdge->contigIndex+1;
                        row[j++] = 1;
                        colno[j] = i+1; 
                        row[j++] = -1; 
                        if(token == false){
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }
                        
                        constraintNumber = constraintNumber+2;
                        
                    }
                    
                    if(token == false){
                        j = 0;
                        colno[j] = contigCount + p;
                        row[j++] = 1;
                        add_constraintex(lp, j, row, colno, LE, 1);
                        j = 0;
                        colno[j] = contigCount + p;
                        row[j++] = 1;
                        add_constraintex(lp, j, row, colno, GE, 0);
                        constraintNumber = constraintNumber+2;
                        
                        weight[p-1] = tempEdge->fitNumber;
                        edgeLeftNode[p-1] = i;
                        edgeRightNode[p-1] = tempEdge->contigIndex;
                        p++;
                    }
                    
                    index[i][tempEdge->contigIndex] = true;
                    index[tempEdge->contigIndex][i] = true;

                }
                tempEdge = tempEdge->next;
            }    
        }
    }
    
    for(i=0;i<contigCount;i++){
        set_binary(lp, i+1, TRUE);
    }

    p--;
    
    set_add_rowmode(lp, FALSE);
    j=0;
    for(i=0;i<p;i++){
        colno[j] = contigCount + i + 1; 
        row[j] = weight[j];
        j++;
    }
    if(!set_obj_fnex(lp, j, row, colno)){
        printf("couldn't set_obj_fnex");
        exit(0);
    }

    set_bb_depthlimit(lp, 10);
    set_maxim(lp);
    set_scaling(lp, 128);  
    ret = solve(lp);

    if(!(ret==0 || ret ==1)){

    }
    
    REAL * pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
    get_primal_solution(lp,pv);

    double temp = 1;
    long int result = 0;
    for(i=contigCount+constraintNumber+1;i<constraintNumber+contigCount+p+1;i++){
        if(pv[i] < 1){
            if(minScore<0.11){
                DeleteSpecialScaffoldEdge(scaffoldGraph, edgeLeftNode[i-contigCount-constraintNumber-1], edgeRightNode[i-contigCount-constraintNumber-1]);
                DeleteSpecialScaffoldEdge(scaffoldGraph, edgeRightNode[i-contigCount-constraintNumber-1], edgeLeftNode[i-contigCount-constraintNumber-1]);
                result++;
            }
        }else{
            fixIndex[edgeRightNode[i-contigCount-constraintNumber-1]][edgeLeftNode[i-contigCount-constraintNumber-1]] = true;
            fixIndex[edgeLeftNode[i-contigCount-constraintNumber-1]][edgeRightNode[i-contigCount-constraintNumber-1]] = true;
        }
    }
    
    for(i=constraintNumber+1;i<contigCount+constraintNumber+1 && minScore < 0.11;i++){
        contigOrientation[i-constraintNumber-1] = pv[i];
    }
    
    delete [] contigIndex;
    for(i=0;i<contigCount;i++){
        delete [] index[i];
    }
    delete [] index;
    delete [] weight;
    delete [] edgeLeftNode;
    delete [] edgeRightNode;
    delete [] pv;
    
    delete_lp(lp);
    
    return result;
    
}


void SetGapDistance(ScaffoldGraph * scaffoldGraph, long int first, long int second, long int gapDistance){
    ScaffoldEdge * tempEdge = NULL;
    tempEdge = scaffoldGraph[first].outLink;
    while(tempEdge != NULL){
        if(tempEdge->contigIndex==second){
            tempEdge->gapDistance = gapDistance;
            return;
        }
        tempEdge = tempEdge->next;
    }
    tempEdge = scaffoldGraph[first].inLink;
    while(tempEdge != NULL){
        if(tempEdge->contigIndex==second){
            tempEdge->gapDistance = gapDistance;
            return;
        }
        tempEdge = tempEdge->next;
    }
                
}


long int * IterativeDetermineOrderOfContigs(ContigSet * contigSet, ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation, long int * tempContigOrder, long int * contigPosition, bool ** fixIndex, long int insertsize, long int std, long int lambda, double minScore){
    
    
    long int i = 0;
    long int j = 0;
    long int p = 1;
    long int c = allContigLength;
    
    long int contigCutOff = insertsize + lambda*std;
    
    bool * contigVisited = new bool[contigCount];
    
    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
        contigVisited[i] = false;
    }
    long int edgeNumber = 0;
    ScaffoldEdge * tempEdge = NULL;
    
    for(i=0;i<contigCount;i++){
        tempEdge = scaffoldGraph[i].outLink;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
        tempEdge = scaffoldGraph[i].inLink;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
    }
    
    edgeNumber = edgeNumber/2;
    
    lprec *lp;
    int Ncol, *colno = NULL, ret = 0;
    REAL *row = NULL;
    
    Ncol = contigCount + edgeNumber;
    lp = make_lp(0, Ncol);
    if(lp == NULL){
        printf("couldn't construct a new linear programming model");
        exit(0);
    }
    double * weight = new double[edgeNumber]; 
    long int * edgeLeftNode = new long int[edgeNumber];
    long int * edgeRightNode = new long int[edgeNumber];
    long int * gapDistance = new long int[edgeNumber];
    
    
    colno = (int *) malloc(Ncol * sizeof(*colno));
    row = (REAL *) malloc(Ncol * sizeof(*row));
    if((colno == NULL) || (row == NULL)){
        printf("couldn't new colno and row");
        exit(0);
    }
    
    set_add_rowmode(lp, TRUE);
    
    long int constraintNumber = 0;
    
    for(i=0;i<contigCount;i++){
        for(int q=0;q<2;q++){
            if(q==0){
                tempEdge = scaffoldGraph[i].outLink;
            }else{
                tempEdge = scaffoldGraph[i].inLink;
            }
            long int start = 0;
            while(tempEdge!=NULL){
                
                bool longIndex = false;
                if(tempEdge->fitNumber>=minScore && fixIndex[i][tempEdge->contigIndex]==true){
                    longIndex = true;
                }
                if(tempEdge->fitNumber<minScore){
                    tempEdge = tempEdge->next;
                    continue;
                }
                
                if(index[i][tempEdge->contigIndex] == false && index[tempEdge->contigIndex][i] == false){
                    
                    if((contigOrientation[i] ==0 && q ==0) || (contigOrientation[i]==1 && q==1)){
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = -1;
                        colno[j] = tempEdge->contigIndex+1; 
                        row[j++] = 1;
                        
                        if(longIndex == false){
                           
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c+scaffoldGraph[i].length+tempEdge->gapDistance)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }  
                        }else{
                            
                            if(!add_constraintex(lp, j, row, colno, LE, scaffoldGraph[i].length+insertsize+lambda*std)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                            
                            if(!add_constraintex(lp, j, row, colno, GE, scaffoldGraph[i].length+1)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }
                        
                        if(longIndex == false){
                            j = 0;
                            colno[j] = i+1;
                            row[j++] = 1;
                            colno[j] = tempEdge->contigIndex+1; 
                            row[j++] = -1;
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c-scaffoldGraph[i].length-tempEdge->gapDistance)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                            edgeLeftNode[p-1] = i;
                            edgeRightNode[p-1] = tempEdge->contigIndex;  
                        }
                        contigVisited[i] = true;
                        contigVisited[tempEdge->contigIndex] = true;
                    }
                    if((contigOrientation[i] ==1 && q ==0) || (contigOrientation[i]==0 && q==1)){
                        j = 0;
                        colno[j] = i+1;
                        row[j++] = 1;
                        colno[j] = tempEdge->contigIndex+1; 
                        row[j++] = -1;
                        
                        if(longIndex == false){
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c+scaffoldGraph[tempEdge->contigIndex].length+tempEdge->gapDistance)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                        }else if(longIndex == true){
                            
                            if(!add_constraintex(lp, j, row, colno, LE, scaffoldGraph[tempEdge->contigIndex].length+insertsize+lambda*std)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                            
                            if(!add_constraintex(lp, j, row, colno, GE, scaffoldGraph[tempEdge->contigIndex].length+1)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                            
                        }
                        
                        if(longIndex == false){
                            j = 0;
                            colno[j] = i+1;
                            row[j++] = -1;
                            colno[j] = tempEdge->contigIndex+1; 
                            row[j++] = 1;
                        
                            colno[j] = contigCount + p; 
                            row[j++] = c;
                            if(!add_constraintex(lp, j, row, colno, LE, c-scaffoldGraph[tempEdge->contigIndex].length-tempEdge->gapDistance)){
                                printf("couldn't add_constraintex");
                                exit(0);
                            }
                            edgeLeftNode[p-1] = tempEdge->contigIndex;
                            edgeRightNode[p-1] = i;
                        }
                        contigVisited[i] = true;
                        contigVisited[tempEdge->contigIndex] = true;
                    }
                    if(longIndex == false){
                        j = 0;
                        colno[j] = contigCount + p;
                        row[j++] = 1;
                        add_constraintex(lp, j, row, colno, LE, 1);
                        j = 0;
                        colno[j] = contigCount + p;
                        row[j++] = 1;
                        add_constraintex(lp, j, row, colno, GE, 0);
                        weight[p-1] = tempEdge->fitNumber;
                        gapDistance[p-1] = tempEdge->gapDistance;
                        p++;
                        constraintNumber = constraintNumber + 4; 
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }else{
                        constraintNumber = constraintNumber + 2; 
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }

                }
                tempEdge = tempEdge->next;
            }
                                     
        }
    }
 
    
    for(i=0;i<contigCount;i++){
        set_int(lp, i+1, TRUE);
    }

    p--;
    
    set_add_rowmode(lp, FALSE);
    
    j=0;
    for(i=0;i<p;i++){
        colno[j] = contigCount + i + 1; 
        row[j] = weight[j];
        j++;
    }
    
    if(!set_obj_fnex(lp, j, row, colno)){
        printf("couldn't set_obj_fnex");
        exit(0);
    }

    set_maxim(lp);
    
    set_bb_depthlimit(lp, 10);

    ret = solve(lp);

    if(!(ret==0 || ret ==1)){

        set_break_at_first(lp, true);
        ret = solve(lp);
        if(!(ret==0 || ret ==1)){
            return NULL;
        }
    }
    
    REAL * pv = new REAL[constraintNumber + contigCount + edgeNumber + 1];
    get_primal_solution(lp,pv);
    
    for(i=0;i<contigCount;i++){

        contigPosition[i] = pv[1+constraintNumber+i];
    }
    
    long int trueNumber = 0;
    long int realTrueNumber = 0;
    
    for(i=contigCount+constraintNumber+1;i<contigCount+constraintNumber+p+1;i++){
        if(pv[i] < 1){
            long int d = pv[1+constraintNumber+edgeRightNode[i-contigCount-constraintNumber-1]] - pv[1+constraintNumber+edgeLeftNode[i-contigCount-constraintNumber-1]] - scaffoldGraph[edgeLeftNode[i-contigCount-constraintNumber-1]].length;
            double varD = double(labs(d - gapDistance[i-contigCount-constraintNumber-1]))/(double)std;
            if(varD>3 || d < 0){
                continue;
            }
        }
        fixIndex[edgeRightNode[i-contigCount-constraintNumber-1]][edgeLeftNode[i-contigCount-constraintNumber-1]] = true;
        fixIndex[edgeLeftNode[i-contigCount-constraintNumber-1]][edgeRightNode[i-contigCount-constraintNumber-1]] = true;

    }
    
    for(i=0;i<contigCount;i++){
        delete [] index[i];
    }
    
    
    delete [] contigVisited;
    delete [] index;
    delete [] weight;
    delete [] edgeLeftNode;
    delete [] edgeRightNode;
    delete [] gapDistance;
    delete [] pv;
    
    delete_lp(lp);
    
}



SubGraph * ReachableNode(ScaffoldGraph * scaffoldGraph, bool * contigOrientation, long int contigCount, int inOrOutIndex){
    long int i = 0;
    long int j = 0;
    
    SubGraph * reachableNodeList = new SubGraph[contigCount];
    SubGraphNode * tempNode = NULL;
    SubGraphNode * lastNode = NULL;
    
    ScaffoldEdge * tempEdge = NULL;
    bool * visited = new bool[contigCount];
    
    for(i = 0; i<contigCount; i++){
        
        if(contigOrientation[i] == inOrOutIndex){
            tempEdge = scaffoldGraph[i].outLink;
        }else{
            tempEdge = scaffoldGraph[i].inLink;
        }
        
        for(long int t = 0; t<contigCount; t++){
            visited[t] = false;
        }
        
        j = 0;
        tempNode = NULL;
        while(tempNode!=NULL || (tempEdge != NULL && j==0) ){
            
            if(tempNode != NULL){
                if(contigOrientation[tempNode->index] == inOrOutIndex){
                    tempEdge = scaffoldGraph[tempNode->index].outLink;
                }else{
                    tempEdge = scaffoldGraph[tempNode->index].inLink;
                }
                
                if(tempEdge == NULL){
                    tempNode = tempNode->next;
                    j++;
                    continue;
                }
            }
            
            
            ScaffoldEdge * temp = tempEdge;
            long int adjacentNodeCount = 0;
            while(temp!=NULL){
                adjacentNodeCount++;
                temp = temp->next;
            }

            long int * gapDistance = new long int[adjacentNodeCount];
            long int * nodeIndex = new long int[adjacentNodeCount];
            temp = tempEdge;
            long int t = 0;
            while(temp!=NULL){
                gapDistance[t] = temp->gapDistance;
                nodeIndex[t] = temp->contigIndex;
                temp = temp->next;
                t++;
            }
            long int p = 0;
            for(t=0;t<adjacentNodeCount-1;t++){
                for(p=t+1;p<adjacentNodeCount;p++){
                    if(gapDistance[t]>gapDistance[p]){
                        
                        long int tempGap = gapDistance[t];
                        gapDistance[t] = gapDistance[p];
                        gapDistance[p] = tempGap;
                        
                        tempGap = nodeIndex[t];
                        nodeIndex[t] = nodeIndex[p];
                        nodeIndex[p] = tempGap;
                        
                    }
                }
            }
            
            t = 0;
            if(tempNode == NULL){

                if(adjacentNodeCount>0){
                    reachableNodeList[i].startNode = new SubGraphNode;
                    reachableNodeList[i].startNode->index = nodeIndex[0];
                    visited[nodeIndex[0]] = true;
                    lastNode = reachableNodeList[i].startNode;
                    tempNode = reachableNodeList[i].startNode;

                }
                t = 1;
            }
            if(t == 0){
                tempNode = tempNode->next;
            }
            for(t = t;t<adjacentNodeCount;t++){
                if(visited[nodeIndex[t]] == true){
                    continue;
                }
                lastNode->next = new SubGraphNode;
                lastNode->next->index = nodeIndex[t];
                lastNode = lastNode->next;
                visited[nodeIndex[t]] = true;
            }
            j++;
            
            delete [] gapDistance;
            delete [] nodeIndex;
        }
        
        
    }
    
    return reachableNodeList;
    

}

SubGraphNode * DFS(ScaffoldGraph * scaffoldGraph, long int contigIndex, SubGraphNode * subGraphNode, bool * contigOrientation, bool * visited){
    visited[contigIndex] = true;
    ScaffoldEdge * tempEdge = NULL;
    for(int i = 0; i<2;i++){
        if(i==0){
            tempEdge = scaffoldGraph[contigIndex].inLink;
        }else{
            tempEdge = scaffoldGraph[contigIndex].outLink;
        }
        while(tempEdge!=NULL){
            if(visited[tempEdge->contigIndex] == false){
                subGraphNode->next = new SubGraphNode;
                subGraphNode = subGraphNode->next;
                subGraphNode->index = tempEdge->contigIndex;
                subGraphNode = DFS(scaffoldGraph, tempEdge->contigIndex, subGraphNode, contigOrientation, visited);
            }
            tempEdge = tempEdge->next;
        }
    }
    return subGraphNode;
}

SubGraph * DFSTranverseScaffoldGraph(ScaffoldGraph * scaffoldGraph, long int contigCount, bool * contigOrientation){
    
    long int i = 0;
    long int j = 0;
    
    bool * visited = new bool[contigCount];
    for(i = 0; i<contigCount; i++){
        visited[i] = false;
    }
    
    SubGraph * subGraph = NULL;
    SubGraph * firstSubGraph = NULL;
    
    
    long int connectedSubGraph = 0;
    
    for(i = 0; i<contigCount; i++){
        if(visited[i] == false){
            if(subGraph == NULL){
                subGraph = new SubGraph;
                firstSubGraph = subGraph;
            }else{
                subGraph->next =new SubGraph;
                subGraph = subGraph->next;
            }
            subGraph->startNode = new SubGraphNode;
            subGraph->startNode->index = i;
            DFS(scaffoldGraph, i, subGraph->startNode, contigOrientation, visited);
            connectedSubGraph++;
        }
    }
    long int subGraphCount = 0;
    SubGraph * first = firstSubGraph;
    while(firstSubGraph!=NULL){
        subGraphCount++;
        firstSubGraph = firstSubGraph->next;
    }
    return first;
}

long int RemovePositionOverlap(ScaffoldGraph * scaffoldGraph, SubGraph * subGraph, SubGraph * reachableNode, long int contigCount, long int * contigPosition){
    
    
    long int i = 0;
    long int j = 0;
    long int p = 1;
    long int c = 100000;
    
    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
    }
    
    bool ** overlapIndex = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        overlapIndex[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            overlapIndex[i][j] = false;
        }   
    }
    
    
    long int edgeNumber = 0;
    long int constraintNumber = 0;
    ScaffoldEdge * tempEdge = NULL;
    
    for(i=0;i<contigCount;i++){
        tempEdge = scaffoldGraph[i].outLink;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
        tempEdge = scaffoldGraph[i].inLink;
        while(tempEdge!=NULL){
            edgeNumber++;
            tempEdge = tempEdge->next;
        }

    }
    
    
    edgeNumber = edgeNumber/2;
    
    constraintNumber = 0;
    
    long int positionOverlapCount = 0;

    while(subGraph!=NULL){
        SubGraphNode * tempNode = subGraph->startNode;
        SubGraphNode * secondNode = subGraph->startNode;
        while(tempNode!=NULL){
            secondNode = subGraph->startNode;
            while(secondNode!=NULL){
                if(secondNode->index == tempNode->index){
                    secondNode = secondNode->next;
                    continue;
                }
                i = tempNode->index;
                j = secondNode->index;     
                
                if(contigPosition[i] + scaffoldGraph[i].length > contigPosition[j] 
                    && contigPosition[i] < contigPosition[j] && contigPosition[j]!=-1){
                    long int overlapLength = 0;
                    if(contigPosition[i] + scaffoldGraph[i].length < contigPosition[j] + scaffoldGraph[j].length){
                        overlapLength = contigPosition[i] + scaffoldGraph[i].length - contigPosition[j];
                    }else{
                        overlapLength = scaffoldGraph[j].length;
                    }
                    
                    double a = (double)(overlapLength)/(double)(scaffoldGraph[i].length);
                    double b = (double)(overlapLength)/(double)(scaffoldGraph[j].length);

                    if(a>0.3 || b>0.3){
                        if(overlapIndex[i][j] == false){
                            positionOverlapCount++;
                            overlapIndex[i][j] = true;
                            
                            SubGraphNode * tempReachableNode = reachableNode[i].startNode;
                            SubGraphNode * tempReachableNode1 = reachableNode[j].startNode;
                            long int previusReachableNodeIndex = i;
                            long int previusReachableNodeIndex1 = j;
                            bool reachIndex = false;
                            
                            while(tempReachableNode != NULL){
                                tempReachableNode1 = reachableNode[j].startNode;
                                previusReachableNodeIndex1 = j;
                                while(tempReachableNode1 != NULL){
                                    if(tempReachableNode->index == tempReachableNode1->index){
                                        reachIndex = true;
                                        break;
                                    }
                                    previusReachableNodeIndex1 = tempReachableNode1->index;
                                    tempReachableNode1 = tempReachableNode1->next;
                                }
                                if(reachIndex == true){
                                    break;
                                }
                                previusReachableNodeIndex = tempReachableNode->index;
                                tempReachableNode = tempReachableNode->next;
                                
                            }
                            
                            double score = 0;
                            double score1 = 0;
                            
                            if(reachIndex == true){
                                
                                long int currentReachableNodeIndex = tempReachableNode->index;
                                tempEdge = scaffoldGraph[currentReachableNodeIndex].outLink;
                                while(tempEdge!=NULL){
                                    if(tempEdge->contigIndex == previusReachableNodeIndex){
                                        score = tempEdge->fitNumber;
                                    }
                                    if(tempEdge->contigIndex == previusReachableNodeIndex1){
                                        score1 = tempEdge->fitNumber;
                                    }
                                    tempEdge = tempEdge->next;
                                }
                                tempEdge = scaffoldGraph[currentReachableNodeIndex].inLink;
                                while(tempEdge!=NULL){
                                    if(tempEdge->contigIndex == previusReachableNodeIndex){
                                        score = tempEdge->fitNumber;
                                    }
                                    if(tempEdge->contigIndex == previusReachableNodeIndex1){
                                        score1 = tempEdge->fitNumber;
                                    }
                                    tempEdge = tempEdge->next;
                                }
                                
                                
                                
                                if(score > score1){
                                    DeleteSpecialScaffoldEdge(scaffoldGraph, currentReachableNodeIndex, previusReachableNodeIndex1);
                                    DeleteSpecialScaffoldEdge(scaffoldGraph, previusReachableNodeIndex1, currentReachableNodeIndex);
                                }else{
                                    DeleteSpecialScaffoldEdge(scaffoldGraph, currentReachableNodeIndex, previusReachableNodeIndex);
                                    DeleteSpecialScaffoldEdge(scaffoldGraph, previusReachableNodeIndex, currentReachableNodeIndex);
                                }
                            }
                            
                        }
                    }
                }    

                secondNode = secondNode->next;
            }
            tempNode = tempNode->next;
        }
        subGraph = subGraph->next;
    }
    
}


ScaffoldSet * GetScaffoldSet(ScaffoldGraph * scaffoldGraph, long int contigCount, long int contigCutOff){
    
    long int i = 0;
    long int j = 0;
    
    int * printIndex = new int[contigCount];
    
    for(i=0;i<contigCount;i++){
        printIndex[i] = 0;
    }
    
    ScaffoldSet * scaffoldSet = new ScaffoldSet;
    ScaffoldSet * scaffoldSetHead = scaffoldSet;
    
    int orientation = 1;
    
    for(i=0;i<contigCount;i++){
        
        if(printIndex[i]==1){
            continue;
        }
        
        if(scaffoldGraph[i].length<contigCutOff){
            continue;
        }
        
        ScaffoldEdge * temp = scaffoldGraph[i].outLink;
        while(temp!=NULL){
            if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                break;
            }
            temp=temp->next;
        }
              
        if(temp==NULL){
            continue;
        }
        
        
        ScaffoldSet * tempScaffoldSet = new ScaffoldSet;
        scaffoldSet->next = tempScaffoldSet;
        scaffoldSet = tempScaffoldSet;
        
        ContigSequence * tempContigSequence = new ContigSequence;
        tempContigSequence->index = i;
        tempContigSequence->orientation = 0;
        scaffoldSet->contigSequence = tempContigSequence;
        
        orientation = temp->mapPosition->mateOrientation;
        j = i;
        while(temp!=NULL && printIndex[j]!=1){
                
            printIndex[j]=1;  
            j = temp->contigIndex; 
                        
            ContigSequence * tempContigSequence1 = new ContigSequence;
            tempContigSequence1->index = temp->contigIndex;
            tempContigSequence->gapDistance = temp->gapDistance;
            tempContigSequence->next = tempContigSequence1;     
            tempContigSequence = tempContigSequence1;
            tempContigSequence1 = NULL;       
            if(orientation==0){

                tempContigSequence->orientation = 0;
                temp = scaffoldGraph[temp->contigIndex].outLink;
                
                while(temp!=NULL){
                    if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                        break;
                    }
                    temp=temp->next;
                }
                
                if(temp==NULL){
                    printIndex[j]=1;
                    break;
                }
                
                int cc = SearchScaffoldEdge(temp->contigIndex, scaffoldSet->contigSequence);
                if(cc==0){
                    printIndex[j]=1;
                    break;
                }
                if(temp->mapPosition->mateOrientation == 0){
                    orientation = 0;
                }else{
                    orientation = 1;
                }
            }else if(orientation==1){
                
                tempContigSequence->orientation = 1;
                temp = scaffoldGraph[temp->contigIndex].inLink;
                
                while(temp!=NULL){
                    if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                        break;
                    }
                    temp=temp->next;
                }
                
                if(temp==NULL){
                    printIndex[j]=1;
                    break;
                }
                int cc = SearchScaffoldEdge(temp->contigIndex, scaffoldSet->contigSequence);
                if(cc==0){
                    printIndex[j]=1;
                    break;
                }
                if(temp->mapPosition->mateOrientation == 0){
                    orientation = 0;
                }else{
                    orientation = 1;
                }
            } 
                                                               
        }
        
        temp = scaffoldGraph[i].inLink;
        
         while(temp!=NULL){
            if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                break;
            }
            temp=temp->next;
        }
        
        if(temp==NULL){
            continue;
        }
        
        orientation = temp->mapPosition->mateOrientation;
        j = i;
        printIndex[j] = 0;
        while(temp!=NULL && printIndex[j]!=1){
            
            int cc = SearchScaffoldEdge(temp->contigIndex, scaffoldSet->contigSequence);
            if(cc==0){
                printIndex[j]=1;
                break;
            }
            
               
            printIndex[j]=1;
            j = temp->contigIndex;
            
            ContigSequence * tempContigSequence1 = new ContigSequence;
            tempContigSequence1->index = temp->contigIndex;
            tempContigSequence1->gapDistance = temp->gapDistance;
            tempContigSequence1->next = scaffoldSet->contigSequence;
            scaffoldSet->contigSequence = tempContigSequence1;
                   
            if(orientation==0){
                scaffoldSet->contigSequence->orientation = 1;
                temp = scaffoldGraph[temp->contigIndex].outLink;
                
                while(temp!=NULL){
                    if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                        break;
                    }
                    temp=temp->next;
                }
                
                if(temp==NULL){
                    printIndex[j]=1;
                    break;
                }
                
                if(temp->mapPosition->mateOrientation == 0){
                    orientation = 0;
                }else{
                    orientation = 1;
                }
            }else if(orientation==1){
                scaffoldSet->contigSequence->orientation = 0;
                temp = scaffoldGraph[temp->contigIndex].inLink;
                
                while(temp!=NULL){
                    if(scaffoldGraph[temp->contigIndex].length>=contigCutOff){
                        break;
                    }
                    temp=temp->next;
                }
                
                if(temp==NULL){
                    printIndex[j]=1;
                    break;
                }
                if(temp->mapPosition->mateOrientation == 1){
                    orientation = 1;
                }else{
                    orientation = 0;
                }
            } 
                                                     
        }
                
    }
    
    for(i = 0; i<contigCount; i++){
        if(printIndex[i]!=1 && scaffoldGraph[i].length>=contigCutOff){
            ScaffoldSet * tempScaffoldSet = new ScaffoldSet;
            scaffoldSet->next = tempScaffoldSet;
            scaffoldSet = tempScaffoldSet;
            
            ContigSequence * tempContigSequence = new ContigSequence;
            tempContigSequence->index = i;
            tempContigSequence->orientation = 0;
            scaffoldSet->contigSequence = tempContigSequence;  
            
        }
    }
    
    
    scaffoldSet = scaffoldSetHead->next;
    
    return scaffoldSetHead->next;
    
}

void OutPutScaffoldTag(ScaffoldSet * scaffoldSet, char * result){
    
    ofstream ocout;
    ocout.open(result);
    
    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    long int i = 0;
    long int j = 0;
    while(tempScaffoldSet!=NULL){
        
        ContigSequence * temp = tempScaffoldSet->contigSequence;
        while(temp!=NULL){
            
            ocout<<temp->index<<",";
            temp = temp->next;
        }
        if(tempScaffoldSet->contigSequence!=NULL){
            ocout<<endl;
        }
        tempScaffoldSet = tempScaffoldSet->next;
    }
    
    ocout.close();
    
}

void OutPutScaffoldSet(ScaffoldSet * scaffoldSet, ContigSet * contigSet, long int contigCount, char * result){    
    long int i = 0;
    long int j = 0;
    
    bool * printContigIndex = new bool[contigCount];
    for(i=0;i<contigCount;i++){
        printContigIndex[i] = false;
    }

    ScaffoldSet * tempScaffoldSet = scaffoldSet;
    scaffoldSet = tempScaffoldSet;
    
    char * scaffoldTagFileName = new char[100];
    strcpy(scaffoldTagFileName, result);
    strcat(scaffoldTagFileName, "_Scaffold_Tag.fa");
    ofstream ocoutTag;
    ocoutTag.open(scaffoldTagFileName);
    
    char * scaffoldSetFileName = new char[100];   
    strcpy(scaffoldSetFileName, result);
    strcat(scaffoldSetFileName, "_ScaffoldSet.fa");
    ofstream ocout1;
    ocout1.open(scaffoldSetFileName);
    j = 0;
    while(tempScaffoldSet!=NULL){
        ContigSequence * tempContigSequence = tempScaffoldSet->contigSequence;
        if(tempContigSequence == NULL){
            tempScaffoldSet = tempScaffoldSet->next;
            continue;
        }
        ocout1<<">"<<j<<endl;
        
        long int allLength = 0;
        
        while(tempContigSequence!=NULL){
            
            if(printContigIndex[tempContigSequence->index] == true){
                ocoutTag<<"--";
            }
            
            printContigIndex[tempContigSequence->index] = true;
            ocoutTag<<tempContigSequence->index<<"("<<tempContigSequence->gapDistance<<"--"<<tempContigSequence->orientation<<"),";
            
            if(tempContigSequence->orientation==1){
                char * tempChar1 = new char[strlen(contigSet[tempContigSequence->index].contig)+1];
                ReverseComplement(contigSet[tempContigSequence->index].contig, tempChar1);
                ocout1<<tempChar1;                    
                allLength = allLength + strlen(contigSet[tempContigSequence->index].contig);     
            }else{
                ocout1<<contigSet[tempContigSequence->index].contig;
                allLength = allLength + strlen(contigSet[tempContigSequence->index].contig);
            }
            
            if(tempContigSequence->gapDistance>0 && tempContigSequence->next!=NULL){
                
                for(int tt = 0; tt<tempContigSequence->gapDistance; tt++){
                    ocout1<<"N";
                }
                
                allLength = allLength + tempContigSequence->gapDistance;
            }
            
            if(tempContigSequence->gapDistance<=0 && tempContigSequence->next!=NULL){
                
                for(int tt = 0; tt<5; tt++){
                    ocout1<<"N";
                }
                allLength = allLength + tempContigSequence->gapDistance;
                
            }
           
            tempContigSequence = tempContigSequence->next;
           
        }
        
        ocoutTag<<endl;
        ocout1<<endl;
        j++;
        tempScaffoldSet = tempScaffoldSet->next;
    }
    
    
    for(i = 0; i<contigCount; i++){
        if(printContigIndex[i]==false && contigSet[i].contig != NULL){
            ocout1<<">"<<j<<endl;
            ocout1<<contigSet[i].contig<<endl;
            ocoutTag<<i<<","<<endl;
            j++;
        }
    }

}

long int GetScaffoldGraphEdgeNumber(ScaffoldGraph * scaffoldGraph, long int contigCount){
    
    long int i = 0;
    long int edgeNumber = 0;
    
    char * outputFile = new char[20];
    strcpy(outputFile, "edgeNumberCount.fa");
    ofstream ocout;
    ocout.open(outputFile);
    
    ScaffoldEdge * tempEdge = NULL;
    for(i =0; i<contigCount; i++){
        
        tempEdge = scaffoldGraph[i].outLink;
        
        while(tempEdge!=NULL){
            ocout<<i<<"--"<<tempEdge->contigIndex<<"--"<<tempEdge->fitNumber<<endl;
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
        
        tempEdge = scaffoldGraph[i].inLink;
        
        while(tempEdge!=NULL){
            ocout<<i<<"--"<<tempEdge->contigIndex<<"--"<<tempEdge->fitNumber<<endl;
            edgeNumber++;
            tempEdge = tempEdge->next;
        }
    
    }
    ocout<<"-------------------------------------"<<endl;
    ocout<<"allEdgeNumber:"<<edgeNumber<<"--"<<edgeNumber/2<<endl;
    
    double cutoff = 1;
    
    while(cutoff>=0){
        edgeNumber = 0;
        for(i =0; i<contigCount; i++){
        
            tempEdge = scaffoldGraph[i].outLink;
        
            while(tempEdge!=NULL){
                if(tempEdge->fitNumber>=cutoff && tempEdge->fitNumber<cutoff+0.1){
                    ocout<<cutoff<<"--"<<i<<"--"<<tempEdge->contigIndex<<"--"<<tempEdge->fitNumber<<endl;
                    edgeNumber++;
                }
                tempEdge = tempEdge->next;
            }
            
            tempEdge = scaffoldGraph[i].inLink;
        
            while(tempEdge!=NULL){
                if(tempEdge->fitNumber>=cutoff && tempEdge->fitNumber<cutoff+0.1){
                    ocout<<cutoff<<"--"<<i<<"--"<<tempEdge->contigIndex<<"--"<<tempEdge->fitNumber<<endl;
                    edgeNumber++;
                }
                tempEdge = tempEdge->next;
            }
            
        }
        ocout<<cutoff<<"--EdgeNumber:"<<edgeNumber<<"--"<<edgeNumber/2<<endl;
        cutoff = cutoff - 0.1;
    }    
    
    
    
    
}


void OutPutScaffoldSet(ScaffoldGraph * scaffoldGraph, ContigSet * contigSet, double cutoff, long int contigCount, char * suffix){
    long int i = 0;
    long int j = 0;
    
    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
    }
    
    long int edgeNumber = 0;
    
    char * outputFile = new char[30];
    strcpy(outputFile, "contigWeight_");
    strcat(outputFile, suffix);
    char cutoffFile[6];
    sprintf(cutoffFile, "%04.2f", cutoff);
    strcat(outputFile, cutoffFile);
    char cutoffFile1[6];
    strcpy(cutoffFile1, ".fa");
    strcat(outputFile, cutoffFile1);
    
    ofstream ocout;
    ocout.open(outputFile);
  
    ScaffoldEdge * tempEdge = NULL;
    for(i =0; i<contigCount; i++){
        
        long int a = 0;
        while(a<=1){
            if(a==0){
                tempEdge = scaffoldGraph[i].outLink;
            }else{
                tempEdge = scaffoldGraph[i].inLink;
            }
            
            while(tempEdge != NULL){
                if(tempEdge->fitNumber >= cutoff - 0.000001 && tempEdge->fitNumber < cutoff + 0.1){
                    if(tempEdge->mapPosition->orientation != tempEdge->mapPosition->mateOrientation && index[i][tempEdge->contigIndex] == false){
                        ocout<<">"<<edgeNumber<<endl;
                        ocout<<contigSet[i].contig;
                        for(j =0;j<tempEdge->gapDistance;j++){
                            ocout<<"N";
                        }
                        ocout<<contigSet[tempEdge->contigIndex].contig<<endl;
                        edgeNumber++;
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }
                    if(tempEdge->mapPosition->orientation == tempEdge->mapPosition->mateOrientation && index[i][tempEdge->contigIndex] == false){
                        ocout<<">"<<edgeNumber<<endl;
                        ocout<<contigSet[i].contig;
                        for(j =0;j<tempEdge->gapDistance;j++){
                              ocout<<"N";
                        }
                        char * tempChar1 = new char[strlen(contigSet[tempEdge->contigIndex].contig)+1];
                        ReverseComplement(contigSet[tempEdge->contigIndex].contig, tempChar1);
                        ocout<<tempChar1<<endl;
                        edgeNumber++;
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }
                }
                tempEdge = tempEdge->next;
            }   
            a++;
        }
        
    }
    
} 

void OutPutScaffoldSetAll(ScaffoldGraph * scaffoldGraph, ContigSet * contigSet, double min, long int contigCount, char * outputFile){
    long int i = 0;
    long int j = 0;
    
    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
    }
    
    long int edgeNumber = 0;
    
    ofstream ocout;
    ocout.open(outputFile);
  
    ScaffoldEdge * tempEdge = NULL;
    for(i =0; i<contigCount; i++){
        long int a = 0;
        while(a<=1){
            if(a==0){
                tempEdge = scaffoldGraph[i].outLink;
            }else{
                tempEdge = scaffoldGraph[i].inLink;
            }
            
            while(tempEdge != NULL){
                if(tempEdge->fitNumber >= min - 0.000001){
                    if(tempEdge->mapPosition->orientation != tempEdge->mapPosition->mateOrientation && index[i][tempEdge->contigIndex] == false){
                        ocout<<">"<<edgeNumber<<endl;
                        ocout<<contigSet[i].contig;
                        for(j =0;j<tempEdge->gapDistance;j++){
                            ocout<<"N";
                        }
                        ocout<<contigSet[tempEdge->contigIndex].contig<<endl;
                        edgeNumber++;
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }
                    if(tempEdge->mapPosition->orientation == tempEdge->mapPosition->mateOrientation && index[i][tempEdge->contigIndex] == false){
                        ocout<<">"<<edgeNumber<<endl;
                        ocout<<contigSet[i].contig;
                        for(j =0;j<tempEdge->gapDistance;j++){
                              ocout<<"N";
                        }
                        char * tempChar1 = new char[strlen(contigSet[tempEdge->contigIndex].contig)+1];
                        ReverseComplement(contigSet[tempEdge->contigIndex].contig, tempChar1);
                        ocout<<tempChar1<<endl;
                        edgeNumber++;
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }
                }
                tempEdge = tempEdge->next;
            }   
            a++;
        }
    }
    
} 

void OutPutScaffoldSetAllNumber(ScaffoldGraph * scaffoldGraph, ContigSet * contigSet, int cutoff, long int contigCount){
    long int i = 0;
    long int j = 0;
    
    bool ** index = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        index[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            index[i][j] = false;
        }   
    }
    
    long int edgeNumber = 0;
    
    char * outputFile = new char[20];
    strcpy(outputFile, "G1.fa");
    
    ofstream ocout;
    ocout.open(outputFile);
  
    ScaffoldEdge * tempEdge = NULL;
    for(i =0; i<contigCount; i++){
        long int a = 0;
        while(a<=1){
            if(a==0){
                tempEdge = scaffoldGraph[i].outLink;
            }else{
                tempEdge = scaffoldGraph[i].inLink;
            }
            
            while(tempEdge != NULL){
                if(tempEdge->edgeWeight >= cutoff){
                    if(tempEdge->mapPosition->orientation != tempEdge->mapPosition->mateOrientation && index[i][tempEdge->contigIndex] == false){
                        ocout<<">"<<edgeNumber<<endl;
                        ocout<<contigSet[i].contig;
                        for(j =0;j<tempEdge->gapDistance;j++){
                            ocout<<"N";
                        }
                        ocout<<contigSet[tempEdge->contigIndex].contig<<endl;
                        edgeNumber++;
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }
                    if(tempEdge->mapPosition->orientation == tempEdge->mapPosition->mateOrientation && index[i][tempEdge->contigIndex] == false){
                        ocout<<">"<<edgeNumber<<endl;
                        ocout<<contigSet[i].contig;
                        for(j =0;j<tempEdge->gapDistance;j++){
                              ocout<<"N";
                        }
                        char * tempChar1 = new char[strlen(contigSet[tempEdge->contigIndex].contig)+1];
                        ReverseComplement(contigSet[tempEdge->contigIndex].contig, tempChar1);
                        ocout<<tempChar1<<endl;
                        edgeNumber++;
                        index[i][tempEdge->contigIndex] = true;
                        index[tempEdge->contigIndex][i] = true;
                    }
                }
                tempEdge = tempEdge->next;
            }   
            a++;
        }
        
    }
    
} 


int AddShortContigToScaffoldSet(ScaffoldSet * scaffoldSet, long int contigCount, bool * index){
    
    long int i = 0;
    
    ScaffoldSet * temp = scaffoldSet;
    ScaffoldSet * pre = NULL;
    ContigSequence * tempContigSequence = NULL;
    
    while(temp!=NULL){
        tempContigSequence = temp->contigSequence;
        while(tempContigSequence != NULL){
            index[tempContigSequence->index] = 1;
            tempContigSequence = tempContigSequence->next;
        }
        pre = temp;
        temp = temp->next;
    }
    
    temp = pre;
    for(i = 0; i<contigCount; i++){
        if(index[i]!=1){
            ScaffoldSet * tempScaffoldSet = new ScaffoldSet;
            temp->next = tempScaffoldSet;
            temp = tempScaffoldSet;
            
            ContigSequence * tempContigSequence = new ContigSequence;
            tempContigSequence->index = i;
            tempContigSequence->orientation = 0;
            temp->contigSequence = tempContigSequence;  
        }
    }
    
}

void WriteScaffoldGraph(ScaffoldGraph * scaffoldGraph, bool * contigOrientation, long int contigCount){
    
    long int i = 0;
    long int j = 0;
    
    char * outPutFileName = new char[20];
    strcpy(outPutFileName, "scaffoldGraph.fa");
    ofstream ocout;
    ocout.open(outPutFileName);
    
    ScaffoldEdge * tempEdge = NULL;
    
    for(i = 0; i<contigCount; i++){
        if(contigOrientation[i]==0){
            tempEdge = scaffoldGraph[i].outLink;
        }else{
            tempEdge = scaffoldGraph[i].inLink;
        }
        while(tempEdge!=NULL){
            ocout<<i<<" "<<tempEdge->contigIndex<<" "<<tempEdge->fitNumber<<endl;
            tempEdge = tempEdge->next;
        }
    }
    
}

double * weightOrder(ScaffoldGraph * scaffoldGraph, long int contigCount){
    long int i = 0;
    long int j = 0;
    
    long int edgeNumber = 0;
    
    ScaffoldEdge * tempEdge = NULL;
        
    for(long int t = 0; t<contigCount; t++){
        for(long int p = 0; p<2;p++){
            if(p==0){
                tempEdge = scaffoldGraph[t].outLink;
            }else{
                tempEdge = scaffoldGraph[t].inLink;
            }
            
            while(tempEdge!=NULL){
                edgeNumber++;
                tempEdge = tempEdge->next;
            }
                
        }    
    } 
    
    double * weight = new double[edgeNumber];
    for(long int t = 0; t<contigCount; t++){
        for(long int p = 0; p<2;p++){
            if(p==0){
                tempEdge = scaffoldGraph[t].outLink;
            }else{
                tempEdge = scaffoldGraph[t].inLink;
            }
            
            while(tempEdge!=NULL){
                weight[i] = tempEdge->fitNumber;
               
                i++;
                tempEdge = tempEdge->next;
            }
                
        }    
    } 
    
    
    
    for(i=0;i<edgeNumber-1;i++){
        for(j=i;j<edgeNumber;j++){
            if(weight[i]<weight[j]){
                double temp = weight[i];
                weight[i] = weight[j];
                weight[j] = temp;
            }
        }
    }
    
    double * order = new double[10];
    long int interval = edgeNumber/10;
    //cout<<"interval--"<<interval<<endl;
    j=1;
    for(i=0;i<edgeNumber;i++){
        if(i==j*interval){
            order[j-1] = weight[i-1];
            //cout<<"order--"<<order[j-1]<<endl;
            j++;    
        }
    }
    order[9] = weight[edgeNumber-1];
    //cout<<"order--"<<order[9]<<endl;
    
    //cout<<"25%=="<<weight[(long int)(edgeNumber*0.75)]<<endl;
    //cout<<"30%=="<<weight[(long int)(edgeNumber*0.7)]<<endl;
    //cout<<"20%=="<<weight[(long int)(edgeNumber*0.8)]<<endl;
    //cout<<"21%=="<<weight[(long int)(edgeNumber*0.79)]<<endl;
    //cout<<"22%=="<<weight[(long int)(edgeNumber*0.78)]<<endl;
    //cout<<"23%=="<<weight[(long int)(edgeNumber*0.77)]<<endl;
    //cout<<"24%=="<<weight[(long int)(edgeNumber*0.76)]<<endl;
    
    delete [] weight;
    weight = NULL;
    return order;
    
    
    
}

ScaffoldSet * OptimizeScaffoldSet(ContigSet * contigSet, ScaffoldSet * scaffoldSet, ScaffoldGraph * scaffoldGraph, long int & contigCount, long int realContigCount, long int * contigLength, long int insertsize, long int std, long int lambda){
    
    long int i = 0;
    long int j = 0;
    long int contigLengthCutOff = insertsize + lambda*std;
    
    
    allContigLength = 0;
    for(i=0;i<contigCount;i++){
        allContigLength = allContigLength + contigLength[i];
    }
    allContigLength = 2*allContigLength; 
    
    bool * contigOrientation = new bool[contigCount];
    BFSGraph * bfsGraph = new BFSGraph[contigCount];
    bool ** visitedIndex = new bool*[contigCount];
    for(i=0;i<contigCount;i++){
        visitedIndex[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            visitedIndex[i][j] = false;
        }   
    }
    long int * contigOrder = NULL;
    long int * contigPosition = new long int[contigCount];
    
    long int conflictOrientation = 1;
    i = 0;
    bool last = false;
    
    double minScore = 0.9;
    time_t timep;
    time(&timep);
    while(minScore > 0.0001){
        if(minScore<0.11){
            last = true;
        }
        conflictOrientation = DetermineOrientationOfContigs(scaffoldGraph, contigCount, contigOrientation, visitedIndex, minScore, last);
        if(conflictOrientation == -1){
            break;
        }
        minScore = minScore - 0.1;
    }
    
    for(i=0;i<contigCount;i++){
        visitedIndex[i] = new bool[contigCount];
        for(j=0;j<contigCount;j++){
            visitedIndex[i][j] = false;
        }   
    }
    
    minScore = 0.9;
    long int * noResult = NULL;
    while(minScore > 0.0001 && conflictOrientation != -1){
        noResult = IterativeDetermineOrderOfContigs(contigSet, scaffoldGraph,contigCount,contigOrientation, contigOrder, contigPosition, visitedIndex, insertsize, std, lambda,minScore);
        if(noResult == NULL){
            break;
        }
        minScore = minScore - 0.1;
    }
    
    
    for(i=0;i<contigCount && noResult!=NULL && conflictOrientation != -1;i++){
        
        ScaffoldEdge * tempEdge = scaffoldGraph[i].outLink;
        ScaffoldEdge * pre = NULL;
        long int tempContigIndex = 0;
        while(tempEdge!=NULL){
            pre = tempEdge->next;
            
            if(visitedIndex[i][tempEdge->contigIndex] == false && tempEdge->fitNumber >=0.1){
                
                tempContigIndex = tempEdge->contigIndex;
                DeleteSpecialScaffoldEdge(scaffoldGraph, i, tempContigIndex);
                
                DeleteSpecialScaffoldEdge(scaffoldGraph, tempContigIndex, i);
                
                tempEdge = pre;
                continue;
            }
            tempEdge = tempEdge->next;
        }
        tempEdge = scaffoldGraph[i].inLink;
        while(tempEdge!=NULL){
            pre = tempEdge->next;
            if(visitedIndex[i][tempEdge->contigIndex] == false && tempEdge->fitNumber >=0.1){
                
                tempContigIndex = tempEdge->contigIndex;
                DeleteSpecialScaffoldEdge(scaffoldGraph, i, tempContigIndex);
                DeleteSpecialScaffoldEdge(scaffoldGraph, tempContigIndex, i);
                
                tempEdge = pre;
                continue;
            }
            tempEdge = tempEdge->next;
        }
    }
    
    SubGraph * subGraph = DFSTranverseScaffoldGraph(scaffoldGraph,contigCount,contigOrientation);
    SubGraph * reachableNode = ReachableNode(scaffoldGraph, contigOrientation,contigCount,0);
    RemovePositionOverlap(scaffoldGraph,subGraph,reachableNode,contigCount,contigPosition);
    
    subGraph = DFSTranverseScaffoldGraph(scaffoldGraph,contigCount,contigOrientation);
    reachableNode = ReachableNode(scaffoldGraph, contigOrientation,contigCount,0);
    RemovePositionOverlap(scaffoldGraph,subGraph,reachableNode,contigCount,contigPosition);
    
    subGraph = DFSTranverseScaffoldGraph(scaffoldGraph,contigCount,contigOrientation);
    reachableNode = ReachableNode(scaffoldGraph, contigOrientation,contigCount,1);
    RemovePositionOverlap(scaffoldGraph,subGraph,reachableNode,contigCount,contigPosition);
    
    subGraph = DFSTranverseScaffoldGraph(scaffoldGraph,contigCount,contigOrientation);
    reachableNode = ReachableNode(scaffoldGraph, contigOrientation,contigCount,1);
    RemovePositionOverlap(scaffoldGraph,subGraph,reachableNode,contigCount,contigPosition);
    
    
    bool * index = new bool[contigCount];
    long int tempContigCount = contigCount;
    for(i = 0; i< contigCount; i++){
        index[i] = 0;
    } 
    
    char * tempFileName = new char[100];
    
    ScaffoldSet * tempScaffoldSet = GetScaffoldSet(scaffoldGraph,contigCount,contigLengthCutOff);     
    
    InsertShortContigInScaffoldSet(scaffoldGraph, tempScaffoldSet, contigCount, contigLengthCutOff);
    
    for(j=0;j<1;j++){        
        
        BFSScaffolding(scaffoldGraph, tempScaffoldSet, contigCount, contigLengthCutOff);
        
    }
    
    ScaffoldSet * ss = tempScaffoldSet;
    ScaffoldSet * tt = NULL;
    while(ss!=NULL){
        if(ss->contigSequence == NULL){
            if(tt == NULL){
                tempScaffoldSet = ss->next;
                ss->next = NULL;
                ss = tempScaffoldSet;
                continue;
            }else{
                tt->next = ss->next;
                ss->next = NULL;
                ss = tt->next;
                continue;
            }
        }
        tt = ss;
        ss = ss->next;
    }
    
    MergeTwoScaffoldSet(scaffoldSet, contigCount, tempScaffoldSet);

    GetScaffoldSetLength(tempScaffoldSet, contigCount, contigLength);
    
  
    return tempScaffoldSet;  
}

#endif
