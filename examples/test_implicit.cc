// $Id$

/*
 Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
 Copyright (C) 2011, Iowa State University Research Foundation, Inc.
 
 This library is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published
 by the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with this library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <string.h>
#include <fstream>

#define _MEDDLY_WITHOUT_IOSTREAM_

#include "meddly.h"
#include "meddly_expert.h"
#include "simple_model.h"
#include "timer.h"
#include "loggers.h"

// #define DUMP_NSF
// #define DUMP_REACHABLE


char** sample_model;
int p1_position;
int PLACES = 9;
int TRANS = 7;
int BOUNDS = 1000;
int N = -1;
int* nxtArray = (int*) malloc(TRANS * PLACES * sizeof(int));


const char* kanban[] = {
  "X-+......-",  // Tin1 TA
  "X.-+....-.",  // Tin1 TA
  "X..-+....+",  // Tin1 TA
  "X...-+...-",  // Tin1 TA
  "X....-+.+.",  // Tin1 TA
  "X.....-+.+",  // Tin1 TA
  "X+.....-..",  // Tin1 TA
  
};

using namespace MEDDLY;

FILE_output meddlyout(stdout);

void buildModel(const char* order)
{
  const char* sample_model1[] = {
    "X-++"  // T1
  };
  
  p1_position = 1;
  for (int i=0; i<PLACES; i++)
    {
    if (order[i]=='1')
      p1_position = i+1;
    }
  
  sample_model = (char**) malloc(TRANS * sizeof(char*));
  
  for(int i=0;i<TRANS;i++)
    {
    sample_model[i] = (char*) malloc((PLACES+2) * sizeof(char));
    for(int j=1;j<(PLACES+1);j++)
      {
      sample_model[i][j]=sample_model1[i][order[j-1]-'0'];
      std::cout<< order[j-1]-'0' << "\n";
      std::cout<< sample_model1[i][order[j-1]-'0'] << "\n";
      }
    sample_model[i][PLACES+1] = '\0';
    std::cout<<sample_model[i];
    }
  
}


int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  printf("\nUsage: %s nnnn \n\n", name);
  printf("\tnnnn: number of initial tokens\n");
  return 1;
}

void printStats(const char* who, const forest* f)
{
  printf("%s stats:\n", who);
  const expert_forest* ef = (expert_forest*) f;
  ef->reportStats(meddlyout, "\t",
                  expert_forest::HUMAN_READABLE_MEMORY  |
                  expert_forest::BASIC_STATS | expert_forest::EXTRA_STATS |
                  expert_forest::STORAGE_STATS | expert_forest::HOLE_MANAGER_STATS
                  );
}


int main(int argc, const char** argv)
{
  char method = 'i';
  int batchsize = 256;
  const char* lfile = 0;
  
  for (int i=1; i<argc; i++) {
    N = atoi(argv[i]);
  }
  
  
  if (N<0) return usage(argv[0]);
  
  domain* d = 0;
  try {
    
    MEDDLY::initialize();
    
    timer start;
    
    printf("+----------------------------------------------------+\n");
    printf("|         Initializing sample_model with %-4d        |\n", N);
    printf("+----------------------------------------------------+\n");
    fflush(stdout);
    
    // Initialize domain
    int* sizes = new int[PLACES];
    for (int i=PLACES-1; i>=0; i--) sizes[i] = BOUNDS;
    d = createDomainBottomUp(sizes, PLACES);
    delete[] sizes;
    forest::policies pr(true);
    pr.setPessimistic();
    forest::policies p(false);
    p.setPessimistic();
    
    
    
    if('i' == method)
      {
      
      
      /*KANBAN
      satimpl_opname::relation_node* TA1 = new satimpl_opname::relation_node(10,1,1); //- 2
      satimpl_opname::relation_node* TA2 = new satimpl_opname::relation_node(11,2,2); //+ 3
      
      satimpl_opname::relation_node* TB1 = new satimpl_opname::relation_node(20,2,1); //- 4
      satimpl_opname::relation_node* TB2 = new satimpl_opname::relation_node(21,3,4); //+ 5
      
      satimpl_opname::relation_node* TC1 = new satimpl_opname::relation_node(30,2,1); //+ 6
      satimpl_opname::relation_node* TC2 = new satimpl_opname::relation_node(31,3,6); //- 7
      
      satimpl_opname::relation_node* TD2 = new satimpl_opname::relation_node(41,4,4); //+ 8
      
      satimpl_opname::relation_node* TE1 = new satimpl_opname::relation_node(50,6,1); //- 9
      satimpl_opname::relation_node* TE2 = new satimpl_opname::relation_node(51,7,9); //+ 10
      
      satimpl_opname::relation_node* TF1 = new satimpl_opname::relation_node(60,6,1); //+ 11
      satimpl_opname::relation_node* TF2 = new satimpl_opname::relation_node(61,7,11); //- 12
      
      satimpl_opname::relation_node* TG2 = new satimpl_opname::relation_node(71,8,9); //+ 13
      
      satimpl_opname::relation_node* TH1 = new satimpl_opname::relation_node(80,1,1); //+ 14
      satimpl_opname::relation_node* TH2 = new satimpl_opname::relation_node(81,4,14);//- 15
      satimpl_opname::relation_node* TH3 = new satimpl_opname::relation_node(82,5,15);//- 16
      satimpl_opname::relation_node* TH4 = new satimpl_opname::relation_node(83,6,16); //+ 17
      satimpl_opname::relation_node* TH5 = new satimpl_opname::relation_node(84,9,17);//- 18
      satimpl_opname::relation_node* TH6 = new satimpl_opname::relation_node(85,10,18);//+ 19
      
      satimpl_opname::relation_node* TI1 = new satimpl_opname::relation_node(90,10,1);//- 20
      satimpl_opname::relation_node* TI2 = new satimpl_opname::relation_node(91,11,20);//+ 21
      
      satimpl_opname::relation_node* TJ1 = new satimpl_opname::relation_node(100,10,1);//+ 22
      satimpl_opname::relation_node* TJ2 = new satimpl_opname::relation_node(101,11,22);//- 23
      
      satimpl_opname::relation_node* TK2 = new satimpl_opname::relation_node(111,12,20);//+ 24
      
      satimpl_opname::relation_node* TL1 = new satimpl_opname::relation_node(120,5,1);//+ 25
      satimpl_opname::relation_node* TL2 = new satimpl_opname::relation_node(121,8,25);//- 26
      satimpl_opname::relation_node* TL2m = new satimpl_opname::relation_node(122,9,26);//+ 27
      satimpl_opname::relation_node* TL3 = new satimpl_opname::relation_node(123,12,27);//- 28
      satimpl_opname::relation_node* TL4 = new satimpl_opname::relation_node(124,13,28);//- 29
      satimpl_opname::relation_node* TL5 = new satimpl_opname::relation_node(125,14,29);//+ 30
      
      satimpl_opname::relation_node* TM1 = new satimpl_opname::relation_node(130,14,1);//- 31
      satimpl_opname::relation_node* TM2 = new satimpl_opname::relation_node(131,15,31);//+ 32
      
      satimpl_opname::relation_node* TN1 = new satimpl_opname::relation_node(140,14,1);//+ 33
      satimpl_opname::relation_node* TN2 = new satimpl_opname::relation_node(141,15,33);//- 34
      
      satimpl_opname::relation_node* TO1 = new satimpl_opname::relation_node(150,13,1);//+ 35
      satimpl_opname::relation_node* TO2 = new satimpl_opname::relation_node(151,16,35);//- 36
    
      satimpl_opname::relation_node* TP2 = new satimpl_opname::relation_node(161,16,31);//+ 37
      */
    
      
      /*ERATOSHENES
      satimpl_opname::relation_node* T11 = new satimpl_opname::relation_node(10,4,1); //0 2
      satimpl_opname::relation_node* T12 = new satimpl_opname::relation_node(11,9,2); //-1 3
      
      satimpl_opname::relation_node* T21 = new satimpl_opname::relation_node(12,1,1); //0 4
      satimpl_opname::relation_node* T22 = new satimpl_opname::relation_node(13,9,4);//-1 5
      
      satimpl_opname::relation_node* T31 = new satimpl_opname::relation_node(14,5,4);//-1 6
      
      satimpl_opname::relation_node* T41 = new satimpl_opname::relation_node(15,2,1);//0 7
      satimpl_opname::relation_node* T42 = new satimpl_opname::relation_node(16,5,7);//-1 8
      
      satimpl_opname::relation_node* T51 = new satimpl_opname::relation_node(17,8,7);//-1 9
      
      satimpl_opname::relation_node* T61 = new satimpl_opname::relation_node(18,7,4);//-1 10
      
      satimpl_opname::relation_node* T71 = new satimpl_opname::relation_node(19,3,1);//-1 11
      satimpl_opname::relation_node* T72 = new satimpl_opname::relation_node(20,7,11);//0 12
      
      satimpl_opname::relation_node* T81 = new satimpl_opname::relation_node(21,3,4);//-1 13*/
      
      
      
      //CREATE FORESTS
      forest* inmdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL,p);
      forest* outmdd = inmdd;
      
      //ADD INITIAL STATE
      int* initialState;
      initialState = new int[PLACES + 1];
      for(int g = 1;g <= PLACES;g++) initialState[g] = 0;
      initialState[7]=N*4; initialState[8]=N*3; initialState[9]=N*2;
      dd_edge first(inmdd);
      dd_edge reachable(outmdd);
      inmdd->createEdge(&initialState, 1, first);
      outmdd->createEdge(&initialState, 1, reachable);
      
      
      //CREATE RELATION
      satimpl_opname::implicit_relation* T = new satimpl_opname::implicit_relation(inmdd,outmdd);
     
      nxtArray = buildImplicitRelation(kanban, TRANS, PLACES, T);
      
      /*KANBAN
      T->registerNode(false,TA1);T->registerNode(true,TA2);
      T->registerNode(false,TB1);T->registerNode(true,TB2);
      T->registerNode(false,TC1);T->registerNode(true,TC2);
      T->registerNode(true,TD2);
      T->registerNode(false,TE1);T->registerNode(true,TE2);
      T->registerNode(false,TF1);T->registerNode(true,TF2);
      T->registerNode(true,TG2);
      T->registerNode(false,TH1);T->registerNode(false,TH2);T->registerNode(false,TH3);T->registerNode(false,TH4);T->registerNode(false,TH5);T->registerNode(true,TH6);
      T->registerNode(false,TI1);T->registerNode(true,TI2);
      T->registerNode(false,TJ1);T->registerNode(true,TJ2);
      T->registerNode(true,TK2);
      T->registerNode(false,TL1);T->registerNode(false,TL2);T->registerNode(false,TL2m);T->registerNode(false,TL3);T->registerNode(false,TL4);T->registerNode(true,TL5);
      T->registerNode(false,TM1);T->registerNode(true,TM2);
      T->registerNode(false,TN1);T->registerNode(true,TN2);
      T->registerNode(false,TO1);T->registerNode(true,TO2);
      T->registerNode(true,TP2);*/
      
      /*ERATOSTHENES
       T->registerNode(false,T11);T->registerNode(true,T12);
       T->registerNode(false,T12);T->registerNode(true,T22);
       T->registerNode(true,T31);
       T->registerNode(false,T41);T->registerNode(true,T42);
       T->registerNode(true,T51);
       T->registerNode(true,T61);
       T->registerNode(false,T71);T->registerNode(true,T72);
       T->registerNode(true,T81);
      */
      
 
      specialized_operation* sat = 0;
      
      
      printf("\nBuilding reachability set using saturation implicit relation");
      if (0==SATURATION_IMPL_FORWARD) {
        throw error(error::UNKNOWN_OPERATION);
      }
      sat = SATURATION_IMPL_FORWARD->buildOperation(T);
      
      if (0==sat) {
        throw error(error::INVALID_OPERATION);
      }
      sat->compute(first, reachable);
      
      start.note_time();
      printf("\nReachability set construction took %.4e seconds\n",
             start.get_last_interval() / 1000000.0);
      fflush(stdout);
      
      #ifdef DUMP_REACHABLE
      printf("Reachable states:\n");
      reachable.show(meddlyout, 2);
      #endif
      
      printStats("MDD", outmdd);
      fflush(stdout);
      
      double c;
      apply(CARDINALITY, reachable, c);
      operation::showAllComputeTables(meddlyout, 3);
      
      printf("Approx. %g reachable states\n", c);
      
      
      }
    return 0;
  }
  catch (MEDDLY::error e) {
    printf("Caught MEDDLY error: %s\n", e.getName());
    return 1;
  }
}

long MEDDLY::satimpl_opname::relation_node::nextOf(long i)
{
  long ans = i+nxtArray[this->getID()];
  if((ans>=0)&&(ans<BOUNDS))
    return  ans;
  else
    return -1;
}

/*KANBAN
long MEDDLY::satimpl_opname::relation_node::nextOf(long i)
{
  
  switch(this->getSignature())
  {
    case 10: if(1<=i) return i-1; else return -1; break;
    case 11: if(1+i<N+1) return i+1; else return -1; break;
    case 20: if(1<=i) return i-1; else return -1; break;
    case 21: if(1+i<N+1) return i+1; else return -1; break;
    case 30: if(1+i<N+1) return i+1; else return -1; break;
    case 31: if(1<=i) return i-1; else return -1; break;
    case 41: if(1+i<N+1) return i+1; else return -1; break;
    case 50: if(1<=i) return i-1; else return -1; break;
    case 51: if(1+i<N+1) return i+1; else return -1; break;
    case 60: if(1+i<N+1) return i+1; else return -1; break;
    case 61: if(1<=i) return i-1; else return -1; break;
    case 71: if(1+i<N+1) return i+1; else return -1; break;
    case 80: if(1+i<N+1) return i+1; else return -1; break;
    case 81: if(1<=i) return i-1; else return -1; break;
    case 82: if(1<=i) return i-1; else return -1; break;
    case 83: if(1+i<N+1) return i+1; else return -1; break;
    case 84: if(1<=i) return i-1; else return -1; break;
    case 85: if(1+i<N+1) return i+1; else return -1; break;
    case 90: if(1<=i) return i-1; else return -1; break;
    case 91:  if(1+i<N+1) return i+1; else return -1; break;
    case 100: if(1+i<N+1) return i+1; else return -1; break;
    case 101: if(1<=i) return i-1; else return -1; break;
    case 111: if(1+i<N+1) return i+1; else return -1; break;
    case 120: if(1+i<N+1) return i+1; else return -1; break;
    case 121: if(1<=i) return i-1; else return -1; break;
    case 122: if(1+i<N+1) return i+1; else return -1; break;
    case 123: if(1<=i) return i-1; else return -1; break;
    case 124: if(1<=i) return i-1; else return -1; break;
    case 125: if(1+i<N+1) return i+1; else return -1; break;
    case 130: if(1<=i) return i-1; else return -1; break;
    case 131: if(1+i<N+1) return i+1; else return -1; break;
    case 140: if(1+i<N+1) return i+1; else return -1; break;
    case 141: if(1<=i) return i-1; else return -1; break;
    case 150: if(1+i<N+1) return i+1; else return -1; break;
    case 151: if(1<=i) return i-1; else return -1; break;
    case 161: if(1+i<N+1) return i+1; else return -1; break;
    default: return i;
  }
  
}*/



/*ERATOSTHENES
long MEDDLY::satimpl_opname::relation_node::nextOf(long i)
{
  
 
  switch(this->getSignature())
  {
    case 10: return i; break;
    case 11: if(1<=i) return i-1; else return -1; break;

    case 12: return i; break;
    case 13: if(1<=i) return i-1; else return -1; break;

    case 14: if(1<=i) return i-1; else return -1; break;
    case 15: return i; break;
    case 16: if(1<=i) return i-1; else return -1; break;

    case 17: if(1<=i) return i-1; else return -1; break;

    case 18: if(1<=i) return i-1; else return -1; break;

    case 19: return i; break;
    case 20: if(1<=i) return i-1; else return -1; break;

    case 21: if(1<=i) return i-1; else return -1; break;
    default: return i;
  }
}*/


