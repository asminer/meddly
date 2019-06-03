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

#include "../src/meddly.h"
#include "../src/meddly_expert.h"
#include "simple_model.h"
#include "../src/timer.h"
#include "../src/loggers.h"

// #define DUMP_NSF
// #define DUMP_REACHABLE


int** model;
int p1_position=1,p13_position=13,p5_position=5,p9_position=9;
int N = -1;

/*Kanban*/
 const int PLACES = 16;
 const int TRANS = 16;
 int BOUNDS = -1;//N+1



using namespace MEDDLY;

FILE_output meddlyout(stdout);

void buildModel(const char* order)
{
  
  int const modelTest[TRANS][PLACES+1] = {
    {0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},  // Tin1 TA
    {0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},  // Tr1 TB
    {0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0},  // Tb1 TC
    {0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0},  // Tg1 TD
    {0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0},  // Tr2 TE
    {0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0},  // Tb2 TF
    {0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0},  // Tg2 TG
    {0,1,0,0,-1,-1,1,0,0,-1,1,0,0,0,0,0,0},  // Ts1_23 TH
    {0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0},  // Tr3 TI
    {0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0},  // Tb3 TJ
    {0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0},  // Tg3 TK
    {0,0,0,0,0,1,0,0,-1,1,0,0,-1,-1,1,0,0},  // Ts23_4 TL
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0},  // Tr4 TM
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0},  // Tb4 TN
    {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1},  // Tout4 TO
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1}   // Tg4 TP
  };

  
  p1_position = 1;
  p13_position = 13;
  p5_position = 5;
  p9_position = 9;
  for (int i=0; i<PLACES; i++)
    {
    if (order[i]=='1')
      p1_position = i+1;
    if (order[i]=='D')
      p13_position = i+1;
    if (order[i]=='5')
      p5_position = i+1;
    if (order[i]=='9')
      p9_position = i+1;
    }
  
  model = (int**) malloc(TRANS * sizeof(int*));
  
  for(int i=0;i<TRANS;i++)
    {
    model[i] = (int*) malloc((PLACES+2) * sizeof(int));
    for(int j=1;j<(PLACES+1);j++)
      {
      if((order[j-1]>'0')&&(order[j-1]<='9'))
        model[i][j]=modelTest[i][order[j-1]-'0'];
      else
        model[i][j]=modelTest[i][order[j-1]-'7'];
        
      }
    }
  
}


int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  printf("\nUsage: %s nnnn <-O> order\n\n", name);
  printf("\tnnnn: number of initial tokens\n");
  printf("\torder: the order of variables:123456789ABCDEFG \n");
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
  char method ;
  
  for (int i=1; i<argc; i++)
    {
    if(i==1)
      N = atoi(argv[i]);

    
    if (strcmp("-O", argv[i])==0) {
      if(i+1 < argc)
        {
        buildModel(argv[i+1]);
        i++;
        }
    }
    }
  BOUNDS = 2;
  
  if (argc<4) return usage(argv[0]);
  if (N<0) return usage(argv[0]);
  
  domain* d = 0;
  try {
    
    MEDDLY::initialize();
    
    timer start;
    
    
    printf("+----------------------------------------------------+\n");
    printf("|         Initializing Kanban with %-4d        |\n", N);
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
    
    //INITIAL STATE
    int* initialState;
    initialState = new int[PLACES + 1];
    for(int g = 1;g <= PLACES;g++) initialState[g] = 0;
    initialState[p1_position]=initialState[p13_position]=initialState[p5_position]=initialState[p9_position]=N;
    
    method ='i';
    std::cout<<"\n********************";
    std::cout<<"\n     Implicit";
    std::cout<<"\n********************";
    
    if('i' == method)
      {
      
      //CREATE FORESTS
      forest* inmdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL,p);
      forest* relmxd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL,p);
      
      expert_domain* dm = static_cast<expert_domain*>(inmdd->useDomain());

      dm->enlargeVariableBound(p1_position, false, N+1);
      dm->enlargeVariableBound(p13_position, false, N+1);
      dm->enlargeVariableBound(p5_position, false, N+1);
      dm->enlargeVariableBound(p9_position, false, N+1);
    
      
      
      //ADD INITIAL STATE
      dd_edge first(inmdd);
      dd_edge reachable(inmdd);
      inmdd->createEdge(&initialState, 1, first);
      
      
      //CREATE RELATION
      satimpl_opname::implicit_relation* T = new satimpl_opname::implicit_relation(inmdd,relmxd,inmdd);
      
      start.note_time();
      buildImplicitRelation(model, TRANS, PLACES, BOUNDS, T);
      printf("\nNext-state function construction took %.4e seconds\n",
             start.get_last_seconds());
      specialized_operation* sat = 0;
      
      
      printf("\nBuilding reachability set using saturation implicit relation");
      if (0==SATURATION_IMPL_FORWARD) {
        throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
      }
      sat = SATURATION_IMPL_FORWARD->buildOperation(T);
      
      if (0==sat) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
      }
      sat->compute(first, reachable);
      
      start.note_time();
      printf("\nReachability set construction took %.4e seconds\n",
             start.get_last_seconds());
      fflush(stdout);
      
      #ifdef DUMP_REACHABLE
      printf("Reachable states:\n");
      reachable.show(meddlyout, 2);
      #endif
      
      printStats("MDD", inmdd);
      fflush(stdout);
      
      double c;
      apply(CARDINALITY, reachable, c);
      operation::showAllComputeTables(meddlyout, 3);
      printf("Approx. %g reachable states\n", c);
      
      /* Building Mxd From Implicit */
      
      
      /*dd_edge mxd_edge_all = T->buildMxdForest();
      
      mxd_edge_all.show(meddlyout,2);*/

      
      }
    
    return 0;
  }
  catch (MEDDLY::error e) {
    printf("Caught MEDDLY error: %s\n", e.getName());
    return 1;
  }
}



