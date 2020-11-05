
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



int p7_position=7,p8_position=8,p9_position=9;
int N = -1;


 /*SwimmingPool*/
 const int PLACES = 2;
 const int TRANS = 1;
 int BOUNDS = -1;
 std::vector<std::map<int, std::map<int, int>>> model;

using namespace MEDDLY;

FILE_output meddlyout(stdout);


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


int main()
{
  BOUNDS = 2;
  char method ;
  std::map<int, std::map<int, int>> one_event;
  std::map<int, int> arc;
  arc.insert({1,-1}); //-p1;
  arc.insert({2,1}); //+p2;
  one_event.insert({1,arc}); // p1' = p1 - p1 + p2
  arc.clear();
  arc.insert({2,-1}); //-p2;
  arc.insert({1,1}); //+p1;
  one_event.insert({2,arc}); // p2' = p2 - p2 + p1
  arc.clear();
  model.push_back(one_event);
  
  
  std::map<int, std::map<int, int>> two_event;
  arc.insert({0,-1}); //-1;
  two_event.insert({2,arc}); // p2' = p2 - 1
  arc.clear();
  arc.insert({2,-1}); //-p2;
  two_event.insert({1,arc}); // p1' = p1 - p2
  model.push_back(two_event);
  
  
  domain* d = 0;
  
  try {
    
    MEDDLY::initialize();
    
    timer start;
    
    
    printf("+----------------------------------------------------+\n");
    printf("|        Initializing swimming pool with %-4d        |\n", N);
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
    // p.setQuasiReduced();
    
    //INITIAL STATE
    int* initialState;
    initialState = new int[PLACES + 1];
    for(int g = 1;g <= PLACES;g++) initialState[g] = 0;
    initialState[1]=4; initialState[2]=0;
    
    method ='i';
    std::cout<<"\n********************";
    std::cout<<"\n     Implicit";
    std::cout<<"\n********************";
    if('i' == method)
      {
      
      //CREATE FORESTS
      forest* inmdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL,p);
      forest* relmxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL,pr);
      
      expert_domain* dm = static_cast<expert_domain*>(inmdd->useDomain());
      
      dm->enlargeVariableBound(1, false, 5);
      dm->enlargeVariableBound(2, false, 5);
      //dm->enlargeVariableBound(3, false, 1);
      
      
      //ADD INITIAL STATE
      dd_edge first(inmdd);
      dd_edge reachable(inmdd);
      inmdd->createEdge(&initialState, 1, first);
      //outmdd->createEdge(&initialState, 1, reachable);
        printf("\nAdded initial state\n");
      
      
      //CREATE RELATION
      satmdimpl_opname::md_implicit_relation* T = new satmdimpl_opname::md_implicit_relation(inmdd, relmxd, inmdd);
      printf("\nAdded the T state\n");
      
      start.note_time();
      buildGenericImplicitRelation(model, TRANS, relmxd, T);
      printf("\nNext-state function construction took %.4e seconds\n",
             start.get_last_seconds());
      specialized_operation* sat = 0;
      T->setConfirmedStates(first);
      
      printf("\nBuilding reachability set using saturation implicit relation");
      if (0==SATURATION_MDIMPL_FORWARD) {
        throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
      }
      sat = SATURATION_MDIMPL_FORWARD->buildOperation(T);
      
      if (0==sat) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
      }
      sat->compute(first, reachable);
      
      start.note_time();
      printf("\nReachability set construction took %.4e seconds\n",
             start.get_last_seconds());
      fflush(stdout);
      
//#ifdef DUMP_REACHABLE
      printf("Reachable states:\n");
      reachable.show(meddlyout, 2);
//#endif
      
      printStats("MDD", inmdd);
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










