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


int p1_position=1,p3_position=3,p5_position=5,p7_position=7;
int MT = -1;
int DC = -1;
int N = -1;



/*SmallOS*/
const int PLACES = 9;
const int TRANS = 8;
int BOUNDS = 7000;
int** model;

using namespace MEDDLY;

FILE_output meddlyout(stdout);

void buildModel(const char* order)
{
  const int modelTest[TRANS][PLACES+1] = {
    {0,1,0,0,0,0,0,0,0,-1},  // freeMem
    {0,-1,1,-1,0,-1,0,0,0,0},  // loadMem
    {0,0,0,1,-1,1,0,0,0,1},  // endUnload
    {0,0,0,-1,1,-1,-1,0,0,0},  // startUnload
    {0,0,-1,1,0,1,0,0,0,1},  // endLoad
    {0,0,0,0,0,0,-1,-1,1,0},  // startNext
    {0,0,0,0,0,0,1,1,-1,0},  // suspend
    {0,0,0,0,0,0,0,-1,1,-1},  // starFirst
  };
  
  p1_position = 1;
  p3_position = 3;
  p5_position = 5;
  p7_position = 7;
  
  for (int i=0; i<PLACES; i++)
    {
    if (order[i]=='1')
      p1_position = i+1;
    if (order[i]=='3')
      p3_position = i+1;
    if (order[i]=='5')
      p5_position = i+1;
    if (order[i]=='7')
      p7_position = i+1;
    }
  
  model = (int**) malloc(TRANS * sizeof(int*));
  
  for(int i=0;i<TRANS;i++)
    {
    model[i] = (int*) malloc((PLACES+2) * sizeof(int));
    for(int j=1;j<(PLACES+1);j++)
      {
      model[i][j]=modelTest[i][order[j-1]-'0'];
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
  printf("\torder: the order of variables:123456789\n");
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
int batchsize = 256;
const char* lfile = 0;

for (int i=1; i<argc; i++)
{
  if(i==1)
    MT = atoi(argv[i]);
  else if(i==2)
    DC = atoi(argv[i]);
  
  if (strcmp("-O", argv[i])==0) {
    if(i+1 < argc)
      {
      buildModel(argv[i+1]);
      i++;
      }
  }
}
  
if (argc<5) return usage(argv[0]);
if (MT<0) return usage(argv[0]);
if (DC<0) return usage(argv[0]);
BOUNDS = 2;
domain* d = 0;

try {

MEDDLY::initialize();

timer start;


printf("+----------------------------------------------------+\n");
printf("|         Initializing smallOS with %-4d %-4d        |\n", MT,DC);
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
initialState[p1_position]=initialState[p3_position]=MT;initialState[p5_position]=DC;initialState[p7_position]=2*DC;

N = 2*DC>MT?(2*DC):MT;
  
method ='i';
  std::cout<<"\n********************";
  std::cout<<"\n     implicit";
  std::cout<<"\n********************";
if('i' == method)
{

//CREATE FORESTS
  forest* inmdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL,p);

  expert_domain* dm = static_cast<expert_domain*>(inmdd->useDomain());
  dm->enlargeVariableBound(p1_position, false, MT+1);
  dm->enlargeVariableBound(p3_position, false, MT+1);
  dm->enlargeVariableBound(p5_position, false, DC+1);
  dm->enlargeVariableBound(p7_position, false, 2*DC+1);
  
    
  //ADD INITIAL STATE
  dd_edge first(inmdd);
  dd_edge reachable(inmdd);
  inmdd->createEdge(&initialState, 1, first);
    //outmdd->createEdge(&initialState, 1, reachable);


  //CREATE RELATION
  satimpl_opname::implicit_relation* T = new satimpl_opname::implicit_relation(inmdd,inmdd);

  start.note_time();
  buildImplicitRelation(model, TRANS, PLACES, BOUNDS, T);
  printf("\nNext-state function construction took %.4e seconds\n",
  start.get_last_interval() / 1000000.0);
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
  start.get_last_interval() / 1000000.0);
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
  dd_edge mxd_edge_all = T->buildMxdForest();
  
  mxd_edge_all.show(meddlyout,2);
  
}

return 0;
}
catch (MEDDLY::error e) {
printf("Caught MEDDLY error: %s\n", e.getName());
return 1;
}
}






