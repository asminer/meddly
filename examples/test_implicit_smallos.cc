
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
#include <iostream>

#define _MEDDLY_WITHOUT_IOSTREAM_

#include "../src/meddly.h"
#include "simple_model.h"
#include "../timing/timer.h"

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
  printf("\nUsage: %s mt dc <-O> order\n\n", name);
  printf("\tnnnn: number of initial tokens\n");
  printf("\torder: the order of variables:123456789\n");
return 1;
}

void printStats(const char* who, const forest* f)
{
printf("%s stats:\n", who);
f->reportStats(meddlyout, "\t",
HUMAN_READABLE_MEMORY  |
BASIC_STATS | EXTRA_STATS |
STORAGE_STATS | HOLE_MANAGER_STATS
);
}


int main(int argc, const char** argv)
{
char method ;

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
domain* dm = 0;

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
dm = domain::createBottomUp(sizes, PLACES);
delete[] sizes;
policies pr(true);
  pr.setPessimistic();
policies p(false);
  p.setPessimistic();
  // p.setQuasiReduced();

N = 2*DC>MT?(2*DC):MT;

  std::cout<<"\n********************";
  std::cout<<"\n     implicit";
  std::cout<<"\n********************";

//CREATE FORESTS
  forest* inmdd = forest::create(dm, 0, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL,p);
  forest* relmxd = forest::create(dm, 1, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL,pr);

  dm->enlargeVariableBound(p1_position, false, MT+1);
  dm->enlargeVariableBound(p3_position, false, MT+1);
  dm->enlargeVariableBound(p5_position, false, DC+1);
  dm->enlargeVariableBound(p7_position, false, 2*DC+1);


  //ADD INITIAL STATE
  dd_edge first(inmdd);
  dd_edge reachable(inmdd);

  minterm initState(inmdd);
  for(unsigned g = 1; g <= PLACES; g++) initState.setVar(g, 0);
  initState.setVar(p1_position, MT);
  initState.setVar(p3_position, MT);
  initState.setVar(p5_position, DC);
  initState.setVar(p7_position, 2*DC);
  initState.buildFunction(first);

    //outmdd->createEdge(&initialState, 1, reachable);


  //CREATE RELATION
  implicit_relation* T = new implicit_relation(inmdd,relmxd,inmdd);

  start.note_time();
  buildImplicitRelation(model, TRANS, PLACES, BOUNDS, inmdd, relmxd, T);
  printf("\nNext-state function construction took %.4e seconds\n",
  start.get_last_seconds());
  saturation_operation* sat = 0;


  printf("\nBuilding reachability set using saturation implicit relation");
  sat = SATURATION_IMPL_FORWARD(inmdd, T, inmdd);

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
  compute_table::showAll(meddlyout, 3);

  printf("Approx. %g reachable states\n", c);


  /* Building Mxd From Implicit */
 /* dd_edge mxd_edge_all = T->buildMxdForest();

  mxd_edge_all.show(meddlyout,2);*/


return 0;
}
catch (MEDDLY::error e) {
printf("Caught MEDDLY error: %s\n", e.getName());
return 1;
}
}






