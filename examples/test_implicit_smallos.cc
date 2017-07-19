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
int MT = -1;
int DC = -1;



/*SmallOS*/
int PLACES = 9;
int TRANS = 8;
int BOUNDS = 100;
const char* model[] = {
"X+.......-",  // freeMem
"X-+-.-....",  // loadMem
"X..+-+...+",  // endUnload
"X..-+--...",  // startUnload
"X.-+.+...+",  // endLoad
"X.....--+.",  // startNext
"X.....++-.",  // suspend
"X......-+-",  // starFirst
};


using namespace MEDDLY;

FILE_output meddlyout(stdout);


int usage(const char* who)
{
/* Strip leading directory, if any: */
const char* name = who;
for (const char* ptr=who; *ptr; ptr++) {
if ('/' == *ptr) name = ptr+1;
}
printf("\nUsage: %s MT DC \n\n", name);
printf("\tMT: number of MT tokens\n");
printf("\tDC: number of DC tokens\n");
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

if(argc<3) return usage(argv[0]);
MT = atoi(argv[1]);
DC = atoi(argv[2]);

if (MT<0) return usage(argv[0]);
if (DC<0) return usage(argv[0]);

domain* d = 0;
try {

MEDDLY::initialize();

timer start;


printf("+----------------------------------------------------+\n");
printf("|     Initializing sample_model with %-4d %-4d        |\n", MT,DC);
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
initialState[1]=initialState[3]=MT;initialState[5]=DC;initialState[7]=2*DC;


method = 'k';
  std::cout<<"\n********************";
  std::cout<<"\n       ksat";
  std::cout<<"\n********************";
if('k' == method)
{
forest* mdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
forest* mxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);


dd_edge init_state(mdd);
mdd->createEdge(&initialState, 1, init_state);


dd_edge nsf(mxd);
satpregen_opname::pregen_relation* ensf = 0;
specialized_operation* sat = 0;

ensf = new satpregen_opname::pregen_relation(mdd, mxd, mdd);
if (ensf) {
start.note_time();
buildNextStateFunction(model, TRANS, ensf, 4);
start.note_time();
} else {
start.note_time();
buildNextStateFunction(model, TRANS, mxd, nsf, 4);
start.note_time();
}
printf("Next-state function construction took %.4e seconds\n",
start.get_last_interval() / 1000000.0);
printStats("MxD", mxd);
dd_edge reachable(mdd);
start.note_time();
printf("\nBuilding reachability set using saturation, relation");

fflush(stdout);
if (0==SATURATION_FORWARD) {
throw error(error::UNKNOWN_OPERATION);
}
sat = SATURATION_FORWARD->buildOperation(ensf);
if (0==sat) {
throw error(error::INVALID_OPERATION);
}
sat->compute(init_state, reachable);
start.note_time();
printf("\nReachability set construction took %.4e seconds\n",
start.get_last_interval() / 1000000.0);
printStats("MDD", mdd);
fflush(stdout);
double c;
apply(CARDINALITY, reachable, c);
operation::showAllComputeTables(meddlyout, 3);

printf("Approx. %g reachable states\n", c);

}
method ='i';
  std::cout<<"\n********************";
  std::cout<<"\n     implicit";
  std::cout<<"\n********************";
if('i' == method)
{

//CREATE FORESTS
forest* inmdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL,p);
forest* outmdd = inmdd;

//ADD INITIAL STATE
dd_edge first(inmdd);
dd_edge reachable(outmdd);
inmdd->createEdge(&initialState, 1, first);
outmdd->createEdge(&initialState, 1, reachable);


//CREATE RELATION
satimpl_opname::implicit_relation* T = new satimpl_opname::implicit_relation(inmdd,outmdd);

start.note_time();
buildImplicitRelation(model, TRANS, PLACES, BOUNDS, T);
printf("\nNext-state function construction took %.4e seconds\n",
start.get_last_interval() / 1000000.0);
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






