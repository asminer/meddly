
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "meddly.h"
#include "meddly_expert.h"
#include "simple_model.h"

const char* kanban[] = {
  "X-+..............",  // Tin1
  "X.-+.............",  // Tr1
  "X.+-.............",  // Tb1
  "X.-.+............",  // Tg1
  "X.....-+.........",  // Tr2
  "X.....+-.........",  // Tb2
  "X.....-.+........",  // Tg2
  "X+..--+..-+......",  // Ts1_23
  "X.........-+.....",  // Tr3
  "X.........+-.....",  // Tb3
  "X.........-.+....",  // Tg3
  "X....+..-+..--+..",  // Ts23_4
  "X.............-+.",  // Tr4
  "X.............+-.",  // Tb4
  "X............+..-",  // Tout4
  "X.............-.+"   // Tg4
};

using namespace MEDDLY;

int usage(const char* name)
{
  printf("\nUsage: %s nnnn (-bfs) (-dfs)\n\n", name);
  printf("\tnnnn: number of parts\n");
  printf("\t-bfs: use traditional iterations\n");
  printf("\t-dfs: use saturation\n\n");
  return 1;
}

void printmem(long m)
{
  if (m<1024) {
    printf("%ld bytes", m);
    return;
  }
  double approx = m;
  approx /= 1024;
  if (approx < 1024) {
    printf("%3.2lf Kbytes", approx);
    return;
  }
  approx /= 1024;
  if (approx < 1024) {
    printf("%3.2lf Mbytes", approx);
    return;
  }
  approx /= 1024;
  if (approx < 1024) {
    printf("%3.2lf Gbytes", approx);
    return;
  }
  approx /= 1024;
  printf("%3.2lf Tbytes", approx);
}


int main(int argc, const char** argv)
{
  int N = -1;
  bool useSaturation = true;

  for (int i=1; i<argc; i++) {
    if (strcmp("-bfs", argv[i])==0) {
      useSaturation = false;
      continue;
    }
    if (strcmp("-dfs", argv[i])==0) {
      useSaturation = true;
      continue;
    }
    N = atoi(argv[i]);
  }

  if (N<0) return usage(argv[0]);

  MEDDLY::initialize();

  // Initialize domain
  int* sizes = new int[16];
  for (int i=15; i>=0; i--) sizes[i] = N+1;
  domain* d = createDomainBottomUp(sizes, 16);

  // Build initial state
  int* initial = new int[17];
  for (int i=16; i; i--) initial[i] = 0;
  initial[1] = initial[5] = initial[9] = initial[13] = N;
  forest* mdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge init_state(mdd);
  mdd->createEdge(&initial, 1, init_state);

  // Build next-state function
  forest* mxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  buildNextStateFunction(kanban, 16, mxd, nsf, 4);

  printf("MxD stats:\n");
  printf("\t%ld current nodes\n", mxd->getCurrentNumNodes());
  printf("\t%ld peak nodes\n", mxd->getPeakNumNodes());
  printf("\t");
  printmem(mxd->getCurrentMemoryUsed());
  printf(" current memory used\n\t");
  printmem(mxd->getPeakMemoryUsed());
  printf(" peak memory used\n\t");
  printmem(mxd->getCurrentMemoryAllocated());
  printf(" current memory allocated\n\t");
  printmem(mxd->getPeakMemoryAllocated());
  printf(" peak memory allocated\n");
  
  dd_edge reachable(mdd);
  if (useSaturation) {
    printf("Building reachability set using saturation\n");
    fflush(stdout);
    apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
  } else {
    printf("Building reachability set using traditional algorithm\n");
    fflush(stdout);
    apply(REACHABLE_STATES_BFS, init_state, nsf, reachable);
  }
  printf("Done\n");

  printf("MDD stats:\n");
  printf("\t%ld current nodes\n", mdd->getCurrentNumNodes());
  printf("\t%ld peak nodes\n", mdd->getPeakNumNodes());
  printf("\t");
  printmem(mdd->getCurrentMemoryUsed());
  printf(" current memory used\n\t");
  printmem(mdd->getPeakMemoryUsed());
  printf(" peak memory used\n\t");
  printmem(mdd->getCurrentMemoryAllocated());
  printf(" current memory allocated\n\t");
  printmem(mdd->getPeakMemoryAllocated());
  printf(" peak memory allocated\n");
  fflush(stdout);

  double c;
  apply(CARDINALITY, reachable, c);
  operation::showAllComputeTables(stdout, 1);

  printf("Approx. %g reachable states\n", c);
  
  // cleanup
  MEDDLY::cleanup();
  return 0;
}

