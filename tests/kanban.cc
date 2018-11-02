
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

#include "../src/meddly.h"
#include "simple_model.h"

#define PROGRESS

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

long expected[] = { 
  1, 160, 4600, 58400, 454475, 2546432, 11261376, 
  41644800, 133865325, 384392800, 1005927208 
};

const int nstart = 1;
const int nstop = 10;

using namespace MEDDLY;

long buildReachset(int N, bool useSat)
{
  int sizes[16];

  for (int i=15; i>=0; i--) sizes[i] = N+1;
  domain* d = createDomainBottomUp(sizes, 16);

  // Build initial state
  int* initial = new int[17];
  for (int i=16; i; i--) initial[i] = 0;
  initial[1] = initial[5] = initial[9] = initial[13] = N;
  forest* mdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge init_state(mdd);
  mdd->createEdge(&initial, 1, init_state);
  delete[] initial;

#ifdef PROGRESS
  fputc('i', stdout);
  fflush(stdout);
#endif

  // Build next-state function
  forest* mxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  buildNextStateFunction(kanban, 16, mxd, nsf); 

#ifdef PROGRESS
  fputc('n', stdout);
  fflush(stdout);
#endif

  dd_edge reachable(mdd);
  if (useSat)
    apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
  else
    apply(REACHABLE_STATES_BFS, init_state, nsf, reachable);

#ifdef PROGRESS
  fputc('r', stdout);
  fflush(stdout);
#endif

  long c;
  apply(CARDINALITY, reachable, c);

#ifdef PROGRESS
  fputc('c', stdout);
  fflush(stdout);
#endif

  destroyDomain(d);
  
#ifdef PROGRESS
  fputc(' ', stdout);
  fflush(stdout);
#endif
  return c;
}

int main()
{
  MEDDLY::initialize();

  printf("Building Kanban reachability sets, using saturation\n");
  for (int n=nstart; n<=nstop; n++) {
    printf("N=%2d:  ", n);
    fflush(stdout);
    long c = buildReachset(n, true);
    printf("%12ld states\n", c);
    if (c != expected[n]) {
      printf("Wrong number of states!\n");
      return 1;
    }
  }

  printf("Building Kanban reachability sets, using traditional iteration\n");
  for (int n=nstart; n<=nstop; n++) {
    printf("N=%2d:  ", n);
    fflush(stdout);
    long c = buildReachset(n, false);
    printf("%12ld states\n", c);
    if (c != expected[n]) {
      printf("Wrong number of states!\n");
      return 1;
    }
  }

  MEDDLY::cleanup();
  printf("Done\n");
  return 0;
}

