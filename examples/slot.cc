
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

using namespace MEDDLY;

int usage(const char* name)
{
  printf("\nUsage: %s nnnn (-bfs) (-dfs)\n\n", name);
  printf("\tnnnn: number of network nodes\n");
  printf("\t-bfs: use traditional iterations\n");
  printf("\t-dfs: use saturation\n\n");
  return 1;
}

inline char* newEvent(int N)
{
  char* ev = new char[N*8+1];
  for (int i=N*8; i; i--) ev[i] = '.';
  ev[0] = '_';
  return ev;
}

char* Get(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  tloc[6] = '-';
  tloc[8] = '-';
  tloc[3] = '+';
  tloc[5] = '+';
  return t;
}

char* Free(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  char* t_rt = i ? t+(8*i-8) : t+(8*N-8);
  tloc[5] = '-';
  tloc[6] = '+';
  t_rt[3] = '-';
  t_rt[2] = '+';
  return t;
}

char* Put(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  tloc[3] = '+';
  tloc[7] = '+';
  tloc[4] = '-';
  tloc[6] = '-';
  return t;
}

char* Used(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  char* t_rt = i ? t+(8*i-8) : t+(8*N-8);
  tloc[7] = '-';
  tloc[6] = '+';
  t_rt[3] = '-';
  t_rt[1] = '+';
  return t;
}

char* Other(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  tloc[1] = '-';
  tloc[4] = '+';
  return t;
}

char* Owner(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  tloc[1] = '-';
  tloc[2] = '+';
  return t;
}

char* Write(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  tloc[2] = '-';
  tloc[4] = '+';
  return t;
}

char* Go(int i, int N)
{
  char* t = newEvent(N);
  char* tloc = t+8*i;
  tloc[2] = '-';
  tloc[8] = '+';
  return t;
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

  printf("Building model for %d-node network\n", N);

  char** events = new char*[8*N];
  char** fill = events;
  for (int i=0; i<N; i++) {
    fill[0] = Other(i, N);
    fill[1] = Owner(i, N);
    fill[2] = Write(i, N);
    fill[3] = Go(i, N);
    fill[4] = Get(i, N);
    fill[5] = Put(i, N);
    fill[6] = Used(i, N);
    fill[7] = Free(i, N);
    fill += 8;
  }

  MEDDLY::initialize();

  // Initialize domain
  int* sizes = new int[N*8];
  for (int i=N*8-1; i>=0; i--) sizes[i] = 2;
  domain* d = createDomainBottomUp(sizes, N*8);

  // Build initial state
  int* initial = new int[1+N*8];
  for (int i=N*8; i; i--) initial[i] = 0;
  int* initLocal = initial;
  for (int i=0; i<N; i++) {
    initLocal[3] = initLocal[5] = 1;
    initLocal += 8;
  }
  forest* mdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge init_state(mdd);
  mdd->createEdge(&initial, 1, init_state);

  // Build next-state function
  forest* mxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  buildNextStateFunction(events, 8*N, mxd, nsf, 2);

  printf("Building reachable states\n");
  fflush(stdout);

  dd_edge reachable(mdd);
  if (useSaturation)
    apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
  else
    apply(REACHABLE_STATES_BFS, init_state, nsf, reachable);

#ifdef SHOW_STATES
  int count = 0;
  for (dd_edge::const_iterator i = reachable.begin(); i; ++i, ++count) {
    const int* element = i.getAssignments();
    printf("State %4d: [%d", count, element[1]);
    for (int j=2; j<=8*N; j++) {
      printf(", %d", element[j]);
    } // for j
    printf("]\n");
  }  // for i
#endif


  printf("Done\n");
  double c;
  apply(CARDINALITY, reachable, c);
  printf("Approx. %g reachable states\n", c);
  
  // cleanup
  operation::showAllComputeTables(stdout, 1);
  MEDDLY::cleanup();
  return 0;
}

