
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

// #define DUMP_NSF
// #define DUMP_REACHABLE

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

int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  printf("\nUsage: %s nnnn [options]\n\n", name);
  printf("\tnnnn: number of parts\n\n");
  printf("\t-bfs: use traditional iterations\n\n");
  printf("\t-dfs: use saturation (monolithic)\n");
  printf("\t-esat: use saturation by events\n");
  printf("\t-ksat: use saturation by levels\n\n");
  printf("\t-exp: use explicit (very slow)\n");
  printf("\t--batch b: specify explicit batch size\n\n");
  return 1;
}

void printStats(const char* who, const forest* f)
{
  printf("%s stats:\n", who);
  const expert_forest* ef = (expert_forest*) f;
  ef->reportStats(stdout, "\t",
    expert_forest::HUMAN_READABLE_MEMORY  |
    expert_forest::BASIC_STATS | expert_forest::EXTRA_STATS |
    expert_forest::STORAGE_STATS | expert_forest::HOLE_MANAGER_STATS
  );
}



int main(int argc, const char** argv)
{
  int N = -1;
  char method = 'd';
  int batchsize = 256;

  for (int i=1; i<argc; i++) {
    if (strcmp("-bfs", argv[i])==0) {
      method = 'b';
      continue;
    }
    if (strcmp("-dfs", argv[i])==0) {
      method = 'd';
      continue;
    }
    if (strcmp("-esat", argv[i])==0) {
      method = 's';
      continue;
    }
    if (strcmp("-ksat", argv[i])==0) {
      method = 'k';
      continue;
    }
    if (strcmp("-exp", argv[i])==0) {
      method = 'e';
      continue;
    }
    if (strcmp("--batch", argv[i])==0) {
      i++;
      if (argv[i]) batchsize = atoi(argv[i]);
      continue;
    }
    N = atoi(argv[i]);
  }

  if (N<0) return usage(argv[0]);

  MEDDLY::initialize();

  printf("+-------------------------------------------+\n");
  printf("|   Initializing Kanban model for N = %-4d  |\n", N);
  printf("+-------------------------------------------+\n");
  fflush(stdout);

  // Initialize domain
  int* sizes = new int[16];
  for (int i=15; i>=0; i--) sizes[i] = N+1;
  domain* d = createDomainBottomUp(sizes, 16);
  delete[] sizes;

  // Build initial state
  int* initial = new int[17];
  for (int i=16; i; i--) initial[i] = 0;
  initial[1] = initial[5] = initial[9] = initial[13] = N;
  forest* mdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge init_state(mdd);
  mdd->createEdge(&initial, 1, init_state);
  delete[] initial;

  // Build next-state function
  forest* mxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  satpregen_opname::pregen_relation ensf(mxd, 16);
  dd_edge* knsf = 0;

  if ('e' != method) {

    switch (method) {
      case 's':
        try {
          buildByEventsNSF(kanban, 16, ensf, 4);
          ensf.finalize();
        }
        catch (MEDDLY::error e) {
          printf("Caught error %s\n", e.getName());
        }
        break;

      case 'k':
        /*
        knsf = new dd_edge[17];
        for (int i=0; i<=16; i++) knsf[i].setForest(mxd);
        buildByLevelsNSF(kanban, 16, mxd, knsf, 4);
        */
        break;

      case 'b':
      case 'd':
        buildNextStateFunction(kanban, 16, mxd, nsf, 4);
#ifdef DUMP_NSF
        printf("Next-state function:\n");
        nsf.show(stdout, 2);
#endif
        break;
    }
    printStats("MxD", mxd);
  }

  dd_edge reachable(mdd);
  switch (method) {
    case 'b':
        printf("Building reachability set using traditional algorithm\n");
        fflush(stdout);
        apply(REACHABLE_STATES_BFS, init_state, nsf, reachable);
        break;

    case 'd':
        printf("Building reachability set using saturation\n");
        fflush(stdout);
        apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
        break;

    case 'e': 
        printf("Building reachability set using explicit search\n");
        printf("Using batch size: %d\n", batchsize);
        fflush(stdout);
        explicitReachset(kanban, 16, mdd, init_state, reachable, batchsize);
        break;

    case 'k':
        // for now - debugging - check the relation
        for (int k=16; k; k--) {
          printf("Relation at level %d:  %ld\n", k, long(knsf[k].getNode()));
        }
        break;

    case 's':
        // for now - debugging - check the relation
        for (int k=16; k; k--) {
          int len = ensf.lengthForLevel(k);
          if (0==len) continue;
          printf("Events at level %d:\n\t", k);
          node_handle* List = ensf.arrayForLevel(k);
          for (int i=0; i<len; i++)
            printf("%d ", List[i]);
          printf("\n");
        }
        break;

    default:
        printf("Error - unknown method\n");
        exit(2);
  };
  printf("Done\n");
  fflush(stdout);

#ifdef DUMP_REACHABLE
  printf("Reachable states:\n");
  reachable.show(stdout, 2);
#endif

  printStats("MDD", mdd);
  fflush(stdout);

  double c;
  apply(CARDINALITY, reachable, c);
  operation::showAllComputeTables(stdout, 2);

  printf("Approx. %g reachable states\n", c);
  
  // cleanup
  MEDDLY::cleanup();
  return 0;
}

