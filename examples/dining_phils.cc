
// $Id$

/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

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



/*
  State space generation for Dining Philosphers (N=2).

  The model has 2 philosophers and 2 forks.

  Each philosopher can be in state {I, W, L, R, E} where
  
  I:    idle philosopher
  WB:   philosopher is waiting for both forks
  HL:   philosopher has left fork
  HR:   philosopher has right fork
  E:    philosopher is eating

  Each fork can be in state {A, NA} where
  
  A:    fork is available
  NA:   fork is not available

  Philosphers can move from one state to another as:
  
  I -> WB
 
  The synchronization between philosopher 1 and the forks:

  WB1 ->  HR1
  A1  ->  NA1

  WB1 ->  HL1
  A2  ->  NA2

  HR1 ->  E1
  A2  ->  NA2

  HL1 ->  E1
  A1  ->  NA1

  E1  ->  I1
  NA1 ->  A1
  NA2 ->  A2

  The synchronization between philosopher 2 and the forks:

  WB2 ->  HR2
  A2  ->  NA2

  WB2 ->  HL2
  A1  ->  NA1

  HR2 ->  E2
  A1  ->  NA1

  HL2 ->  E2
  A2  ->  NA2

  E2  ->  I2
  NA1 ->  A1
  NA2 ->  A2

  Initially, all forks are in state "A" and all philosophers are in state "I".
  How many reachable states?
  Exceptions?
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../config.h"
#if HAVE_LIBGMP
  #include <gmp.h>
#endif

#include "meddly.h"
#include "meddly_expert.h"
#include "timer.h"


using namespace MEDDLY;

// #define NAME_VARIABLES
// #define SHOW_MXD

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

variable** initializeVariables(int nLevels)
{
  // set bounds for each variable
  // forks: 2 states, philosophers: 5 states
  variable** vars = (variable**) malloc((1+nLevels) * sizeof(void*));
  vars[0] = 0;
#ifdef NAME_VARIABLES
  char buffer[32];
  for (int i = nLevels; i; i -= 2) {
    // forks
    buffer[0] = 0;
    snprintf(buffer, 32, "fork%d", i/2);
    vars[i] = createVariable(2, strdup(buffer));
    // philosophers one level below corresponding forks
    buffer[0] = 0;
    snprintf(buffer, 32, "phil%d", i/2);
    vars[i-1] = createVariable(5, strdup(buffer));
  }
#else
  for (int i = nLevels; i; i -= 2) {
    // forks
    vars[i] = createVariable(2, 0);
    // philosophers one level below corresponding forks
    vars[i-1] = createVariable(5, 0);
  }
#endif
  return vars;
}


int* initializeInitialState(int nLevels)
{
  // initial state -- all levels at 0
  int *initialState = (int *) malloc((nLevels + 1) * sizeof(int));
  memset(initialState, 0, (nLevels + 1) * sizeof(int));
  return initialState;
}


void SetIntArray(int *p, int p_size, int c)
{
  for (int i = 0; i < p_size; ++i) p[i] = c;
}


dd_edge MakeSynchP_Forks(int philosopher, int nPhilosophers, forest* mxd)
{
  MEDDLY_CHECK_RANGE(0, philosopher, nPhilosophers);

  dd_edge nsf(mxd);
  dd_edge temp(mxd);

  int sz = nPhilosophers * 2 + 1;

  static int* from = NULL;
  static int* to = NULL;
  static int** addrFrom = NULL;
  static int** addrTo = NULL;

  if (from == NULL) {
    assert(to == NULL);
    assert(addrFrom == NULL);
    assert(addrTo == NULL);

    from = (int *) malloc(sz * sizeof(int));
    to = (int *) malloc(sz * sizeof(int));
    assert(from != NULL);
    assert(to != NULL);

    addrFrom = (int **) malloc(sizeof(int *));
    addrTo = (int **) malloc(sizeof(int *));
    assert(addrFrom != NULL);
    assert(addrTo != NULL);
    addrFrom[0] = from;
    addrTo[0] = to;
  }

  assert(from != NULL);
  assert(to != NULL);
  assert(addrFrom != NULL);
  assert(addrTo != NULL);

  from[0] = 0;
  to[0] = 0;

  /* Right Fork level is above Philosopher level */
  int ph = philosopher * 2 + 1;  // 0th philosopher at (0*2+1) = 1st level
  int rf = ph + 1;
  int lf = (philosopher > 0)? ph - 1: nPhilosophers * 2;

#if 0
  printf("%s: philosopher = %d, ph = %d, rf = %d, lf = %d\n",
      __func__, philosopher, ph, rf, lf);
#endif

  /* I(ph) -> WB(ph) */
  SetIntArray(from+1, sz-1, DONT_CARE);
  SetIntArray(to+1, sz-1, DONT_CHANGE);
  from[ph] = 0;
  to[ph] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* WB(ph) -> HR(ph), A(rf) -> NA(rf) */
  SetIntArray(from+1, sz-1, DONT_CARE);
  SetIntArray(to+1, sz-1, DONT_CHANGE);
  from[ph] = 1;
  to[ph] = 3;
  from[rf] = 0;
  to[rf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* WB(ph) -> HL(ph), A(lf) -> NA(lf) */
  SetIntArray(from+1, sz-1, DONT_CARE);
  SetIntArray(to+1, sz-1, DONT_CHANGE);
  from[ph] = 1;
  to[ph] = 2;
  from[lf] = 0;
  to[lf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* HR(ph) -> E(ph), A(lf) -> NA(lf) */
  SetIntArray(from+1, sz-1, DONT_CARE);
  SetIntArray(to+1, sz-1, DONT_CHANGE);
  from[ph] = 3;
  to[ph] = 4;
  from[lf] = 0;
  to[lf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* HL(ph) -> E(ph), A(rf) -> NA(rf) */
  SetIntArray(from+1, sz-1, DONT_CARE);
  SetIntArray(to+1, sz-1, DONT_CHANGE);
  from[ph] = 2;
  to[ph] = 4;
  from[rf] = 0;
  to[rf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* E(ph) -> I(ph), NA(rf) -> A(rf), NA(lf) -> A(lf) */
  SetIntArray(from+1, sz-1, DONT_CARE);
  SetIntArray(to+1, sz-1, DONT_CHANGE);
  from[ph] = 4;
  to[ph] = 0;
  from[rf] = 1;
  to[rf] = 0;
  from[lf] = 1;
  to[lf] = 0;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  return nsf;
}


void usage()
{
  printf("Usage: dining_phils_batch [-n<#phils>|-m<#MB>|-dfs|-p|-pgif]\n");
  printf("-n   : \
      number of philosophers\n");
  printf("-dfs : \
      use depth-first algorithm to compute reachable states\n");
  printf("-exact : \
      display the exact number of states\n");
  printf("-cs   : \
      set cache size (default is 262144).\n");
  printf("-nc   : \
      use cache with no-chaining instead of default.\n");
  printf("-pess: \
      use pessimistic node deletion (lower mem usage)\n");
  printf("\n");
}


// Test Index Set
void testIndexSet(const dd_edge& mdd, dd_edge& indexSet)
{
  apply(CONVERT_TO_INDEX_SET, mdd, indexSet);

#if 1
  indexSet.show(stdout, 3);
#else
  indexSet.show(stdout, 1);
#endif
}

int main(int argc, char *argv[])
{
  timer start;
  int nPhilosophers = 0; // number of philosophers
  bool pessimistic = false;
  bool exact = false;
  bool dfs = false;
  int cacheSize = 0;
  bool chaining = true;
  bool printReachableStates = false;
  if (argc > 1) {
    assert(argc > 1 && argc < 6);
    if (argc > 1) {
      assert(argc < 6);
      for (int i=1; i<argc; i++) {
        char *cmd = argv[i];
        if (strncmp(cmd, "-pess", 6) == 0) pessimistic = true;
        else if (strncmp(cmd, "-nc", 4) == 0) {
          chaining = false;
        }
        else if (strncmp(cmd, "-cs", 3) == 0) {
          cacheSize = strtol(&cmd[3], NULL, 10);
          if (cacheSize < 1) {
            usage();
            exit(1);
          }
        }
        else if (strncmp(cmd, "-dfs", 5) == 0) dfs = true;
        else if (strncmp(cmd, "-exact", 7) == 0) exact = true;
        else if (strncmp(cmd, "-print", 7) == 0) printReachableStates = true;
        else if (strncmp(cmd, "-n", 2) == 0) {
          nPhilosophers = strtol(&cmd[2], NULL, 10);
        } else {
          usage();
          exit(1);
        }
      }
    }
  }
  while (nPhilosophers < 2) {
    printf("Enter the number of philosophers (at least 2): ");
    scanf("%d", &nPhilosophers);
  }

  // Initialize MEDDLY

  MEDDLY::settings s;
  // TBD
  // s.doComputeTablesUseChaining = chaining;
  if (cacheSize > 0) {
    s.computeTable.maxSize = cacheSize;
  }
  MEDDLY::initialize(s);

  // Number of levels in domain (excluding terminals)
  int nLevels = nPhilosophers * 2;

  // Set up arrays bounds based on nPhilosophers
  variable** vars = initializeVariables(nLevels);

  printf("Initiailzing forests\n");

  // Create a domain and set up the state variables.
  domain *d = createDomain(vars, nLevels);
  assert(d != NULL);

  // Set up MDD options
  forest::policies pmdd(false);
  if (pessimistic)  pmdd.setPessimistic();
  else              pmdd.setOptimistic();

  // Create an MDD forest in this domain (to store states)
  forest* mdd =
    d->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL, pmdd);
  assert(mdd != NULL);

  // Set up MXD options
  forest::policies pmxd(true);
  if (pessimistic)  pmdd.setPessimistic();
  else              pmdd.setOptimistic();

  // Create a MXD forest in domain (to store transition diagrams)
  forest* mxd = 
    d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL, pmxd);
  assert(mxd != NULL);

  // Set up initial state array based on nPhilosophers
  int *initSt = initializeInitialState(nLevels);
  int *addrInitSt[1] = { initSt };
  dd_edge initialStates(mdd);
  mdd->createEdge(reinterpret_cast<int**>(addrInitSt), 1, initialStates);

  printf("Building next-state function for %d dining philosophers\n", 
          nPhilosophers);
  fflush(stdout);
  start.note_time();

  // Create a matrix diagram to represent the next-state function
  // Next-State function is computed by performing a union of next-state
  // functions that each represent how a philosopher's state can change.
  dd_edge nsf(mxd);
  for (int i = 0; i < nPhilosophers; i++) {
    nsf += MakeSynchP_Forks(i, nPhilosophers, mxd);
  }
  start.note_time();
  printf("Next-state function construction took %.4e seconds\n",
          start.get_last_interval()/1000000.0);

  // Show stats for nsf construction
  printStats("MxD", mxd);

#ifdef SHOW_MXD
  printf("Next-State Function:\n");
  nsf.show(stdout, 2);
#endif

  dd_edge reachableStates(initialStates);

  printf("Building reachability set using %s\n", 
    dfs
    ? "saturation"
    : "traditional iteration"
  );
  fflush(stdout);
  start.note_time();
  apply(
      dfs?
      REACHABLE_STATES_DFS:
      REACHABLE_STATES_BFS,
      reachableStates, nsf, reachableStates);
  start.note_time();
  printf("Reachability set construction took %.4e seconds\n",
          start.get_last_interval()/1000000.0);
  fflush(stdout);
  printf("#Nodes: %d\n", reachableStates.getNodeCount());
  printf("#Edges: %d\n", reachableStates.getEdgeCount());


  // Show stats for rs construction
  printStats("MDD", mdd);
  
  operation::showAllComputeTables(stdout, 1);

  double c;
  apply(CARDINALITY, reachableStates, c);
  printf("Approximately %e reachable states\n", c);
  fflush(stdout);

#if HAVE_LIBGMP
  if (exact) {
    mpz_t nrs;
    mpz_init(nrs);
    apply(CARDINALITY, reachableStates, nrs);
    printf("Exactly ");
    mpz_out_str(0, 10, nrs);
    printf(" reachable states\n");
    fflush(stdout);
    mpz_clear(nrs);
  }
#endif

  if (printReachableStates) {
    // Create a EV+MDD forest in this domain (to store index set)
    forest* evplusmdd =
      d->createForest(false, forest::INTEGER, forest::EVPLUS);
    assert(evplusmdd != NULL);

    // Test Convert MDD to Index Set EV+MDD
    dd_edge indexSet(evplusmdd);
    testIndexSet(reachableStates, indexSet);
    int* element = (int *) malloc((nLevels + 1) * sizeof(int));

    double cardinality = indexSet.getCardinality();
    for (int index = 0; index < int(cardinality); index++)
    {
      evplusmdd->getElement(indexSet, index, element);
      printf("Element at index %d: [ ", index);
      for (int i = nLevels; i > 0; i--)
      {
        printf("%d ", element[i]);
      }
      printf("]\n");
    }
  }

  if (false) {
    start.note_time();
    unsigned counter = 0;
    for (enumerator iter(reachableStates); 
        iter; ++iter, ++counter)
    {
      const int* element = iter.getAssignments();
      assert(element != 0);
#if 1
      printf("%d: [", counter);
      for (int i = 2*nPhilosophers; i > 0; --i)
      {
        printf("%d ", element[i]);
      }
      printf("]\n");
#endif
    }
    start.note_time();
    printf("Iterator traversal time (%0.4e elements): %0.4e seconds\n",
        double(counter), start.get_last_interval()/double(1000000.0));
  }

  // Cleanup
  MEDDLY::destroyDomain(d);
  MEDDLY::cleanup();

  printf("\n\nDONE\n");
  return 0;
}

