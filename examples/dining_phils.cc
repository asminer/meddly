
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

struct switches {
  bool pessimistic;
  bool exact;
  char method;
  // bool chaining;
  bool printReachableStates;

public:
  switches() {
    pessimistic = false;
    exact = false;
    method = 'b';
    // chaining = true;
    printReachableStates = false;
  }
};


class philsModel {
  public:
    philsModel(int nPhils, forest* mxd);
    ~philsModel();

    // event builders
    void Idle2WaitBoth(int phil, dd_edge &e);
    void WaitBoth2HaveRight(int phil, dd_edge &e);
    void WaitBoth2HaveLeft(int phil, dd_edge &e);
    void HaveRight2Eat(int phil, dd_edge &e);
    void HaveLeft2Eat(int phil, dd_edge &e);
    void Eat2Idle(int phil, dd_edge &e);

    // build everything for a given phil
    void eventsForPhil(int phil, dd_edge &e);

  private:
    inline void setMinterm(int* m, int c) {
      for (int i = 1; i<sz; i++) m[i] = c;
    }

    inline int philVar(int p) const {
      return 2*p +1;  // 0th philosopher at 1st variable
    }
    inline int leftVar(int p) const {
      // left fork for phil p
      return (p > 0) ? 2*p : nPhils*2;
    }
    inline int rightVar(int p) const {
      // right fork for phil p
      return 2*p+2;
    }

  private:
    int* from;
    int* to;
    int nPhils;
    int sz;
    forest* mxd;
};


philsModel::philsModel(int nP, forest* _mxd)
{
  nPhils = nP;
  mxd = _mxd;

  sz = nPhils * 2 + 1;

  from = new int[sz];
  to = new int[sz];

  from[0] = to[0] = 0;
}

philsModel::~philsModel()
{
  delete[] from;
  delete[] to;
}

void philsModel::Idle2WaitBoth(int phil, dd_edge &e)
{
  const int ph = philVar(phil);

  /* I(ph) -> WB(ph) */
  setMinterm(from, DONT_CARE);
  setMinterm(to, DONT_CHANGE);
  from[ph] = 0;
  to[ph] = 1;
  mxd->createEdge(&from, &to, 1, e);
}

void philsModel::WaitBoth2HaveRight(int phil, dd_edge &e)
{
  const int ph = philVar(phil);
  const int rf = rightVar(phil);

  /* WB(ph) -> HR(ph), A(rf) -> NA(rf) */
  setMinterm(from, DONT_CARE);
  setMinterm(to, DONT_CHANGE);
  from[ph] = 1;
  to[ph] = 3;
  from[rf] = 0;
  to[rf] = 1;
  mxd->createEdge(&from, &to, 1, e);
}

void philsModel::WaitBoth2HaveLeft(int phil, dd_edge &e)
{
  const int ph = philVar(phil);
  const int lf = leftVar(phil);

  /* WB(ph) -> HR(ph), A(lf) -> NA(lf) */
  setMinterm(from, DONT_CARE);
  setMinterm(to, DONT_CHANGE);
  from[ph] = 1;
  to[ph] = 3;
  from[lf] = 0;
  to[lf] = 1;
  mxd->createEdge(&from, &to, 1, e);
}

void philsModel::HaveRight2Eat(int phil, dd_edge &e)
{
  const int ph = philVar(phil);
  const int lf = leftVar(phil);

  /* HR(ph) -> E(ph), A(lf) -> NA(lf) */
  setMinterm(from, DONT_CARE);
  setMinterm(to, DONT_CHANGE);
  from[ph] = 3;
  to[ph] = 4;
  from[lf] = 0;
  to[lf] = 1;
  mxd->createEdge(&from, &to, 1, e);
}

void philsModel::HaveLeft2Eat(int phil, dd_edge &e)
{
  const int ph = philVar(phil);
  const int rf = rightVar(phil);

  /* HL(ph) -> E(ph), A(rf) -> NA(rf) */
  setMinterm(from, DONT_CARE);
  setMinterm(to, DONT_CHANGE);
  from[ph] = 2;
  to[ph] = 4;
  from[rf] = 0;
  to[rf] = 1;
  mxd->createEdge(&from, &to, 1, e);
}

void philsModel::Eat2Idle(int phil, dd_edge &e)
{
  const int ph = philVar(phil);
  const int rf = rightVar(phil);
  const int lf = leftVar(phil);

  /* E(ph) -> I(ph), NA(rf) -> A(rf), NA(lf) -> A(lf) */
  setMinterm(from, DONT_CARE);
  setMinterm(to, DONT_CHANGE);
  from[ph] = 4;
  to[ph] = 0;
  from[rf] = 1;
  to[rf] = 0;
  from[lf] = 1;
  to[lf] = 0;
  mxd->createEdge(&from, &to, 1, e);
}

void philsModel::eventsForPhil(int phil, dd_edge &e)
{
  dd_edge temp(mxd);
  Idle2WaitBoth(phil, e);
  WaitBoth2HaveRight(phil, temp);
  e += temp;
  WaitBoth2HaveLeft(phil, temp);
  e += temp;
  HaveRight2Eat(phil, temp);
  e += temp;
  HaveLeft2Eat(phil, temp);
  e += temp;
  Eat2Idle(phil, temp);
  e += temp;
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

void runWithOptions(int nPhilosophers, const switches &sw)
{
  timer start;

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
  if (sw.pessimistic) pmdd.setPessimistic();
  else                pmdd.setOptimistic();

  // Create an MDD forest in this domain (to store states)
  forest* mdd =
    d->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL, pmdd);
  assert(mdd != NULL);

  // Set up MXD options
  forest::policies pmxd(true);
  if (sw.pessimistic) pmdd.setPessimistic();
  else                pmdd.setOptimistic();

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
  philsModel model(nPhilosophers, mxd);
  dd_edge nsf(mxd);
  satpregen_opname::pregen_relation* ensf = 0;
  specialized_operation* sat = 0;

  if ('s' == sw.method) {
    ensf = new satpregen_opname::pregen_relation(mdd, mxd, mdd, 6*nPhilosophers);
  }
  if ('k' == sw.method) {
    ensf = new satpregen_opname::pregen_relation(mdd, mxd, mdd);
  }

  if (ensf) {
    dd_edge temp(mxd);
    for (int i = 0; i < nPhilosophers; i++) {
      model.Idle2WaitBoth(i, temp);
      ensf->addToRelation(temp);
      model.WaitBoth2HaveRight(i, temp);
      ensf->addToRelation(temp);
      model.WaitBoth2HaveLeft(i, temp);
      ensf->addToRelation(temp);
      model.HaveRight2Eat(i, temp);
      ensf->addToRelation(temp);
      model.HaveLeft2Eat(i, temp);
      ensf->addToRelation(temp);
      model.Eat2Idle(i, temp);
      ensf->addToRelation(temp);
    }
  } else {
    dd_edge phil(mxd);
    for (int i = 0; i < nPhilosophers; i++) {
      model.eventsForPhil(i, phil);
      nsf += phil;
    }
  }
  start.note_time();
  printf("Next-state function construction took %.4e seconds\n",
          start.get_last_interval()/1000000.0);

  //
  // Show stats for nsf construction
  //
  printStats("MxD", mxd);

#ifdef SHOW_MXD
  printf("Next-State Function:\n");
  nsf.show(stdout, 2);
#endif


  //
  // Build reachable states
  //
  dd_edge reachableStates(initialStates);
  start.note_time();

  switch (sw.method) {
    case 'b':
        printf("Building reachability set using traditional algorithm\n");
        fflush(stdout);
        apply(REACHABLE_STATES_BFS, initialStates, nsf, reachableStates);
        break;

    case 'm':
        printf("Building reachability set using saturation, monolithic relation\n");
        fflush(stdout);
        apply(REACHABLE_STATES_DFS, initialStates, nsf, reachableStates);
        break;

    case 'k':
    case 's':
        printf("Building reachability set using saturation, relation");
        if ('k'==sw.method) printf(" by levels\n");
        else                printf(" by events\n");
        fflush(stdout);
        if (0==SATURATION_FORWARD) {
          throw error(error::UNKNOWN_OPERATION);
        }
        sat = SATURATION_FORWARD->buildOperation(ensf);
        if (0==sat) {
          throw error(error::INVALID_OPERATION);
        }
        sat->compute(initialStates, reachableStates);
        break;

    default:
        printf("Error - unknown method\n");
        exit(2);
  };
  start.note_time();
  printf("Done\n");

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
  if (sw.exact) {
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

  if (sw.printReachableStates) {
    // Create a EV+MDD forest in this domain (to store index set)
    forest* evplusmdd =
      d->createForest(false, forest::INTEGER, forest::INDEX_SET);
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

}



int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }

  printf("\nUsage: %s [options]\n\n", name);

  printf("\t-n<#phils>: set number of philosophers\n\n");

  printf("\t-exact:     display the exact number of states\n");
  printf("\t-cs<cache>: set cache size (0 for library default)\n");
  printf("\t-pess:      use pessimistic node deletion (lower mem usage)\n\n");

  printf("\t-bfs:   use traditional iterations (default)\n\n");
  printf("\t-dfs:   use fastest saturation (currently, -msat)\n");
  printf("\t-esat:  use saturation by events\n");
  printf("\t-ksat:  use saturation by levels\n");
  printf("\t-msat:  use monolithic saturation\n");
  printf("\n");
  return 0;
}


int main(int argc, char *argv[])
{
  int nPhilosophers = 0; // number of philosophers
  switches sw;
  int cacheSize = 0;

  for (int i=1; i<argc; i++) {
    const char* cmd = argv[i];
    if (strcmp(cmd, "-pess") == 0) {
      sw.pessimistic = true;
      continue;
    }
    if (strncmp(cmd, "-cs", 3) == 0) {
      cacheSize = strtol(&cmd[3], NULL, 10);
      if (cacheSize < 1) {
        return 1+usage(argv[0]);
      }
      continue;
    }
    if (strcmp(cmd, "-exact") == 0) {
      sw.exact = true;
      continue;
    }
    if (strcmp(cmd, "-print") == 0) {
      sw.printReachableStates = true;
      continue;
    }
    if (strncmp(cmd, "-n", 2) == 0) {
      nPhilosophers = strtol(cmd+2, NULL, 10);
      if (nPhilosophers < 1) {
        return 1+usage(argv[0]);
      }
      continue;
    }

    if (strcmp(cmd, "-bfs") == 0) {
      sw.method = 'b';
      continue;
    }

    if (strcmp(cmd, "-dfs") == 0) {
      sw.method = 'm';
      continue;
    }

    if (strcmp(cmd, "-msat") == 0) {
      sw.method = 'm';
      continue;
    }
    if (strcmp(cmd, "-esat") == 0) {
      sw.method = 's';
      continue;
    }
    if (strcmp(cmd, "-ksat") == 0) {
      sw.method = 'k';
      continue;
    }

    return 1+usage(argv[0]);
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

  try {
    runWithOptions(nPhilosophers, sw);
    MEDDLY::cleanup();
    printf("\n\nDONE\n");
    return 0;
  }
  catch (error e) {
    printf("Caught MEDDLY error: %s\n", e.getName());
    return 1;
  }
}

