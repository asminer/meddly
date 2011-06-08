
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

#include "meddly.h"
#include "meddly_expert.h"
#include "timer.h"

using namespace MEDDLY;

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

int* initializeLevelBounds(int nLevels)
{
  // set bounds (node size) for each level
  // forks: 2 states, philosophers: 5 states
  int *bounds = (int *) malloc(nLevels * sizeof(int));
  for (int i = nLevels - 1; i >= 0; i -= 2) {
    // forks
    bounds[i] = 2;
    // philosophers one level below corresponding forks
    bounds[i - 1] = 5;
  }
  return bounds;
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
  CHECK_RANGE(0, philosopher, nPhilosophers);

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
  SetIntArray(from+1, sz-1, -2);
  SetIntArray(to+1, sz-1, -2);
  from[ph] = 0;
  to[ph] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* WB(ph) -> HR(ph), A(rf) -> NA(rf) */
  SetIntArray(from+1, sz-1, -2);
  SetIntArray(to+1, sz-1, -2);
  from[ph] = 1;
  to[ph] = 3;
  from[rf] = 0;
  to[rf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* WB(ph) -> HL(ph), A(lf) -> NA(lf) */
  SetIntArray(from+1, sz-1, -2);
  SetIntArray(to+1, sz-1, -2);
  from[ph] = 1;
  to[ph] = 2;
  from[lf] = 0;
  to[lf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* HR(ph) -> E(ph), A(lf) -> NA(lf) */
  SetIntArray(from+1, sz-1, -2);
  SetIntArray(to+1, sz-1, -2);
  from[ph] = 3;
  to[ph] = 4;
  from[lf] = 0;
  to[lf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* HL(ph) -> E(ph), A(rf) -> NA(rf) */
  SetIntArray(from+1, sz-1, -2);
  SetIntArray(to+1, sz-1, -2);
  from[ph] = 2;
  to[ph] = 4;
  from[rf] = 0;
  to[rf] = 1;
  temp.clear();
  mxd->createEdge(reinterpret_cast<int**>(addrFrom),
      reinterpret_cast<int**>(addrTo), 1, temp);
  nsf += temp;

  /* E(ph) -> I(ph), NA(rf) -> A(rf), NA(lf) -> A(lf) */
  SetIntArray(from+1, sz-1, -2);
  SetIntArray(to+1, sz-1, -2);
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
  compute_manager* cm = getComputeManager();
  cm->apply(compute_manager::CONVERT_TO_INDEX_SET, mdd, indexSet);

#if 1
  indexSet.show(stdout, 3);
#else
  indexSet.show(stdout, 1);
#endif
}


// Test Pre-Image
dd_edge testPreImage(const dd_edge& mdd, const dd_edge& mxd)
{
  dd_edge preImage(mdd.getForest());
  compute_manager* cm = getComputeManager();
  cm->apply(compute_manager::PRE_IMAGE, mdd, mxd, preImage);

#if 1
  preImage.show(stdout, 3);
#else
  preImage.show(stdout, 1);
#endif

  return preImage;
}


// Test Post-Image
dd_edge testPostImage(const dd_edge& mdd, const dd_edge& mxd)
{
  dd_edge postImage(mdd.getForest());
  compute_manager* cm = getComputeManager();
  cm->apply(compute_manager::POST_IMAGE, mdd, mxd, postImage);

#if 1
  postImage.show(stdout, 3);
#else
  postImage.show(stdout, 1);
#endif

  return postImage;
}


// Test SubMatrix
dd_edge testSubMatrix(int* bounds, int nLevels, const dd_edge& nsf)
{
  bool** vlist = (bool **) malloc((nLevels + 1) * sizeof(bool *));
  bool** vplist = (bool **) malloc((nLevels + 1) * sizeof(bool *));
  for (int i = 0; i <= nLevels; i++)
  {
    int levelSize = (i == 0)? 2: bounds[i - 1];
    unsigned arraySize = levelSize * sizeof(bool);
    vlist[i] = (bool *) malloc(arraySize);
    vplist[i] = (bool *) malloc(arraySize);
    memset(vlist[i], 0, arraySize);
    memset(vplist[i], 0, arraySize);
  }
  for (int i = 1; i <= nLevels; i++)
  {
    int levelSize = bounds[i - 1];
    for (int j = 0; j < levelSize - 1; j++)
    {
      vlist[i][j] = true;
      vplist[i][j] = true;
    }
  }

#if 1
  nsf.show(stdout, 3);
#else
  nsf.show(stdout, 1);
#endif

  forest* mxd = nsf.getForest();
  dd_edge subMatrix(mxd);
  mxd->createSubMatrix(vlist, vplist, nsf, subMatrix);

#if 1
  subMatrix.show(stdout, 3);
#else
  subMatrix.show(stdout, 1);
#endif

  return subMatrix;
}


int main(int argc, char *argv[])
{
  timer start;
  int nPhilosophers = 0; // number of philosophers
  bool pessimistic = false;
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

  // Number of levels in domain (excluding terminals)
  int nLevels = nPhilosophers * 2;

  // Set up arrays bounds based on nPhilosophers
  int *bounds = initializeLevelBounds(nLevels);

  expert_compute_manager* ecm = 
    static_cast<expert_compute_manager*>(getComputeManager());
  assert(ecm != 0);

  if (cacheSize > 0) {
    ecm->setHashTablePolicy(chaining, cacheSize);
  } else {
    ecm->setHashTablePolicy(chaining);
  }

  printf("Initiailzing forests\n");

  // Create a domain
  domain *d = createDomain();
  assert(d != NULL);

  // Set up the state variables.
  // Use one per "machine", with 4 values each: {W = 0, M, B, G}
  d->createVariablesBottomUp(bounds, nLevels);

  // Create an MDD forest in this domain (to store states)
  forest* mdd =
    d->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL);
  assert(mdd != NULL);

  // Set up MDD options
#if 1
  mdd->setReductionRule(forest::FULLY_REDUCED);
  mdd->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE);
#else
  mdd->setReductionRule(forest::QUASI_REDUCED);
  mdd->setNodeStorage(forest::FULL_STORAGE);
#endif
  mdd->setNodeDeletion(
    pessimistic
    ? forest::PESSIMISTIC_DELETION
    : forest::OPTIMISTIC_DELETION
  );

  // Create a MXD forest in domain (to store transition diagrams)
  forest* mxd = d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL);
  assert(mxd != NULL);

  // Set up MXD options
  mxd->setReductionRule(forest::IDENTITY_REDUCED);
#if 1
  mxd->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE);
#else
  mxd->setNodeStorage(forest::FULL_STORAGE);
#endif
  mxd->setNodeDeletion(
    pessimistic
      ? forest::PESSIMISTIC_DELETION
      : forest::OPTIMISTIC_DELETION
  );

  // Set up initial state array based on nPhilosophers
  int *initSt = initializeInitialState(nLevels);
  int *addrInitSt[1] = { initSt };
  dd_edge initialStates(mdd);
  mdd->createEdge(reinterpret_cast<int**>(addrInitSt), 1, initialStates);

  printf("Building next-state function for %d dining philosophers\n", 
          nPhilosophers);
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

  printf("MxD stats:\n");
  printf("\t%ld current nodes\n", mxd->getCurrentNumNodes());
  printf("\t%ld peak nodes\n", mxd->getPeakNumNodes());
  printf("\t");
  printmem(mxd->getCurrentMemoryUsed());
  printf(" current memory\n\t");
  printmem(mxd->getPeakMemoryUsed());
  printf(" peak memory\n");
  
  fflush(stdout);

#if 0
  printf("Initial states:\n");
  initialStates.show(stdout, 2);

  printf("\nNext-State Function:\n");
  nsf.show(stdout, 2);
#endif

  dd_edge reachableStates(initialStates);


#if 1

  printf("Building reachability set using %s\n", 
    dfs
    ? "saturation"
    : "traditional iteration"
  );
  start.note_time();
  ecm->apply(
      dfs?
      compute_manager::REACHABLE_STATES_DFS:
      compute_manager::REACHABLE_STATES_BFS,
      reachableStates, nsf, reachableStates);
  start.note_time();
  printf("Reachability set construction took %.4e seconds\n",
          start.get_last_interval()/1000000.0);
  fflush(stdout);
  printf("#Nodes: %d\n", reachableStates.getNodeCount());
  printf("#Edges: %d\n", reachableStates.getEdgeCount());

#else

  printf("Building reachability set\n");
  start.note_time();

  // set up aliases
  // traditional reachability analysis:
  // reachableStates = initialStates
  // postImage = getPostImage(reachableStates, nextStateFunction)
  // while (postImage - reachableStates != empty)
  // do
  //   reachableStates += postImage
  //   postImage = getPostImage(postImage)
  //

  const int nOperands = 3;
  forest* forests[nOperands] = {mdd, mdd, mdd};
  op_info* unionOp =
    ecm->getOpInfo(compute_manager::UNION, forests, nOperands);
  assert(unionOp != 0);
  forests[1] = mxd;
  op_info* postImageOp =
    ecm->getOpInfo(compute_manager::POST_IMAGE, forests, nOperands);
  assert(postImageOp != 0);

  dd_edge prevReachableStates(mdd);
  dd_edge postImage(mdd);

  ecm->apply(postImageOp, reachableStates, nsf, postImage);
  prevReachableStates = reachableStates;
  ecm->apply(unionOp, reachableStates, postImage, reachableStates);

  while(prevReachableStates != reachableStates)
  {
    ecm->apply(postImageOp, postImage, nsf, postImage);
    prevReachableStates = reachableStates;
    ecm->apply(unionOp, reachableStates, postImage, reachableStates);
  }

  start.note_time();
  printf("Time interval: %.4e seconds\n",
      start.get_last_interval()/1000000.0);

#endif

  // Show stats for rs construction

  printf("MDD stats:\n");
  printf("\t%ld current nodes\n", mdd->getCurrentNumNodes());
  printf("\t%ld peak nodes\n", mdd->getPeakNumNodes());
  printf("\t");
  printmem(mdd->getCurrentMemoryUsed());
  printf(" current memory\n\t");
  printmem(mdd->getPeakMemoryUsed());
  printf(" peak memory\n");
  

#if 0   // determine verbosity of output
  reachableStates.show(stdout, 2);
#else
  reachableStates.show(stdout, 1);
#endif

#if 0
  printf("\nCompute table info:\n");
  ecm->showComputeTable(stdout);
#endif

#if 0

  printf("\nMDD forest info:\n");
  mdd->showInfo(stdout);
  printf("\nMXD forest info:\n");
  mxd->showInfo(stdout);
  printf("\nCompute table info:\n");
  ecm->showComputeTable(stdout);

  start.note_time();


#if 1
  printf("\nClearing compute table... ");
  ecm->clearComputeTable();
  printf("done.\n");
#else
  printf("\nForcing garbage collection on MDD forest... ");
  mdd->garbageCollect();
  mxd->garbageCollect();
  printf("done.\n");
#endif

  start.note_time();
  printf("Time interval: %.4e seconds\n",
      start.get_last_interval()/1000000.0);

  printf("\nMDD forest info:\n");
  mdd->showInfo(stdout);
  printf("\nMXD forest info:\n");
  mxd->showInfo(stdout);
  printf("\nCompute table info:\n");
  ecm->showComputeTable(stdout);

#endif

#if 0
  // Test Pre-Image
  dd_edge preImage = testPreImage(initialStates, nsf);

  // Test Post-Image
  dd_edge postImage = testPostImage(initialStates, nsf);

  // Test SubMatrix
  dd_edge subMatrix = testSubMatrix(bounds, nLevels, nsf);
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
    for (dd_edge::const_iterator iter = reachableStates.begin();
        iter; ++iter, ++counter)
    {
      const int* element = iter.getAssignments();
      assert(element != 0);
#if 0
      printf("%d: [", counter);
      for (int i = N - 1; i > 0; --i)
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

  if (false) {
    // Test findFirstElement()
    int* element = (int *) malloc((nLevels + 1) * sizeof(int));
    int** elements = &element;
    dd_edge rsCopy(reachableStates);
    start.note_time();
    double cardinality = rsCopy.getCardinality();
    for (int index = 0; index < cardinality; index++)
    {
      // memset(element, 0, (nLevels + 1) * sizeof(int));
      mdd->findFirstElement(rsCopy, element);
#if 0
      printf("Element at index %d: [ ", index);
      for (int i = N - 1; i > 0; i--)
      {
        printf("%d ", element[i]);
      }
      printf("]\n");
#endif
      dd_edge temp(mdd);
      mdd->createEdge(elements, 1, temp);
      rsCopy -= temp;
    }
    start.note_time();
    printf("findFirstElement traversal time (%0.4e elements): %0.4e seconds\n",
        cardinality, start.get_last_interval()/double(1000000.0));
    free(element);
  }

  if (false) {
    // Create a EV+MDD forest in this domain (to store index set)
    forest* evplusmdd = d->createForest(false, forest::INTEGER, forest::EVPLUS);
    assert(evplusmdd != NULL);

    // Convert MDD to Index Set EV+MDD and print the states
    int* element = (int *) malloc((nLevels + 1) * sizeof(int));
    dd_edge indexSet(evplusmdd);
    start.note_time();
    compute_manager* cm = getComputeManager();
    cm->apply(compute_manager::CONVERT_TO_INDEX_SET, reachableStates, indexSet);
    double cardinality = indexSet.getCardinality();
    for (int index = 0, card = int(cardinality); index < card; index++)
    {
      evplusmdd->getElement(indexSet, index, element);
#if 0
      printf("Element at index %d: [ ", index);
      for (int i = N - 1; i > 0; i--)
      {
        printf("%d ", element[i]);
      }
      printf("]\n");
#endif
    }
    start.note_time();
    printf("Index Set traversal time (%0.4e elements): %0.4e seconds\n",
        cardinality, start.get_last_interval()/double(1000000.0));
    free(element);
  }

  // Cleanup
  delete d;

  printf("\n\nDONE\n");
  return 0;
}

