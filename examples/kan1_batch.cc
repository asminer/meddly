
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
    State space generation for a small, simple model (Kanban, N=1).
    The model is described in the article

    A. Miner and G. Ciardo.  "Efficient reachability set generation 
    and storage using decision diagrams", in Proc. ICATPN 1999, LNCS 1639,
    pp. 6-25.

    The model has 4 components.  Each component can be in state {W, M, B, G}
    where
    W:	the machine is waiting for a part
    M:	the machine is processing a part
    B:	a bad part was produced
    G:	a good part was produced

    The 4 machines are connected, so that output of machine 1 is used for
    machines 2 and 3, and output of machines 2 and 3 is used for machine 4.

    Machine 1 can change state locally as:
    
    W -> M
    M -> B
    B -> M
    M -> G

    Machines 2 and 3 can change state locally as:
    
    M -> B
    B -> M
    M -> G

    Machine 4 can change state locally as:
    
    M -> B
    B -> M
    M -> G
    G -> W

    The synchronization between machines 1,2,3 is:
    
    G1 -> W1
    W2 -> M2
    W3 -> M3

    The synchronization between machines 2,3,4 is:
    
    G2 -> W2
    G3 -> W3
    W4 -> M4

    Initially, all machines are in local state "W".
    There are 160 reachable states.  All combinations are possible, EXCEPT
    that machines 2 and 3 are either BOTH in state W, or BOTH NOT in state W.


    TODO in the library:

   	Write the function 
	AddMatrixElementAtLevel(dd_tempedge, level lh, int from, int to, bool term)

*/


/* Probably we need a makefile and a mechanism to specify
   the location of these header files...
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../include/meddly.h"
#include "../include/meddly_expert.h"


const int N = 5; // number of machines + 1
int sizes[N-1] = { 4, 4, 4, 4 };


dd_edge MakeLocalTransitions(int machine, forest* mxd)
{
  dd_edge nsf(mxd);

  int from[N];
  int to[N];
  int* addrFrom[1] = { from };
  int* addrTo[1] = { to };

  for (int i = 1; i < N; ++i)
  {
    from[i] = -2;
    to[i] = -2;
  }

  // W -> M
  if (1 == machine) {
    // adjust the vector
    from[machine] = 0; to[machine] = 1;

    // add it to nsf
    dd_edge temp(mxd);
    mxd->createEdge(reinterpret_cast<int**>(&addrFrom),
        reinterpret_cast<int**>(&addrTo), 1, temp);
    nsf += temp;

    // revert the vector
    from[machine] = -2; to[machine] = -2;
  }

  /* M -> B */
  {
    // adjust the vector
    from[machine] = 1; to[machine] = 2;

    // add it to nsf
    dd_edge temp(mxd);
    mxd->createEdge(reinterpret_cast<int**>(&addrFrom),
        reinterpret_cast<int**>(&addrTo), 1, temp);
    nsf += temp;

    // revert the vector
    from[machine] = -2; to[machine] = -2;
  }

  /* B -> M */
  {
    // adjust the vector
    from[machine] = 2; to[machine] = 1;

    // add it to nsf
    dd_edge temp(mxd);
    mxd->createEdge(reinterpret_cast<int**>(&addrFrom),
        reinterpret_cast<int**>(&addrTo), 1, temp);
    nsf += temp;

    // revert the vector
    from[machine] = -2; to[machine] = -2;
  }

  /* M -> G */
  {
    // adjust the vector
    from[machine] = 1; to[machine] = 3;

    // add it to nsf
    dd_edge temp(mxd);
    mxd->createEdge(reinterpret_cast<int**>(&addrFrom),
        reinterpret_cast<int**>(&addrTo), 1, temp);
    nsf += temp;

    // revert the vector
    from[machine] = -2; to[machine] = -2;
  }

  /* G -> W */
  if (4 == machine) {
    // adjust the vector
    from[machine] = 3; to[machine] = 0;

    // add it to nsf
    dd_edge temp(mxd);
    mxd->createEdge(reinterpret_cast<int**>(&addrFrom),
        reinterpret_cast<int**>(&addrTo), 1, temp);
    nsf += temp;

    // revert the vector
    from[machine] = -2; to[machine] = -2;
  }

  return nsf;
}


dd_edge MakeSynch1_23(forest* mxd)
{
  dd_edge nsf(mxd);

  int from[N] = {0, 3, 0, 0, -2};
  int to[N]   = {0, 0, 1, 1, -2};
  int* addrFrom[1] = { from };
  int* addrTo[1] = { to };

  mxd->createEdge(reinterpret_cast<int**>(&addrFrom),
      reinterpret_cast<int**>(&addrTo), 1, nsf);

  return nsf;
}

dd_edge MakeSynch23_4(forest* mxd)
{
  dd_edge nsf(mxd);

  int from[N] = {0, -2, 3, 3, 0};
  int to[N]   = {0, -2, 0, 0, 1};
  int* addrFrom[1] = { from };
  int* addrTo[1] = { to };

  mxd->createEdge(reinterpret_cast<int**>(&addrFrom),
      reinterpret_cast<int**>(&addrTo), 1, nsf);

  return nsf;
}

void getIndexSet(const dd_edge& mdd, dd_edge& indexSet)
{
  compute_manager* cm = MEDDLY_getComputeManager();
  assert(compute_manager::SUCCESS ==
      cm->apply(compute_manager::CONVERT_TO_INDEX_SET, mdd, indexSet));
}

int main(int argc, char* argv[])
{
  bool dfs = false;
  bool make_gifs = false;
  if (argc > 1) {
    assert(argc <= 3);
    char *cmd = NULL;
    for (int i = 1; i < argc; i++) {
      cmd = argv[i];
      if (strncmp(cmd, "-gif", 6) == 0) make_gifs = true;
      else if (strncmp(cmd, "-dfs", 5) == 0) dfs = true;
      else {
        printf("Usage: $ kan1_batch [-gif|-dfs]\n");
        printf("-gif : create gif representing the reachable states\n");
        printf("-dfs : use depth-first algorithm to compute reachable states");
        printf("\n\n");
        // return 1;
        exit(1);
      }
    }
  }
  if (argc > 1) dfs = true;
  int i;

  // Create a domain
  domain *d = MEDDLY_createDomain();
  assert(d != NULL);

  // Set up the state variables.
  // Use one per "machine", with 4 values each: {W = 0, M, B, G}
  assert(domain::SUCCESS == d->createVariablesBottomUp(sizes, N-1));

  // Create an MDD forest in this domain (to store states)
  forest* mdd = d->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL);
  assert(mdd != NULL);

  // Set up MDD options
  assert(forest::SUCCESS == mdd->setReductionRule(forest::FULLY_REDUCED));
  assert(forest::SUCCESS ==
      mdd->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
  assert(forest::SUCCESS ==
      mdd->setNodeDeletion(forest::OPTIMISTIC_DELETION));

  // Set up initial set of states
  int initst[N] = { 0, 0, 0, 0, 0 };
  int* addrInitst[1] = { initst };
  dd_edge initialStates(mdd);
  assert(forest::SUCCESS == mdd->createEdge(
        reinterpret_cast<int**>(&addrInitst), 1, initialStates));

  // Create a MXD forest in domain (to store transition diagrams)
  forest* mxd = d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL);
  assert(mxd != NULL);

  // Set up MDD options
  assert(forest::SUCCESS == mxd->setReductionRule(forest::IDENTITY_REDUCED));
  assert(forest::SUCCESS ==
      mxd->setNodeStorage(forest::FULL_OR_SPARSE_STORAGE));
  assert(forest::SUCCESS ==
      mxd->setNodeDeletion(forest::OPTIMISTIC_DELETION));

  // Initialize Next-State Function (nsf)
  dd_edge nsf(mxd);

  // Build and add local transitions to nsf
  for (i=1; i<=4; i++) {
    nsf += MakeLocalTransitions(i, mxd);
  }

  // Build and add synchronizing transitions to nsf
  nsf += MakeSynch1_23(mxd);
  nsf += MakeSynch23_4(mxd);

#if 0
  printf("Initial states:\n");
  initialStates.show(stdout, 2);

  printf("\nNext-State Function:\n");
  nsf.show(stdout, 2);
#endif

  dd_edge reachableStates(mdd);
  assert(compute_manager::SUCCESS ==
      MEDDLY_getComputeManager()->apply(compute_manager::REACHABLE_STATES_BFS,
      initialStates, nsf, reachableStates));
  // reachableStates.show(stdout, 3);

  // Create a EV+MDD forest in this domain (to store index set)
  forest* evplusmdd = d->createForest(false, forest::INTEGER, forest::EVPLUS);
  assert(evplusmdd != NULL);

  // Convert MDD to Index Set EV+MDD and print the states
  dd_edge indexSet(evplusmdd);
  getIndexSet(reachableStates, indexSet);
  int* element = (int *) malloc(N * sizeof(int));

  double cardinality = indexSet.getCardinality();
  for (int index = 0; index < int(cardinality); index++)
  {
    assert(forest::SUCCESS == evplusmdd->getElement(indexSet, index, element));
    printf("Element at index %d: [ ", index);
    for (int i = N - 1; i > 0; i--)
    {
      printf("%d ", element[i]);
    }
    printf("]\n");
  }

  // Test findFirstElement()
  if (true) {
    int** elements = &element;
    dd_edge rsCopy(reachableStates);
    cardinality = rsCopy.getCardinality();
    for (int index = 0; index < cardinality; index++)
    {
      memset(element, 0, N * sizeof(int));
      assert(forest::SUCCESS ==
          mdd->findFirstElement(rsCopy, element));
      printf("Element at index %d: [ ", index);
      for (int i = N - 1; i > 0; i--)
      {
        printf("%d ", element[i]);
      }
      printf("]\n");
      dd_edge temp(mdd);
      mdd->createEdge(elements, 1, temp);
      rsCopy -= temp;
    }
  }
  free(element);

  if (true) {
    unsigned counter = 0;
    for (dd_edge::const_iterator iter = reachableStates.begin(),
        endIter = reachableStates.end(); iter != endIter; ++iter, ++counter)
    {
      const int* element = iter.getAssignments();
      const int* curr = element + N - 1;
      const int* end = element - 1;
      printf("%d: [%d", counter, *curr--);
      while (curr != end) { printf(" %d", *curr--); }
      printf("]\n");
    }
    printf("Iterator traversal: %0.4e elements\n", double(counter));
    printf("Cardinality: %0.4e\n", reachableStates.getCardinality());
    counter = 0;
    for (dd_edge::const_iterator iter = nsf.begin(),
        endIter = nsf.end(); iter != endIter; ++iter, ++counter)
    {
      const int* element = iter.getAssignments();
      const int* pelement = iter.getPrimedAssignments();
      assert(element != 0 && pelement != 0);

      const int* curr = element + N - 1;
      const int* end = element - 1;
      printf("%d: [%d", counter, *curr--);
      while (curr != end) { printf(" %d", *curr--); }

      curr = pelement + N - 1;
      end = pelement - 1;
      printf("] --> [%d", *curr--);
      while (curr != end) { printf(" %d", *curr--); }
      printf("]\n");
    }
    printf("Iterator traversal: %0.4e elements\n", double(counter));
    printf("Cardinality: %0.4e\n", nsf.getCardinality());
  }

  // Cleanup
  delete d;

  printf("\nDone\n");
  return 0;
}

