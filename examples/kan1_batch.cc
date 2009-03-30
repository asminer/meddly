
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

#include "../include/mddlib.h"


const int N = 5; // number of machines + 1
int sizes[N-1] = { 4, 4, 4, 4 };


dd_edge MakeLocalTransitions(int machine, forest* mxd)
{
  dd_edge nsf(mxd);
  dd_edge temp(mxd);

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
    temp.clear();
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
    temp.clear();
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
    temp.clear();
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
    temp.clear();
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
    temp.clear();
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
  domain *d = MDDLIB_createDomain();
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

  compute_manager* cm = MDDLIB_getComputeManager();
  dd_edge reachableStates(initialStates);
  dd_edge prevReachableStates(mdd);
  dd_edge postImage(mdd);

  while(prevReachableStates != reachableStates)
  {
    prevReachableStates = reachableStates;
    // printf("\nPost-Image (mdd:%d, mxd:%d): ",
    //    reachableStates.getNode(), nsf.getNode());
    cm->apply(compute_manager::POST_IMAGE, reachableStates, nsf, postImage);
    // printf("%d\n", postImage.getNode());
    // postImage.show(stdout, 2);
    // printf("\nUnion (mdd:%d, mdd:%d): ",
    //    reachableStates.getNode(), postImage.getNode());
    cm->apply(compute_manager::UNION, reachableStates, postImage,
        reachableStates);
    // printf("%d\n", reachableStates.getNode());
  }
  reachableStates.show(stdout, 3);

#if 0
  // Image operations
  printf("\nInitial State: ");
  ShowDDEdge(stdout, initial);
  dd_edge *curr = NULL;
  if (dfs) {
    vector<dd_edge *> *xd = NULL;
    assert(SUCCESS == SplitMxd(nsf, xd));
    assert(SUCCESS == Saturate(initial, xd, curr));
  } else {
    Saturate(initial, nsf, curr);
  }

  printf("\nStates reachable from Initial State: ");
  ShowDDEdge(stdout, curr);

  // Number of reachable states
  printf("\nCounting reachable nodes... ");
  fflush(stdout);
  double card = Cardinality(curr);
  printf("done\n");
  printf("# of reachable states: %1.3e\n", card);
  fflush(stdout);

  if (make_gifs) {
    const char gif[] = "gif";
    const char filename[] = "kan1_batch_image";
    CreateDDEdgePic(filename, gif, curr);
    printf("Wrote reachable states to %s.%s\n", filename, gif);

    const char filename2[] = "kan1_batch_nsf";
    CreateDDEdgePic(filename2, gif, nsf);
    printf("Wrote reachable states to %s.%s\n", filename2, gif);
  }

  /*
  // Check if the reachable set is correct
  bool found = false;
  error_code ec = SUCCESS;
  int v[N];

  v[variables[0]]=0;
  v[variables[4]]=1; v[variables[3]]=0; v[variables[2]]=0; v[variables[1]]=1;
  printf("\nLooking for 1,0,0,1... ");
  ec = EvaluateVectorBool(curr, v, N, found);
  DCASSERT(ec == SUCCESS);
  if (found) { printf("found\n"); } else { printf("not found\n"); }

  v[variables[4]]=1; v[variables[3]]=1; v[variables[2]]=0; v[variables[1]]=1;
  printf("\nLooking for 1,1,0,1... ");
  ec = EvaluateVectorBool(curr, v, N, found);
  DCASSERT(ec == SUCCESS);
  if (found) { printf("found\n"); } else { printf("not found\n"); }

  v[variables[4]]=1; v[variables[3]]=1; v[variables[2]]=1; v[variables[1]]=1;
  printf("\nLooking for 1,1,1,1... ");
  ec = EvaluateVectorBool(curr, v, N, found);
  DCASSERT(ec == SUCCESS);
  if (found) { printf("found\n"); } else { printf("not found\n"); }
  */

  DestroyForest(states);
  if (INVALID_FOREST != states) {
    fprintf(stderr, "Couldn't destroy forest of states\n");
    return 1;
  } else {
    fprintf(stderr, "Destroyed forest of states\n");
  }

  DestroyForest(relation);
  if (INVALID_FOREST != relation) {
    fprintf(stderr, "Couldn't destroy forest of relations\n");
    return 1;
  } else {
    fprintf(stderr, "Destroyed forest of relations\n");
  }

  DestroyDomain(d);
#endif

  // Cleanup
  delete d;

  printf("\nDone\n");
  return 0;
}

