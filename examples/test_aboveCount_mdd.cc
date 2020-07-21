
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



/**
 * testmddlib.cc
 *
 * This is stub to test the simple (or high-level) MDD interface.
 */


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "../src/meddly.h"
#include "../src/meddly_expert.h"

using namespace MEDDLY;

/**
 * Model: A simple service counter
 *
 * The queue can be atmost of length 2.
 * The service counter is either occupied or not.
 * Once a job enters the queue, it can only exit via the service counter.
 *
 * Two variables: Q {0, 1, 2}, C {0, 1}
 * State: (Q, C)
 * Initial state: (0, 0)
 *
 * Transitions:
 *
 * Adding to queue:
 * (0, x) -> (1, x)
 * (1, x) -> (2, x)
 *
 * Moving from queue to service counter:
 * (1, 0) -> (0, 1)
 * (2, 0) -> (1, 1)
 *
 * Exiting from service counter:
 * (x, 1) -> (x, 0)
 */

int main(int argc, char *argv[])
{
  FILE_output meddlyout(stdout);

  initialize();
  // Create a domain
  domain *d = createDomain();
  const int N = 2;
  const int bounds[N] = {4, 2};
  // d->createVariablesTopDown(bounds, N);
  d->createVariablesBottomUp(bounds, N);

  // Create an MDD forest in this domain (to store states)
 forest* states = d->createForest(false, forest::BOOLEAN,
      forest::MULTI_TERMINAL, forest::policies(false));

#if 1
  printf("Constructing initial set of states\n");
#if 1
  // Create an edge in MDD forest
  int** v = (int **) malloc(2 * sizeof(int*));
  v[0] = (int *) malloc((N+1) * sizeof(int));
  v[0][0] = 0; v[0][1] = 0; v[0][2] = 0;
  v[1] = (int *) malloc((N+1) * sizeof(int));
  v[1][0] = 0; v[1][1] = 1; v[1][2] = 0;
  v[2] = (int *) malloc((N+1) * sizeof(int));
  v[2][0] = 0; v[2][1] = 1; v[2][2] = 1;
  dd_edge initial_state(states);
  states->createEdge(v, 3, initial_state);
  initial_state.show(meddlyout, 2);
  /*expert_forest* f=static_cast<expert_forest*>(states);
  node_handle nn=initial_state.getNode();
  printf("nn= %d\n",nn);
  long l=f->getIncomingCount(nn);
  node_handle nn1=f->getDownPtr(nn,0);
  printf("nn1=%d\n",nn1);
  long l1=f->getIncomingCount(nn1);
  printf("l=%d l1=%d\n",l,l1);*/

 double c;

 // node_headers nnh=f->getIncomingCount(nn);
  //l=nnh.getIncomingCount(nn);
//  long l=initial_state.getCardinality();

// node_handle nn=initial_state.getNode();
// expert_forest* f=static_cast<expert_forest*>(states);
// node_headers nnh=node_headers(f);
//
// static_cast<expert_forest*>(states)->getEdgeCount()
// node_headers nnh=node_headers((MEDDLY::expert_forest)states);
// node_headers nh=node_headers(&states);
//
// node_handle nhh=initial_state.getNode();
// long l=nnh.getIncomingCount(nn);
 //long l=0;
// unpacked_node* A = MEDDLY::unpacked_node::newFromNode( states,initial_state, false);
 //printf("initial_state %d\n", initial_state.getIncomingCount());
   apply(IEC, initial_state, c);
   apply(AC, initial_state, c);
   for(int i=0;i<4;i++){
 	  printf("%d \t %d\n", i, abovecount[i]);
    }
//
//  printf(" reachable states\n");
//  fflush(stdout);
  // states->showInfo(meddlyout);
  // initial_state.clear();
  // initial_state.show(meddlyout, 2);
  // states->showInfo(meddlyout);
#else
  // Create an edge in MDD forest
  int** v = (int **) malloc(1 * sizeof(int*));
  v[0] = (int *) malloc((N+1) * sizeof(int));

  v[0][0] = 0; v[0][1] = 0; v[0][2] = 0;
  dd_edge stateA(states);
  states->createEdge(v, 1, stateA);
  stateA.show(meddlyout, 2);

  dd_edge stateD = stateA;
  stateD.show(meddlyout, 2);
  states->showInfo(meddlyout);

  v[0][0] = 0; v[0][1] = 1; v[0][2] = 0;
  dd_edge stateB(states);
  states->createEdge(v, 1, stateB);
  stateB.show(meddlyout, 2);

  stateD += stateB;
  stateD.show(meddlyout, 2);
  states->showInfo(meddlyout);

  v[0][0] = 0; v[0][1] = 2; v[0][2] = 0;
  dd_edge stateC(states);
  states->createEdge(v, 1, stateC);
  stateC.show(meddlyout, 2);

  stateD += stateC;
  stateD.show(meddlyout, 2);
  states->showInfo(meddlyout);

  stateA.clear();
  stateB.clear();
  stateC.clear();
  stateD.clear();
  states->showInfo(meddlyout);
#endif
#endif


  // Cleanup; in this case simply delete the domain
  destroyDomain(d);
  cleanup();

  return 0;
}
