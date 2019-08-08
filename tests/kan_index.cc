
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

 #define SHOW_INDEXES

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

bool equal(const int* a, const int* b, int N)
{
  for (int i=1; i<=N; i++)
    if (a[i] != b[i]) return false;
  return true;
}

bool checkReachset(int N)
{
  printf("Running test for N=%d...\n", N);
  FILE_output meddlyout(stdout);

  int sizes[16];

  for (int i=15; i>=0; i--) sizes[i] = N+1;
  domain* dom = createDomainBottomUp(variable::variableTypes::boundedClass,sizes, 16);

  // Build initial state
  int_extra* initial = new int_extra[17];
  for (int i=16; i; i--) initial[i] = 0;
  initial[1] = initial[5] = initial[9] = initial[13] = N;
//  int_extra* ginitial = new int_extra[17];
//  for (int i=16; i; i--) ginitial[i] = initial[i];
  forest* mdd = dom->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge init_state(mdd);
  mdd->createEdge(&initial, 1, init_state);
  init_state.show(meddlyout,2);
  delete[] initial;
//  delete[] ginitial;
  printf("\tbuilt initial state\n");
  fflush(stdout);

  // Build next-state function
  forest* mxd = dom->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  buildNextStateFunction(kanban, 16, mxd, nsf); 
  printf("\tbuilt next-state function\n");
  fflush(stdout);

  // Build reachable states
  dd_edge reachable(mdd);
  apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
  printf("\tbuilt reachable states\n");
  reachable.show(meddlyout,2);
  fflush(stdout);

  // Build index set for reachable states
  forest* evmdd = dom->createForest(0, forest::INTEGER, forest::INDEX_SET);
  dd_edge reach_index(evmdd);
  apply(CONVERT_TO_INDEX_SET, reachable, reach_index);
#ifdef SHOW_INDEXES
  printf("\tbuilt index set:\n");
  reach_index.show(meddlyout, 2);

//  reachable.show(meddlyout, 2);
  // reach_index.show(stdout, 2);
#else
  printf("\tbuilt index set\n");
#endif
  fflush(stdout);

  // Verify indexes
  int c = 0;
  /*for (enumerator s(reachable); s; ++s) {
    const int* state = s.getAssignments();
//    printf("state SIZE IS   %d %d\n", state[0],state[1]);
//    int size=sizeof(state)/sizeof( int);
//
//    int_extra gstate[size];
//    for(int i=0;i<size;i++)
//    {  gstate[i]=state[i];
//
//    }
//    printf("gstate SIZE IS  %d %d\n", gstate[0],gstate[1]);

//     const int_extra* ggstate=s.getGIAssignments();
//    printf("ggstate SIZE IS  %d %d\n", ggstate[0].getInteger(),ggstate[1].getInteger());
//    fflush(stdout);
    long index;
    evmdd->evaluate(reach_index, state, index);
    if (c != index) {
      printf("\nState number %d has index %ld\n", c, index);
      return false;
    }
    c++;
  } // for s
  printf("\tverified `forward'\n");
  fflush(stdout);*/

  // verify the other way
  int d = 0;
  for (enumerator s(reach_index); s; ++s) {
    const int* state = s.getAssignments();
    int size=sizeof(state)/sizeof(const int);
    printf("size %d %d\n",state[0],state[1]);
    const int_extra* gstate=s.getGIAssignments();
    printf("gsize %d %d\n",gstate[0].getInteger(),gstate[1].getInteger());
//        int_extra* gstate= new int_extra[size];
//        for(int i=0;i<size;i++)
//          gstate[i]=state[i];
//       const int_extra* ggstate=gstate;
    bool ok;
    mdd->evaluate(reachable,gstate, ok);
    if (!ok) {
      printf("\nIndex number %d does not appear in reachability set\n", d);
      return false;
    }
    d++;
  }
  printf("\tverified `backward'\n");
  fflush(stdout);

  // Verify index search
  c = 0;
  int elem[17];
  for (enumerator s(reachable); s; ++s) {
    const int* state = s.getAssignments();
    evmdd->getElement(reach_index, c, elem); 
    if (!equal(state, elem, 16)) {
      printf("\nFetch index %d got wrong state\n", c);
      return false;
    }
    c++;
  } // for s
  printf("\tverified `getElement'\n");
  fflush(stdout);


  if (c!=d) {
    printf("\nCardinality mismatch\n");
    return false;
  }

  destroyDomain(dom);
  
  return true;
}

int main()
{
  MEDDLY::initialize();

  for (int n=1; n<4; n++) {
    if (!checkReachset(n)) return 1;
  }

  MEDDLY::cleanup();
  printf("Done\n");
  return 0;
}

