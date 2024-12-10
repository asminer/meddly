
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

// #define SHOW_INDEXES

#define NEW_EVAL

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

    int sizes[16];

    for (int i=15; i>=0; i--) sizes[i] = N+1;
    domain* dom = domain::createBottomUp(sizes, 16);
    forest* mdd = forest::create(dom, SET, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(dom, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);
    forest* evmdd = forest::create(dom, SET, range_type::INTEGER,
                        edge_labeling::INDEX_SET);

    //
    // Build initial state
    //
    minterm initial(mdd);
    dd_edge init_state(mdd);
    for (int i=16; i; i--) initial.setVar(i, 0);
    initial.setVar(1, N);
    initial.setVar(5, N);
    initial.setVar(9, N);
    initial.setVar(13, N);
    initial.buildFunction(init_state);
    printf("\tbuilt initial state\n");
    fflush(stdout);

    //
    // Build next-state function
    //
    dd_edge nsf(mxd);
    buildNextStateFunction(kanban, 16, mxd, nsf);
    printf("\tbuilt next-state function\n");
    fflush(stdout);

    //
    // Build reachable states
    //
    dd_edge reachable(mdd);
    apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
    printf("\tbuilt reachable states\n");
    fflush(stdout);

    //
    // Build index set for reachable states
    //
    dd_edge reach_index(evmdd);
    apply(CONVERT_TO_INDEX_SET, reachable, reach_index);
#ifdef SHOW_INDEXES
    FILE_output out(stdout);
    printf("\tbuilt index set:\n");
    reach_index.showGraph(out);
#else
    printf("\tbuilt index set\n");
#endif
    fflush(stdout);


  // Verify indexes
  int c = 0;
#ifdef NEW_EVAL
  minterm assign(mdd);
#endif
  for (enumerator s(reachable); s; ++s) {
    long index;
    const int* state = s.getAssignments();
#ifdef NEW_EVAL
    assign.setAll(state, true);   // UGH, copy, for now
    reach_index.evaluate(assign, index);
#else
    evmdd->evaluate(reach_index, state, index);
#endif
    if (c != index) {
      printf("\nState number %d has index %ld\n", c, index);
      return false;
    }
    c++;
  } // for s
  printf("\tverified `forward'\n");
  fflush(stdout);

  // verify the other way
  int d = 0;
  for (enumerator s(reach_index); s; ++s) {
    const int* state = s.getAssignments();
    bool ok;
#ifdef NEW_EVAL
    assign.setAll(state, true);
    reachable.evaluate(assign, ok);
#else
    mdd->evaluate(reachable, state, ok);
#endif
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

  domain::destroy(dom);

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

