
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
  Testing Saturation.

  Test1: Test fully reduced MDDs with Saturation.
*/

#include "../src/meddly.h"
#include <cstdlib>
#include <cassert>

#define VERBOSE

using namespace MEDDLY;

int main(int argc, char *argv[])
{
    // Initialize MEDDLY
    MEDDLY::initialize();

    // Number of levels in domain (excluding terminals)
    const int nLevels = 2;

    // Set up arrays bounds
    variable** vars = new variable*[1+nLevels];
    vars[0] = nullptr;
    for (unsigned i = nLevels; i; i--) {
        vars[i] = createVariable(2, "");
    }

    // Create a domain and set up the state variables.
    domain *d = domain::create(vars, nLevels);
    assert(d);

    // Create an MDD forest in this domain (to store states)
    policies pmdd(false);
    pmdd.setFullyReduced();
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
                            edge_labeling::MULTI_TERMINAL, pmdd);
    assert(mdd);

    // Create a MXD forest in domain (to store transition diagrams)
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                            edge_labeling::MULTI_TERMINAL);
    assert(mxd);

    // Set up initial state
    minterm initial(mdd);
    dd_edge initialStates(mdd);
    initial.setVar(1, 0);
    initial.setVar(2, DONT_CARE);
    initial.buildFunction(initialStates);

    // Create a matrix diagram to represent the next-state function
    //
    minterm_coll nsf_coll(2, mxd);
    dd_edge nsf(mxd);
    nsf_coll.unused().setVars(1, DONT_CARE, 1);
    nsf_coll.unused().setVars(2, 0, 0);
    nsf_coll.pushUnused();
    nsf_coll.unused().setVars(1, DONT_CARE, 1);
    nsf_coll.unused().setVars(2, 0, 1);
    nsf_coll.pushUnused();

    nsf_coll.buildFunction(nsf);

    /*
  int from1[nLevels+1] = {0, DONT_CARE, 0};
  int from2[nLevels+1] = {0, DONT_CARE, 0};
  int to1[nLevels+1] = {0, 1, 0};
  int to2[nLevels+1] = {0, 1, 1};
  int *from[2] = { from1, from2 };
  int *to[2] = { to1, to2 };
  mxd->createEdge(reinterpret_cast<int**>(from),
      reinterpret_cast<int**>(to), 2, nsf);
      */

    dd_edge reachBFS(initialStates);
    dd_edge reachDFS(initialStates);

    apply(REACHABLE_STATES_BFS, reachBFS, nsf, reachBFS);
    apply(REACHABLE_STATES_DFS, reachDFS, nsf, reachDFS);

    int retval = (reachBFS == reachDFS)? 0: 1;

#ifdef VERBOSE
    FILE_output meddlyout(stdout);

    printf("Initial States:\n");
    initialStates.showGraph(meddlyout);

    printf("Next-State Minterms:\n");
    nsf_coll.show(meddlyout, nullptr, "\n");

    printf("Next-State Function:\n");
    nsf.showGraph(meddlyout);

    printf("BFS states\n");
    reachBFS.showGraph(meddlyout);

    printf("DFS states\n");
    reachDFS.showGraph(meddlyout);

    if (retval) {
        printf("\nReachable states DO NOT match\n\n");
    } else {
        printf("\nReachable states match\n\n");
    }

#endif

    // Cleanup
    MEDDLY::cleanup();
    return retval;
}
