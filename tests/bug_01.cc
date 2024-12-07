
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

    // Create a domain
    const int varsizes[] = { 2, 2 };
    domain *d = domain::createBottomUp(varsizes, 2);
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

#ifdef VERBOSE
    FILE_output out(stdout);
#endif

    // Set up initial state
#ifdef VERBOSE
    out << "Building initial state\n";
    out.flush();
#endif
    minterm initial(mdd);
    dd_edge initialStates(mdd);
    initial.setVar(1, 0);
    initial.setVar(2, DONT_CARE);
    initial.buildFunction(initialStates);
#ifdef VERBOSE
    printf("Initial States:\n");
    initialStates.showGraph(out);
#endif

    // Create a matrix diagram to represent the next-state function
    //
#ifdef VERBOSE
    out << "Building next-state function\n";
    out.flush();
#endif
    minterm_coll nsf_coll(2, mxd);
    dd_edge nsf(mxd);
    nsf_coll.unused().setVars(1, DONT_CARE, 1);
    nsf_coll.unused().setVars(2, 0, 0);
    nsf_coll.pushUnused();
    nsf_coll.unused().setVars(1, DONT_CARE, 1);
    nsf_coll.unused().setVars(2, 0, 1);
    nsf_coll.pushUnused();
    nsf_coll.buildFunction(nsf);
#ifdef VERBOSE
    out << "Next-State Minterms:\n";
    nsf_coll.show(out);
    out << "Next-State Function:\n";
    nsf.showGraph(out);
#endif


    // Generate reachable states
    //
    dd_edge reachBFS(initialStates);
    dd_edge reachDFS(initialStates);

#ifdef VERBOSE
    out << "Building reachable using BFS\n";
    out.flush();
#endif
    apply(REACHABLE_STATES_BFS, reachBFS, nsf, reachBFS);
#ifdef VERBOSE
    out << "BFS states\n";
    reachBFS.showGraph(out);
#endif

#ifdef VERBOSE
    out << "Building reachable using DFS\n";
    out.flush();
#endif
    apply(REACHABLE_STATES_DFS, reachDFS, nsf, reachDFS);
#ifdef VERBOSE
    out << "DFS states\n";
    reachBFS.showGraph(out);
#endif


    int retval = (reachBFS == reachDFS)? 0: 1;
    if (retval) {
        printf("\nReachable states DO NOT match\n\n");
    } else {
        printf("\nReachable states match\n\n");
    }

    // Cleanup
    MEDDLY::cleanup();
    return retval;
}
