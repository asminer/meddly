
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
    domain *d = domain::create();
    const int N = 2;
    const int bounds[N] = {4, 2};
    // d->createVariablesTopDown(bounds, N);
    d->createVariablesBottomUp(bounds, N);

    // Create an MDD forest in this domain (to store states)
    forest* states = forest::create(d, SET, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL);

    printf("Constructing initial set of states\n");
    dd_edge initial_state(states);
    minterm init(states);
    init.setAllVars(0);
    init.buildFunction(initial_state);
    initial_state.showGraph(meddlyout);


    // Create a MXD forest in domain (to store transition diagrams)
    forest* transitions = forest::create(d, RELATION, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL);

    // Construct a transition diagram in the MXD forest (using +, *)
    // Note: x here denotes "value does not change"
    // (0, x) -> (1, x)
    // (1, x) -> (2, x)
    // (1, 0) -> (0, 1)
    // (2, 0) -> (1, 1)
    // (x, 1) -> (x, 0)
    const unsigned num_of_transitions = 5;
    minterm_coll tlist(num_of_transitions, transitions);

    tlist.unused().setVars(1, 0, 1);
    tlist.unused().setVars(2, DONT_CARE, DONT_CHANGE);
    tlist.pushUnused();

    // vlist[0][0] = 0; vlist[0][1] = 0; vlist[0][2] = -1;
    // vplist[0][0] = 0; vplist[0][1] = 1; vplist[0][2] = -2;

    tlist.unused().setVars(1, 1, 2);
    tlist.unused().setVars(2, DONT_CARE, DONT_CHANGE);
    tlist.pushUnused();

    // vlist[0][0] = 0; vlist[1][1] = 1; vlist[1][2] = -1;
    // vplist[0][0] = 0; vplist[1][1] = 2; vplist[1][2] = -2;

    tlist.unused().setVars(1, 1, 0);
    tlist.unused().setVars(2, 0, 1);
    tlist.pushUnused();

    // vlist[0][0] = 0; vlist[2][1] = 1; vlist[2][2] = 0;
    // vplist[0][0] = 0; vplist[2][1] = 0; vplist[2][2] = 1;

    tlist.unused().setVars(1, 2, 1);
    tlist.unused().setVars(2, 0, 1);
    tlist.pushUnused();

    // vlist[0][0] = 0; vlist[3][1] = 2; vlist[3][2] = 0;
    // vplist[0][0] = 0; vplist[3][1] = 1; vplist[3][2] = 1;

    tlist.unused().setVars(1, DONT_CARE, DONT_CHANGE);
    tlist.unused().setVars(2, 1, 0);
    tlist.pushUnused();

    // vlist[0][0] = 0; vlist[4][1] = -1; vlist[4][2] = 1;
    // vplist[0][0] = 0; vplist[4][1] = -2; vplist[4][2] = 0;

    // Create a edge representing the model's transition diagram
    dd_edge xd(transitions);
    printf("Constructing transition diagram\n");
    tlist.buildFunction(xd);
    xd.showGraph(meddlyout);
    // transitions->showInfo(meddlyout);

    printf("\nCompute Table:\n");
    compute_table::showMonolithicComputeTable(meddlyout, true);

    dd_edge reachableStates(initial_state);
    dd_edge prevReachableStates(states);
    dd_edge postImage(states);

    while(prevReachableStates != reachableStates)
    {
        prevReachableStates = reachableStates;
        printf("\nPost-Image (mdd:%ld, mxd:%ld): ",
                long(reachableStates.getNode()), long(xd.getNode()));
        apply(POST_IMAGE, reachableStates, xd, postImage);
        printf("%ld\n", long(postImage.getNode()));
        // postImage.show(meddlyout, 2);
        printf("\nUnion (mdd:%ld, mdd:%ld): ",
                long(reachableStates.getNode()), long(postImage.getNode()));
        apply(UNION, reachableStates, postImage,
                reachableStates);
        printf("%ld\n", long(reachableStates.getNode()));
    }
    reachableStates.showGraph(meddlyout);

    // Cleanup; in this case simply delete the domain
    domain::destroy(d);
    cleanup();

    return 0;
}
