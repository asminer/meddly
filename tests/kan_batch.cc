
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

#include "kan_rs1.h"
#include "kan_rs2.h"
#include "kan_rs3.h"

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

// 160 states for N=1
long expected[] = {
  1, 160, 4600, 58400, 454475, 2546432, 11261376,
  41644800, 133865325, 384392800, 1005927208
};

using namespace MEDDLY;

dd_edge buildReachsetSAT(forest* mdd, forest* mxd, int N)
{
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

    //
    // Build next-state function
    //
    dd_edge nsf(mxd);
    buildNextStateFunction(kanban, 16, mxd, nsf);

    //
    // Build reachable states
    //
    dd_edge reachable(mdd);
    apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);

    return reachable;
}

void mark2minterm(const char* mark, minterm &m)
{
    for (unsigned i=m.getNumVars(); i; i--) {
        m.setVar(i, mark[i]-'0');
    }
}

dd_edge buildReachsetBatch(int b, forest* mdd, const char* rs[], long states)
{
    minterm_coll batch(b, mdd);
    dd_edge reachable(mdd);

    dd_edge tmp(mdd);
    for (long s=0; s<states; s++) {
        if (batch.size() >= batch.maxsize()) {
            batch.buildFunction(tmp);
            reachable += tmp;
            batch.clear();
        }
        mark2minterm(rs[s]-1, batch.unused());
        batch.pushUnused();
    }
    batch.buildFunction(tmp);
    reachable += tmp;

    return reachable;
}

bool checkRS(int N, const char* rs[])
{
    bool ok = true;
    const int batches[] = { 1000, 100, 10, 1, 0 };
    int sizes[16];

    for (int i=15; i>=0; i--) sizes[i] = N+1;
    domain* d = domain::createBottomUp(sizes, 16);
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL);

    dd_edge reachable = buildReachsetSAT(mdd, mxd, N);

    for (int b=0; batches[b]; b++) {
        printf("\tBuilding %ld markings by hand, %d at a time\n",
                expected[N], batches[b]);
        fflush(stdout);

        dd_edge brs = buildReachsetBatch(batches[b], mdd, rs, expected[N]);

        if (brs != reachable) {
            printf("\t\tError, didn't match\n");
            ok = false;
            break;
        }
    }

    domain::destroy(d);
    return ok;
}


int main()
{
    MEDDLY::initialize();

    printf("Checking Kanban reachability set, N=1\n");
    if (!checkRS(1, kanban_rs1)) return 1;

    printf("Checking Kanban reachability set, N=2\n");
    if (!checkRS(2, kanban_rs2)) return 1;

    printf("Checking Kanban reachability set, N=3\n");
    if (!checkRS(3, kanban_rs3)) return 1;

    MEDDLY::cleanup();
    printf("Done\n");
    return 0;
}

