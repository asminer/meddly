
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

// #define OLD_ITERATOR

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

dd_edge buildReachset(domain* d, int N)
{
    //
    // Build forests
    //
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);

    //
    // Build initial state
    //
    minterm initial(mdd);
    dd_edge init_state(mdd);
    initial.setAllVars(0);
    initial.setVar(1, N);
    initial.setVar(5, N);
    initial.setVar(9, N);
    initial.setVar(13, N);
    initial.buildFunction(false, init_state);
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

    return reachable;
}

#ifdef OLD_ITERATOR

bool matches(const char* mark, const int* minterm, int np)
{
    for (int i=0; i<np; i++)
        if (mark[i]-48 != minterm[i]) return false;
    return true;
}

#else

bool matches(const char* mark, const minterm &m)
{
    for (unsigned i=m.getNumVars(); i; --i) {
        if (mark[i-1]-48 != m.from(i)) return false;
    }
    return true;
}

#endif

void checkRS(int N, const char* rs[])
{
    int sizes[16];

    for (int i=15; i>=0; i--) sizes[i] = N+1;
    domain* d = domain::createBottomUp(sizes, 16);

    dd_edge reachable = buildReachset(d, N);

    // enumerate states
    long c = 0;

#ifdef OLD_ITERATOR
    for (enumerator i(reachable); i; ++i)
#else
    for (dd_edge::iterator i = reachable.begin(); i; ++i)
#endif
    {
        if (c>=expected[N]) {
            throw "too many markings";
        }
#ifdef OLD_ITERATOR
        if (!matches(rs[c], i.getAssignments()+1, 16))
#else
        if (!matches(rs[c], *i))
#endif
        {
            fprintf(stderr, "Marking %ld mismatched\n", c);
            break;
        }
        c++;
    }
    if (c<expected[N]) {
        throw "not enough markings";
    }

    domain::destroy(d);
}


/*
void showRS(int N)
{
    int sizes[16];

    for (int i=15; i>=0; i--) sizes[i] = N+1;
    domain* d = domain::createBottomUp(sizes, 16);

    dd_edge reachable = buildReachset(d, N);

    // enumerate states
    long c = 0;
    for (enumerator i(reachable); i; ++i) {
        const int* minterm = i.getAssignments();
        printf("%5ld: %d", c, minterm[1]);
        for (int l=2; l<=16; l++) printf(", %d", minterm[l]);
        printf("\n");
        c++;
    }

    domain::destroy(d);
}
*/

int main()
{
    try {
        MEDDLY::initialize();

        printf("Checking Kanban reachability set, N=1\n");
        checkRS(1, kanban_rs1);
        printf("\t%ld markings checked out\n", expected[1]);

        printf("Checking Kanban reachability set, N=2\n");
        checkRS(2, kanban_rs2);
        printf("\t%ld markings checked out\n", expected[2]);

        printf("Checking Kanban reachability set, N=3\n");
        checkRS(3, kanban_rs3);
        printf("\t%ld markings checked out\n", expected[3]);

        MEDDLY::cleanup();
        printf("Done\n");
        return 0;
    }
    catch (MEDDLY::error e) {
        fprintf(stderr, "\nCaught meddly error %s\n", e.getName());
        fprintf(stderr, "    thrown in %s line %u\n",
                e.getFile(), e.getLine());
        return 1;
    }
    catch (const char* e) {
        fprintf(stderr, "\nCaught our own error: %s\n", e);
        return 2;
    }
    return 4;
}

