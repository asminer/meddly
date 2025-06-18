
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

#define HAVE_ORACLE

#ifdef HAVE_ORACLE
#include "kan_dist1.h"
#endif

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

void show_states(const bool fwd, unsigned dist, const dd_edge &ftr)
{
    using namespace std;
    if (fwd) {
        cout << "//\n// States at distance " << dist << "\n//\n";
        cout << "\nconst char* kanban_fwd_" << dist << "[] = {\n";
    } else {
        cout << "//\n// States at reverse distance " << dist << "\n//\n";
        cout << "\nconst char* kanban_bwd_" << dist << "[] = {\n";
    }
    for (dd_edge::iterator i = ftr.begin(); i; ++i)
    {
        cout << "    \"b";
        for (unsigned k=1; k <= 16; ++k) {
            cout << (*i).from(k);
        }
        cout << "\",\n";
    }
    cout << "    nullptr\n};\n\n";
}

void show_collection(const bool fwd, unsigned dist)
{
    using namespace std;
    const char* name = fwd ? "fwd" : "bwd";
    cout << "\nconst char** kanban_" << name << "[] = {\n";
    for (unsigned i=0; i<dist; i++) {
        cout << "    kanban_" << name << "_" << i << ",\n";
    }
    cout << "    nullptr\n};\n\n";
}


bool matches(const char* mark, const minterm &m)
{
    for (unsigned i=m.getNumVars(); i; --i) {
        if (mark[i]-48 != m.from(i)) return false;
    }
    return true;
}

void compare_states(const bool fwd, unsigned dist, const dd_edge &ftr,
        const char** oracle)
{
    using namespace std;

    cout << "    checking distance " << dist << "\n";
    unsigned card = 0;

    ostream_output mout(cout);

    for (dd_edge::iterator i = ftr.begin(); i; ++i)
    {
        if (!oracle[card]) {
            cout << "        more states than expected\n";
            throw "too large set";
        }
        if (!matches(oracle[card], *i))
        {
            cout << "        marking " << card << " mismatched\n";
            throw "mismatched set";
        }
        ++card;
    }
    if (oracle[card]) {
        cout << "        fewer states than expected";
        throw "too small set";
    }

    cout << "        " << card << " markings matched\n";
}

//
// Build next forward or backward step.
//
// On input:
//     reachable: all states reachable up to current distance d
//     frontier:  ignored
//
// On output:
//      reachable: all states reachable up to distance d+1
//      frontier:  states at distance d+1
//
void next_step(const bool fwd, dd_edge &reachable,
        const dd_edge &nsf, dd_edge &frontier)
{
    dd_edge s_x_R(reachable);

    if (fwd) {
        apply(POST_IMAGE, reachable, nsf, s_x_R);
    } else {
        apply(PRE_IMAGE, reachable, nsf, s_x_R);
    }
    apply(DIFFERENCE, s_x_R, reachable, frontier);
    apply(UNION, reachable, s_x_R, reachable);
}

//
// Generate forward/backward distance sets
//
void genStates(const bool fwd)
{
    int sizes[16];

    for (int i=15; i>=0; i--) sizes[i] = 2;

    domain* d = domain::createBottomUp(sizes, 16);
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL);

    // Build next-state function
    dd_edge nsf(mxd);
    buildNextStateFunction(kanban, 16, mxd, nsf);

    // Build initial state
    dd_edge reachset(mdd);
    minterm initial(mdd);
    initial.setAllVars(0);
    initial.setVar(1, 1);
    initial.setVar(5, 1);
    initial.setVar(9, 1);
    initial.setVar(13, 1);
    initial.buildFunction(false, reachset);
    dd_edge frontier(reachset);

    // Classic BFS-style iteration
    unsigned dist = 0;
    do {
        show_states(fwd, dist, frontier);
        next_step(fwd, reachset, nsf, frontier);
        ++dist;
    } while (frontier.getNode());
    show_collection(fwd, dist);
}

//
// Generate and Check forward/backward distance sets
//
void checkStates(const bool fwd)
{
    std::cout << "Checking " << (fwd ? "forward" : "backward") << " states\n";
    int sizes[16];

    for (int i=15; i>=0; i--) sizes[i] = 2;

    domain* d = domain::createBottomUp(sizes, 16);
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL);

    // Build next-state function
    dd_edge nsf(mxd);
    buildNextStateFunction(kanban, 16, mxd, nsf);

    // Build initial state
    dd_edge reachset(mdd);
    minterm initial(mdd);
    initial.setAllVars(0);
    initial.setVar(1, 1);
    initial.setVar(5, 1);
    initial.setVar(9, 1);
    initial.setVar(13, 1);
    initial.buildFunction(false, reachset);
    dd_edge frontier(reachset);

    // Classic BFS-style iteration
#ifdef HAVE_ORACLE
    const char*** oracle = fwd ? kanban_fwd : kanban_bwd;
#endif
    unsigned dist = 0;
    do {
#ifdef HAVE_ORACLE
        if (!oracle[dist]) {
            throw "distance mismatch";
        }
        compare_states(fwd, dist, frontier, oracle[dist]);
#endif
        next_step(fwd, reachset, nsf, frontier);
        ++dist;
    } while (frontier.getNode());
}

// ============================================

int Usage(const char* exe)
{
    using namespace std;

    cerr << "\n";
    cerr << "Usage: " << exe << " [ switches ]\n";
    cerr << "\n";
    cerr << "Check pre- and post-image computations on Kanban N=1.\n";
    cerr << "\n";
    cerr << "Switches:\n";
    cerr << "    -c: check against expected states (default)\n";
    cerr << "    -g: generate header for checking states; do not perform check\n";
    cerr << "\n";

    return 1;
}

// ============================================


int main(int argc, const char** argv)
{
    bool generate = false;
    for (int i=1; i<argc; i++) {

        if (0==strcmp("-c", argv[i])) {
            generate = false;
            continue;
        }

        if (0==strcmp("-g", argv[i])) {
            generate = true;
            continue;
        }

        return Usage(argv[0]);
    }

    using namespace std;

    try {
        MEDDLY::initialize();

        if (generate) {
            genStates(true);
            genStates(false);
        } else {
            checkStates(true);
            checkStates(false);
        }

        MEDDLY::cleanup();
        return 0;
    }

    catch (MEDDLY::error e) {
        cerr << "\nCaught Meddly error '" << e.getName() << "'\n"
             << "    thrown in " << e.getFile()
             << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    cerr << "\nSome other error?\n";
    return 4;
}

