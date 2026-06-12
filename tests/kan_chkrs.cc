
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

struct statedist {
    const char* state;
    unsigned distance;
};

#define HAVE_ORACLE

#ifdef HAVE_ORACLE
#include "kan_fwd_df1.h"
#include "kan_fwd_df2.h"
#include "kan_fwd_df3.h"
#include "kan_bwd_df1.h"
#include "kan_bwd_df2.h"
#include "kan_bwd_df3.h"
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

bool FORWD;

using namespace MEDDLY;

// ============================================

domain* buildDomain(int N)
{
    static int sizes[16];
    for (unsigned i=0; i<16; i++) {
        sizes[i] = N+1;
    }
    return domain::createBottomUp(sizes, 16);
}

// ============================================

// Adjust multi-terminal, integer distance function
// to be distance plus one, for reachable states,
// zero for unreachable states.

void distance_adjust(const rangeval &in, rangeval &out)
{
    long x = long(in);
    if (x<0) {
        out = 0;
    } else {
        out = x+1;
    }
}

// ============================================

dd_edge buildReachset(domain* d, int N, range_type R, edge_labeling E,
        char method)
{
    //
    // Build forests
    //
    forest* mdd = forest::create(d, SET, R, E);
    if (range_type::BOOLEAN == R) {
        std::cout << "\tusing MTMDD for reachable states\n";
    } else {
        if (edge_labeling::EVPLUS == E) {
            std::cout << "\tusing EV+MDD for distance\n";
        } else {
            std::cout << "\tusing MTMDD for distance\n";
        }
    }
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


    if (range_type::BOOLEAN == R) {
        initial.setValue(true);
        initial.buildFunction(false, init_state);
    } else {
        initial.setValue(0);
        if (edge_labeling::EVPLUS == E) {
            rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);
            initial.buildFunction(infty, init_state);
        } else {
            initial.buildFunction(-1, init_state);
        }
    }

    std::cout << "\t\tbuilt initial state\n";

    //
    // Build next-state function
    //
    // TBD: monolithic vs implicit vs ...
    //
    dd_edge nsf(mxd);
    buildNextStateFunction(kanban, 16, mxd, nsf);
    std::cout << "\t\tbuilt next-state function\n";

    //
    // Build reachable states
    //
    const char* fwd = FORWD ? " forward " : " backward ";

    dd_edge reachable(mdd);
    switch (method) {
        case 'F':
            std::cout << "\t\tUsing traditional" << fwd << "generation, with frontier\n";
            apply(REACHABLE_TRAD_FS(FORWD), init_state, nsf, reachable);
            break;

        case 'T':
            std::cout << "\t\tUsing traditional" << fwd << "generation, without frontier\n";
            apply(REACHABLE_TRAD_NOFS(FORWD), init_state, nsf, reachable);
            break;

        case '1':
            std::cout << "\t\tUsing" << fwd << "saturation v1\n";
            apply(REACHABLE_SATUR(FORWD, 1), init_state, nsf, reachable);
            break;

        case '2':
            std::cout << "\t\tUsing" << fwd << "saturation v2\n";
            apply(REACHABLE_SATUR(FORWD, 2), init_state, nsf, reachable);
            break;

        default:
#ifdef ALLOW_DEPRECATED_0_18_1
            std::cout << "\t\tusing original" << fwd << "saturation\n";
            if (FORWD) {
                apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
            } else {
                apply(REVERSE_REACHABLE_DFS, init_state, nsf, reachable);
            }
#else
            throw "Unexpected method";
#endif
    };

    if ( (range_type::INTEGER == R) && (edge_labeling::EVPLUS != E) )
    {
        // Adjust distances if needed
        user_unary_factory adjust("distadjust", distance_adjust);
        apply(adjust, reachable, reachable);
    }

    return reachable;
}

// ============================================

bool matches(const long adjust, const statedist &sd, const minterm &m)
{
    if (!sd.state) {
        throw "not enough reachable markings";
    }

    for (unsigned k=1; k <= 16; ++k) {
        if (sd.state[k]-48 != m.from(k)) return false;
    }

    const rangeval& oracle = m.getValue();
    if (oracle.isBoolean()) return true;

    if (oracle.isInteger()) {
        return long(oracle) == sd.distance + adjust;
    }

    throw "Unexpected state value type";
}

// ============================================

void checkRS(int N, range_type R, edge_labeling E, char method)
{
    std::cout << "Checking Kanban reachability set, N=" << N << "\n";

    domain* d = buildDomain(N);
    dd_edge reachable = buildReachset(d, N, R, E, method);

    const statedist* krs = nullptr;
#ifdef HAVE_ORACLE
    switch (N) {
        case 1:
            krs = FORWD ? kanban_fwd_rs1 : kanban_bwd_rs1;
            break;

        case 2:
            krs = FORWD ? kanban_fwd_rs2 : kanban_bwd_rs2;
            break;

        case 3:
            krs = FORWD ? kanban_fwd_rs3 : kanban_bwd_rs3;
            break;

        default:
            krs = nullptr;
    }
#endif

    if (krs) {
        std::cout << "\tverifying against oracle\n";

        // enumerate states
        long c = 0;

        const long adjust =
            (range_type::INTEGER == R) && (edge_labeling::MULTI_TERMINAL == E)
            ? 1
            : 0;

        for (dd_edge::iterator i = reachable.begin(); i; ++i)
        {
            if (!matches( adjust, krs[c], *i) ) {
                std::cerr << "Marking " << c << " mismatched\n";

                ostream_output merr(std::cerr);
                merr << "Marking was ";
                (*i).show(merr);
                merr << "\n";

                throw "mismatch";
            }
            c++;
        }
        if (krs[c].state) {
            throw "not enough markings";
        }

        std::cout << "\tperfect match!\n";
    }

    domain::destroy(d);
}

// ============================================

void genRS(int N)
{
    using namespace std;

    domain* d = buildDomain(N);
    dd_edge reachable = buildReachset(d, N, range_type::INTEGER,
                            edge_labeling::EVPLUS, 'T');

    cout << "//\n";
    cout << "// Finite "
         << (FORWD ? "forward" : "backward")
         << " distance function for Kanban N=" << N << "\n";
    cout << "//\n";
    cout << "\n";
    cout << "const statedist kanban_" << (FORWD ? "fwd" : "bwd")
         << "_rs" << N << "[] = {\n";

    // enumerate states
    for (dd_edge::iterator i = reachable.begin(); i; ++i)
    {
        cout << "    { \"b";
        for (unsigned k=1; k <= 16; ++k) {
            cout << (*i).from(k);
        }
        cout << "\", " << long((*i).getValue()) << "},\n";
    }
    cout << "    { nullptr, 0 }\n";
    cout << "};\n\n";

    domain::destroy(d);

}

// ============================================

int Usage(const char* exe)
{
    using namespace std;

    cerr << "\n";
    cerr << "Usage: " << exe << " [ switches ]\n";
    cerr << "\n";
    cerr << "Check reachable states (with distances as appropriate) on small\n";
    cerr << "Kanban instances.\n";
    cerr << "\n";
    cerr << "Switches:\n";
    cerr << "    -b     : backward generation\n";
    cerr << "    -f     : forward generation (default)\n";
    cerr << "\n";
    cerr << "    -c     : check against expected states and distances (default)\n";
    cerr << "    -g n   : generate header for N=n instance. Use n>0. Does not\n";
    cerr << "             perform any checking. Ignores all -- switches.\n";
    cerr << "\n";
    cerr << "    --mtb  : use boolean MTMDDs (reachable states)\n";
    cerr << "    --mti  : use integer MTMDDs (distances)\n";
    cerr << "    --ev   : use EV+MDDs (distances)\n";
    cerr << "\n";
#ifdef ALLOW_DEPRECATED_0_18_1
    cerr << "    --dfs  : use original saturation\n";
#endif
    cerr << "    --sat1 : use saturation v1 (default)\n";
    cerr << "    --sat2 : use saturation v2\n";
    cerr << "    --front: use traditional iteration, with frontier set\n";
    cerr << "    --trad : use traditional iteration, without frontier set\n";

    return 1;
}

// ============================================



int main(int argc, const char** argv)
{
    long geninstance = 0;
    range_type rtype = range_type::BOOLEAN;
    edge_labeling elmdd = edge_labeling::MULTI_TERMINAL;
    char method = '1';
    FORWD = true;

    //
    // Process switches
    //

    for (int i=1; i<argc; i++) {
        if (0==strcmp("-g", argv[i])) {
            if (argv[i+1]) {
                geninstance = atol(argv[i+1]);
                if (geninstance<1) {
                    return Usage(argv[0]);
                }
            }
            ++i;
            continue;
        }
        if (0==strcmp("-c", argv[i])) {
            geninstance = 0;
            continue;
        }

        //
        // FORWD or not
        //
        if (0==strcmp("-b", argv[i])) {
            FORWD = false;
            continue;
        }
        if (0==strcmp("-f", argv[i])) {
            FORWD = true;
            continue;
        }


        //
        // mt vs ev
        //

        if (0==strcmp("--mtb", argv[i])) {
            rtype = range_type::BOOLEAN;
            elmdd = edge_labeling::MULTI_TERMINAL;
            continue;
        }
        if (0==strcmp("--mti", argv[i])) {
            rtype = range_type::INTEGER;
            elmdd = edge_labeling::MULTI_TERMINAL;
            continue;
        }
        if (0==strcmp("--ev", argv[i])) {
            rtype = range_type::INTEGER;
            elmdd = edge_labeling::EVPLUS;
            continue;
        }

        //
        // saturation vs traditional
        //

#ifdef ALLOW_DEPRECATED_0_18_1
        if (0==strcmp("--dfs", argv[i])) {
            method = 'S';
            continue;
        }
#endif
        if (0==strcmp("--sat1", argv[i])) {
            method = '1';
            continue;
        }
        if (0==strcmp("--sat2", argv[i])) {
            method = '2';
            continue;
        }
        if (0==strcmp("--trad", argv[i])) {
            method = 'T';
            continue;
        }
        if (0==strcmp("--front", argv[i])) {
            method = 'F';
            continue;
        }

        return Usage(argv[0]);
    }

    //
    // Run
    //

    try {
        MEDDLY::initialize();

        if (geninstance) {
            genRS(geninstance);
            MEDDLY::cleanup();
            return 0;
        }

        checkRS(1, rtype, elmdd, method);
        checkRS(2, rtype, elmdd, method);
        checkRS(3, rtype, elmdd, method);

        MEDDLY::cleanup();
        std::cout << "Done\n";
        return 0;
    }
    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error " << e.getName() << "\n"
                    << "    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr   << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    return 4;
}

