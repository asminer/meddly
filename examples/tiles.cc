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

#include <iostream>
#include <vector>
#include <cassert>

#include "../src/meddly.h"
#include "../timing/timer.h"
#include "simple_model.h"

// #define SHOW_EVENTS

// Number of rows.
int R;

// Number of columns.
int C;

// Approximate cardinality.
bool approx_count;

// Verbose output: show iteration details
bool verbose;

// #define DEBUG_EVENTS

// There are two encodings.
// If this variable is true, then the state encoding is,
//     "for each position, which tile is currently there."
// Otherwise, the state encoding is,
//     "for each tile, where is it."
//
bool for_each_position_which_tile;


void my_progress(unsigned iter, char st)
{
    if (' ' == st) {
        std::cerr << "    Iteration " << iter << ": ";
        return;
    }
    if (';' == st) {
        std::cerr << std::endl;
        return;
    }
    std::cerr << st << " ";
}


// **********************************************************************
// build the initial state, as a minterm
// **********************************************************************

void buildInitial(MEDDLY::minterm &init)
{
    //
    // The initial state:
    //      starting from the upper left, reading from left to right,
    //      the tiles are 1, 2, 3, ...
    //      and the bottom right corner is empty.
    //

    if (for_each_position_which_tile) {
        // p0: 1, p1: 2, ..., plast-1: R*C-1, plast: 0
        for (int p=0; p<R*C-1; p++) {
            init.setVar(p+1, p+1);
        }
        init.setVar(R*C, 0);
    } else {
        // t0: R*C-1, t1: 0, t2: 1, ...
        for (int t=1; t<R*C; t++) {
            init.setVar(t+1, t-1);
        }
        init.setVar(1, R*C-1);
    }
}


// **********************************************************************
// build a minterm to move a tile
//
//    @param T              Tile number, from 1 to R*C-1.
//    @param p              Position currently containing tile T,
//                              from 0 to R*C-1.
//    @param q              Neighboring position containing nothing (0),
//                              from 0 to R*C-1.
//    @param m              Minterm to fill
//
//    @return true, iff this is a valid move
//
// **********************************************************************

bool moveTile(const int T, const int p, const int q, MEDDLY::minterm &m)
{
    if (p == q)     return false;
    if (T<=0)       return false;
    if (T>= R*C)    return false;
    if (p<0)        return false;
    if (p>= R*C)    return false;
    if (q<0)        return false;
    if (q>= R*C)    return false;

    m.setAllVars(MEDDLY::DONT_CARE, MEDDLY::DONT_CHANGE);

    if (for_each_position_which_tile) {

        // x_p: T -> 0
        //
        m.setVars(p+1, T, 0);

        // x_q: 0 -> T
        //
        m.setVars(q+1, 0, T);

    } else {

        // x_T: p -> q
        //
        m.setVars(T+1, p, q);

        // x_0: q -> p
        //
        m.setVars(1, q, p);
    }

    m.setValue(true);
#ifdef DEBUG_EVENTS
    std::cerr << "   move tile " << T << " from " << p << " to " << q << "\n";
#endif
    return true;
}

// **********************************************************************
// build overall, monolithic, transition relation
//    @param trel           Transition relation, on output
// **********************************************************************

void buildRelation(MEDDLY::dd_edge& trel)
{
    using namespace MEDDLY;

    timer stopwatch;
    std::cout << "Building minterm collection..." << std::endl;
    stopwatch.note_time();

    // Lazy bound: number of tiles, squared, times 4 events
    minterm_coll mtc(R*C*R*C*4, trel.getForest());

    for (int p=0; p<R*C; p++) {
        for (int t=1; t<R*C; t++) {
            // Move to the left, if we can
            if (0 != p%C) {
                if (moveTile(t, p, p-1, mtc.unused())) {
                    mtc.pushUnused();
                }
            }
            // Move to the right, if we can
            if (C-1 != p%C) {
                if (moveTile(t, p, p+1, mtc.unused())) {
                    mtc.pushUnused();
                }
            }
            // Move up
            if (moveTile(t, p, p-C, mtc.unused())) {
                mtc.pushUnused();
            }

            // Move down
            if (moveTile(t, p, p+C, mtc.unused())) {
                mtc.pushUnused();
            }
        }
    }

    stopwatch.note_time();
    std::cout << "    minterm collection construction took "
              << stopwatch.get_last_seconds() << " seconds" << std::endl;

    std::cout << "    there are " << mtc.size() << " events\n";

    std::cout << "Building relation from minterms..." << std::endl;
    stopwatch.note_time();

    mtc.buildFunctionMax(false, trel);

    stopwatch.note_time();
    std::cout << "    relation construction took "
              << stopwatch.get_last_seconds() << " seconds" << std::endl;
}

// **********************************************************************
// build reachable states or whatever was asked for
//    @param method         What to do:
//                              'm': saturation, monolithic
//
//                              'f': BFS with frontier
//                              'b': BFS without frontier
//
//    @param prel           Partitioned pre-generated transition relation
//    @param initial        Initial state
//    @param reachable      (output) reachable states
// **********************************************************************

void buildReachable(bool dist, char method, const MEDDLY::dd_edge &relation,
        const MEDDLY::dd_edge &initial, MEDDLY::dd_edge &reachable)
{
    using namespace MEDDLY;
    timer watch;

    watch.note_time();

    if (dist) {
        std::cout << "Building distance function using ";
    } else {
        std::cout << "Building reachable states using ";
    }

    binary_operation* bfs = nullptr;
    forest* initF = initial.getForest();
    forest* relF = relation.getForest();

    switch (method) {
        case 'm':
                    std::cout << "saturation..." << std::endl;
                    apply(REACHABLE_STATES_DFS, initial, relation, reachable);
                    break;

        case 'f':
                    std::cout << "traditional with frontier..." << std::endl;
                    bfs = build(REACHABLE_TRAD_FS(true), initF, relF, initF);
                    if (!bfs) throw "null bfs";
                    if (verbose) bfs->setProgressNotifier(my_progress);
                    bfs->compute(initial, relation, reachable);
                    break;

        case 't':
                    std::cout << "traditional without frontier..." << std::endl;
                    bfs = build(REACHABLE_TRAD_NOFS(true), initF, relF, initF);
                    if (!bfs) throw "null bfs";
                    if (verbose) bfs->setProgressNotifier(my_progress);
                    bfs->compute(initial, relation, reachable);
                    break;

        default:
                    throw "unknown method";
    }

    watch.note_time();
    std::cout << "    reachable states construction took "
              << watch.get_last_seconds() << " seconds" << std::endl;

    if (dist) {
        std::cout << "EV+MDD for distance function requires ";
    } else {
        std::cout << "MDD for reachable states requires ";
    }
    std::cout << reachable.getNodeCount() << " nodes" << std::endl;
}

// **********************************************************************
// Add comma separators to a large integer
// **********************************************************************
#ifdef HAVE_LIBGMP
void showNumber(mpz_t x)
{
    unsigned digits = mpz_sizeinbase(x, 10)+2;
    char* str = new char[digits];
    mpz_get_str(str, 10, x);
    unsigned nextcomma = strlen(str) % 3;
    if (0==nextcomma) {
        nextcomma = 3;
    }
    for (unsigned i=0; str[i]; i++) {
        if (i == nextcomma) {
            std::cout << ",";
            nextcomma += 3;
        }
        std::cout << str[i];
    }
    delete[] str;
}
#endif

// **********************************************************************
// compute and show cardinality
//    @param reachable      reachable states
// **********************************************************************
void showCardinality(const MEDDLY::dd_edge &reachable)
{
    using namespace std;
    using namespace MEDDLY;

#ifdef HAVE_LIBGMP
    mpz_t ns_mpz;
#endif

    cout << "Counting states..." << endl;

    oper_item numstates;
    if (approx_count) {
        numstates.init(0.0);
    } else {
#ifdef HAVE_LIBGMP
        mpz_init(ns_mpz);
        numstates.init(ns_mpz);
#endif
    }
    timer watch;
    watch.note_time();
    apply(CARDINALITY, reachable, numstates);
    watch.note_time();
    cout << "    counting took " << watch.get_last_seconds()
         << " seconds" << endl;

    if (approx_count) {
        cout << "Approx. " << numstates.getReal()
             << " reachable states" << endl;

        double theory = 1.0;
        for (unsigned i=2; i<=R*C; i++) {
            theory *= i;
        }
        theory /= 2.0;
        cout << "Theory: " << theory << " reachable states" << endl;
    } else {
#ifdef HAVE_LIBGMP
        cout << "Exactly ";
        showNumber(numstates.getHugeint());
        cout << " reachable states" << endl;

        mpz_t theory;
        mpz_init_set_ui(theory, 1);
        for (unsigned i=2; i<=R*C; i++) {
            mpz_mul_ui(theory, theory, i);
        }
        mpz_div_ui(theory, theory, 2);
        cout << "Theory: ";
        showNumber(theory);
        cout << " reachable states" << endl;
#endif
    }

}

// **********************************************************************
// print usage instructions
// **********************************************************************

int usage(const char* exe)
{
    using namespace std;
    const char* base = exe;
    for (; *exe; ++exe)
    {
        if ('/' == *exe)
        {
            base = exe+1;
        }
    }

    cerr << "\n";
    cerr << "Usage: " << base << " [switches] #rows #cols\n";
    cerr << "\n";
    cerr << "Generate reachable states for a 2-d sliding puzzle.\n";
    cerr << "The number of rows and columns given should be natural numbers.\n";
    cerr << "The standard \"16 puzzle\" uses 4 rows and 4 columns.\n";
    cerr << "\n";
    cerr << "Switches:\n";
    cerr << "    -h:    This help\n\n";

#ifdef HAVE_LIBGMP
    cerr << "    -a:    Approximate count for number of states (default).\n";
    cerr << "    -e:    Exact count for number of states.\n\n";
#endif

    cerr << "    -q:    Quiet. Do not show BFS iteration details (default).\n";
    cerr << "    -v:    Verbose. Show BFS iteration details.\n";
    cerr << "\n";

    cerr << "    --pos:     Use a position-based encoding: for each position,\n";
    cerr << "               which tile is there.\n";
    cerr << "    --tile:    Use a tile-based encoding: for each tile,\n";
    cerr << "               where is it (default)\n";
    cerr << "\n";

    cerr << "    --reach:   Build the reachability set (default).\n";
    cerr << "    --dist:    Build the distance function.\n";
    cerr << "\n";

    cerr << "    --msat     Saturation, monolithic relation\n";
    cerr << "    --trad     Traditional BFS without frontier set\n";
    cerr << "    --front    Traditional BFS with frontier set (RS only)\n";
    cerr << "\n";
    return 1;
}

inline int naturalArg(const char* arg)
{
    if (0==arg) return -1;
    return atol(arg);
}

// **********************************************************************
// main
// **********************************************************************

int main(int argc, const char** argv)
{
    using namespace std;
    using namespace MEDDLY;

    R = -1;
    C = -1;
    approx_count = true;
    verbose = false;
    char satmethod = 'm';
    for_each_position_which_tile = false;
    bool distances = false;

    //
    // Process command line
    //
    unsigned i;
    for (i=1; i<argc; i++)
    {
        const char *arg = argv[i];
        if ('-' == arg[0])
        {
            // single character switches
            if (arg[1] && 0==arg[2]) {
                switch (arg[1])
                {
#ifdef HAVE_LIBGMP
                    case 'a':   approx_count = true;
                                continue;

                    case 'e':   approx_count = false;
                                continue;
#endif
                    case 'q':
                                verbose = false;
                                continue;

                    case 'v':
                                verbose = true;
                                continue;

                    default:
                                if (arg[1] != 'h') {
                                    cerr << "Unknown switch \"-" << arg[1] << "\"\n";
                                }
                                return usage(argv[0]);
                }
            }
            // double dash switches
            if ('-' == arg[1]) {
                if (0==strcmp("--pos", arg)) {
                    for_each_position_which_tile = true;
                    continue;
                }
                if (0==strcmp("--tile", arg)) {
                    for_each_position_which_tile = false;
                    continue;
                }

                if (0==strcmp("--reach", arg)) {
                    distances = false;
                    continue;
                }
                if (0==strcmp("--dist", arg)) {
                    distances = true;
                    continue;
                }

                if (0==strcmp("--msat", arg)) {
                    satmethod = 'm';
                    continue;
                }

                if (0==strcmp("--trad", arg)) {
                    satmethod = 't';
                    continue;
                }
                if (0==strcmp("--front", arg)) {
                    satmethod = 'f';
                    continue;
                }
            }
            cerr << "Unknown switch \"" << arg << "\"\n";
            return usage(argv[0]);
        }
        break;
    }
    if (i < argc) {
        R = naturalArg(argv[i]);
        i++;
    }
    if (i < argc) {
        C = naturalArg(argv[i]);
        i++;
    }
    if ((i != argc) || (R<0) || (C<0)) {
        return usage(argv[0]);
    }

    if (0==R) return 0;
    if (0==C) return 0;

    try {
        //
        // Bring up MDD library
        //
        MEDDLY::initialize();

        //
        // Build domain, forests
        //
        int* sizes = new int[R*C];
        for (int i=R*C-1; i>=0; i--) sizes[i] = R*C;
        domain* d = domain::createBottomUp(sizes, R*C);
        delete[] sizes;

        policies pmdd(false), pmxd(true);
        // any policy changes?
        forest* mdd = nullptr;

        if (distances) {
            mdd = forest::create(d, false, range_type::INTEGER, edge_labeling::EVPLUS, pmdd);
        } else {
            mdd = forest::create(d, false, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, pmdd);
        }

        forest* mxd = forest::create(d, true,  range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, pmxd);

        cout << "Sliding tile puzzle, board has " << R << " rows and "
             << C << " columns.\n";
        if (for_each_position_which_tile) {
            cout << "Encoding using 'for each position, which tile is there'\n";
        } else {
            cout << "Encoding using 'for each tile, where is it'\n";
        }

        //
        // Build initial state
        //
        dd_edge initial(mdd);
        minterm init_minterm(mdd);
        buildInitial(init_minterm);

        if (distances) {
            init_minterm.setValue(0);
            rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);
            init_minterm.buildFunction(infty, initial);
        } else {
            init_minterm.setValue(true);
            init_minterm.buildFunction(false, initial);
        }

        //
        // Build transition relation
        //
        dd_edge rel(mxd);
        buildRelation(rel);

        //
        // Build reachable states
        //
        dd_edge reachable(mdd);
        buildReachable(distances, satmethod, rel, initial, reachable);
        showCardinality(reachable);

        //
        // Reporting
        //
        ostream_output meddlyout(cout);
        meddlyout << "\n======================================================================\n\n";
        meddlyout << "MxD stats:\n";
        mxd->reportStats(meddlyout, "    ",
                HUMAN_READABLE_MEMORY | BASIC_STATS
        );
        if (distances) {
            meddlyout << "EV+MDD stats:\n";
        } else {
            meddlyout << "MDD stats:\n";
        }
        mdd->reportStats(meddlyout, "    ",
            HUMAN_READABLE_MEMORY | BASIC_STATS | EXTRA_STATS
        );


        cout << "Done!\n";
        return 0;
    }
    catch (MEDDLY::error e) {
        cerr    << "\nCaught meddly error '" << e.getName()
                << "'\n    thrown in " << e.getFile()
                << " line " << e.getLine() << "\n";
        return 2;
    }
    catch (const char* e) {
        cerr    << "\nCaught our own error: " << e << "\n";
        return 4;
    }
    cerr << "\nSome other error?\n";
    return 6;
}
