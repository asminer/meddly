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

#include "../src/meddly.h"
#include "../timing/timer.h"

// #define SHOW_EVENTS

// Number of rings.
int N;

// Approximate cardinality.
bool approx_count;

// **********************************************************************
// build an event to move a ring
//    @param ring2level     Mapping from rings to level
//    @param ring           Which ring to move
//    @param src            Source peg (0, 1, 2)
//    @param dest           Destination peg (0, 1, 2)
//    @param aux            Unused peg (0, 1, 2)
//
//    @param E              in:  blank dd_edge attached to relation forest
//                          out: dd_edge holding the event
// **********************************************************************

void moveRing(std::vector <int> &ring2level, int ring,
        int src, int dest, int aux, MEDDLY::dd_edge &E)
{
    assert(src != dest);
    assert(src != aux);
    assert(dest != aux);

    MEDDLY::minterm delta(E.getForest());

    for (int i=N; i>0; i--) {
        const int k = ring2level[i];
        if (i > ring) {
            delta.setVars(k, MEDDLY::DONT_CARE, MEDDLY::DONT_CHANGE);
            continue;
        }
        if (i < ring) {
            delta.setVars(k, aux, aux);
            continue;
        }
        delta.setVars(k, src, dest);
    }

    delta.setValue(true);
    /*
    MEDDLY::ostream_output merr(std::cerr);
    delta.show(merr);
    */
    delta.buildFunction(false, E);
}

// **********************************************************************
// build overall (partitioned) transition relation
//    @param ring2level     Mapping from rings to level
//    @param prel           Partitioned pre-generated transition relation
// **********************************************************************

void buildRelation(std::vector <int> &ring2level, MEDDLY::pregen_relation& prel)
{
    using namespace MEDDLY;

    timer stopwatch;
    std::cout << "Building transition relation (events grouped by top level)..." << std::endl;
    stopwatch.note_time();
    dd_edge event(prel.getRelForest());

    for (int ring=1; ring<=N; ring++) {
        moveRing(ring2level, ring, 0, 1, 2, event);
        prel.addToRelation(event);
        moveRing(ring2level, ring, 0, 2, 1, event);
        prel.addToRelation(event);
        moveRing(ring2level, ring, 1, 0, 2, event);
        prel.addToRelation(event);
        moveRing(ring2level, ring, 1, 2, 0, event);
        prel.addToRelation(event);
        moveRing(ring2level, ring, 2, 0, 1, event);
        prel.addToRelation(event);
        moveRing(ring2level, ring, 2, 1, 0, event);
        prel.addToRelation(event);
    } // for ring
    prel.finalize();

    stopwatch.note_time();
    std::cout << "    transition relation construction took "
              << stopwatch.get_last_seconds() << " seconds" << std::endl;

#ifdef SHOW_EVENTS
    ostream_output merr(std::cerr);

    for (int k=1; k<=N; k++) {
        merr.indent_more();
        merr << "Events with top level=" << k << ":\n";
        for (unsigned i=0; i<prel.lengthForLevel(k); i++) {
            merr.indent_more();
            merr << "event " << i << ":\n";
            prel.arrayForLevel(k)[i].showGraph(merr);
            merr.indent_less();
        }
        merr.indent_less();
        merr << "\n";
    }
#endif
}

// **********************************************************************
// build reachable states or whatever was asked for
//    @param prel           Partitioned pre-generated transition relation
//    @param initial        Initial state
//    @param reachable      (output) reachable states
// **********************************************************************

void buildReachable(MEDDLY::pregen_relation* prel, MEDDLY::dd_edge &initial,
        MEDDLY::dd_edge &reachable)
{
    using namespace MEDDLY;

    //
    // TBD: do different things based on command line switches
    //

    const char* what = "reachability set";

    std::cout << "Building " << what << " using saturation..." << std::endl;

    timer watch;
    watch.note_time();

    saturation_operation* sat =
        SATURATION_FORWARD(initial.getForest(), prel, reachable.getForest());
    if (!sat) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    sat->compute(initial, reachable);
    watch.note_time();
    std::cout << "    " << what << " construction took "
              << watch.get_last_seconds() << " seconds" << std::endl;


    std::cout << "MDD for " << what << " requires "
              << reachable.getNodeCount() << " nodes" << std::endl;
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

    cout << "Counting states..." << endl;

    oper_item numstates;
    if (approx_count) {
        numstates.init(0.0);
    } else {
#ifdef HAVE_LIBGMP
        mpz_t ns_mpz;
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
        for (unsigned i=0; i<N; i++) {
            theory *= 3.0;
        }
        cout << "Theory: " << theory << " reachable states" << endl;
    } else {
#ifdef HAVE_LIBGMP
        cout << "Exactly ";
        showNumber(numstates.getHugeint());
        cout << " reachable states" << endl;

        mpz_t theory;
        mpz_init_set_ui(theory, 1);
        for (unsigned i=0; i<N; i++) {
            mpz_mul_ui(theory, theory, 3);
        }
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
    cerr << "Usage: " << base << " [switches] #rings\n";
    cerr << "\n";
    cerr << "Generate reachable states for the Towers of Hanoi puzzle.\n";
    cerr << "The number of rings given should be a natural number.\n";
    cerr << "\n";
    cerr << "Switches:\n";
    cerr << "    -h:    This help\n\n";

#ifdef HAVE_LIBGMP
    cerr << "    -a:    Approximate count for number of states (default).\n";
    cerr << "    -e:    Exact count for number of states.\n\n";
#endif

    cerr << "    -o:    Order the variables so that the smallest ring is the bottom-most\n";
    cerr << "           variable in the MDD. (Default)\n";
    cerr << "    -O:    Order the variables so that the largest ring is the bottom-most\n";
    cerr << "           variable in the MDD.\n";

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

    N = -1;
    approx_count = true;
    bool reverse_order = false;
    //
    // Process command line
    //
    for (unsigned i=1; i<argc; i++)
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
                    case 'o':
                                reverse_order = false;
                                continue;

                    case 'O':
                                reverse_order = true;
                                continue;

                    default:
                                if (arg[1] != 'h') {
                                    cerr << "Unknown switch \"-" << arg[1] << "\"\n";
                                }
                                return usage(argv[0]);
                }
            }
            cerr << "Unknown switch \"" << arg << "\"\n";
            return usage(argv[0]);
        }
        if (N>=0) {
            cerr << "More than one #rings specified?\n";
            return usage(argv[0]);
        }
        N = naturalArg(arg);
        if (N<0) {
            return usage(argv[0]);
        }
    }
    if (N<0) {
        cerr << "#rings not specified\n";
        return usage(argv[0]);
    }

    //
    // Build ring to level mappings
    //
    cout << "Towers of Hanoi with " << N << " rings.\n";
    if (0==N) return 0;

    vector <int> ring2level(N+1);
    vector <int> level2ring(N+1);
    ring2level[0] = 0;
    level2ring[0] = 0;
    if (reverse_order) {
        for (unsigned i=N; i>0; i--) {
            ring2level[i] = N-i+1;
            level2ring[N-i+1] = i;
        }
    } else {
        for (unsigned i=N; i>0; i--) {
            ring2level[i] =i;
            level2ring[i] = i;
        }
    }

    cout << "Variable order is TOP, ";
    for (unsigned i=N; i>0; i--) {
        cout << "ring " << level2ring[i] << ", ";
        if (i>3) {
            cout << "..., ";
            i = 3;
        }
    }
    cout << "BOTTOM\n";

    try {
        //
        // Bring up MDD library
        //
        MEDDLY::initialize();

        //
        // Build domain, forests
        //
        int* sizes = new int[N];
        for (int i=N-1; i>=0; i--) sizes[i] = 3;
        domain* d = domain::createBottomUp(sizes, N);
        delete[] sizes;

        policies pmdd(false), pmxd(true);
        // any policy changes?
        forest* mdd = forest::create(d, false, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, pmdd);
        forest* mxd = forest::create(d, true,  range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, pmxd);


        //
        // Build initial state
        //
        minterm init_minterm(mdd);
        init_minterm.setAllVars(0);
        dd_edge initial(mdd);
        init_minterm.buildFunction(false, initial);


        //
        // Build transition relation
        //
        pregen_relation* lnsf = new pregen_relation(mxd);
        buildRelation(ring2level, *lnsf);

        //
        // Build reachable states
        //
        dd_edge reachable(mdd);
        buildReachable(lnsf, initial, reachable);
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
        meddlyout << "MDD stats:\n";
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
