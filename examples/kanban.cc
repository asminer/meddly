
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
#include <iostream>
#include <iomanip>
#include <string.h>
#include <gmp.h>

#include "../src/meddly.h"
#include "simple_model.h"
#include "../timing/timer.h"
#include "../src/log_simple.h"

// #define DUMP_NSF
// #define DUMP_REACHABLE

const char* kanban[] = {
    "X-+..............",    // Tin1 TA
    "X.-+.............",    // Tr1 TB
    "X.+-.............",    // Tb1 TC
    "X.-.+............",    // Tg1 TD
    "X.....-+.........",    // Tr2 TE
    "X.....+-.........",    // Tb2 TF
    "X.....-.+........",    // Tg2 TG
    "X+..--+..-+......",    // Ts1_23 TH
    "X.........-+.....",    // Tr3 TI
    "X.........+-.....",    // Tb3 TJ
    "X.........-.+....",    // Tg3 TK
    "X....+..-+..--+..",    // Ts23_4 TL
    "X.............-+.",    // Tr4 TM
    "X.............+-.",    // Tb4 TN
    "X............+..-",    // Tout4 TO
    "X.............-.+"     // Tg4 TP
};

using namespace MEDDLY;
using namespace std;

ostream_output meddlyout(cout);

bool verbose;

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


int usage(const char* who)
{
    /* Strip leading directory, if any: */
    const char* name = who;
    for (const char* ptr=who; *ptr; ptr++) {
        if ('/' == *ptr) name = ptr+1;
    }
    cout << "\nUsage: " << name << " [options] nnnn\n\n";
    cout << "\tnnnn  : number of parts\n";
    cout << "\t-trad : use traditional iterations, without frontier\n";
    cout << "\t-front: use traditional iterations, with frontier\n";
    cout << "\t-dfs  : use saturation\n";
    cout << "\t-esat : use saturation by events\n";
    cout << "\t-ksat : use saturation by levels\n";
    cout << "\t-msat : use monolithic saturation (default)\n\n";

    cout << "\t-edges:  count number of (actual) edges in reachability graph\n\n";
    cout << "\t-approx: approximate the number of states (default)\n";
    cout << "\t-exact:  determine the exact number of states (requires gmp)\n\n";
    cout << "\t-exp: use explicit (very slow)\n\n";
    cout << "\t--batch b: specify explicit batch size\n\n";
    cout << "\t -l lfile: Write logging information to specified file\n\n";
    cout << "\t-pdf: write the MDD representing the reachable states to Kanban.pdf\n\n";
    cout << "\t-ms: mark and sweep\n";
    cout << "\t-opt: optimistic reference counts\n";
    cout << "\t-pess: pessimistic reference counts\n\n";

    cout << "\t-v   : verbose\n";
    return 1;
}

void printStats(const char* who, const forest* f)
{
    cout << who << " stats:\n";
    f->reportStats(meddlyout, "\t",
        HUMAN_READABLE_MEMORY  |
        BASIC_STATS | EXTRA_STATS |
        STORAGE_STATS | STORAGE_DETAILED |
        HOLE_MANAGER_STATS | HOLE_MANAGER_DETAILED
    );
}

void showCard(const oper_item &x, const char* what)
{
    if (x.hasType(opnd_type::REAL)) {
        cout << "Approx. " << setprecision(0) << x.getReal()
             << ' ' << what << endl;
        return;
    }
#ifdef HAVE_LIBGMP
    unsigned digits = mpz_sizeinbase(x.getHugeint(), 10) + 2;
    char* str = new char[digits];
    mpz_get_str(str, 10, x.getHugeint());
    cout << "Exactly " << str << ' ' << what << endl;
    delete[] str;
    return;
#endif
    throw "unknown cardinality type";
}

int main(int argc, const char** argv)
{
    int N = -1;
    char method = 'm';
    unsigned batchsize = 256;
    const char* lfile = 0;
    bool build_pdf = false;
    bool approx_count = true;
    bool count_edges = false;
    char node_del = 'o';
    verbose = false;

    //
    // Process command line
    //
    for (int i=1; i<argc; i++) {
        if (strcmp("-trad", argv[i])==0) {
            method = 't';
            continue;
        }
        if (strcmp("-front", argv[i])==0) {
            method = 'f';
            continue;
        }
        if (strcmp("-dfs", argv[i])==0) {
            method = 'm';
            continue;
        }
        if (strcmp("-esat", argv[i])==0) {
            method = 's';
            continue;
        }
        if (strcmp("-ksat", argv[i])==0) {
            method = 'k';
            continue;
        }
        if (strcmp("-msat", argv[i])==0) {
            method = 'm';
            continue;
        }
        if (strcmp("-exp", argv[i])==0) {
            method = 'e';
            continue;
        }
        if (strcmp("-edges", argv[i])==0) {
            count_edges = true;
            continue;
        }
        if (strcmp("-approx", argv[i])==0) {
            approx_count = true;
            continue;
        }
        if (strcmp("-exact", argv[i])==0) {
            approx_count = false;
            continue;
        }
        if (strcmp("-l", argv[i])==0) {
            lfile = argv[i+1];
            i++;
            continue;
        }
        if (strcmp("--batch", argv[i])==0) {
            i++;
            if (argv[i]) {
                int bs = atoi(argv[i]);
                if (bs>0) batchsize = unsigned(bs);
            }
            continue;
        }
        if (strcmp("-pdf", argv[i])==0) {
            build_pdf = true;
            continue;
        }
        if (strcmp("-ms", argv[i])==0) {
            node_del='m';
            continue;
        }
        if (strcmp("-opt", argv[i])==0) {
            node_del='o';
            continue;
        }
        if (strcmp("-pess", argv[i])==0) {
            node_del='p';
            continue;
        }
        if (strcmp("-v", argv[i]) == 0) {
            verbose = true;
            continue;
        }

        if (argv[i][0] == '-') {
            return usage(argv[0]);
        }

        N = atoi(argv[i]);
    }

    if (N<0) return usage(argv[0]);

    try {
        MEDDLY::initialize();

        timer start;

        cout << "+-------------------------------------------+\n";
        cout << "|   Initializing Kanban model for N = "
             << left << setw(4) << N << "  |\n";
        cout << "+-------------------------------------------+" << endl;

        // Initialize domain
        int* sizes = new int[16];
        for (int i=15; i>=0; i--) sizes[i] = N+1;
        domain* d = domain::createBottomUp(sizes, 16);
        delete[] sizes;

        // Initialize forests
        policies pmdd(false), pmxd(true);
        switch (node_del) {
            case 'm':
                cout << "Using mark and sweep\n";
                pmdd.useReferenceCounts = false;
                pmxd.useReferenceCounts = false;
                break;

            case 'p':
                cout << "Using pessimistic node deletion\n";
                pmdd.useReferenceCounts = true;
                pmxd.useReferenceCounts = true;
                pmdd.setPessimistic();
                pmxd.setPessimistic();
                break;

            default:
                cout << "Using optimistic node deletion\n";
                pmdd.useReferenceCounts = true;
                pmxd.useReferenceCounts = true;
                pmdd.setOptimistic();
                pmxd.setOptimistic();
                break;
        }
        forest* mdd = forest::create(d, 0, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, pmdd);
        forest* mxd = forest::create(d, 1, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL, pmxd);

        // associate loggers
        FILE_output log;
        logger* LOG = nullptr;
        if (lfile) {
            FILE* fout = fopen(lfile, "w");
            if (!fout) {
                cout << "Couldn't open " << lfile
                     << " for writing, no logging\n";
            } else {
                log.setFILE(fout);
                LOG = new simple_logger(log);
                LOG->recordNodeCounts();
                LOG->recordTimeStamps();
                char comment[80];
                snprintf(comment, 80, "Automatically generated by kanban (N=%d)", N);
                LOG->addComment(comment);
                mdd->setLogger(LOG, "MDD");
                mxd->setLogger(LOG, "MxD");
            }
        }

        //
        // Build initial state
        //
        if (LOG) LOG->newPhase(mdd, "Building initial state");

        minterm initial(mdd);
        dd_edge init_state(mdd);
        initial.setAllVars(0);
        initial.setVar(1, N);
        initial.setVar(5, N);
        initial.setVar(9, N);
        initial.setVar(13, N);
        initial.buildFunction(false, init_state);

        //
        // Build next-state function
        //
        if (LOG) LOG->newPhase(mxd, "Building next-state function");
        dd_edge nsf(mxd);
        pregen_relation* ensf = nullptr;
        saturation_operation* sat = nullptr;

        if ('s' == method) {
            ensf = new pregen_relation(mxd, 16);
        }
        if ('k' == method) {
            ensf = new pregen_relation(mxd);
        }

        if ('e' != method) {

            if (ensf) {
                start.note_time();
                buildNextStateFunction(kanban, 16, ensf, 4);
                start.note_time();
            } else {
                start.note_time();
                buildNextStateFunction(kanban, 16, mxd, nsf, 4);
                start.note_time();
#ifdef DUMP_NSF
                cout << "Next-state function:\n";
                nsf.showGraph(meddlyout);
#endif
            }
            cout << "Next-state function construction took "
                 << start.get_last_seconds() << " seconds\n";
            printStats("MxD", mxd);
        }

        if (LOG) LOG->newPhase(mdd, "Building reachability set");
        dd_edge reachable(mdd);
        start.note_time();
        binary_operation* rsgen = nullptr;
        switch (method) {
            case 't':
                cout << "Building reachability set using traditional algorithm (no frontier)" << endl;
                rsgen = build(REACHABLE_TRAD_NOFS(true), mdd, mxd, mdd);
                if (verbose) rsgen->setProgressNotifier(my_progress);
                rsgen->compute(init_state, nsf, reachable);
                break;

            case 'f':
                cout << "Building reachability set using traditional algorithm (with frontier)" << endl;
                rsgen = build(REACHABLE_TRAD_FS(true), mdd, mxd, mdd);
                if (verbose) rsgen->setProgressNotifier(my_progress);
                rsgen->compute(init_state, nsf, reachable);
                break;

            case 'm':
                cout << "Building reachability set using saturation, monolithic relation"
                     << endl;
                apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
                break;

            case 'e':
                cout << "Building reachability set using explicit search\n";
                cout << "Using batch size: " << batchsize << endl;
                explicitReachset(kanban, 16, mdd, init_state, reachable, batchsize);
                break;

            case 'k':
            case 's':
                cout << "Building reachability set using saturation, relation";
                if ('k'==method)  cout << " by levels" << endl;
                else              cout << " by events" << endl;
                sat = SATURATION_FORWARD(mdd, ensf, mdd);
                if (0==sat) {
                    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
                }
                sat->compute(init_state, reachable);
                break;

            default:
                cout << "Error - unknown method\n";
                exit(2);
        };
        start.note_time();
        cout << "Done\n";
        cout << "Reachability set construction took "
             << start.get_last_seconds() << " seconds" << endl;
#ifdef DUMP_REACHABLE
        cout << "Reachable states:\n";
        reachable.showGraph(meddlyout);
#endif

        printStats("MDD", mdd);

        oper_item numstates, numedges;
#ifdef HAVE_LIBGMP
        mpz_t ns_mpz, ne_mpz;
        if (!approx_count) {
            mpz_init(ns_mpz);
            mpz_init(ne_mpz);
            numstates.init(ns_mpz);
            numedges.init(ne_mpz);
        } else {
#endif
            numstates.init(0.0);
            numedges.init(0.0);
#ifdef HAVE_LIBGMP
        }
#endif
        apply(CARDINALITY, reachable, numstates);
        start.note_time();
        cout << "RS Cardinality took " << start.get_last_seconds()
             << " seconds" << endl;

        if (count_edges) {
            dd_edge rsxrs(mxd);
            apply(CROSS, reachable, reachable, rsxrs);
            apply(INTERSECTION, rsxrs, nsf, rsxrs);
            start.note_time();
            cout << "Reachability graph took " << start.get_last_seconds()
                 << " seconds" << endl;
            apply(CARDINALITY, rsxrs, numedges);
            start.note_time();
            cout << "RG Cardinality took " << start.get_last_seconds()
                 << " seconds" << endl;
        }

        compute_table::showAll(meddlyout, 3);

        showCard(numstates, "reachable states");
        if (count_edges) {
            showCard(numedges, "graph edges");
        }

        if (build_pdf) {
            reachable.setLabel("reachable");
            dot_maker mdd_dot(mdd, "kanban");
            mdd_dot.addRootEdge(reachable);
            mdd_dot.doneGraph();
            mdd_dot.runDot("pdf");
            if ('m' == method) {
                nsf.setLabel("next-state");
                dot_maker mxd_dot(mxd, "kanban-nsf");
                mxd_dot.addRootEdge(nsf);
                mxd_dot.doneGraph();
                mxd_dot.runDot("pdf");
            }
        }

        // cleanup
        if (LOG) {
            LOG->newPhase(mdd, "Cleanup");
            LOG->newPhase(mxd, "Cleanup");
            domain::destroy(d);
            delete LOG;
        }
        MEDDLY::cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        cerr    << "\nCaught meddly error '" << e.getName()
                << "'\n    thrown in " << e.getFile()
                << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        cerr    << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    cerr << "\nSome other error?\n";
    return 4;
}

