
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

#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "../src/meddly.h"
#include "simple_model.h"
#include "../timing/timer.h"
#include "../src/log_simple.h"

using namespace MEDDLY;

FILE_output meddlyout(stdout);

int usage(const char* who)
{
    /* Strip leading directory, if any: */
    const char* name = who;
    for (const char* ptr=who; *ptr; ptr++) {
        if ('/' == *ptr) name = ptr+1;
    }
    printf("\nUsage: %s nnnn [options]\n\n", name);
    printf("\tnnnn: number of parts\n\n");
    printf("\t-bfs: use traditional iterations\n\n");
    printf("\t-dfs: use fastest saturation (currently, -msat)\n");
    printf("\t-esat: use saturation by events\n");
    printf("\t-ksat: use saturation by levels\n");
    printf("\t-msat: use monolithic saturation (default)\n\n");
    printf("\t-exp: use explicit (very slow)\n");
    printf("\t-pdf: Write MDD for reachable states to out.pdf\n\n");
    printf("\t--batch b: specify explicit batch size\n\n");
    printf("\t -l lfile: Write logging information to specified file\n\n");
    return 1;
}

void printStats(const char* who, const forest* f)
{
    printf("%s stats:\n", who);
    f->reportStats(meddlyout, "\t",
            HUMAN_READABLE_MEMORY  |
            BASIC_STATS | EXTRA_STATS |
            STORAGE_STATS | HOLE_MANAGER_STATS |
            HOLE_MANAGER_DETAILED
    );
}



inline char* newEvent(int N)
{
    char* ev = new char[N*8+1];
    for (int i=N*8; i; i--) ev[i] = '.';
    ev[0] = '_';
    return ev;
}

char* Get(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[6] = '-';
    tloc[8] = '-';
    tloc[3] = '+';
    tloc[5] = '+';
    return t;
}

char* Free(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    char* t_rt = i ? t+(8*i-8) : t+(8*N-8);
    tloc[5] = '-';
    tloc[6] = '+';
    t_rt[3] = '-';
    t_rt[2] = '+';
    return t;
}

char* Put(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[3] = '+';
    tloc[7] = '+';
    tloc[4] = '-';
    tloc[6] = '-';
    return t;
}

char* Used(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    char* t_rt = i ? t+(8*i-8) : t+(8*N-8);
    tloc[7] = '-';
    tloc[6] = '+';
    t_rt[3] = '-';
    t_rt[1] = '+';
    return t;
}

char* Other(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[1] = '-';
    tloc[4] = '+';
    return t;
}

char* Owner(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[1] = '-';
    tloc[2] = '+';
    return t;
}

char* Write(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[2] = '-';
    tloc[4] = '+';
    return t;
}

char* Go(int i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[2] = '-';
    tloc[8] = '+';
    return t;
}

void runWithArgs(int N, char method, int batchsize, bool build_pdf, logger* LOG)
{
    timer start;

    printf("+-------------------------------------------------+\n");
    printf("|   Initializing Slotted ring model for N = %-4d  |\n", N);
    printf("+-------------------------------------------------+\n");
    fflush(stdout);

    char** events = new char*[8*N];
    char** fill = events;
    for (int i=0; i<N; i++) {
        fill[0] = Other(i, N);
        fill[1] = Owner(i, N);
        fill[2] = Write(i, N);
        fill[3] = Go(i, N);
        fill[4] = Get(i, N);
        fill[5] = Put(i, N);
        fill[6] = Used(i, N);
        fill[7] = Free(i, N);
        fill += 8;
    }

    // Initialize domain
    int* sizes = new int[N*8];
    for (int i=N*8-1; i>=0; i--) sizes[i] = 2;
    domain* d = domain::createBottomUp(sizes, N*8);

    // Initialize forests
    forest* mdd = forest::create(d, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL);
    if (LOG) {
        mdd->setLogger(LOG, "MDD");
        mxd->setLogger(LOG, "MxD");
    }

    //
    // Build initial state
    //
    if (LOG) LOG->newPhase(mdd, "Building initial state");
    minterm initial(mdd);
    dd_edge init_state(mdd);
    initial.setAllVars(0);
    for (int i=0; i<N; i++) {
        initial.setVar(i*8+3, 1);
        initial.setVar(i*8+5, 1);
    }
    initial.buildFunction(false, init_state);

    //
    // Build next-state function
    //
    if (LOG) LOG->newPhase(mxd, "Building next-state function");
    dd_edge nsf(mxd);
    pregen_relation* ensf = nullptr;
    saturation_operation* sat = nullptr;

    if ('s' == method) {
        ensf = new pregen_relation(mxd, 8*N);
    }
    if ('k' == method) {
        ensf = new pregen_relation(mxd);
    }

    if ('e' != method) {

        if (ensf) {
            start.note_time();
            buildNextStateFunction(events, 8*N, ensf, 2);
            start.note_time();
        } else {
            start.note_time();
            buildNextStateFunction(events, 8*N, mxd, nsf, 2);
            start.note_time();
#ifdef DUMP_NSF
            printf("Next-state function:\n");
            nsf.show(meddlyout, 2);
#endif
        }
        printf("Next-state function construction took %.4e seconds\n",
                start.get_last_seconds() );
        printStats("MxD", mxd);
    }

    //
    // Build reachable states
    //
    if (LOG) LOG->newPhase(mdd, "Building reachability set");
    dd_edge reachable(mdd);
    start.note_time();
    switch (method) {
        case 'b':
            printf("Building reachability set using traditional algorithm\n");
            fflush(stdout);
            apply(REACHABLE_STATES_BFS, init_state, nsf, reachable);
            break;

        case 'm':
            printf("Building reachability set using saturation, monolithic relation\n");
            fflush(stdout);
            apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
            break;

        case 'e':
            printf("Building reachability set using explicit search\n");
            printf("Using batch size: %d\n", batchsize);
            fflush(stdout);
            explicitReachset(events, 8*N, mdd, init_state, reachable, batchsize);
            break;

        case 'k':
        case 's':
            printf("Building reachability set using saturation, relation");
            if ('k'==method)  printf(" by levels\n");
            else              printf(" by events\n");
            fflush(stdout);
            sat = SATURATION_FORWARD(mdd, ensf, mdd);
            if (!sat) {
                throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
            }
            sat->compute(init_state, reachable);
            break;

        default:
            printf("Error - unknown method\n");
            exit(2);
    }
    start.note_time();
    printf("Done\n");
    printf("Reachability set construction took %.4e seconds\n",
            start.get_last_seconds());
    fflush(stdout);

#ifdef SHOW_STATES
    int count = 0;
    for (enumerator i(reachable); i; ++i, ++count) {
        const int* element = i.getAssignments();
        printf("State %4d: [%d", count, element[1]);
        for (int j=2; j<=8*N; j++) {
            printf(", %d", element[j]);
        } // for j
        printf("]\n");
    }  // for i
#endif

    printStats("MDD", mdd);
    fflush(stdout);

    double c;
    apply(CARDINALITY, reachable, c);
    compute_table::showAll(meddlyout, 3);

    printf("Approx. %g reachable states\n", c);
    operation::destroy(sat);
    // or, don't, and let cleanup() take care of it?

    if (build_pdf) {
        reachable.setLabel("reachable");
        dot_maker dm(mdd, "out");
        dm.addRootEdge(reachable);
        dm.doneGraph();
        dm.runDot("pdf");
    }

    if (LOG) {
        LOG->newPhase(mdd, "Cleanup");
        LOG->newPhase(mxd, "Cleanup");
        domain::destroy(d);
    }
}

int main(int argc, const char** argv)
{
    int N = -1;
    char method = 'm';
    int batchsize = 256;
    const char* lfile = 0;
    bool build_pdf = false;

    for (int i=1; i<argc; i++) {
        if (strcmp("-bfs", argv[i])==0) {
            method = 'b';
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
        if (strcmp("-l", argv[i])==0) {
            lfile = argv[i+1];
            i++;
            continue;
        }
        if (strcmp("-pdf", argv[i])==0) {
            build_pdf = true;
            continue;
        }
        if (strcmp("--batch", argv[i])==0) {
            i++;
            if (argv[i]) batchsize = atoi(argv[i]);
            continue;
        }
        N = atoi(argv[i]);
    }

    if (N<0) return usage(argv[0]);

    MEDDLY::initialize();

    //
    // Set up logger, if any
    //

    FILE_output log;
    logger* LOG = 0;
    if (lfile) {
        FILE* fout = fopen(lfile, "w");
        if (!fout) {
            printf("Couldn't open %s for writing, no logging\n", lfile);
        } else {
            log.setFILE(fout);
            LOG = new simple_logger(log);
            LOG->recordNodeCounts();
            char comment[80];
            snprintf(comment, 80, "Automatically generated by slot (N=%d)", N);
            LOG->addComment(comment);
        }
    }

    try {
        runWithArgs(N, method, batchsize, build_pdf, LOG);
        delete LOG;
        MEDDLY::cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        printf("Caught MEDDLY error: %s\n", e.getName());
        return 1;
    }

}


