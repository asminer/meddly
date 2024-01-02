
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
#include <iostream>
#include <iomanip>

#include "../src/meddly.h"
#include "simple_model.h"
#include "timer.h"
#include "reporting.h"

using namespace MEDDLY;

// #define DEBUG_EVENTS

inline char* newEvent(unsigned N)
{
    char* ev = new char[N*8+2];
    ev[N*8+1] = 0;
    for (unsigned i=N*8; i; i--) ev[i] = '.';
    ev[0] = '_';
    return ev;
}

char* Get(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    t[loc+6] = '-';
    t[loc+8] = '-';
    t[loc+3] = '+';
    t[loc+5] = '+';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event   Get(" << i << "): " << t << "\n";
#endif
    return t;
}

char* Free(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    const unsigned _rt = i ? (8*i-8) : (8*N-8);
    t[loc+5] = '-';
    t[loc+6] = '+';
    t[_rt+3] = '-';
    t[_rt+2] = '+';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event  Free(" << i << "): " << t << "\n";
#endif
    return t;
}

char* Put(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    t[loc+3] = '+';
    t[loc+7] = '+';
    t[loc+4] = '-';
    t[loc+6] = '-';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event   Put(" << i << "): " << t << "\n";
#endif
    return t;
}

char* Used(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    const unsigned _rt = i ? (8*i-8) : (8*N-8);
    t[loc+7] = '-';
    t[loc+6] = '+';
    t[_rt+3] = '-';
    t[_rt+1] = '+';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event  Used(" << i << "): " << t << "\n";
#endif
    return t;
}

char* Other(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    t[loc+1] = '-';
    t[loc+4] = '+';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event Other(" << i << "): " << t << "\n";
#endif
    return t;
}

char* Owner(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    t[loc+1] = '-';
    t[loc+2] = '+';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event Owner(" << i << "): " << t << "\n";
#endif
    return t;
}

char* Write(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    t[loc+2] = '-';
    t[loc+4] = '+';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event Write(" << i << "): " << t << "\n";
#endif
    return t;
}

char* Go(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    const unsigned loc = 8*i;
    t[loc+2] = '-';
    t[loc+8] = '+';
#ifdef DEBUG_EVENTS
    std::cerr << "Built event    Go(" << i << "): " << t << "\n";
#endif
    return t;
}

void markTest(const char* name, const dd_edge &E, const unsigned marks,
        const unsigned counts)
{
    const unsigned mdots = marks/16;
    const unsigned cdots = counts/16;

    forest* ef = E.getForest();
    if (!ef) throw "null forest";
    node_marker M(ef);

    std::cout << "Testing " << name << " (" << M.getSize() << " total nodes)\n";
    std::cout << "     Marking " << std::setw(6) << marks << "x ";
    std::cout.flush();

    timer T;

    for (unsigned i=marks; i; --i) {
        if (0==i % mdots) {
            std::cout << '.';
            std::cout.flush();
        }
        M.unmarkAll();
        M.mark(E.getNode());
    }

    T.note_time();
    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__)) {
        report  << name << " mark $ "
                << "Marked " << name << marks << " times" << std::endl;
    }

    std::cout << "        " << M.countMarked() << " marked nodes\n";
    std::cout << "        " << E.getNodeCount() << " according to dd_edge\n";

    std::cout << "    Counting " << std::setw(6) << counts << "x ";
    std::cout.flush();

    size_t eco = 0;
    for (unsigned i=counts; i; --i) {
        if (0==i % cdots) {
            std::cout << '.';
            std::cout.flush();
        }
        eco = M.countNonzeroEdges();
    }
    T.note_time();
    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__)) {
        report  << name << " count $ "
                << "Counted nodes in " << name << counts << " times" << std::endl;
    }

    std::cout << "        " << eco << " non-zero edges\n";
    std::cout << "        " << E.getEdgeCount(false) << " according to dd_edge\n";

}


void runWithArgs(unsigned N, unsigned marks, unsigned counts)
{
    std::cout << "Building Slotted ring model for N = " << N << std::endl;

    char** events = new char*[8*N];
    char** fill = events;
    for (unsigned i=0; i<N; i++) {
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

    std::cout << "Building domain and forests" << std::endl;

    // Initialize domain
    int* sizes = new int[N*8];
    for (int i=int(N*8-1); i>=0; i--) sizes[i] = 2;
    domain* d = domain::createBottomUp(sizes, N*8);

    // Initialize forests
    forest* mdd = forest::create(d, 0, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, 1, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);

    //
    // Build initial state
    //
    std::cout << "Building initial state" << std::endl;
    int* initial = new int[1+N*8];
    for (unsigned i=N*8; i; i--) initial[i] = 0;
    int* initLocal = initial;
    for (unsigned i=0; i<N; i++) {
        initLocal[3] = initLocal[5] = 1;
        initLocal += 8;
    }
    dd_edge init_state(mdd);
    mdd->createEdge(&initial, 1, init_state);

    //
    // Build next-state function
    //
    std::cout << "Building next-state function" << std::endl;
    dd_edge nsf(mxd);
    buildNextStateFunction(events, 8*N, mxd, nsf);

    //
    // Build reachable states
    //
    std::cout << "Building reachability set" << std::endl;
    dd_edge reachable(mdd);
    apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);

    //
    // Mark timing tests, finally
    //
    markTest("RS ", reachable, marks, counts);
    markTest("nsf ", nsf, marks*16, counts*16);
}


int main(int argc, const char** argv)
{
    try {
        setReport(argc, argv);

        MEDDLY::initialize();

//          runWithArgs(5, 1, 1);
//          return 0;

#ifdef DEVELOPMENT_CODE
        runWithArgs(50, 256, 64);
#else
        runWithArgs(100, 16384, 4096);
#endif
        MEDDLY::cleanup();
        return 0;
    }
    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error '" << e.getName()
                    << "'\n    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    catch (two_strings p) {
        std::cerr << "\t" << p.first << p.second << "\n";
        return 3;
    }
    std::cerr << "\nSome other error?\n";
    return 4;
}


