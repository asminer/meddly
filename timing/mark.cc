
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
#include <fstream>

#include "../src/meddly.h"
#include "simple_model.h"
#include "timer.h"

using namespace MEDDLY;



inline char* newEvent(unsigned N)
{
    char* ev = new char[N*8+1];
    for (unsigned i=N*8; i; i--) ev[i] = '.';
    ev[0] = '_';
    return ev;
}

char* Get(unsigned i, int N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[6] = '-';
    tloc[8] = '-';
    tloc[3] = '+';
    tloc[5] = '+';
    return t;
}

char* Free(unsigned i, unsigned N)
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

char* Put(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[3] = '+';
    tloc[7] = '+';
    tloc[4] = '-';
    tloc[6] = '-';
    return t;
}

char* Used(unsigned i, unsigned N)
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

char* Other(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[1] = '-';
    tloc[4] = '+';
    return t;
}

char* Owner(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[1] = '-';
    tloc[2] = '+';
    return t;
}

char* Write(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[2] = '-';
    tloc[4] = '+';
    return t;
}

char* Go(unsigned i, unsigned N)
{
    char* t = newEvent(N);
    char* tloc = t+8*i;
    tloc[2] = '-';
    tloc[8] = '+';
    return t;
}

void markTest(const char* name, const dd_edge &E, unsigned marks, unsigned dots)
{
    expert_forest* ef = dynamic_cast<expert_forest*> (E.getForest());
    if (!ef) throw "null expert_forest";
    node_marker* M = ef->makeNodeMarker();
    if (!M) throw "null node_marker";

    std::cout << "Marking " << name;
    std::cout.flush();
    timer T;

    while (marks) {
        if (0==marks % dots) {
            std::cout << '.';
            std::cout.flush();
        }
        --marks;
        M->unmarkAll();
        M->mark(E.getNode());
    }

    T.note_time();
    std::cout << ' ' << T.get_last_seconds() << " seconds\n";
    std::cout << "    " << M->countMarked() << " marked nodes\n";
    std::cout << "    " << ef->getNodeCount(E.getNode()) << " according to expert_forest\n";
    delete M;

}

void runWithArgs(unsigned N, unsigned marks, unsigned dots)
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
    for (int i=N*8-1; i>=0; i--) sizes[i] = 2;
    domain* d = domain::createBottomUp(sizes, N*8);

    // Initialize forests
    forest* mdd = forest::create(d, 0, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
    forest* mxd = forest::create(d, 1, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);

    //
    // Build initial state
    //
    std::cout << "Building initial state" << std::endl;
    int* initial = new int[1+N*8];
    for (int i=N*8; i; i--) initial[i] = 0;
    int* initLocal = initial;
    for (int i=0; i<N; i++) {
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
    markTest("reachability set    ", reachable, marks, dots);
    markTest("transition relation ", nsf, marks, dots);
}

int main()
{
    try {
        MEDDLY::initialize();
        runWithArgs(100, 256*8, 256);
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
    std::cerr << "\nSome other error?\n";
    return 3;
}


