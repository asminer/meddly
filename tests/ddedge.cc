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

#include "../src/meddly.h"

#include <iostream>
#include <iomanip>

// #define REPORTING

using namespace MEDDLY;

const unsigned EDGES = 1024;
const unsigned CHANGES = 1024 * 1024;

long seed = -1;

inline double Random()
{
    static const long MODULUS = 2147483647L;
    static const long MULTIPLIER = 48271L;
    static const long Q = MODULUS / MULTIPLIER;
    static const long R = MODULUS % MULTIPLIER;

    long t = MULTIPLIER * (seed % Q) - R * (seed / Q);
    if (t > 0) {
        seed = t;
    } else {
        seed = t + MODULUS;
    }
    return ((double) seed / MODULUS);
}

inline unsigned Equilikely(unsigned a, unsigned b)
{
    return (a + (unsigned) ((b - a + 1) * Random()));
}


void test_edge_registry(forest* F1, forest* F2)
{
    std::cout << "Creating array of " << EDGES << " edges\n";

    dd_edge* E = new dd_edge[EDGES];

    unsigned nocurr = EDGES;
    unsigned f1curr = 0;
    unsigned f2curr = 0;
    unsigned f1max = 0;
    unsigned f2max = 0;

    std::cout << "Updating array " << CHANGES*16 << " times ";

    for (unsigned i=0; i<16; i++) {
        std::cout << '.';
        std::cout.flush();
        for (unsigned j=0; j<CHANGES; j++) {
            unsigned u = Equilikely(0, EDGES-1);

            if (nullptr == E[u].getForest()) {
                // switch from none to F1
                E[u].attach(F1);
                --nocurr;
                ++f1curr;
                if (f1curr > f1max) {
                    f1max = f1curr;
                }
                continue;
            }

            if (F1 == E[u].getForest()) {
                // switch from F1 to F2
                E[u].attach(F2);
                --f1curr;
                ++f2curr;
                if (f2curr > f2max) {
                    f2max = f2curr;
                }
                continue;
            }

            if (F2 == E[u].getForest()) {
                // switch from F2 to null
                E[u].detach();
                --f2curr;
                ++nocurr;
                continue;
            }

            throw "unknown forest";
        }
    }

    std::cout << "\n\n";
    std::cout << "Max edges in forest 1: " << f1max << "\n";
    std::cout << "Max edges in forest 2: " << f2max << "\n\n";

    std::cout << "Current null edges: " << nocurr << "\n";
    std::cout << "Current  F1  edges: " << f1curr << "\n";
    std::cout << "Current  F2  edges: " << f2curr << "\n";
    std::cout << "                    " << nocurr + f1curr + f2curr << "\n\n";

    const unsigned F1count = F1->countRegisteredEdges();
    const unsigned F2count = F2->countRegisteredEdges();

    std::cout << "F1 reports " << F1count << " edges\n";
    std::cout << "F2 reports " << F2count << " edges\n\n";

    if (F1count != f1curr) throw "edge count mismatch";
    if (F2count != f2curr) throw "edge count mismatch";

    std::cout << "Deleting all edges\n\n";
    delete[] E;
}



int main()
{
    try {
        MEDDLY::initialize();

        seed = 123456789;

        int bounds[3];
        bounds[0] = bounds[1] = bounds[2] = 5;
        domain *D = domain::createBottomUp(bounds, 3);
        forest *F1, *F2;
        policies p;
        p.useDefaults(false);

        F1 = forest::create(D, false, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);
        F2 = forest::create(D, false, range_type::INTEGER,
                        edge_labeling::EVPLUS, p);

        test_edge_registry(F1, F2);

        std::cout << "F1 reports " << F1->countRegisteredEdges() << " edges\n";
        std::cout << "F2 reports " << F2->countRegisteredEdges() << " edges\n\n";

        if (F1->countRegisteredEdges()) throw "F1 count mismatch";
        if (F2->countRegisteredEdges()) throw "F2 count mismatch";

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
