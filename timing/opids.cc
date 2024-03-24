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

#include <iostream>
#include <iomanip>
#include <vector>

#include "../src/meddly.h"
#include "timer.h"
#include "reporting.h"
#include "park_random.h"

using namespace MEDDLY;

const unsigned FORESTS = 16;
const unsigned OPCOMBS = FORESTS * FORESTS * FORESTS;
const unsigned DOTS    = 16;

#ifdef DEVELOPMENT_CODE
const unsigned CREATIONS = 2;
const unsigned FASTSRCH  = 4;
const unsigned SLOWSRCH  = 4;
#else
const unsigned CREATIONS = 2;
const unsigned FASTSRCH  = 4;
const unsigned SLOWSRCH  = 4;
#endif

void createOperations(const std::vector <dd_edge> &Edges,
        std::vector <operation*> &ops)
{
    unsigned u=0;
    binary_opname* Inter = INTERSECTION();

    for (unsigned i=0; i<Edges.size(); i++) {
        for (unsigned j=0; j<Edges.size(); j++) {
            for (unsigned k=0; k<Edges.size(); k++) {
                ops[u++] = Inter->getOperation(Edges[i], Edges[j], Edges[k]);
            } // for k
        } // for j
    } // for i
    if (u != ops.size()) {
        throw "size mismatch in createOperations";
    }
}

void destroyOperations(std::vector <operation*> &ops)
{
    for (unsigned u=0; u<ops.size(); u++) {
        destroyOperation(ops[u]);
    }
}

void searchFast(const std::vector <dd_edge> &Edges)
{
    const unsigned numops = operation::getOpListSize();
    binary_opname* Inter = INTERSECTION();

    for (unsigned i=0; i<Edges.size(); i++) {
        for (unsigned j=0; j<Edges.size(); j++) {
            for (unsigned k=0; k<Edges.size(); k++) {
                for (unsigned s=0; s<FASTSRCH; s++) {
                    Inter->getOperation(Edges[i], Edges[j], Edges[k]);
                }
            } // for k
        } // for j
    } // for i

    if (operation::getOpListSize() != numops) {
        throw "some operation was not found in searchFast";
    }
}

void searchSlow(const std::vector <dd_edge> &Edges)
{
    const unsigned numops = operation::getOpListSize();
    binary_opname* Inter = INTERSECTION();

    for (unsigned s=0; s<SLOWSRCH; s++) {
        for (unsigned i=0; i<Edges.size(); i++) {
            for (unsigned j=0; j<Edges.size(); j++) {
                for (unsigned k=0; k<Edges.size(); k++) {
                    Inter->getOperation(Edges[i], Edges[j], Edges[k]);
                } // for k
            } // for j
        } // for i
    } // for s

    if (operation::getOpListSize() != numops) {
        throw "some operation was not found in searchSlow";
    }
}

int main(int argc, const char** argv)
{
    try {
        setReport(argc, argv);

        MEDDLY::initialize();

        SeedRNG(123456789);

        //
        // Build a small domain
        //
        int bounds[3];
        bounds[0] = bounds[1] = bounds[2] = 5;
        domain *D = domain::createBottomUp(bounds, 3);

        //
        // Build several identical forests,
        // just so we can have different operations
        //
        std::vector <forest*> F(FORESTS);
        std::vector <dd_edge> E(FORESTS);
        std::vector <operation*> Ops(OPCOMBS);
        policies p(SET);

        for (unsigned i=0; i<FORESTS; i++) {
            F[i] = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);
            E[i].attach(F[i]);
        }

        //
        // Create/destroy test
        //
        std::cout << "\nTiming test: operation create/destroy\n\n";
        std::cout << "Building/destroying " << DOTS*CREATIONS*OPCOMBS << " ops ";

        timer T;
        for (unsigned i=0; i<DOTS; i++) {
            std::cout << '.';
            std::cout.flush();
            for (unsigned j=0; j<CREATIONS; j++) {
                createOperations(E, Ops);
                destroyOperations(Ops);
            }
        } // for i
        T.note_time();
        show_sec(std::cout, T, 3, 3);
        if (startReport(T, __FILE__)) {
            report  << "opbuild $ "
                    << "Created and destroyed " << DOTS*CREATIONS*OPCOMBS
                    << " operations"
                    << std::endl;
        }

        createOperations(E, Ops);

        // -------------------------------------------------------------------

        std::cout << "Fast search         " << DOTS*FASTSRCH*OPCOMBS << " ops ";
        T.note_time();
        for (unsigned i=0; i<DOTS; i++) {
            std::cout << '.';
            std::cout.flush();
            searchFast(E);
        }
        T.note_time();
        show_sec(std::cout, T, 3, 3);
        if (startReport(T, __FILE__)) {
            report  << "fastsrch $ "
                    << "Fast search " << DOTS*FASTSRCH*OPCOMBS
                    << " operations"
                    << std::endl;
        }

        // -------------------------------------------------------------------

        // -------------------------------------------------------------------

        std::cout << "Slow search         " << DOTS*SLOWSRCH*OPCOMBS << " ops ";
        T.note_time();
        for (unsigned i=0; i<DOTS; i++) {
            std::cout << '.';
            std::cout.flush();
            searchSlow(E);
        }
        T.note_time();
        show_sec(std::cout, T, 3, 3);
        if (startReport(T, __FILE__)) {
            report  << "slowsrch $ "
                    << "Slow search " << DOTS*SLOWSRCH*OPCOMBS
                    << " operations"
                    << std::endl;
        }

        // -------------------------------------------------------------------

        std::cout << std::endl;

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
