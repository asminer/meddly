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

#include "../src/meddly.h"
#include "timer.h"
#include "reporting.h"
#include "park_random.h"

using namespace MEDDLY;

const unsigned FORESTS = 1024;
const unsigned FTYPES = 9;

#ifdef DEVELOPMENT_CODE
const unsigned CREATIONS = 32 * 1024;
#else
const unsigned CREATIONS = 1024 * 1024;
#endif

int main(int argc, const char** argv)
{
    try {
        setReport(argc, argv);

        MEDDLY::initialize();

        //
        // Build array of forest combinations
        //
        //  0: set, boolean, MT
        //  1: set, integer, MT
        //  2: set, integer, ev+
        //  3: set, real, MT
        //
        //  4: relation, boolean, MT
        //  5: relation, integer, MT
        //  6: relation, integer, ev+
        //  7: relation, real, MT
        //  8: relation, real, ev*
        //
        policies p[FTYPES];
        set_or_rel sr[FTYPES];
        range_type rt[FTYPES];
        edge_labeling el[FTYPES];

        for (unsigned i=0; i<FTYPES; i++) {
            p[i].useDefaults(i>=4);
            sr[i] = i<4 ? SET : RELATION;
            el[i] = edge_labeling::MULTI_TERMINAL;  // overwrite later if needed
        }

        rt[0] = rt[4] = range_type::BOOLEAN;
        rt[1] = rt[2] = rt[5] = rt[6] = range_type::INTEGER;
        rt[3] = rt[7] = rt[8] = range_type::REAL;

        el[2] = el[6] = edge_labeling::EVPLUS;
        el[8] = edge_labeling::EVTIMES;

        SeedRNG(123456789);

        //
        // Build a small domain
        //
        int bounds[3];
        bounds[0] = bounds[1] = bounds[2] = 5;
        domain *D = domain::createBottomUp(bounds, 3);

        //
        // Build a bunch of forests, at random, holding
        // a bunch at a time.
        //
        forest* F[FORESTS];
        for (unsigned i=0; i<FORESTS; i++) {
            F[i] = nullptr;
        }

        std::cout << "\nTiming test: forest creation / destruction\n\n";
        std::cout << "Building " << 16*CREATIONS << " forests ";

        timer T;
        unsigned ft = 0;
        for (unsigned i=0; i<16; i++) {
            std::cout << '.';
            std::cout.flush();
            for (unsigned j=0; j<CREATIONS; j++) {
                unsigned u = Equilikely(0, FORESTS-1);

                forest::destroy(F[u]);
                F[u] = forest::create(D, sr[ft], rt[ft], el[ft], p[ft]);

                ft = (ft + 1 ) % FTYPES;
            } // for j
        } // for i

        for (unsigned u=0; u<FORESTS; u++) {
            forest::destroy(F[u]);
        }

        T.note_time();
        show_sec(std::cout, T, 3, 3);

        if (startReport(T, __FILE__, "fbuild")) {
            report  << "Created and destroyed " << 16*CREATIONS << " forests "
                    << std::endl;
        }

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
