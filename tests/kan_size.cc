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

#include "kan_base.h"

// ============================================

long expected[] = {
    1, 160, 4600, 58400, 454475, 2546432, 11261376,
    41644800, 133865325, 384392800, 1005927208
};

const int nstart = 1;
const int nstop = 10;

// ============================================

class kanban : public kanban_base {
    public:
        kanban();

        inline void setN(int n) {
            N = n;
        }

        virtual void after_generation(MEDDLY::output &out,
                const MEDDLY::dd_edge &e);
};

// ============================================

kanban::kanban()
{
    disable_generation_types();
    disable_output_options();
}

// ============================================

void kanban::after_generation(MEDDLY::output &out,
                const MEDDLY::dd_edge &e)
{
    if (N>nstop) {
        throw "N is too large for an oracle";
    }

    long c = -1;
    MEDDLY::apply(MEDDLY::CARDINALITY, e, c);
    out.put(c, 12);
    out << " states\n";
    if (c != expected[N]) {
        out << "\tWrong number of states!\n";
        throw "size mismatch";
    }
}

// ============================================

int main(int argc, const char** argv)
{
    kanban K;

    unsigned nArgs = (unsigned)argc;
    unsigned i=1;

    if (1 == argc) {
        K.show_usage(std::cerr, argv[0]);
        return 1;
    }

    while (i < nArgs) {
        if (!K.process_arg(i, argv)) {
            K.show_usage(std::cerr, argv[0]);
            return 1;
        }
    }

    try {
        MEDDLY::initialize();

        MEDDLY::ostream_output out(std::cout);

        for (unsigned n=nstart; n<=nstop; n++) {
            out.indent_more();
            K.setN(n);
            K.build_model(out);
            K.generate(out, false, false);
            K.destroy_model();
            out.indent_less();
            out << "\n";
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
        std::cerr   << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "\nSome other error?\n";
    return 4;

}

