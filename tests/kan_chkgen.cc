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

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"


struct statedist {
    const char* state;
    unsigned distance;
};

#define HAVE_ORACLE

#ifdef HAVE_ORACLE
#include "kan_fwd_df1.h"
#include "kan_fwd_df2.h"
#include "kan_fwd_df3.h"
#include "kan_bwd_df1.h"
#include "kan_bwd_df2.h"
#include "kan_bwd_df3.h"
#endif

// ============================================

// Adjust multi-terminal, integer distance function
// to be distance plus one, for reachable states,
// zero for unreachable states.

void distance_adjust(const MEDDLY::rangeval &in, MEDDLY::rangeval &out)
{
    long x = long(in);
    if (x<0) {
        out = 0;
    } else {
        out = x+1;
    }
}

// ============================================

bool matches(const long adjust, const statedist &sd, const MEDDLY::minterm &m)
{
    if (!sd.state) {
        throw "not enough reachable markings";
    }

    for (unsigned k=1; k <= 16; ++k) {
        if (sd.state[k]-48 != m.from(k)) return false;
    }

    const MEDDLY::rangeval& oracle = m.getValue();
    if (oracle.isBoolean()) return true;

    if (oracle.isInteger()) {
        return long(oracle) == sd.distance + adjust;
    }

    throw "Unexpected state value type";
}

// ============================================

class kanban : public rsgen_model {
    public:
        kanban();

        inline void setN(int n) {
            N = n;
        }

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

        virtual void after_generation(MEDDLY::output &out,
                const MEDDLY::dd_edge &e);

    private:
        int N;

};

// ============================================

kanban::kanban()
{
    N = -1;
    disable_output_options();
}

// ============================================

bool kanban::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    // TBD
    return false;
}

// ============================================

void kanban::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options]\n\n";
    }
    // TBD
    rsgen_model::show_usage(out);
}

// ============================================

bool kanban::build_model(MEDDLY::output &out)
{
    if (N<0) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    out << "Checking Kanban " << get_set_name() << ", N=" << N << "\n";

    //
    // Initialize domain
    //
    int* sizes = new int[16];
    for (int i=15; i>=0; i--) sizes[i] = N+1;
    set_domain(MEDDLY::domain::createBottomUp(sizes, 16));
    delete[] sizes;

    //
    // Make the events
    //
    add_event(new simple_event("Tin1",    "X-+.............."));
    add_event(new simple_event("Tr1",     "X.-+............."));
    add_event(new simple_event("Tb1",     "X.+-............."));
    add_event(new simple_event("Tg1",     "X.-.+............"));
    add_event(new simple_event("Tr2",     "X.....-+........."));
    add_event(new simple_event("Tb2",     "X.....+-........."));
    add_event(new simple_event("Tg2",     "X.....-.+........"));
    add_event(new simple_event("Ts1_23",  "X+..--+..-+......"));
    add_event(new simple_event("Tr3",     "X.........-+....."));
    add_event(new simple_event("Tb3",     "X.........+-....."));
    add_event(new simple_event("Tg3",     "X.........-.+...."));
    add_event(new simple_event("Ts23_4",  "X....+..-+..--+.."));
    add_event(new simple_event("Tr4",     "X.............-+."));
    add_event(new simple_event("Tb4",     "X.............+-."));
    add_event(new simple_event("Tout4",   "X............+..-"));
    add_event(new simple_event("Tg4",     "X.............-.+"));

    return true;
}

// ============================================

void kanban::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(0);
    init.setVar(1, N);
    init.setVar(5, N);
    init.setVar(9, N);
    init.setVar(13, N);
}

// ============================================

void kanban::after_generation(MEDDLY::output &out,
                const MEDDLY::dd_edge &e)
{
    const statedist* krs = nullptr;
#ifdef HAVE_ORACLE
    switch (N) {
        case 1:
            krs = is_forward() ? kanban_fwd_rs1 : kanban_bwd_rs1;
            break;

        case 2:
            krs = is_forward() ? kanban_fwd_rs2 : kanban_bwd_rs2;
            break;

        case 3:
            krs = is_forward() ? kanban_fwd_rs3 : kanban_bwd_rs3;
            break;

        default:
            krs = nullptr;
    }
#endif

    if (krs) {
        out << "Verifying against oracle...\n";

        // enumerate states
        long c = 0;

        long adjust = 0;

        MEDDLY::dd_edge reachable(e);
        if ('i' == get_set_type()) {
            // need to add one (for efficient enumeration)
            // for integer multi-terminal distances
            MEDDLY::user_unary_factory plus1("distadjust", distance_adjust);
            apply(plus1, reachable, reachable);
            adjust = 1;
        }

        for (MEDDLY::dd_edge::iterator i = reachable.begin(); i; ++i)
        {
            if (!matches( adjust, krs[c], *i) ) {
                out << "Marking " << c << " mismatched\n";
                out << "Marking was ";
                (*i).show(out);
                out << "\n";

                throw "mismatch";
            }
            c++;
        }
        if (krs[c].state) {
            throw "not enough markings";
        }

        out << "\tperfect match!\n";
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

        out.indent_more();
        K.setN(1);
        K.build_model(out);
        K.generate(out, false, false);
        K.destroy_model();
        out.indent_less();
        out << "\n";

        out.indent_more();
        K.setN(2);
        K.build_model(out);
        K.generate(out, false, false);
        K.destroy_model();
        out.indent_less();
        out << "\n";

        out.indent_more();
        K.setN(3);
        K.build_model(out);
        K.generate(out, false, false);
        K.destroy_model();
        out.indent_less();
        out << "\n";

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

