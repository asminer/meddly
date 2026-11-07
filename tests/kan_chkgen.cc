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

bool matches(const long adjust, const MEDDLY::minterm &morc,
        const MEDDLY::minterm &m)
{
    for (unsigned k=1; k <= m.getNumVars(); ++k) {
        if (morc.from(k) != m.from(k)) return false;
    }
    const MEDDLY::rangeval& mv = m.getValue();
    if (mv.isBoolean()) {
        return true;
    }
    if (mv.isInteger()) {
        return long(morc.getValue()) + adjust == long(mv);
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

    protected:
        virtual void write_oracle_header(MEDDLY::output &out) const;
        virtual void read_oracle_header(MEDDLY::input &inf) const;

    private:
        int N;
        bool generate_oracle;

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

    if (0==strcmp("-g", argv[i])) {
        generate_oracle = true;
        ++i;
        return true;
    }
    if (0==strcmp("-c", argv[i])) {
        generate_oracle = false;
        ++i;
        return true;
    }


    return false;
}

// ============================================

void kanban::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options]\n";
    }
    out << "\n";
    out << "    MAIN OPTIONS\n";
    out << "\n";
    out << "        -c         : Check answer against an oracle (default).\n";
    out << "        -g         : Generate an oracle.\n";

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
    if (N>5) {
        throw "N is too large for an oracle";
    }
    char fname[32];
    strlcpy(fname, "oracle_kanban_x_d_N.txt", sizeof(fname));
    // replace x
    fname[14] = is_forward() ? 'f' : 'b';
    // replace N
    fname[18] = '0' + N;

    if (generate_oracle) {
        if ('+' != get_set_type()) {
            out << "Use -evdist for generating oracles!\n";
            throw "Bad settings for oracle generation";
        }

        out << "Generating oracle " << fname << "\n";
        std::ofstream outf(fname);

        outf << "#\n# Oracle for Kanban N=" << N << ", "
             << (is_forward() ? "forward" : "backward")
             << "\n#\n";

        MEDDLY::ostream_output mout(outf);
        write_oracle(mout, e);
        return;
    }

    char pname[48];
    strlcpy(pname, "../oracles/", sizeof(pname));
    strlcat(pname, fname, sizeof(pname));

    out << "Reading oracle file " << pname << "\n";
    std::ifstream inf(pname);
    if (!inf) {
        out << "    Missing oracle file!\n";
        out << "    Generate using -evdist -g, and move into the oracles directory\n";
        throw "Missing oracle";
    }
    MEDDLY::istream_input minf(inf);
    MEDDLY::minterm_coll MC(16, e.getForest());
    read_oracle(minf, MC);

    out << "Verifying against oracle\n";

    long c = 0; // state counter
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
        if (!matches( adjust, MC.at(c), *i)) {
            out << "\tMarking #" << c << " mismatched\n";
            out << "\tOracle: ";
            MC.at(c).show(out);
            out << "\n";
            out << "\tOurs  : ";
            (*i).show(out);
            out << "\n";
            throw "mismatch";
        }
        c++;
    } // for loop

    if (c < MC.size()) {
        out << "\tNot enough markings!\n";
        throw "mismatch";
    }

    out << "\tperfect match!\n";
}

// ============================================

void kanban::write_oracle_header(MEDDLY::output &out) const
{
    out << "kanban " << N << "\n";
    rsgen_model::write_oracle_header(out);
}

// ============================================

void kanban::read_oracle_header(MEDDLY::input &inf) const
{
    inf.stripWS();
    inf.consumeKeyword("kanban");
    if (inf.get_integer() != N) {
        throw "oracle file mismatch";
    }
    rsgen_model::read_oracle_header(inf);
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

