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
#include <string.h>

#include "kan_base.h"

// ============================================

bool equal(const MEDDLY::minterm &a, const MEDDLY::minterm &b)
{
    if (a.getNumVars() != b.getNumVars()) return false;
    for (unsigned k = a.getNumVars(); k; --k)
    {
        if (a.from(k) != b.from(k)) return false;
    }
    return true;
}

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
    set_gen_method('t');
}

// ============================================

void kanban::after_generation(MEDDLY::output &out,
                const MEDDLY::dd_edge &reachable)
{
    using namespace MEDDLY;

    out << "Running test for N=" << N << "\n";
    //
    // Build index set for reachable states
    //
    forest* evmdd = forest::create(get_domain(), SET, range_type::INTEGER,
                        edge_labeling::INDEX_SET);

    dd_edge reach_index(evmdd);
    apply(MEDDLY::CONVERT_TO_INDEX_SET, reachable, reach_index);
    out << "\tbuilt index set\n";
    out.flush();

    //
    // Iterate through reachable states, and make sure
    // the reach_index matches.
    //
    long c = 0;
    for (dd_edge::iterator s = reachable.begin(); s; ++s)
    {
        long index;
        reach_index.evaluate(*s, index);
        if (c != index) {
            out << "\nState number " << c << " has index " << index << "\n";
            out << "  state: ";
            (*s).show(out);
            out << "\n";
            throw "mismatch";
        }
        c++;
    }
    out << "\tverified set -> index\n";
    out.flush();

    //
    // Iterate through reach_index, and make sure all indexed
    // states are in fact reachable.
    //
    c = 0;
    for (dd_edge::iterator s = reach_index.begin(); s; ++s)
    {
        long index = (*s).getValue();
        if (index != c) {
            out << "\nState number " << c << " has index " << index << "\n";
            throw "mismatch";
        }
        bool ok;
        reachable.evaluate(*s, ok);
        if (!ok) {
            out << "\nIndex number " << c
                << "does not appear in reachability set\n";
            throw "mismatch";
        }
        c++;
    }
    out << "\tverified index -> set\n";
    out.flush();

    //
    // Verify index search
    //
    c = 0;
    minterm assign(get_domain(), MEDDLY::SET);
    for (dd_edge::iterator s = reachable.begin(); s; ++s)
    {
        bool ok = reach_index.getElement(c, assign);
        if (!ok) {
            out << "\nCan't find state with index " << c << "\n";
            throw "not found";
        }
        if (!equal(*s, assign)) {
            out << "\nFetch index " << c << " got wrong state\n";
            throw "mismatch";
        }
        c++;
    }
    out << "\tverified index search\n";
    out << "\tdone verifying " << c << " state indexes\n";
    out.flush();
}

// ============================================

int main()
{
    kanban K;

    try {
        MEDDLY::initialize();

        MEDDLY::ostream_output out(std::cout);

        for (unsigned n=1; n<5; n++) {
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
