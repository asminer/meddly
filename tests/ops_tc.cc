/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2009, Iowa State University Research Foundation, Inc.

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

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"
#include <cassert>

/*
    Simple and adjustable PN model for testing
    transitive closure.

    Petri net:

      p1      t1     p2      t2     p3     t3      p4
      (N) --> [] --> ( ) --> [] --> ( ) --> [] --> ( )

          <-- [] <--                    <-- [] <--
              t5                            t6

    Transitions t5 and t6 may be removed.
*/
class smallPN : public rsgen_model {
        int N;
        bool allow_t5;
        bool allow_t6;
    public:
        smallPN(int N, bool t5, bool t6);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

        void test_tc(MEDDLY::output &out);
};

// =============================================================

smallPN::smallPN(int _N, bool t5, bool t6)
{
    assert(_N>=0 );
    N = _N;
    allow_t5 = t5;
    allow_t6 = t6;
}

// =============================================================

bool smallPN::build_model(MEDDLY::output &out)
{
    out << "Building model:\n";
    out << "    N=" << N << "\n";
    out << "    t5: " << (allow_t5 ? "enabled" : "disabled") << "\n";
    out << "    t6: " << (allow_t6 ? "enabled" : "disabled") << "\n";

    //
    // Initialize domain
    //
    int sizes[4];
    sizes[0] = sizes[1] = sizes[2] = sizes[3] = N+1;
    set_domain(MEDDLY::domain::createBottomUp(sizes, 4));

    //
    // Make the events
    //
    add_event(new simple_event("t1",    "X-+.."));
    add_event(new simple_event("t2",    "X.-+."));
    add_event(new simple_event("t3",    "X..-+"));

    if (allow_t5) {
        add_event(new simple_event("t5",    "X+-.."));
    }
    if (allow_t6) {
        add_event(new simple_event("t3",    "X..+-"));
    }
    return true;
}

// =============================================================

void smallPN::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(0);
    init.setVar(1, N);
}


// =============================================================

void smallPN::test_tc(MEDDLY::output &out)
{
    using namespace MEDDLY;

    forest* mdd = forest::create(get_domain(), SET, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);

    forest* mxd = forest::create(get_domain(), RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL);

    dd_edge constr(mdd);
    mdd->createConstant(true, constr);

    dd_edge nsf(mxd), tc(mxd);
    build_monolithic(nsf);

    // apply(TRANSITIVE_CLOSURE_DFS, constr, nsf, nsf, tc);

    // out << "Transitive closure:\n";
    // tc.showRelation(out);

    // TBD: test it

    out << "Done\n";

    // cleanup
    forest::destroy(mdd);
    forest::destroy(mxd);
}

// =============================================================

int main()
{
    try {
        MEDDLY::initialize();
        MEDDLY::ostream_output out(std::cout);

        smallPN m1(3, false, false);
        m1.build_model(out);
        m1.test_tc(out);

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

