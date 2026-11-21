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
#include "kan_base.h"

// ============================================

class kanban : public kanban_base {
    public:
        kanban();

        inline void fill_initial(MEDDLY::dd_edge &init_state)
        {
            MEDDLY::minterm init(mdd);
            build_initial(init);
            init.buildFunction(false, init_state);
        }

        void check_states(MEDDLY::output &out);

    protected:
        virtual void read_oracle_header(MEDDLY::input &inf) const;

        void make_forests();
        void done_forests();

        void compare_states(MEDDLY::output &out, int dist,
                const MEDDLY::dd_edge &ftr) const;

        //
        // Build next forward or backward step.
        //
        // On input:
        //     reachable: all states reachable up to current distance d
        //     frontier:  ignored
        //
        // On output:
        //      reachable: all states reachable up to distance d+1
        //      frontier:  states at distance d+1
        //
        inline void next_step(MEDDLY::dd_edge &reachable,
            const MEDDLY::dd_edge &nsf, MEDDLY::dd_edge &frontier) const
        {
            using namespace MEDDLY;
            dd_edge s_x_R(reachable);

            if (is_forward()) {
                apply(POST_IMAGE, reachable, nsf, s_x_R);
            } else {
                apply(PRE_IMAGE, reachable, nsf, s_x_R);
            }
            apply(DIFFERENCE, s_x_R, reachable, frontier);
            apply(UNION, reachable, s_x_R, reachable);
        }



    private:
        MEDDLY::forest* mdd;
        MEDDLY::forest* mxd;
        MEDDLY::forest* evmdd;
};

// ============================================

kanban::kanban()
{
    disable_output_options();
    mdd = nullptr;
    mxd = nullptr;
    evmdd = nullptr;
}

// ============================================

void kanban::check_states(MEDDLY::output &out)
{
    using namespace MEDDLY;

    make_forests();
    //
    // Build next-state function
    //
    dd_edge nsf(mxd);
    build_monolithic(nsf);

    //
    // Build initial state, frontier sets
    //
    dd_edge reachset(mdd), frontier(mdd);

    // Classic BFS-style iteration, forward

    for (int iter=1; iter >= 0; --iter) {
        fill_initial(reachset);
        fill_initial(frontier);
        out.indent_more();
        set_forward(iter);
        out << ( is_forward() ? "Forward\n" : "Backward\n" );
        int dist = 0;
        do {
            out << "Verifying distance ";
            out.put(dist, 3);
            out << " ... ";
            out.flush();
            compare_states(out, dist, frontier);
            next_step(reachset, nsf, frontier);
            out << " match!\n";
            ++dist;
        } while (frontier.getNode());
        out.indent_less();
        out << "\n";
    }

    done_forests();
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

void kanban::make_forests()
{
    using namespace MEDDLY;
    mdd = forest::create(get_domain(), SET, range_type::BOOLEAN,
                            edge_labeling::MULTI_TERMINAL);
    mxd = forest::create(get_domain(), RELATION, range_type::BOOLEAN,
                            edge_labeling::MULTI_TERMINAL);
    evmdd = forest::create(get_domain(), SET, range_type::INTEGER,
                            edge_labeling::INDEX_SET);
}

// ============================================

void kanban::done_forests()
{
    using namespace MEDDLY;
    forest::destroy(evmdd);
    forest::destroy(mxd);
    forest::destroy(mdd);
}

// ============================================

void kanban::compare_states(MEDDLY::output &out, int dist,
        const MEDDLY::dd_edge &ftr) const
{
    std::string pname = "../oracles/oracle_kanban_";
    if (is_forward()) {
        pname += "f_d_";
    } else {
        pname += "b_d_";
    }
    pname += '0' + N;
    pname += ".txt";

    std::ifstream inf(pname);
    if (!inf) {
        out << "\n";
        out << "    Missing oracle file!\n";
        out << "    Use kan_chkgen -evdist -g to generate\n";
        throw "Missing oracle";
    }

    //
    // Load oracle but just states with given distance
    //
    MEDDLY::istream_input minf(inf);
    MEDDLY::minterm_coll MC(16, mdd);
    read_oracle(minf, MC, dist);

    //
    // Compare with ftr
    //
    verify_set(out, MC, ftr);
}

// ============================================

/*
    Main program, of course
*/
int main()
{
    MEDDLY::initialize();
    kanban K;

    MEDDLY::ostream_output out(std::cout);

    try {
        for (unsigned n=1; n<4; n++) {
            out.indent_more();
            K.setN(n);
            K.build_model(out);
            K.check_states(out);
            K.destroy_model();
            out.indent_less();
            out << "\n";
        }
        out << "Done\n";
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

