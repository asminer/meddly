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

#include <iostream>
#include <cassert>

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"

class hanoi : public rsgen_model {
        class move_ring : public rsgen_model::event {
                const std::vector <int> &ring2level;
                int ring;
                int src;
                int dest;
                int aux;
            public:
                // **********************************************************
                // build an event to move a ring
                //    @param ring2level     Mapping from rings to level
                //    @param ring           Which ring to move
                //    @param src            Source peg (0, 1, 2)
                //    @param dest           Destination peg (0, 1, 2)
                //    @param aux            Unused peg (0, 1, 2)
                // **********************************************************
                move_ring(const std::vector <int> &ring2level, int ring,
                    int src, int dest, int aux);

                virtual bool fire(const MEDDLY::minterm &curr,
                        MEDDLY::minterm &);

                virtual bool build_nsf(MEDDLY::dd_edge &R);
        };
    public:
        hanoi();

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

        virtual void after_generation(MEDDLY::output &out,
                    const MEDDLY::dd_edge &e);

    private:
        int N;
        bool reverse_order;
        std::vector <int> ring2level;
};

// ******************************************************************
// *                    hanoi::move_ring methods                    *
// ******************************************************************

hanoi::move_ring::move_ring(const std::vector <int> &r2lev,
        int _ring, int _src, int _dest, int _aux)
    : event("move"), ring2level(r2lev)
{
    ring = _ring;
    src = _src;
    dest = _dest;
    aux = _aux;

    assert(src != dest);
    assert(src != aux);
    assert(dest != aux);

    assert(ring > 0);
    assert(src >= 0);
    assert(dest >= 0);
    assert(aux >= 0);
    assert(src < 3);
    assert(dest < 3);
    assert(aux < 3);
}

// ******************************************************************

bool hanoi::move_ring::fire(const MEDDLY::minterm &curr,
        MEDDLY::minterm &next)
{
    for (unsigned i = curr.getNumVars(); i; --i) {
        const int k = ring2level[i];
        if (i > unsigned(ring)) {
            next.setVar(k, curr.getVar(k));
            continue;
        }
        if (i < unsigned(ring)) {
            if (curr.getVar(k) != aux) return false;
            next.setVar(k, aux);
            continue;
        }
        if (curr.getVar(k) != src) return false;
        next.setVar(k, dest);
    }
    return true;
}

// ******************************************************************

bool hanoi::move_ring::build_nsf(MEDDLY::dd_edge &R)
{
    MEDDLY::minterm delta(R.getForest());

    for (unsigned i = delta.getNumVars(); i; --i) {
        const int k = ring2level[i];
        if (i > unsigned(ring)) {
            delta.setVars(k, MEDDLY::DONT_CARE, MEDDLY::DONT_CHANGE);
            continue;
        }
        if (i < unsigned(ring)) {
            delta.setVars(k, aux, aux);
            continue;
        }
        delta.setVars(k, src, dest);
    }

    delta.setValue(true);
    delta.buildFunction(false, R);
    return true;
}

// ******************************************************************
// *                         hanoi  methods                         *
// ******************************************************************

hanoi::hanoi()
{
    N = -1;
    reverse_order = false;
    default_approx();
}

// ******************************************************************

bool hanoi::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    if (0 == strcmp("-o", argv[i])) {
        reverse_order = false;
        ++i;
        return true;
    }
    if (0 == strcmp("-O", argv[i])) {
        reverse_order = true;
        ++i;
        return true;
    }

    if ('-' == argv[i][0])  return false;

    // assume this is the size
    N = atoi(argv[i]);
    ++i;
    return true;
}

// ******************************************************************

void hanoi::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options] rings\n\n";
    }
    out << "\trings  : number of rings\n";
    out << "\n";
    out << "    ORDERING OPTIONS\n";
    out << "\n";
    out << "        -o         : Order the variables so that the smallest ring is the\n";
    out << "                     bottom-most variable in the MDD (default).\n";
    out << "        -O         : Order the variables so that the largest ring is the\n";
    out << "                     bottom-most variable in the MDD.\n";

    rsgen_model::show_usage(out);
}

// ******************************************************************

bool hanoi::build_model(MEDDLY::output &out)
{
    if (N<=0) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }

    out << "+-------------------------------------------+\n";
    out << "|      Towers of Hanoi with ";
    out.put(N, 4);
    out << " rings      |\n";
    out << "+-------------------------------------------+\n";

    //
    // Initialize domain
    //
    int* sizes = new int[N];
    for (int i=N-1; i>=0; i--) sizes[i] = 3;
    set_domain(MEDDLY::domain::createBottomUp(sizes, N));
    delete[] sizes;

    //
    // Initialize level mappings
    //
    ring2level.resize(N+1);
    std::vector <int> level2ring(N+1);
    ring2level[0] = 0;
    level2ring[0] = 0;
    if (reverse_order) {
        for (unsigned i=N; i>0; i--) {
            ring2level[i] = N-i+1;
            level2ring[N-i+1] = i;
        }
    } else {
        for (unsigned i=N; i>0; i--) {
            ring2level[i] =i;
            level2ring[i] = i;
        }
    }

    out << "Variable order is TOP, ";
    for (unsigned i=N; i>0; i--) {
        out << "ring " << level2ring[i] << ", ";
        if (i>3) {
            out << "..., ";
            i = 3;
        }
    }
    out << "BOTTOM\n";


    //
    // Make the events
    //
    for (int ring=1; ring<=N; ring++) {
        add_event(new move_ring(ring2level, ring, 0, 1, 2));
        add_event(new move_ring(ring2level, ring, 0, 2, 1));
        add_event(new move_ring(ring2level, ring, 1, 0, 2));
        add_event(new move_ring(ring2level, ring, 1, 2, 0));
        add_event(new move_ring(ring2level, ring, 2, 0, 1));
        add_event(new move_ring(ring2level, ring, 2, 1, 0));
    } // for ring
    return true;
}

// ******************************************************************

void hanoi::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(0);
}

// ******************************************************************

void hanoi::after_generation(MEDDLY::output &out,
        const MEDDLY::dd_edge &e)
{
    if ('b' != get_set_type()) return;

    show_cardinality(out, e, false);

    if (approx_cardinality()) {
        double theory = 1.0;
        for (int i=0; i<N; i++) {
            theory *= 3.0;
        }
        out << "(Expecting    ";
        out.put(theory, 10, 6);
        out << " states)\n";
    } else {
#ifdef HAVE_LIBGMP
        mpz_t theory;
        mpz_init_set_ui(theory, 1);
        for (int i=0; i<N; i++) {
            mpz_mul_ui(theory, theory, 3);
        }

        out << "(expect ";
        out.put_huge(theory, true);
        out << " states)\n";

        mpz_clear(theory);
#endif
    }
}


// ******************************************************************
// *                          main program                          *
// ******************************************************************

int main(int argc, const char** argv)
{
    hanoi H;

    unsigned nArgs = (unsigned)argc;
    unsigned i=1;

    if (1 == argc) {
        H.show_usage(std::cerr, argv[0]);
        return 1;
    }

    while (i < nArgs) {
        if (!H.process_arg(i, argv)) {
            H.show_usage(std::cerr, argv[0]);
            return 1;
        }
    }

    try {
        MEDDLY::initialize();
        MEDDLY::ostream_output out(std::cout);

        if (!H.build_model(out)) {
            H.show_usage(std::cerr, argv[0]);
            return 1;
        }
        H.generate(out, true, true);

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
