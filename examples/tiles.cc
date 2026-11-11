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

// ******************************************************************

class tilemodel : public rsgen_model {
        template <bool FEPWT>   // for each position which tile?
        class move_tile_event : public rsgen_model::event {
                int RxC;    // number of board squares
                int p_old;  // tile's current position
                int p_new;  // tile's new position
            public:
                /**
                    An event to move the tile from one position to another.
                        @param RC   Rows times columns
                        @param po   Position currently containing tile,
                                    from 0 to RC-1.
                        @param pn   New position for tile, from 0 to RC-1;
                                    must currently contain nothing (0)
                */
                move_tile_event(int RC, int po, int pn);

                virtual bool fire(const MEDDLY::minterm &curr,
                        MEDDLY::minterm &);

                virtual bool build_nsf(MEDDLY::dd_edge &R);
        };
    public:
        tilemodel();

        bool has_dimensions() const {
            return (R>0) && (C>0);
        }

        inline int get_R() const {
            return R;
        }
        inline int get_C() const {
            return C;
        }

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

    protected:

        virtual bool known_cardinality(double &card) const;
#ifdef HAVE_LIBGMP
        virtual bool known_cardinality(mpz_t card) const;
#endif

        rsgen_model::event* move_tile(int po, int pn) {
            const int RxC = R*C;
            if (po == pn)   return nullptr;
            if (po<0)       return nullptr;
            if (po>= RxC)   return nullptr;
            if (pn<0)       return nullptr;
            if (pn>= RxC)   return nullptr;
            if (for_each_position_which_tile) {
                return new move_tile_event<true>(RxC, po, pn);
            } else {
                return new move_tile_event<false>(RxC, po, pn);
            }
        }


    private:
        int R;
        int C;
        bool for_each_position_which_tile;
};

// ******************************************************************
// *                  tilemodel::move_tile methods                  *
// ******************************************************************

template <bool FEPWT>
tilemodel::move_tile_event<FEPWT>::move_tile_event(int RC, int po, int pn)
{
    assert(po != pn);
    assert(po >= 0);
    assert(po < RC);
    assert(pn >= 0);
    assert(pn < RC);

    RxC = RC;
    if (FEPWT) {
        p_old = po+1;
        p_new = pn+1;
    } else {
        p_old = po;
        p_new = pn;
    }

}

// ******************************************************************

template <bool FEPWT>
bool tilemodel::move_tile_event<FEPWT>::fire(const MEDDLY::minterm &curr,
        MEDDLY::minterm &next)
{
    if (FEPWT) {
        /*
            For each position which tile.

            Element p_old:  cannot be zero; becomes zero
            Element p_new:  must be zero, becomes current[p_old].

        */
        if (curr.getVar(p_old) == 0) return false;
        if (curr.getVar(p_new) != 0) return false;
        for (unsigned i=next.getNumVars(); i; i--) {
            next.setVar(i, curr.getVar(i));
        }
        next.setVar(p_old, 0);
        next.setVar(p_new, curr.getVar(p_old));
        return true;

    } else {
        /*
            For each tile, where is it.

            Check all tiles except blank (level 1):
                If in position p_new, disabled.
                If in position p_old, change position to p_new.
                Otherwise, no change.

            Check blank tile (level 1):
                If not in position p_new, disabled.
                Change to position p_old.
        */
        for (unsigned i=next.getNumVars(); i>1; i--) {
            if (curr.getVar(i) == p_new) return false;
            next.setVar(i,
                    (curr.getVar(i) == p_old)
                    ? p_new
                    : curr.getVar(i)
            );
        }
        if (curr.getVar(1) != p_new) return false;
        next.setVar(1, p_old);
        return true;

    }
}

// ******************************************************************

template <bool FEPWT>
bool tilemodel::move_tile_event<FEPWT>::build_nsf(MEDDLY::dd_edge &R)
{
    using namespace MEDDLY;
    minterm_coll mtc(RxC, R.getForest());

    if (FEPWT) {
        //  For each position, which tile
        //  p_old: T -> 0
        //  p_new: 0 -> T
        for (int T=1; T<RxC; T++) {
            minterm &m = mtc.unused();
            m.setAllVars(DONT_CARE, DONT_CHANGE);
            m.setVars(p_old, T, 0);
            m.setVars(p_new, 0, T);
            m.setValue(true);
            mtc.pushUnused();
        }
    } else {
        //  For each tile, where is it
        //  x_T: p_old -> p_new
        //  x_1: p_new -> p_old
        for (int T=2; T<=RxC; T++) {
            minterm &m = mtc.unused();
            m.setAllVars(DONT_CARE, DONT_CHANGE);
            m.setVars(T, p_old, p_new);
            m.setVars(1, p_new, p_old);
            m.setValue(true);
            mtc.pushUnused();
        }
    }
    mtc.buildFunctionMax(false, R);
    return true;
}

// ******************************************************************
// *                       tilemodel  methods                       *
// ******************************************************************

tilemodel::tilemodel()
{
    R = -1;
    C = -1;
    for_each_position_which_tile = false;
}

// ******************************************************************

bool tilemodel::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    if (0==strcmp("-pos", argv[i])) {
        ++i;
        for_each_position_which_tile = true;
        return true;
    }
    if (0==strcmp("-tile", argv[i])) {
        ++i;
        for_each_position_which_tile = false;
        return true;
    }

    //
    // Should be the size in format RxC
    // where R,C are characters '1' .. '9'.
    //

    if (3 != strlen(argv[i])) {
        return false;
    }
    if ('x' != argv[i][1]) {
        return false;
    }
    if ( ('1' > argv[i][0]) || ('9' < argv[i][0]) ||
         ('1' > argv[i][2]) || ('9' < argv[i][2]) )
    {
        return false;
    }

    R = argv[i][0] - '0';
    C = argv[i][2] - '0';
    i++;
    return true;
}

// ******************************************************************

void tilemodel::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options] RxC\n\n";
    }
    out << "\tR          : number of rows, from 1 to 9\n";
    out << "\tC          : number of columns, from 1 to 9\n";

    out << "\n";
    out << "    ENCODING OPTIONS\n";
    out << "\n";

    out << "        -pos       : Use a position-based encoding: for each position,\n";
    out << "                     which tile is there. Position #0 is the bottom DD\n";
    out << "                     variable, and position #rc-1 is the top DD variable.\n";

    out << "        -tile      : Use a tile-based encoding: for each tile,\n";
    out << "                     where is it (default). The empty tile is the bottom\n";
    out << "                     DD variable, tile 1 is next, ..., tile rc-1 is the top\n";
    out << "                     DD variable.\n";

    rsgen_model::show_usage(out);
}

// ******************************************************************

bool tilemodel::build_model(MEDDLY::output &out)
{
    if ((R<0) || (C<0)) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    out << "Sliding tile puzzle, board has " << R << " rows and "
        << C << " columns.\n";
    if (for_each_position_which_tile) {
        out << "Encoding using 'for each position, which tile is there'\n";
    } else {
        out << "Encoding using 'for each tile, where is it'\n";
    }
    out.flush();

    //
    // Initialize domain
    //
    int* sizes = new int[R*C];
    for (int i=R*C-1; i>=0; i--) sizes[i] = R*C;
    set_domain(MEDDLY::domain::createBottomUp(sizes, R*C));
    delete[] sizes;

    //
    // Add events
    //
    for (int p=0; p<R*C; p++) {

        // Move left
        if (0 != p%C) {
            add_event(move_tile(p, p-1));
        }

        // Move right
        if (C-1 != p%C) {
            add_event(move_tile(p, p+1));
        }

        // Move up
        add_event(move_tile(p, p-C));

        // Move down
        add_event(move_tile(p, p+C));
    }

    return true;
}

// ******************************************************************

void tilemodel::build_initial(MEDDLY::minterm &init)
{
    if (for_each_position_which_tile) {
        // p0: 1, p1: 2, ..., plast-1: R*C-1, plast: 0
        for (int p=0; p<R*C-1; p++) {
            init.setVar(p+1, p+1);
        }
        init.setVar(R*C, 0);
    } else {
        // t0: R*C-1, t1: 0, t2: 1, ...
        for (int t=1; t<R*C; t++) {
            init.setVar(t+1, t-1);
        }
        init.setVar(1, R*C-1);
    }
}

// ******************************************************************

bool tilemodel::known_cardinality(double &theory) const
{
    //
    // (R*C) ! / 2
    //
    theory = 1.0;
    for (int i=3; i<=R*C; i++) {
        theory *= i;
    }
    return true;
}

// ******************************************************************

#ifdef HAVE_LIBGMP

bool tilemodel::known_cardinality(mpz_t theory) const
{
    //
    // (R*C) ! / 2
    //
    mpz_set_ui(theory, 1);
    for (int i=3; i<=R*C; i++) {
        mpz_mul_ui(theory, theory, i);
    }
    return true;
}

#endif

// **********************************************************************
// main
// **********************************************************************

int main(int argc, const char** argv)
{
    tilemodel T;

    unsigned nArgs = (unsigned)argc;
    unsigned i=1;

    while (i < nArgs) {
        if (!T.process_arg(i, argv)) {
            T.show_usage(std::cerr, argv[0]);
            return 1;
        }
    }
    if (!T.has_dimensions()) {
        T.show_usage(std::cerr, argv[0]);
        return 1;
    }

    try {
        MEDDLY::initialize();
        MEDDLY::ostream_output out(std::cout);

        if (!T.build_model(out)) {
            T.show_usage(std::cerr, argv[0]);
            return 1;
        }
        T.generate(out, true, true);

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

