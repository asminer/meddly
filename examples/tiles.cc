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

#include "meddly/meddly.h"
#include "utils/rsgen_model.h"

class tilemodel : public rsgen_model {
        template <bool FEPWT>   // for each position which tile?
        class move_tile : public rsgen_model::event {
                int T;      // tile to move
                int p_old;  // tile's current position
                int p_new;  // tile's new position
            public:
                /**
                    Build an event to move a tile.
                        @param T    Tile number, from 1 to R*C-1.
                        @param po   Position currently containing tile T,
                                    from 0 to R*C-1.
                        @param pn   New position for tile, from 0 to R*C-1;
                                    must currently contain nothing (0)
                */
                move_tile(int T, int po, int pn);

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

    private:
        int R;
        int C;
        bool for_each_position_which_tile;
};

// ******************************************************************
// *                  tilemodel::move_tile methods                  *
// ******************************************************************

template <bool FEPWT>
tilemodel::move_tile<FEPWT>::move_tile(int t, int po, int pn)
{
    T = t;
    p_old = po;
    p_new = pn;

    assert(T>0);
    assert(p_old != p_new);
    assert(p_old >= 0);
    assert(p_new >= 0);
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

    // TBD FROM HERE

    /*

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
    */

    out << "Not implemented\n";
    return false;
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


// **********************************************************************
// build a minterm to move a tile
//
//    @param T              Tile number, from 1 to R*C-1.
//    @param p              Position currently containing tile T,
//                              from 0 to R*C-1.
//    @param q              Neighboring position containing nothing (0),
//                              from 0 to R*C-1.
//    @param m              Minterm to fill
//
//    @return true, iff this is a valid move
//
// **********************************************************************

/*

bool moveTile(const int T, const int p, const int q, MEDDLY::minterm &m)
{
    if (p == q)     return false;
    if (T<=0)       return false;
    if (T>= R*C)    return false;
    if (p<0)        return false;
    if (p>= R*C)    return false;
    if (q<0)        return false;
    if (q>= R*C)    return false;

    m.setAllVars(MEDDLY::DONT_CARE, MEDDLY::DONT_CHANGE);

    if (for_each_position_which_tile) {

        // x_p: T -> 0
        //
        m.setVars(p+1, T, 0);

        // x_q: 0 -> T
        //
        m.setVars(q+1, 0, T);

    } else {

        // x_T: p -> q
        //
        m.setVars(T+1, p, q);

        // x_0: q -> p
        //
        m.setVars(1, q, p);
    }

    m.setValue(true);
#ifdef DEBUG_EVENTS
    std::cerr << "   move tile " << T << " from " << p << " to " << q << "\n";
#endif
    return true;
}
*/

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

