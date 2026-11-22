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

/*
    Model for knight's tour on a RxC chess board.
    Model state is, for each square on the board:
        3,  if there's a knight present, and we can return to this square
            (reserved for initial state only)
        2,  if there's a knight present, and we cannot return to this square
        1,  if we can move a knight to this square
        0,  if we cannot move a knight to this square.

    To get a solution, we must generate the distance function,
    and find a backwards path from state [2, 0, 0, ...]
    to initial state [3, 1, 1, ...].
*/
class knight_tour : public rsgen_model {
    public:
        knight_tour();

        inline bool has_dimensions() const {
            return (R>0) && (C>0);
        }
        inline unsigned get_R() const {
            return R;
        }
        inline unsigned get_C() const {
            return C;
        }
        inline unsigned pos2level(unsigned r, unsigned c) const {
            assert(r<R);
            assert(c<C);
            return 1 + r + c*R;
        }

        virtual bool process_arg(unsigned &i, const char** argv);
        virtual void show_usage(std::ostream &out, const char* arg0);

        virtual bool build_model(MEDDLY::output &out);
        virtual void build_initial(MEDDLY::minterm &init);

        virtual void after_generation(MEDDLY::output &out,
                const MEDDLY::dd_edge &e);

    protected:
        kron_event* make_move(unsigned r, unsigned c, int dr, int dc);

    private:
        unsigned R;
        unsigned C;
        bool closed;

        std::vector <level_delta> updates;
};

// ******************************************************************

knight_tour::knight_tour()
{
    R = C = 0;
    closed = true;
}

// ******************************************************************

bool knight_tour::process_arg(unsigned &i, const char** argv)
{
    if (rsgen_model::process_arg(i, argv)) return true;

    //
    // Should be the size in format RxC
    // We allow more than one digit for R and C
    // because we are optimists.
    //

    bool saw_x = false;
    R = 0;
    C = 0;
    for (const char* arg = argv[i]; *arg; ++arg) {

        if (('x' == *arg) || ('X' == *arg)) {
            if (saw_x) return false;
            saw_x = true;
            continue;
        }

        if ((*arg < '0') || (*arg > '9')) {
            return false;
        }

        if (saw_x) {
            // next digit of C
            C *= 10;
            C += *arg - '0';
        } else {
            // next digit of R
            R *= 10;
            R += *arg - '0';
        }
    }

    ++i;
    return true;
}

// ******************************************************************

void knight_tour::show_usage(std::ostream &out, const char* arg0)
{
    show_synopsis(out, arg0, "[options] RxC");

    show_description(out, "MODEL DESCRIPTION", "Knight's tour model. We start with a knight in one corner of the chess board. We can use the standard knight moves in chess except we remember if a square has been visited before or not (with slight exception for the starting position). If reachable states are generated, we show the cardinality, and indicate if a complete tour is possible. If distances are generated, we can construct a tour if one is possible.");

    show_section(out, "SIZE PARAMETERS");
    show_switch(out, "R", "Number of rows in the chess board");
    show_switch(out, "C", "Number of columns in the chess board");

    show_section(out, "TOUR TYPES");
    show_switch(out, "-closed", "Only consider closed tours (end at the starting position). Default.");
    show_switch(out, "-open", "Allow open tours starting from anywhere.");

    rsgen_model::show_usage(out);
}

// ******************************************************************

bool knight_tour::build_model(MEDDLY::output &out)
{
    if ((R<=0) || (C<=0)) {
        std::cerr << "Missing or bad size.\n";
        return false;
    }
    out << "Knight's tour model, chess board has " << R << " rows and "
        << C << " columns.\n";
    if (closed) {
        out << "Considering only closed tours.\n";
    } else {
        out << "Considering open and closed tours.\n";
    }

    //
    // Initialize domain
    //
    int* sizes = new int[R*C];
    for (int i=R*C-1; i>=0; i--) sizes[i] = 3;
    if (closed) {
        sizes[0] = 4;
    }
    set_domain(MEDDLY::domain::createBottomUp(sizes, R*C));
    delete[] sizes;

    //
    // Add possible moves
    //
    for (unsigned r=0; r<R; r++) {
        for (unsigned c=0; c<C; c++) {
            add_event(make_move(r, c, +1, +2));
            add_event(make_move(r, c, +1, -2));
            add_event(make_move(r, c, +2, +1));
            add_event(make_move(r, c, +2, -1));
            add_event(make_move(r, c, -1, +2));
            add_event(make_move(r, c, -1, -2));
            add_event(make_move(r, c, -2, +1));
            add_event(make_move(r, c, -2, -1));
        } // for c
    } // for r

    return true;
}

// ******************************************************************

void knight_tour::build_initial(MEDDLY::minterm &init)
{
    init.setAllVars(1);
    init.setVar(pos2level(0,0), 3);
}

// ******************************************************************

void knight_tour::after_generation(MEDDLY::output &out,
                const MEDDLY::dd_edge &rs)
{
    using namespace MEDDLY;

    forest* mdd = rs.getForest();

    if (!mdd) return;

    minterm finish(mdd);
    finish.setAllVars(0);
    finish.setVar(pos2level(0,0), 2);
    rangeval ans;
    rs.evaluate(finish, ans);

    if ('b' == get_set_type()) {
        show_cardinality(out, rs);

        if (ans) {
            out << "There is at least one closed tour.\n";
        } else {
            out << "There are no closed tours.\n";
        }
    } else {
        if (ans.isNormal()) {
            out << "There is at least one closed tour.\n";

            out << "TBD: show one!\n";
        } else {
            out << "There are no closed tours.\n";
        }
    }
}

// ******************************************************************

rsgen_model::kron_event*
knight_tour::make_move(unsigned r, unsigned c, int dr, int dc)
{
    const int newr = int(r) + dr;
    if (newr<0) return nullptr;
    const unsigned rp = unsigned(newr);
    if (rp >= R) return nullptr;

    const int newc = int(c) + dc;
    if (newc<0) return nullptr;
    const unsigned cp = unsigned(newc);
    if (cp >= C) return nullptr;

    updates.clear();
    level_delta LD;
    if (closed && (0==r) && (0==c)) {
        LD.set(pos2level(r, c), 3, 1);
        updates.push_back(LD);
    }
    LD.set(pos2level(r, c), 2, 0);
    updates.push_back(LD);
    LD.set(pos2level(rp, cp), 1, 2);
    updates.push_back(LD);

    return new kron_event("move", updates);
}

// ******************************************************************

int main(int argc, const char** argv)
{
    knight_tour T;

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

