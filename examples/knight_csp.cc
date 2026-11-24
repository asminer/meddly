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
#include <iomanip>
#include <cassert>

#include "meddly/meddly.h"

// ******************************************************************

struct constraint {
        MEDDLY::dd_edge edge;
        std::vector <bool> depends;
        unsigned dcount;
    public:
        constraint(MEDDLY::forest* f);
        void show_depends(std::ostream &s) const;
        void count_depends();
};

// ******************************************************************

constraint::constraint(MEDDLY::forest* f) : edge(f)
{
    depends.resize(1+f->getNumVariables());
    dcount = 0;
}

// ******************************************************************

void constraint::show_depends(std::ostream &s) const
{
    for (unsigned i=1; i<depends.size(); i++) {
        s << ( depends[i] ? 'X' : '.' );
    }
}

// ******************************************************************

void constraint::count_depends()
{
    dcount = 0;
    for (unsigned i=1; i<depends.size(); i++) {
        if (depends[i]) {
            ++dcount;
        }
    }
}

// ******************************************************************

bool less_than(const constraint* c1, const constraint* c2)
{
    assert(c1->depends.size() == c2->depends.size());
    for (unsigned i=c1->depends.size(); i; ) {
        --i;
        if (c1->depends[i]) {
            if (c2->depends[i]) continue;
            else                return false;
        }
        if (c2->depends[i])     return true;
    }
    // equal
    return false;
}

// ******************************************************************

void sort_by_depends(std::vector <constraint*> &A)
{
    for (unsigned i=0; i+1<A.size(); i++) {
        // put smallest constraint in i..end, at position i
        unsigned min = i;
        for (unsigned j=i+1; j<A.size(); j++) {
            if (less_than(A[j], A[min])) {
                min = j;
            }
        }
        if (min != i) {
            constraint* t = A[i];
            A[i] = A[min];
            A[min] = t;
        }
    }
}

// ******************************************************************

void sort_by_dcount(std::vector <constraint*> &A)
{
    for (unsigned i=0; i<A.size(); i++) {
        A[i]->count_depends();
    }
    for (unsigned i=0; i+1<A.size(); i++) {
        // put smallest constraint in i..end, at position i
        unsigned min = i;
        for (unsigned j=i+1; j<A.size(); j++) {
            if (A[j]->dcount < A[min]->dcount) {
                min = j;
            }
        }
        if (min != i) {
            constraint* t = A[i];
            A[i] = A[min];
            A[min] = t;
        }
    }
}

// ******************************************************************

class knight_csp {
        const int N;
        const int M;

        bool for_each_square_when;

        char verbosity;

        MEDDLY::domain* D;
        MEDDLY::forest* mdd;
        MEDDLY::minterm_coll* MC;

        std::vector <int> level4pos;

    public:
        knight_csp(int _N, int _M, bool square, char v);

        inline int posnum(int n, int m) const
        {
            if (n<0) return -1;
            if (m<0) return -1;
            if (n>=N) return -1;
            if (m>=M) return -1;

            return n + m*N;
        }
        inline void getpos(int num, int &n, int &m) const
        {
            if ((num<0) || (num >= N*M)) {
                n=-1;
                m=-1;
                return;
            }
            n = num % N;
            m = num / N;
        }

        void var_order_by_rows();
        void var_order_by_moves();
        void show_var_order(std::ostream &s) const;

        void fill_movelist(int n, int m, std::vector <int> &moves) const;

        void make_forest();
        void done_forest();

        // for each square, when

        constraint* constraint_for_square(int n, int m);
        constraint* open_tour_squares();
        constraint* closed_tour_squares();

        // for each time, which square

        // generic solver
        void gen_and_solve();
    protected:
        void a_equals_b_plus_one(int va, int vb, MEDDLY::dd_edge &e);
        void a_equals_n(int va, unsigned n, MEDDLY::dd_edge &e);
};

// ******************************************************************

knight_csp::knight_csp(int _N, int _M, bool square, char v)
    : N(_N), M(_M)
{
    assert(_N>0);
    assert(_M>0);

    for_each_square_when = square;
    verbosity = v;

    std::cout << "Building constraint model for knight's tour\n";
    std::cout << "Chess board is " << N << "x" << M << "\n";
    std::cout << "Encoding is ";
    if (for_each_square_when) {
        std::cout << "for each board position, when is it visited\n";
    } else {
        std::cout << "for each move, which position is visited\n";
    }

    //
    // sane initial values
    //
    D = nullptr;
    mdd = nullptr;
    MC = nullptr;
}

// ******************************************************************

void knight_csp::var_order_by_rows()
{
    level4pos.resize(N*M);
    int lvl=1;
    for (int n=0; n<N; n++) {
        for (int m=0; m<M; m++) {
            level4pos[ posnum(n, m) ] = lvl;
            ++lvl;
        }
    } // for n
}

// ******************************************************************

void knight_csp::var_order_by_moves()
{
    level4pos.resize(N*M);
    unsigned i;
    for (i=0; i<level4pos.size(); i++) {
        level4pos[i] = 0;
    }

    std::vector <int> queue;
    std::vector <int> movelist;

    int lvl = 1;
    level4pos[0] = lvl;
    ++lvl;
    queue.push_back(0);
    i=0;
    while (i<queue.size()) {
        int n, m;
        getpos(queue[i], n, m);
        fill_movelist(n, m, movelist);

        for (unsigned j=0; j<movelist.size(); j++) {
            const int mv = movelist[j];
            if (level4pos[mv]) continue;
            queue.push_back(mv);
            level4pos[mv] = lvl;
            ++lvl;
        }

        ++i;
    }
}

// ******************************************************************

void knight_csp::show_var_order(std::ostream &s) const
{
    if (!for_each_square_when) return;
    s << "Using variable order:\n";
    for (int n=0; n<N; n++) {
        for (int m=0; m<M; m++) {
            s << std::setw(3) << level4pos[ posnum(n, m) ] << " ";
        }
        s << "\n";
    } // for n
}

// ******************************************************************

void knight_csp::fill_movelist(int n, int m, std::vector <int> &moves) const
{
    assert(n>=0);
    assert(m>=0);
    assert(n<N);
    assert(m<M);

    int p;
    moves.clear();
    p = posnum(n +1, m +2);
    if (p>=0) moves.push_back(p);
    p = posnum(n +1, m -2);
    if (p>=0) moves.push_back(p);
    p = posnum(n -1, m +2);
    if (p>=0) moves.push_back(p);
    p = posnum(n -1, m -2);
    if (p>=0) moves.push_back(p);

    p = posnum(n +2, m +1);
    if (p>=0) moves.push_back(p);
    p = posnum(n +2, m -1);
    if (p>=0) moves.push_back(p);
    p = posnum(n -2, m +1);
    if (p>=0) moves.push_back(p);
    p = posnum(n -2, m -1);
    if (p>=0) moves.push_back(p);
}

// ******************************************************************

void knight_csp::make_forest()
{
    using namespace MEDDLY;

    assert(!D);

    //
    // Build domain
    //
    const int NxM = N*M;
    int* sizes = new int[NxM];
    for (int i=0; i<NxM; ++i) {
        sizes[i] = NxM;
    }
    D = domain::createBottomUp(sizes, NxM);
    delete[] sizes;

    //
    // Build forest
    //
    mdd = forest::create(
            D, SET, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL
    );

    //
    // Build collection of minterms
    //
    MC = new minterm_coll(10*NxM, mdd);
}

// ******************************************************************

void knight_csp::done_forest()
{
    using namespace MEDDLY;
    if (!mdd) return;
    ostream_output out(std::cout);

    out << "\n";
    out << "======================================================================\n";
    out << "\n";
    out << "MDD stats:\n";
    mdd->reportStats(out, "\t",
            HUMAN_READABLE_MEMORY | BASIC_STATS
            | EXTRA_STATS
            | STORAGE_STATS
            // | STORAGE_DETAILED
            // | HOLE_MANAGER_STATS
            // | HOLE_MANAGER_DETAILED
    );
    out << "\n";
    out.flush();

    forest::destroy(mdd);
}

// ******************************************************************

constraint* knight_csp::constraint_for_square(int n, int m)
{
    const int pa = level4pos[posnum(n, m)];
    assert(pa>=0);

    if (verbosity>0) {
        std::cout << "    Building constraint for square ("
                  << n << ", " << m << ") = #" << pa << "\n";
    }

    //
    // Set up C
    //
    constraint* C = new constraint(mdd);
    for (unsigned i=0; i<C->depends.size(); i++) {
        C->depends[i] = false;
    }

    //
    // What's reachable from n,m
    //
    std::vector <int> moves;
    fill_movelist(n, m, moves);

    //
    // The sequence number for a square must be either
    //      zero, or
    //      one plus the sequence number of some square reachable from us
    //
    if (verbosity > 2) {
        std::cout << "\tconstraints (OR'd):\n";
    }

    MEDDLY::minterm mt(mdd);
    mt.setAllVars(MEDDLY::DONT_CARE);
    mt.setVar(pa, 0);
    mt.setValue(true);
    mt.buildFunction(false, C->edge);
    C->depends[pa] = true;
    if (verbosity>2) {
        std::cout << "\t x_" << pa << " == 0\n";
    }

    MEDDLY::dd_edge ab(mdd);
    for (unsigned i=0; i<moves.size(); i++) {
        const int pb = level4pos[moves[i]];
        C->depends[pb] = true;
        a_equals_b_plus_one(pa, pb, ab);
        MEDDLY::apply(MEDDLY::UNION, C->edge, ab, C->edge);
        if (verbosity>2) {
            std::cout << "\t x_" << pa << " == 1+x_" << pb << "\n";
        }
    }

    //
    // done
    //
    if (verbosity>1) {
        std::cout << "\tDependencies for square ("
                  << n << ", " << m << ") = #" << pa << "\n\t    ";
        C->show_depends(std::cout);
        std::cout << "\n";
    }

    return C;
}

// ******************************************************************

constraint* knight_csp::open_tour_squares()
{
    if (verbosity>0) {
        std::cout << "Building open tour constraint\n";
        std::cout.flush();
    }

    //
    // Set up C
    //
    constraint* C = new constraint(mdd);
    C->depends[0] = false;
    for (unsigned i=1; i<C->depends.size(); i++) {
        C->depends[i] = true;
    }
    const unsigned NxM = unsigned(N*M);
    a_equals_n(1, NxM, C->edge);
    MEDDLY::dd_edge an(mdd);
    for (int va=2; va<=N*M; va++) {
        a_equals_n(va, NxM, an);
        MEDDLY::apply(MEDDLY::UNION, C->edge, an, C->edge);
    }
    if (verbosity>0) {
        std::cout << "   some square has sequence# " << NxM << "\n";
    }
    return C;
}

// ******************************************************************

constraint* knight_csp::closed_tour_squares()
{
    if (verbosity>0) {
        std::cout << "Building closed tour constraint\n";
        std::cout.flush();
    }
    //
    // Set up C
    //
    constraint* C = new constraint(mdd);
    for (unsigned i=0; i<C->depends.size(); i++) {
        C->depends[i] = false;
    }

    //
    // What's reachable from 0,0
    //
    std::vector <int> moves;
    fill_movelist(0, 0, moves);

    if (verbosity > 2) {
        std::cout << "    final position (OR'd):\n";
    }

    const unsigned NxM = unsigned(N*M);
    MEDDLY::dd_edge an(mdd);
    for (unsigned i=0; i<moves.size(); i++) {
        const int pb = level4pos[moves[i]];
        C->depends[pb] = true;
        if (i) {
            a_equals_n(pb, NxM, an);
            MEDDLY::apply(MEDDLY::UNION, C->edge, an, C->edge);
        } else {
            a_equals_n(pb, NxM, C->edge);
        }
        if (verbosity>2) {
            std::cout << "\t x_" << pb << " == " << NxM << "\n";
        }
    }

    const int pa = level4pos[posnum(0, 0)];

    if (verbosity > 2) {
        std::cout << "    starting position:\n";
    }
    C->depends[pa] = true;
    a_equals_n(pa, 0, an);
    MEDDLY::apply(MEDDLY::INTERSECTION, C->edge, an, C->edge);
    if (verbosity > 2) {
        std::cout << "\t x_" << pa << " == 0\n";
    }
    return C;
}

// ******************************************************************

void knight_csp::gen_and_solve()
{
    std::vector <constraint*> the_constraints;

    std::cout << "Building constraints\n";

    if (for_each_square_when) {
        for (int n=0; n<N; n++) {
            for (int m=0; m<M; m++) {
                constraint* C = constraint_for_square(n, m);
                if (C) the_constraints.push_back(C);
            } // for m
        } // for n
    } else {
        std::cout << "Sequence encoding not implemented yet\n";
    }

    std::cout << "Sorting constraints\n";
    // sort_by_dcount(the_constraints);
    sort_by_depends(the_constraints);

    if (verbosity>2) {
        for (unsigned i=0; i<the_constraints.size(); i++) {
            std::cout << "\t" << std::setw(2) << i << " : ";
            the_constraints[i]->show_depends(std::cout);
            std::cout << "\n";
        }
    }

    std::cout << "Accumulating constraints in order\n";
    for (unsigned i=1; i<the_constraints.size(); i++) {
        std::cout << "    0 := 0 AND " << std::setw(2) << i;
        std::cout.flush();

        MEDDLY::apply(MEDDLY::INTERSECTION,
                the_constraints[0]->edge, the_constraints[i]->edge,
                the_constraints[0]->edge);
        delete the_constraints[i];
        the_constraints[i] = nullptr;

        std::cout << "\t(" << mdd->getCurrentNumNodes() << " nodes)" << std::endl;
    }

    /*
    std::cout << "Accumulating constraints in pairs\n";
    for (;;) {
        unsigned fi=0;
        bool have_fi = false;
        for (unsigned i=0; i<the_constraints.size(); i++) {
            if (!the_constraints[i]) continue;
            if (!have_fi) {
                fi = i;
                have_fi = true;
                continue;
            }
            std::cout   << "    " << std::setw(2) << fi << " : "
                        << std::setw(2) << fi << " AND "
                        << std::setw(2) << i;
            std::cout.flush();

            MEDDLY::apply(MEDDLY::INTERSECTION,
                the_constraints[fi]->edge, the_constraints[i]->edge,
                the_constraints[fi]->edge);
            delete the_constraints[i];
            the_constraints[i] = nullptr;

            have_fi = false;
            std::cout   << "\t(" << mdd->getCurrentNumNodes()
                        << " nodes)" << std::endl;
        } // for

        if ((have_fi) && (fi==0)) break;

        std::cout << "    ------------------------------\n";
    }
    */

}

// ******************************************************************

void knight_csp::a_equals_b_plus_one(int va, int vb, MEDDLY::dd_edge &e)
{
    MC->clear();
    const unsigned NxM = unsigned (N * M);
    for (unsigned bv = 0; bv+1 < NxM; ++bv) {
        MC->unused().setAllVars(MEDDLY::DONT_CARE);
        MC->unused().setVar(vb, bv);
        MC->unused().setVar(va, bv+1);
        MC->unused().setValue(true);
        MC->pushUnused();
    }
    MC->buildFunctionMax(false, e);
}

// ******************************************************************

void knight_csp::a_equals_n(int va, unsigned n, MEDDLY::dd_edge &e)
{
    MEDDLY::minterm mt(mdd);
    mt.setAllVars(MEDDLY::DONT_CARE);
    mt.setVar(va, n);
    mt.setValue(true);
    mt.buildFunction(false, e);
}

// ******************************************************************

knight_csp* usage(const char* arg)
{
    const char* a = arg;
    for (; *arg; ++arg) {
        if ('/' == *arg) {
            a = arg+1;
        }
    }

    std::cout
        << "\nUsage: " << a << " [switches] NxM\n"
        << "\n"
        << "Build and count solutions to the Knight's tour on an NxM chessboard.\n"
        << "\nLegal switches:\n"
        << "\t-h:\tThis help\n"
        << "\n"
        << "\t-closed:\tConsider closed tours (start/stop at the same position)\n"
        << "\t-open:\tConsider open and closed tours\n"
        << "\n"
        << "\t-Omove:\tVariable order based on moves\n"
        << "\t-Orow:\tVariable order based on rows\n"
        << "\n"
        << "\t-pos:\tUse encoding: for each position, when to visit\n"
        << "\t-seq:\tUse encoding: for each sequence#, which position\n"
        << "\t-vL:\tSet verbosity level (default 0)\n"
        << "\n";

    return nullptr;
}

// ******************************************************************

knight_csp* process_args(int argc, const char** argv, bool &closed)
{
    bool pos = true;
    char verb = 0;
    unsigned N=0;
    unsigned M=0;
    bool moveorder = true;

    for (int i=1; i<argc; i++) {
        // easy switches
        if ('-' == argv[i][0]) {
            if (0==strcmp("-h", argv[i])) {
                return usage(argv[0]);
            }
            if (0==strcmp("-closed", argv[i])) {
                closed = true;
                continue;
            }
            if (0==strcmp("-open", argv[i])) {
                closed = false;
                continue;
            }
            if (0==strcmp("-pos", argv[i])) {
                pos = true;
                continue;
            }
            if (0==strcmp("-seq", argv[i])) {
                pos = false;
                continue;
            }
            if (0==strcmp("-Omove", argv[i])) {
                moveorder = true;
                continue;
            }
            if (0==strcmp("-Orow", argv[i])) {
                moveorder = false;
                continue;
            }
            if ( ('v' == argv[i][1]) &&
                 ('9' >= argv[i][2]) && ('0' <= argv[i][2]) &&
                 (0==argv[i][3]) )
            {
                verb = argv[i][2] - '0';
                continue;
            }
            std::cout << "Unknown switch: '" << argv[i] << "'\n";
            return usage(argv[0]);
        }

        // must be NxM
        bool sawx = false;
        for (const char* arg = argv[i]; *arg; ++arg) {

            if (('x' == *arg) || ('X' == *arg)) {
                if (sawx) return usage(argv[0]);
                sawx = true;
                continue;
            }

            if ((*arg < '0') || (*arg > '9')) {
                return usage(argv[0]);
            }

            if (sawx) {
                // next digit of M
                M *= 10;
                M += *arg - '0';
            } else {
                // next digit of N
                N *= 10;
                N += *arg - '0';
            }
        } // for arg
    }

    if ((0==M) || (0==N)) {
        std::cout << "Missing board dimension NxM\n";
        return usage(argv[0]);
    }

    knight_csp* K = new knight_csp(N, M, pos, verb);
    if (pos) {
        if (moveorder) {
            K->var_order_by_moves();
        } else {
            K->var_order_by_rows();
        }
    }
    return K;
}

// ******************************************************************

int main(int argc, const char** argv)
{
    bool closed = true;
    knight_csp* K = process_args(argc, argv, closed);

    if (!K) return 2;

    try {
        MEDDLY::initialize();
        K->show_var_order(std::cout);
        K->make_forest();
        K->gen_and_solve();
        K->done_forest();
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
