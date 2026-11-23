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

class knight_csp {
        int N, M;

        bool for_each_square_when;

        char verbosity;

        MEDDLY::domain* D;
        MEDDLY::forest* mdd;
        MEDDLY::minterm_coll* MC;

    public:
        struct constraint {
            MEDDLY::dd_edge edge;
            std::vector <bool> depends;
        };

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

        void fill_movelist(int n, int m, std::vector <int> &moves) const;

        void make_forest();

        void constraint_for_square(int n, int m, constraint &C);

    protected:
        void a_equals_b_plus_one(int va, int vb, MEDDLY::dd_edge &e);
};

// ******************************************************************

knight_csp::knight_csp(int _N, int _M, bool square, char v)
{
    assert(_N>0);
    assert(_M>0);

    N = _N;
    M = _M;
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

void knight_csp::constraint_for_square(int n, int m, constraint &C)
{
    const int pa = posnum(n, m);
    assert(pa>=0);

    if (verbosity>0) {
        std::cout << "Building constraint for square ("
                  << n << ", " << m << ") = #" << pa << "\n";
    }

    //
    // Set up C
    //
    C.edge.attach(mdd);
    C.depends.resize(N*M+1);
    for (unsigned i=0; i<C.depends.size(); i++) {
        C.depends[i] = false;
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
        std::cout << "    constraints (OR'd):\n";
    }

    MEDDLY::minterm mt(mdd);
    mt.setAllVars(MEDDLY::DONT_CARE);
    mt.setVar(1+pa, 0);
    mt.setValue(true);
    mt.buildFunction(false, C.edge);
    C.depends[1+pa] = true;
    if (verbosity>2) {
        std::cout << "\t x_" << (1+pa) << " == 0\n";
    }

    MEDDLY::dd_edge ab(mdd);
    for (unsigned i=0; i<moves.size(); i++) {
        C.depends[1+moves[i]] = true;
        a_equals_b_plus_one(1+pa, 1+moves[i], ab);
        MEDDLY::apply(MEDDLY::UNION, C.edge, ab, C.edge);
        if (verbosity>2) {
            std::cout << "\t x_" << (1+pa) << " == 1+x_" << (1+moves[i]) << "\n";
        }
    }

    //
    // done
    //
    if (verbosity>1) {
        std::cout << "Dependencies for square ("
                  << n << ", " << m << ") = #" << pa << "\n";
        for (unsigned i=1; i<N*M+1; i++) {
            std::cout << ( C.depends[i] ? 'X' : '.' );
        }
        std::cout << "\n";
    }

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

