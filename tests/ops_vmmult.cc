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

/*
    Tests vector-matrix and matrix-vector multiplication.
*/


#include <cstdlib>
#include <time.h>
#include <cassert>
#include <iostream>
#include <iomanip>

#include "../src/meddly.h"

#include "randomize.h"

vectorgen SG(MEDDLY::SET,      5, 3);
vectorgen RG(MEDDLY::RELATION, 5, 3);

const unsigned MIN_SET_CARD = 1;
const unsigned MAX_SET_CARD = 16;
const unsigned MULT_SET_CARD = 4;

const unsigned MIN_REL_CARD = 8;
const unsigned MAX_REL_CARD = 512;
const unsigned MULT_REL_CARD = 8;

using namespace MEDDLY;

void usage(const char* arg0)
{
    /* Strip leading directory, if any: */
    const char* name = arg0;
    for (const char* ptr=arg0; *ptr; ptr++) {
        if ('/' == *ptr) name = ptr+1;
    }
    std::cerr << "\nUsage: " << name << "options seed\n\n";
    std::cerr << "Options:\n";
    std::cerr << "    --vF:     Input/output vector is fully-reduced\n";
    std::cerr << "    --vQ:     Input/output vector is quasi-reduced\n";
    std::cerr << "\n";

    std::cerr << "    --mF:     Matrix is fully-reduced\n";
    std::cerr << "    --mI:     Matrix is identity-reduced\n";
    std::cerr << "    --mQ:     Matrix is quasi-reduced\n";
    std::cerr << "\n";

    std::cerr << "    --int:    Integer-valued vector/matrix\n";
    std::cerr << "    --real:   Real-valued vector/matrix\n";
    std::cerr << "\n";

    exit(1);
}


struct matrlist {
    unsigned index;
    MEDDLY::rangeval value;
    matrlist* next;

    matrlist(unsigned i, MEDDLY::rangeval v, matrlist* n) {
        index = i;
        value = v;
        next = n;
    }
};

/*
void showList(intlist* L)
{
    while (L) {
        std::cout << L->index;
        std::cout << " -> ";
        L = L->next;
    }
    std::cout << " null\n";
}
*/

void makeExplicitGraph(const dd_edge &E, std::vector <matrlist*> &elist)
{
    elist.resize(SG.potential());
    for (unsigned i=0; i<elist.size(); i++) {
        elist[i] = nullptr;
    }
#ifdef SHOW_GRAPH
    ostream_output out(std::cout);
    out << "Explicit graph:\n";
#endif
    for (dd_edge::iterator I = E.begin(); I != E.end(); I++) {
        unsigned f, t;
#ifdef SHOW_GRAPH
        out << "    ";
        (*I).show(out);
        out << "\n";
#endif
        RG.minterm2indexes(*I, f, t);
        // std::cout << "    " << f << " -> " << t << "\n";

        elist[f] = new matrlist(t, (*I).getValue(), elist[f]);
    }
}

/*
void showExplicitGraph(const std::vector <intlist*> &elist)
{
    for (unsigned i=0; i<elist.size(); i++) {
        if (elist[i]) {
            std::cout << i << ": ";
            showList(elist[i]);
        }
    }
}
*/

inline void dist_update(MEDDLY::rangeval &x, long d)
{
    if (x.isPlusInfinity()) {
        x = d;
        return;
    }
    if (long(x) < 0) {
        x = d;
        return;
    }
    if (d < long(x)) {
        x = d;
    }
}

// y = A x
template <typename TYPE>
void matr_vect(const std::vector <matrlist*> &A,
        const std::vector <MEDDLY::rangeval> &x,
        std::vector <MEDDLY::rangeval> &y)
{
    if (x.size() != A.size()) {
        throw "vect_matr size mismatch x, A";
    }
    if (y.size() != A.size()) {
        throw "vect_matr size mismatch y, A";
    }
    for (unsigned i=0; i<y.size(); i++) {
        y[i] = TYPE(0);
    }
    for (unsigned i=0; i<A.size(); i++) {
        const matrlist* curr = A[i];
        while (curr) {
            y[i] = TYPE(y[i]) + TYPE(x[curr->index]) * TYPE(curr->value);
            curr = curr->next;
        }
    }

}

// y = x A
template <typename TYPE>
void vect_matr(const std::vector <MEDDLY::rangeval> &x,
        const std::vector <matrlist*> &A, std::vector <MEDDLY::rangeval> &y)
{
    if (x.size() != A.size()) {
        throw "vect_matr size mismatch x, A";
    }
    if (y.size() != A.size()) {
        throw "vect_matr size mismatch y, A";
    }
    for (unsigned i=0; i<y.size(); i++) {
        y[i] = TYPE(0);
    }
    for (unsigned i=0; i<A.size(); i++) {
        if (0 == TYPE(x[i])) continue;
        const matrlist* curr = A[i];
        while (curr) {
            y[curr->index]
                = TYPE(y[curr->index]) + TYPE(x[i]) * TYPE(curr->value);
            curr = curr->next;
        }
    }

}

void checkEqual(const char* what, const dd_edge &e1, const dd_edge &e2)
{
    if (e1 == e2) return;

    ostream_output out(std::cout);

    out << "\nMismatch on " << what << "\n";
    out << "Expected DD:\n";
    e2.showGraph(out);
    out << "Obtained DD:\n";
    e1.showGraph(out);

    out << "\nExpected minterms:\n";
    for (dd_edge::iterator s = e2.begin(); s; ++s)
    {
        out << "    ";
        (*s).show(out);
        out << "\n";
    }

    out << "Obtained minterms:\n";
    for (dd_edge::iterator s = e1.begin(); s; ++s)
    {
        out << "    ";
        (*s).show(out);
        out << "\n";
    }


    throw "mismatch";
}

template <typename TYPE>
void compare(const std::vector <MEDDLY::rangeval> &x,
        const dd_edge &ddA, forest* Fset)
{
    std::vector <MEDDLY::rangeval> xA(x.size());
    std::vector <MEDDLY::rangeval> Ax(x.size());

    std::vector <matrlist*> A;
    makeExplicitGraph(ddA, A);

    vect_matr<TYPE>(x, A, xA);
    matr_vect<TYPE>(A, x, Ax);

    dd_edge ddx(Fset);
    SG.explicit2edgeMax(x, ddx, TYPE(0));
#ifdef SHOW_INITIAL
    /*
    std::cout << "Initial minterm set:\n";
    SG.showMinterms(std::cout, Fset->getDomain(), S0, unreachable);
    */
    ostream_output out(std::cout);
    out << "Initial minterms:\n";
    for (dd_edge::iterator s = dd0.begin(); s; ++s)
    {
        out << "    ";
        (*s).show(out);
        out << "\n";
    }
#endif

    //
    // Check vect-matr
    //
    dd_edge ddxA(Fset), symxA(Fset);
    SG.explicit2edgeMax(xA, ddxA, TYPE(0));
    apply(VM_MULTIPLY, ddx, ddA, symxA);
    checkEqual("vector-matrix multiply", symxA, ddxA);

    //
    // Check matr-vect
    //
    dd_edge ddAx(Fset), symAx(Fset);
    SG.explicit2edgeMax(Ax, ddAx, TYPE(0));
    apply(MV_MULTIPLY, ddA, ddx, symAx);
    checkEqual("matrix-vector multiply", symAx, ddAx);
}

template <typename TYPE>
void testOnForests(unsigned scard, forest* Fset,
        unsigned rcard, forest* Frel)
{
    if (!Fset) throw "null Fset";
    if (!Frel) throw "null frel";

    std::vector <MEDDLY::rangeval> expset(SG.potential());
    std::vector <MEDDLY::rangeval> exprel(RG.potential());
    dd_edge Rel(Frel);

    std::vector <MEDDLY::rangeval> vals(5);
    vals[0] = TYPE(0);
    vals[1] = TYPE(1);
    vals[2] = TYPE(-1);
    vals[3] = TYPE(2);
    vals[4] = TYPE(-2);

    const unsigned N = 8;

    std::cerr << "    " << std::setw(5) << scard << " x "
              << std::setw(5) << rcard << " ";

    for (unsigned i=0; i<N; i++) {
        std::cerr << '.';
        RG.randomizeVector(exprel, rcard, vals);
        RG.explicit2edgeMax(exprel, Rel, false);

        SG.randomizeVector(expset, scard, vals);
        compare<TYPE>(expset, Rel, Fset);

        SG.randomizeFully(expset, scard, vals);
        compare<TYPE>(expset, Rel, Fset);
    }
    for (unsigned i=0; i<N; i++) {
        std::cerr << "x";
        RG.randomizeFully(exprel, rcard, vals);
        RG.explicit2edgeMax(exprel, Rel, false);

        SG.randomizeVector(expset, scard, vals);
        compare<TYPE>(expset, Rel, Fset);

        SG.randomizeFully(expset, scard, vals);
        compare<TYPE>(expset, Rel, Fset);
    }
    for (unsigned i=0; i<N; i++) {
        std::cerr << "i";
        RG.randomizeIdentity(exprel, rcard, vals);
        RG.explicit2edgeMax(exprel, Rel, false);

        SG.randomizeVector(expset, scard, vals);
        compare<TYPE>(expset, Rel, Fset);

        SG.randomizeFully(expset, scard, vals);
        compare<TYPE>(expset, Rel, Fset);
    }
    std::cerr << std::endl;
}


void runTests(domain* D, range_type RT, const policies &sp, const policies &rp)
{
    forest* Fset = forest::create(D, SET, RT,
                        edge_labeling::MULTI_TERMINAL, sp);

    forest* Frel = forest::create(D, RELATION, RT,
                        edge_labeling::MULTI_TERMINAL, rp);

    std::cerr << "Using " << nameOf(RT) << " values\n";

    std::cerr << "Vectors:  " << nameOf(Fset->getReductionRule()) << "\n";
    std::cerr << "Matrices: " << nameOf(Frel->getReductionRule()) << "\n";


    std::cerr << "Running tests for #nonzeroes in vector x matrix:\n";

    for (unsigned rc=MIN_REL_CARD; rc<=MAX_REL_CARD; rc*=MULT_REL_CARD)
    {
        for (unsigned sc=MIN_SET_CARD; sc<=MAX_SET_CARD; sc*=MULT_SET_CARD)
        {
            if (range_type::INTEGER == RT) {
                testOnForests<int>(sc, Fset, rc, Frel);
            } else {
                testOnForests<float>(sc, Fset, rc, Frel);
            }
        }
    }

}

int main(int argc, const char** argv)
{
    /*
     * Command-line options
     */
    long seed = 0;
    char vred = 'F';
    char mred = 'I';

    range_type RT = range_type::INTEGER;

    for (int i=1; i<argc; i++) {

        if (0==strcmp("--int", argv[i])) {
            RT = range_type::INTEGER;
            continue;
        }
        if (0==strcmp("--real", argv[i])) {
            RT = range_type::REAL;
            continue;
        }

        if (0==strcmp("--vF", argv[i])) {
            vred = 'F';
            continue;
        }
        if (0==strcmp("--vQ", argv[i])) {
            vred = 'Q';
            continue;
        }

        if (0==strcmp("--mF", argv[i])) {
            mred = 'F';
            continue;
        }
        if (0==strcmp("--mQ", argv[i])) {
            mred = 'Q';
            continue;
        }
        if (0==strcmp("--mI", argv[i])) {
            mred = 'I';
            continue;
        }

        if ((argv[i][0] < '0') || (argv[i][0] > '9')) {
            usage(argv[0]);
        }

        seed = atol(argv[i]);
    }
    vectorgen::setSeed(seed);

    /*
     * Run requested test
     */

    try {
        MEDDLY::initialize();
        domain* D = RG.makeDomain();

        policies vp;
        policies mp;

        vp.useDefaults(SET);
        mp.useDefaults(RELATION);

        switch (vred) {
            case 'Q':
                vp.setQuasiReduced();
                break;

            case 'F':
                vp.setFullyReduced();
                break;

            default:
                throw "internal error 1";
        }
        switch (mred) {
            case 'Q':
                mp.setQuasiReduced();
                break;

            case 'F':
                mp.setFullyReduced();
                break;

            case 'I':
                mp.setIdentityReduced();
                break;

            default:
                throw "internal error 2";
        }

        runTests(D, RT, vp, mp);

        domain::destroy(D);
        MEDDLY::cleanup();
        return 0;
    }

    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error " << e.getName()
                    << "\n    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    std::cerr << "\nSome other error?\n";
    return 4;
}


