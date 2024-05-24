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
    Tests union, intersection, set difference operations.
*/


#include <cstdlib>
#include <time.h>
#include <cassert>

#include "../src/meddly.h"

const unsigned VARS = 5;
const unsigned DOM = 4;
const unsigned POTENTIAL = 1024;    // DOM ^ VARS
const unsigned MINCARD = 16;
const unsigned MAXCARD = 512;

using namespace MEDDLY;

// #define DEBUG_MINTERMS

double Random(long newseed=0)
{
    static long seed = 1;

    if (newseed) {
        seed = newseed;
    }

    const long MODULUS = 2147483647L;
    const long MULTIPLIER = 48271L;
    const long Q = MODULUS / MULTIPLIER;
    const long R = MODULUS % MULTIPLIER;

    long t = MULTIPLIER * (seed % Q) - R * (seed / Q);
    if (t > 0) {
        seed = t;
    } else {
        seed = t + MODULUS;
    }
    return ((double) seed / MODULUS);
}

unsigned Equilikely(unsigned a, unsigned b)
{
    return (a + (unsigned) ((b - a + 1) * Random()));
}

void showSet(std::ostream &out, const std::vector <bool> &elems)
{
    out << "{ ";
    bool printed = false;
    for (unsigned i=0; i<elems.size(); i++) {
        if (!elems[i]) continue;
        if (printed) out << ", ";
        out << i;
        printed = true;
    }
    out << " }";
}

void randomizeSet(std::vector <bool> &elems, unsigned card)
{
    //
    // Fill array with card 1s, the rest 0s
    //
    for (unsigned i=0; i<elems.size(); i++) {
        elems[i] = (i < card);
    }

    //
    // Shuffle the array
    //
    for (unsigned i=0; i<elems.size()-1; i++) {
        unsigned j = Equilikely(i, elems.size()-1);
        if (elems[i] != elems[j]) {
            // swap
            bool t = elems[i];
            elems[i] = elems[j];
            elems[j] = t;
        }
    }

}

unsigned set_cardinality(const std::vector <bool> &A)
{
    unsigned card = 0;
    for (unsigned i=0; i<A.size(); i++) {
        if (A[i]) ++card;
    }
    return card;
}

void set_intersection(const std::vector <bool> &A,
        const std::vector <bool> &B, std::vector <bool> &C)
{
    if (A.size() != B.size()) {
        throw "set intersection size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "set intersection size mismatch A,C";
    }

    for (unsigned i=0; i<C.size(); i++) {
        C[i] = A[i] && B[i];
    }
}

void set_union(const std::vector <bool> &A,
        const std::vector <bool> &B, std::vector <bool> &C)
{
    if (A.size() != B.size()) {
        throw "set union size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "set union size mismatch A,C";
    }

    for (unsigned i=0; i<C.size(); i++) {
        C[i] = A[i] || B[i];
    }
}

void set_difference(const std::vector <bool> &A,
        const std::vector <bool> &B, std::vector <bool> &C)
{
    if (A.size() != B.size()) {
        throw "set difference size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "set difference size mismatch A,C";
    }

    for (unsigned i=0; i<C.size(); i++) {
        C[i] = B[i] ? false : A[i];
    }
}


inline void fillMinterm(unsigned x, int* mt)
{
#ifdef DEBUG_MINTERMS
    std::cout << std::hex << x << ": [";
#endif
    for (unsigned j=1; j<=VARS; j++) {
        mt[j] = x % DOM;
        x /= DOM;
    }
#ifdef DEBUG_MINTERMS
    for (unsigned j=VARS; j; j--) {
        if (j<VARS) std::cout << ", ";
        std::cout << mt[j];
    }
    std::cout << "]\n";
#endif
}

void set2mdd(const std::vector<bool> &S, forest *F, dd_edge &s)
{
    if (!F) throw "null forest";

    // Determine S cardinality
    unsigned card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (S[i]) ++card;
    }

    // Special case - empty set
    if (0==card) {
        F->createEdge(false, s);
        return;
    }

    // Convert set S to list of minterms
    int** mtlist = new int* [card];
    card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (!S[i]) continue;
        mtlist[card] = new int[VARS+1];
        fillMinterm(i, mtlist[card]);
        ++card;
    }

    F->createEdge(mtlist, card, s);
}

inline char getReductionType(forest* f)
{
    if (f->isFullyReduced())    return 'f';
    if (f->isQuasiReduced())    return 'q';
    if (f->isIdentityReduced()) return 'i';
    throw "Unknown reduction type";
}

void test_pairs(unsigned scard, forest* f1, forest* f2, forest* fres)
{
    if (!f1) throw "null f1";
    if (!f2) throw "null f2";
    if (!fres) throw "null fres";

    std::cout << "    " << getReductionType(f1) << getReductionType(f2)
              << ':' << getReductionType(fres) << ' ';

    std::vector <bool> Aset(POTENTIAL);
    std::vector <bool> Bset(POTENTIAL);
    std::vector <bool> AiBset(POTENTIAL);
    std::vector <bool> AuBset(POTENTIAL);
    std::vector <bool> AmBset(POTENTIAL);

    for (unsigned i=0; i<16; i++) {
        std::cout << '.';
        randomizeSet(Aset, scard);
        randomizeSet(Bset, scard);
        set_intersection(Aset, Bset, AiBset);
        set_union(Aset, Bset, AuBset);
        set_difference(Aset, Bset, AmBset);

        using namespace MEDDLY;

        dd_edge Add(f1), Bdd(f2), AiBdd(fres), AuBdd(fres), AmBdd(fres);

        set2mdd(Aset, f1, Add);
        set2mdd(Bset, f2, Bdd);
        set2mdd(AiBset, fres, AiBdd);
        set2mdd(AuBset, fres, AuBdd);
        set2mdd(AmBset, fres, AmBdd);

        dd_edge AiBsym(fres), AuBsym(fres), AmBsym(fres);

        apply(INTERSECTION, Add, Bdd, AiBsym);
        apply(UNION, Add, Bdd, AuBsym);
        apply(DIFFERENCE, Add, Bdd, AmBsym);

        if (AiBsym != AiBdd) {
            throw "intersection mismatch";
        }
        if (AuBsym != AuBdd) {
            throw "union mismatch";
        }
        if (AmBsym != AmBdd) {
            throw "difference mismatch";
        }
    }
    std::cout << std::endl;
}

int main()
{
    using namespace MEDDLY;

    Random(11235);

    try {
        MEDDLY::initialize();
        int bounds[VARS];
        for (unsigned i=0; i<VARS; i++) {
            bounds[i] = DOM;
        }
        domain* D = domain::createBottomUp(bounds, VARS);

        policies p;
        p.useDefaults(false);

        forest* F1 = forest::create(D, false, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

        p.setQuasiReduced();

        forest* F2 = forest::create(D, false, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

        for (unsigned i=MINCARD; i<=MAXCARD; i*=2) {
            std::cout << "Testing sets of size " << i << " out of " << POTENTIAL << "\n";

            test_pairs(i, F1, F1, F1);
            test_pairs(i, F1, F1, F2);
            test_pairs(i, F1, F2, F1);
            test_pairs(i, F1, F2, F2);
            test_pairs(i, F2, F1, F1);
            test_pairs(i, F2, F1, F2);
            test_pairs(i, F2, F2, F1);
            test_pairs(i, F2, F2, F2);
        }

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
