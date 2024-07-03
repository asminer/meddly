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
const unsigned RELDOM = 2;
const unsigned SETDOM = RELDOM * RELDOM;
const unsigned SETBITS = 2;
const unsigned POTENTIAL = 1024;    // SETDOM ^ VARS
const unsigned MAX_SET_CARD = 128;
const unsigned MAX_REL_CARD = 16;

using namespace MEDDLY;

// #define DEBUG_MXDOPS

#define TEST_SETS
#define TEST_RELATIONS

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

void showRelMinterms(std::ostream &out, const std::vector <bool> &elems)
{
    out << "{ ";
    bool printed = false;
    for (unsigned i=0; i<elems.size(); i++) {
        if (!elems[i]) continue;
        if (printed) out << ",\n      ";
        out << "[";
        printed = true;

        unsigned x = i;
        for (unsigned j=1; j<=VARS; j++) {
            unsigned un = x % RELDOM;
            x /= RELDOM;
            unsigned pr = x % RELDOM;
            x /= RELDOM;
            if (j>1) out << ", ";
            out << un << "->" << pr;
        }
        out << "]";
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

bool has_common_element(const std::vector <bool> &A,
        const std::vector <bool> &B)
{
    if (A.size() != B.size()) {
        throw "set intersection size mismatch A,B";
    }

    for (unsigned i=0; i<A.size(); i++) {
        if (A[i] && B[i]) return true;
    }
    return false;
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

void set_complement(const std::vector <bool> &A, std::vector <bool> &C)
{
    if (A.size() != C.size()) {
        throw "set difference size mismatch A,C";
    }

    for (unsigned i=0; i<C.size(); i++) {
        C[i] = !A[i];
    }
}


inline void fillMinterm(unsigned x, int* mt)
{
    for (unsigned j=1; j<=VARS; j++) {
        mt[j] = x % SETDOM;
        x /= SETDOM;
    }
}

inline void fillMinterm(unsigned x, int* un, int* pr)
{
    for (unsigned j=1; j<=VARS; j++) {
        un[j] = x % RELDOM;
        x /= RELDOM;
        pr[j] = x % RELDOM;
        x /= RELDOM;
    }
}

inline unsigned whichMinterm(const int* mt)
{
    unsigned x=0;
    for (unsigned j=VARS; j; j--) {
        x <<= SETBITS;
        x |= mt[j];
    }
    return x;
}

inline unsigned whichMinterm(const int* un, const int* pr)
{
    unsigned x=0;
    for (unsigned j=VARS; j; j--) {
        x <<= SETBITS;
        x |= (un[j] << (SETBITS/2)) | pr[j];
    }
    return x;
}

void flipFullyElements(std::vector <bool> &elems, unsigned seed, unsigned k,
        bool on)
{
    // Convert 'seed' to minterm
    int tmp[VARS+1];
    fillMinterm(seed, tmp);
    // Replace variable k with all possible values in minterm,
    // and convert back to unsigned
    for (tmp[k]=0; tmp[k]<SETDOM; tmp[k]++) {
        unsigned x = whichMinterm(tmp);
        elems[x] = on;
    }
}

void flipFullyElements(std::vector <bool> &elems, unsigned seed, unsigned k1,
        unsigned k2, bool on)
{
    // Convert 'seed' to minterm
    int tmp[VARS+1];
    fillMinterm(seed, tmp);
    // Replace variable k with all possible values in minterm,
    // and convert back to unsigned
    for (tmp[k1]=0; tmp[k1]<SETDOM; tmp[k1]++) {
        for (tmp[k2]=0; tmp[k2]<SETDOM; tmp[k2]++) {
            unsigned x = whichMinterm(tmp);
            elems[x] = on;
        }
    }
}


void flipIdentityElements(std::vector <bool> &elems, unsigned seed,
        unsigned k, bool on)
{
    // Convert 'seed' to minterm
    int untmp[VARS+1];
    int prtmp[VARS+1];
    fillMinterm(seed, untmp, prtmp);
    for (unsigned v=0; v<RELDOM; v++) {
        untmp[k] = v;
        prtmp[k] = v;
        unsigned x = whichMinterm(untmp, prtmp);
        elems[x] = on;
    }
}


void flipIdentityElements(std::vector <bool> &elems, unsigned seed,
        unsigned k1, unsigned k2, bool on)
{
    // Convert 'seed' to minterm
    int untmp[VARS+1];
    int prtmp[VARS+1];
    fillMinterm(seed, untmp, prtmp);
    for (unsigned v1=0; v1<RELDOM; v1++) {
        untmp[k1] = v1;
        prtmp[k1] = v1;
        for (unsigned v2=0; v2<RELDOM; v2++) {
            untmp[k2] = v2;
            prtmp[k2] = v2;
            unsigned x = whichMinterm(untmp, prtmp);
            elems[x] = on;
        }
    }
}


void randomizeFully(std::vector <bool> &elems, unsigned card)
{
    for (unsigned i=0; i<elems.size(); i++) {
        elems[i] = false;
    }
    for (unsigned i=0; i<card; i++) {
        unsigned x = Equilikely(0, elems.size()-1);
        unsigned k1 = Equilikely(1, VARS);
        unsigned k2 = Equilikely(1, VARS);

        if (k1 != k2) {
            flipFullyElements(elems, x, k1, k2, true);
        } else {
            flipFullyElements(elems, x, k1, true);
        }
    }
}


void randomizeIdentity(std::vector <bool> &elems, unsigned card)
{
    for (unsigned i=0; i<elems.size(); i++) {
        elems[i] = false;
    }
    for (unsigned i=0; i<card; i++) {
        unsigned x = Equilikely(0, elems.size()-1);
        unsigned k1 = Equilikely(1, VARS);
        unsigned k2 = Equilikely(1, VARS);

        if (k1 != k2) {
            flipIdentityElements(elems, x, k1, k2, true);
        } else {
            flipIdentityElements(elems, x, k1, true);
        }
    }
}




void set2mdd(const std::vector<bool> &S1, const std::vector<bool> &S2,
        bool s1prio, forest *F, dd_edge &s)
{
    if (!F) throw "null forest";
    if (S1.size() != S2.size()) throw "size mismatch";

    // Determine S cardinality
    unsigned card = 0;

    for (unsigned i=0; i<S1.size(); i++) {
        if (S1[i] || S2[i]) ++card;
    }

    // Special case - empty set
    if (0==card) {
        F->createEdge(false, s);
        return;
    }

    // Convert set S to list of minterms
    int** mtlist = new int* [card];
    long* vals = new long[card];
    card = 0;
    for (unsigned i=0; i<S1.size(); i++) {
        if (S1[i] || S2[i]) {
            mtlist[card] = new int[VARS+1];
            fillMinterm(i, mtlist[card]);

            if (s1prio) {
                vals[card] = S1[i] ? 1 : 2;
            } else {
                vals[card] = S2[i] ? 2 : 1;
            }

            ++card;
        }
    }

    F->createEdge(mtlist, vals, card, s);

    // Cleanup
    for (unsigned i=0; i<card; i++) {
        delete[] mtlist[i];
    }
    delete[] mtlist;
    delete[] vals;
}

/*
void set2mxd(const std::vector<bool> &S, forest *F, dd_edge &s)
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
    int** unlist = new int* [card];
    int** prlist = new int* [card];
    card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (!S[i]) continue;
        unlist[card] = new int[VARS+1];
        prlist[card] = new int[VARS+1];
        fillMinterm(i, unlist[card], prlist[card]);
        ++card;
    }

    F->createEdge(unlist, prlist, card, s);

    // Cleanup
    for (unsigned i=0; i<card; i++) {
        delete[] unlist[i];
        delete[] prlist[i];
    }
    delete[] unlist;
    delete[] prlist;
}
*/
inline char getReductionType(forest* f)
{
    if (f->isFullyReduced())    return 'f';
    if (f->isQuasiReduced())    return 'q';
    if (f->isIdentityReduced()) return 'i';
    throw "Unknown reduction type";
}

void test_sets_over(unsigned scard, forest* F)
{
    //
    // FOR NOW: just test random 'sets'
    //
    ostream_output out(std::cout);

    out << "\nGenerating random sets of size " << scard << "\n";

    std::vector <bool> ones(POTENTIAL);
    std::vector <bool> twos(POTENTIAL);

    randomizeSet(ones, scard);
    randomizeSet(twos, scard);

    dd_edge edge(F);

    set2mdd(ones, twos, true, F, edge);

    out << "random set:\n";
    edge.showGraph(out);

    if (has_common_element(ones, twos)) {
        set2mdd(ones, twos, false, F, edge);
        out << "alternate random set:\n";
        edge.showGraph(out);
    }

    randomizeFully(ones, scard);
    randomizeFully(twos, scard);

    set2mdd(ones, twos, true, F, edge);

    out << "random fully:\n";
    edge.showGraph(out);

    if (has_common_element(ones, twos)) {
        set2mdd(ones, twos, false, F, edge);
        out << "alternate random fully:\n";
        edge.showGraph(out);
    }
}

void test_sets(domain* D)
{
    using namespace MEDDLY;
    policies p;

    p.useDefaults(SET);

    forest* fully = forest::create(D, SET, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* quasi = forest::create(D, SET, range_type::INTEGER,
                        edge_labeling::MULTI_TERMINAL, p);

    test_sets_over(16, fully);

    test_sets_over(16, quasi);
}

int main()
{
    using namespace MEDDLY;

    Random(112358);

    try {
        MEDDLY::initialize();

        int bs[VARS];
        for (unsigned i=0; i<VARS; i++) {
            bs[i] = SETDOM;
        }
        domain* Ds = domain::createBottomUp(bs, VARS);
        test_sets(Ds);

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
