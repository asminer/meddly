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

void showSet(std::ostream &out, const std::vector <char> &elems)
{
    out << "{ ";
    bool printed = false;
    for (unsigned i=0; i<elems.size(); i++) {
        if (!elems[i]) continue;
        if (printed) out << ", ";
        out << i << ":" << int(elems[i]);
        printed = true;
    }
    out << " }";
}

void showRelMinterms(std::ostream &out, const std::vector <char> &elems)
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
        out << "]:" << int(elems[i]);
    }
    out << " }";
}


void EQ(const std::vector <char> &A, const std::vector <char> &B,
        std::vector<bool> &C)
{
    if (A.size() != B.size()) {
        throw "EQ size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "EQ size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = (A[i] == B[i]);
    }
}

void NE(const std::vector <char> &A, const std::vector <char> &B,
        std::vector<bool> &C)
{
    if (A.size() != B.size()) {
        throw "NE size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "NE size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = (A[i] != B[i]);
    }
}

void GE(const std::vector <char> &A, const std::vector <char> &B,
        std::vector<bool> &C)
{
    if (A.size() != B.size()) {
        throw "GE size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "GE size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = (A[i] >= B[i]);
    }
}

void GT(const std::vector <char> &A, const std::vector <char> &B,
        std::vector<bool> &C)
{
    if (A.size() != B.size()) {
        throw "GT size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "GT size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = (A[i] > B[i]);
    }
}

void LE(const std::vector <char> &A, const std::vector <char> &B,
        std::vector<bool> &C)
{
    if (A.size() != B.size()) {
        throw "LE size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "LE size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = (A[i] <= B[i]);
    }
}

void LT(const std::vector <char> &A, const std::vector <char> &B,
        std::vector<bool> &C)
{
    if (A.size() != B.size()) {
        throw "LT size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "LT size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = (A[i] < B[i]);
    }
}

void randomizeSet(std::vector <char> &elems, unsigned card, char vals)
{
    //
    // Fill array with card 1s, card 2s, etc.
    //
    unsigned i=0;
    for (unsigned v=1; v<=vals; v++) {
        for (unsigned j=0; j<card; j++) {
            if (i>=elems.size()) break;
            elems[i++] = v;
        }
    }
    for (; i<elems.size(); i++) {
        elems[i] = 0;
    }

    //
    // Shuffle the array
    //
    for (unsigned i=0; i<elems.size()-1; i++) {
        unsigned j = Equilikely(i, elems.size()-1);
        if (elems[i] != elems[j]) {
            // swap
            char t = elems[i];
            elems[i] = elems[j];
            elems[j] = t;
        }
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

void flipFullyElements(std::vector <char> &elems, unsigned seed, unsigned k,
        char val)
{
    // Convert 'seed' to minterm
    int tmp[VARS+1];
    fillMinterm(seed, tmp);
    // Replace variable k with all possible values in minterm,
    // and convert back to unsigned
    for (tmp[k]=0; tmp[k]<SETDOM; tmp[k]++) {
        unsigned x = whichMinterm(tmp);
        elems[x] = val;
    }
}

void flipFullyElements(std::vector <char> &elems, unsigned seed, unsigned k1,
        unsigned k2, char val)
{
    // Convert 'seed' to minterm
    int tmp[VARS+1];
    fillMinterm(seed, tmp);
    // Replace variable k with all possible values in minterm,
    // and convert back to unsigned
    for (tmp[k1]=0; tmp[k1]<SETDOM; tmp[k1]++) {
        for (tmp[k2]=0; tmp[k2]<SETDOM; tmp[k2]++) {
            unsigned x = whichMinterm(tmp);
            elems[x] = val;
        }
    }
}


void flipIdentityElements(std::vector <char> &elems, unsigned seed,
        unsigned k, char val)
{
    // Convert 'seed' to minterm
    int untmp[VARS+1];
    int prtmp[VARS+1];
    fillMinterm(seed, untmp, prtmp);
    for (unsigned v=0; v<RELDOM; v++) {
        untmp[k] = v;
        prtmp[k] = v;
        unsigned x = whichMinterm(untmp, prtmp);
        elems[x] = val;
    }
}


void flipIdentityElements(std::vector <char> &elems, unsigned seed,
        unsigned k1, unsigned k2, char val)
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
            elems[x] = val;
        }
    }
}


void randomizeFully(std::vector <char> &elems, unsigned card, unsigned vals)
{
    for (unsigned i=0; i<elems.size(); i++) {
        elems[i] = false;
    }
    unsigned v=0;
    for (unsigned i=0; i<card; i++) {
        v = (v % vals)+1;
        unsigned x = Equilikely(0, elems.size()-1);
        unsigned k1 = Equilikely(1, VARS);
        unsigned k2 = Equilikely(1, VARS);

        if (k1 != k2) {
            flipFullyElements(elems, x, k1, k2, v);
        } else {
            flipFullyElements(elems, x, k1, v);
        }
    }
}


void randomizeIdentity(std::vector <char> &elems, unsigned card, unsigned vals)
{
    for (unsigned i=0; i<elems.size(); i++) {
        elems[i] = false;
    }
    unsigned v=0;
    for (unsigned i=0; i<card; i++) {
        v = (v % vals)+1;
        unsigned x = Equilikely(0, elems.size()-1);
        unsigned k1 = Equilikely(1, VARS);
        unsigned k2 = Equilikely(1, VARS);

        if (k1 != k2) {
            flipIdentityElements(elems, x, k1, k2, v);
        } else {
            flipIdentityElements(elems, x, k1, v);
        }
    }
}




void set2mdd(const std::vector<char> &S, forest *F, dd_edge &s)
{
    if (!F) throw "null forest";

    // Determine S cardinality
    unsigned card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (S[i]) ++card;
    }

    // Special case - zero function
    if (0==card) {
        F->createEdge(0L, s);
        return;
    }

    // Convert set S to list of minterms
    int** mtlist = new int* [card];
    long* vals = new long[card];
    card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (!S[i]) continue;
        vals[card] = S[i];
        mtlist[card] = new int[VARS+1];
        fillMinterm(i, mtlist[card]);

        ++card;
    }

    F->createEdge(mtlist, vals, card, s);

    // Cleanup
    for (unsigned i=0; i<card; i++) {
        delete[] mtlist[i];
    }
    delete[] vals;
    delete[] mtlist;
}

void set2mxd(const std::vector<char> &S, forest *F, dd_edge &s)
{
    if (!F) throw "null forest";

    // Determine S cardinality
    unsigned card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (S[i]) ++card;
    }

    // Special case - zero function
    if (0==card) {
        F->createEdge(0L, s);
        return;
    }

    // Convert set S to list of minterms
    int** unlist = new int* [card];
    int** prlist = new int* [card];
    long* vals = new long[card];
    card = 0;
    for (unsigned i=0; i<S.size(); i++) {
        if (!S[i]) continue;
        vals[card] = S[i];
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
    delete[] vals;
    delete[] unlist;
    delete[] prlist;
}

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

    std::vector <char> Avec(POTENTIAL);
    dd_edge Add(F);

    randomizeSet(Avec, scard, 3);

    set2mdd(Avec, F, Add);

    out << "random set:\n";
    Add.showGraph(out);

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

    test_sets_over(1, fully);
    test_sets_over(2, fully);
    test_sets_over(4, fully);
    test_sets_over(8, fully);

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
