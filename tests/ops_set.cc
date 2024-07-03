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
const unsigned MAX_SET_CARD = 512;
const unsigned MAX_REL_CARD = 64;

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

    // Cleanup
    for (unsigned i=0; i<card; i++) {
        delete[] mtlist[i];
    }
    delete[] mtlist;
}

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


inline char getReductionType(forest* f)
{
    if (f->isFullyReduced())    return 'f';
    if (f->isQuasiReduced())    return 'q';
    if (f->isIdentityReduced()) return 'i';
    throw "Unknown reduction type";
}

void checkEqual(const char* what, const dd_edge &e1, const dd_edge &e2,
        const std::vector<bool> &set)
{
    if (e1 == e2) return;

    ostream_output out(std::cout);

    out << "\nMismatch on " << what << "\n";
    out << "Expected DD:\n";
    e2.showGraph(out);
    out << "Obtained DD:\n";
    e1.showGraph(out);

    /*
    std::cout << "Encoding set ";
    showSet(std::cout, set);
    std::cout << "\nMinterms\n    ";
    showRelMinterms(std::cout, set);
    */

    throw "mismatch";
}

void compare_sets(const std::vector <bool> &Aset,
        const std::vector <bool> &Bset, forest* f1, forest* f2, forest* fres)
{
    std::vector <bool> AiBset(POTENTIAL);
    std::vector <bool> AuBset(POTENTIAL);
    std::vector <bool> AmBset(POTENTIAL);
    std::vector <bool> cAset(POTENTIAL);

    dd_edge Add(f1), Bdd(f2), AiBdd(fres), AuBdd(fres), AmBdd(fres),
                cABdd(fres);

    set_intersection(Aset, Bset, AiBset);
    set_union(Aset, Bset, AuBset);
    set_difference(Aset, Bset, AmBset);
    set_complement(Aset, cAset);

    set2mdd(Aset, f1, Add);
    set2mdd(Bset, f2, Bdd);
    set2mdd(AiBset, fres, AiBdd);
    set2mdd(AuBset, fres, AuBdd);
    set2mdd(AmBset, fres, AmBdd);
    set2mdd(cAset, fres, cABdd);

    dd_edge AiBsym(fres), AuBsym(fres), AmBsym(fres), cAsym(fres);

    apply(INTERSECTION, Add, Bdd, AiBsym);
    apply(UNION, Add, Bdd, AuBsym);
    apply(DIFFERENCE, Add, Bdd, AmBsym);
    apply(COMPLEMENT, Add, cAsym);

    checkEqual("intersection", AiBsym, AiBdd, AiBset);
    checkEqual("union", AuBsym, AuBdd, AuBset);
    checkEqual("difference", AmBsym, AmBdd, AmBset);
    checkEqual("complement", cAsym, cABdd, cAset);
}

void compare_rels(const std::vector <bool> &Aset,
        const std::vector <bool> &Bset, forest* f1, forest* f2, forest* fres)
{
    std::vector <bool> AiBset(POTENTIAL);
    std::vector <bool> AuBset(POTENTIAL);
    std::vector <bool> AmBset(POTENTIAL);
    std::vector <bool> cAset(POTENTIAL);

    dd_edge Add(f1), Bdd(f2), AiBdd(fres), AuBdd(fres), AmBdd(fres),
                cABdd(fres);

    set_intersection(Aset, Bset, AiBset);
    set_union(Aset, Bset, AuBset);
    set_difference(Aset, Bset, AmBset);
    set_complement(Aset, cAset);

    set2mxd(Aset, f1, Add);
    set2mxd(Bset, f2, Bdd);
    set2mxd(AiBset, fres, AiBdd);
    set2mxd(AuBset, fres, AuBdd);
    set2mxd(AmBset, fres, AmBdd);
    set2mxd(cAset, fres, cABdd);

    dd_edge AiBsym(fres), AuBsym(fres), AmBsym(fres), cAsym(fres);

#ifdef DEBUG_MXDOPS
    ostream_output out(std::cout);
    std::cout << "==================================================================\n";
    std::cout << "Sets:\n";
    std::cout << "    A: ";
    showSet(std::cout, Aset);
    std::cout << "\n  = ";
    showRelMinterms(std::cout, Aset);
    std::cout << "\n    B: ";
    showSet(std::cout, Bset);
    std::cout << "\n  = ";
    showRelMinterms(std::cout, Bset);
    std::cout << "\nMXDs:\n";
    std::cout << "    A:\n";
    Add.showGraph(out);
    std::cout << "    B:\n";
    Bdd.showGraph(out);
#endif

    apply(INTERSECTION, Add, Bdd, AiBsym);
    apply(UNION, Add, Bdd, AuBsym);
    apply(DIFFERENCE, Add, Bdd, AmBsym);
    apply(COMPLEMENT, Add, cAsym);

    checkEqual("intersection", AiBsym, AiBdd, AiBset);
    checkEqual("union", AuBsym, AuBdd, AuBset);
    checkEqual("difference", AmBsym, AmBdd, AmBset);
    checkEqual("complement", cAsym, cABdd, cAset);
}


void test_sets_over(unsigned scard, forest* f1, forest* f2, forest* fres)
{
    if (!f1) throw "null f1";
    if (!f2) throw "null f2";
    if (!fres) throw "null fres";

    std::cerr << "    " << getReductionType(f1) << getReductionType(f2)
              << ':' << getReductionType(fres) << ' ';

    std::vector <bool> Aset(POTENTIAL);
    std::vector <bool> Bset(POTENTIAL);

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        randomizeSet(Aset, scard);
        randomizeSet(Bset, scard);

        compare_sets(Aset, Bset, f1, f2, fres);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        randomizeFully(Aset, scard);
        randomizeFully(Bset, scard);

        compare_sets(Aset, Bset, f1, f2, fres);
    }

    std::cerr << std::endl;
}

void test_sets(domain* D)
{
    policies p;
    p.useDefaults(SET);

    forest* F1 = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* F2 = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    for (unsigned i=1; i<=MAX_SET_CARD; i*=2) {
        std::cout << "Testing sets of size " << i << " out of " << POTENTIAL << "\n";

        test_sets_over(i, F1, F1, F1);
        test_sets_over(i, F1, F1, F2);
        test_sets_over(i, F1, F2, F1);
        test_sets_over(i, F1, F2, F2);
        test_sets_over(i, F2, F1, F1);
        test_sets_over(i, F2, F1, F2);
        test_sets_over(i, F2, F2, F1);
        test_sets_over(i, F2, F2, F2);
    }
}

void test_rels_over(unsigned scard, forest* f1, forest* f2, forest* fres)
{
    if (!f1) throw "null f1";
    if (!f2) throw "null f2";
    if (!fres) throw "null fres";

    std::cerr << "    " << getReductionType(f1) << getReductionType(f2)
              << ':' << getReductionType(fres) << ' ';

    std::vector <bool> Aset(POTENTIAL);
    std::vector <bool> Bset(POTENTIAL);

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        randomizeSet(Aset, scard);
        randomizeSet(Bset, scard);

        compare_rels(Aset, Bset, f1, f2, fres);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        randomizeFully(Aset, scard);
        randomizeFully(Bset, scard);

        compare_rels(Aset, Bset, f1, f2, fres);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "i";
        randomizeIdentity(Aset, scard);
        randomizeIdentity(Bset, scard);

        compare_rels(Aset, Bset, f1, f2, fres);
    }
    std::cerr << std::endl;
}

void test_rels(domain* D)
{
    policies p;
    p.useDefaults(RELATION);

    p.setQuasiReduced();

    forest* R1 = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

    p.setFullyReduced();

    forest* R2 = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

    p.setIdentityReduced();

    forest* R3 = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);


    for (unsigned i=1; i<=MAX_REL_CARD; i*=2) {
        std::cout << "Testing relations of size " << i << " out of " << POTENTIAL << "\n";

        test_rels_over(i, R1, R1, R1);
        test_rels_over(i, R1, R1, R2);
        test_rels_over(i, R1, R1, R3);

        test_rels_over(i, R1, R2, R1);
        test_rels_over(i, R1, R2, R2);
        test_rels_over(i, R1, R2, R3);

        test_rels_over(i, R1, R3, R1);
        test_rels_over(i, R1, R3, R2);
        test_rels_over(i, R1, R3, R3);
    // ---
        test_rels_over(i, R2, R1, R1);
        test_rels_over(i, R2, R1, R2);
        test_rels_over(i, R2, R1, R3);

        test_rels_over(i, R2, R2, R1);
        test_rels_over(i, R2, R2, R2);
        test_rels_over(i, R2, R2, R3);

        test_rels_over(i, R2, R3, R1);
        test_rels_over(i, R2, R3, R2);
        test_rels_over(i, R2, R3, R3);
    // ---
        test_rels_over(i, R3, R1, R1);
        test_rels_over(i, R3, R1, R2);
        test_rels_over(i, R3, R1, R3);

        test_rels_over(i, R3, R2, R1);
        test_rels_over(i, R3, R2, R2);
        test_rels_over(i, R3, R2, R3);

        test_rels_over(i, R3, R3, R1);
        test_rels_over(i, R3, R3, R2);
        test_rels_over(i, R3, R3, R3);
    }
}

int main()
{
    Random(11235);

    try {
        MEDDLY::initialize();
        policies p;

        //
        // Test sets
        //

#ifdef TEST_SETS
        int bs[VARS];
        for (unsigned i=0; i<VARS; i++) {
            bs[i] = SETDOM;
        }
        domain* Ds = domain::createBottomUp(bs, VARS);

        test_sets(Ds);

        domain::destroy(Ds);
#endif

        //
        // Test relations
        //

#ifdef TEST_RELATIONS
        int br[VARS];
        for (unsigned i=0; i<VARS; i++) {
            br[i] = RELDOM;
        }
        domain* Dr = domain::createBottomUp(br, VARS);

        test_rels(Dr);

        domain::destroy(Dr);

#endif

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
