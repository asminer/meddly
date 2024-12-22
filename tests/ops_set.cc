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

#include "randomize.h"

vectorgen<bool> SG(MEDDLY::SET, 5, 4, 1);
vectorgen<bool> RG(MEDDLY::RELATION, 5, 2, 1);

const unsigned MAX_SET_CARD = 512;
const unsigned MAX_REL_CARD = 64;

using namespace MEDDLY;

// #define DEBUG_OPS

#define TEST_SETS
#define TEST_RELATIONS

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

void showMinterms(std::ostream &out, const std::vector <bool> &elems,
        forest* f)
{
    minterm mt(f);
    out << "{ ";
    bool printed = false;
    for (unsigned i=0; i<elems.size(); i++) {
        if (!elems[i]) continue;

        if (printed) out << ",\n      ";

        if (f->isForRelations()) {
            RG.index2minterm(i, mt);
        } else {
            SG.index2minterm(i, mt);
        }
        ostream_output mout(out);
        mt.show(mout);
    }
    out << " }";
}

void set2edge(const std::vector<bool> &S, forest *F, dd_edge &s)
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
    minterm_coll mtlist(card, F);
    for (unsigned i=0; i<S.size(); i++) {
        if (!S[i]) continue;
        if (F->isForRelations()) {
            RG.index2minterm(i, mtlist.unused());
        } else {
            SG.index2minterm(i, mtlist.unused());
        }
        mtlist.pushUnused();
    }
    mtlist.buildFunction(s);
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

    throw "mismatch";
}

void checkCardinality(const dd_edge &e, const std::vector<bool> &set)
{
    long ecard;
    apply(CARDINALITY, e, ecard);

    long scard = set_cardinality(set);

    if (ecard == scard) return;

    ostream_output out(std::cout);

    out << "\nCardinality mismatch\n";
    out << "  Explicit set: " << scard << "\n";
    out << "DD cardinality: " << ecard << "\n";

    showMinterms(std::cout, set, e.getForest());
    e.showGraph(out);

    throw "mismatch";
}

void compare(const std::vector <bool> &Aset, const std::vector <bool> &Bset,
        forest* f1, forest* f2, forest* fres)
{
    const unsigned POTENTIAL = Aset.size();

    std::vector <bool> AiBset(POTENTIAL);
    std::vector <bool> AuBset(POTENTIAL);
    std::vector <bool> AmBset(POTENTIAL);
    std::vector <bool> cAset(POTENTIAL);

    dd_edge Add(f1), Bdd(f2), AiBdd(fres), AuBdd(fres), AmBdd(fres),
                cABdd(fres), ccABdd(fres);

    set_intersection(Aset, Bset, AiBset);
    set_union(Aset, Bset, AuBset);
    set_difference(Aset, Bset, AmBset);
    set_complement(Aset, cAset);

    set2edge(Aset, f1, Add);
    set2edge(Bset, f2, Bdd);
    set2edge(AiBset, fres, AiBdd);
    set2edge(AuBset, fres, AuBdd);
    set2edge(AmBset, fres, AmBdd);
    set2edge(cAset, fres, cABdd);
    set2edge(Aset, fres, ccABdd);

    checkCardinality(Add, Aset);
    checkCardinality(Bdd, Bset);

    dd_edge AiBsym(fres), AuBsym(fres), AmBsym(fres),
            cAsym(fres), ccAsym(fres);

#ifdef DEBUG_OPS
    ostream_output out(std::cout);
    std::cout << "==================================================================\n";
    std::cout << "Sets:\n";
    std::cout << "    A: ";
    showSet(std::cout, Aset);
    std::cout << "\n  = ";
    showMinterms(std::cout, Aset, f1);
    std::cout << "\n    B: ";
    showSet(std::cout, Bset);
    std::cout << "\n  = ";
    showMinterms(std::cout, Bset, f2);
    std::cout << "\nMDD/MXDs:\n";
    std::cout << "    A:\n";
    Add.showGraph(out);
    std::cout << "    B:\n";
    Bdd.showGraph(out);
#endif

    apply(INTERSECTION, Add, Bdd, AiBsym);
    apply(UNION, Add, Bdd, AuBsym);
    apply(DIFFERENCE, Add, Bdd, AmBsym);
    apply(COMPLEMENT, Add, cAsym);
    apply(COMPLEMENT, cAsym, ccAsym);

    checkEqual("intersection", AiBsym, AiBdd, AiBset);
    checkEqual("union", AuBsym, AuBdd, AuBset);
    checkEqual("difference", AmBsym, AmBdd, AmBset);
    checkEqual("complement", cAsym, cABdd, cAset);
    checkEqual("complement^2", ccAsym, ccABdd, Aset);
}


void test_on_forests(unsigned scard, forest* f1, forest* f2, forest* fres)
{
    if (!f1) throw "null f1";
    if (!f2) throw "null f2";
    if (!fres) throw "null fres";

    vectorgen <bool> &Gen = f1->isForRelations() ? RG : SG;

    std::cerr << "    " << getReductionType(f1) << getReductionType(f2)
              << ':' << getReductionType(fres) << ' ';

    std::vector <bool> Aset(Gen.potential());
    std::vector <bool> Bset(Gen.potential());

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        Gen.randomizeVector(Aset, scard);
        Gen.randomizeVector(Bset, scard);

        compare(Aset, Bset, f1, f2, fres);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        Gen.randomizeFully(Aset, scard);
        Gen.randomizeFully(Bset, scard);

        compare(Aset, Bset, f1, f2, fres);
    }
    if (fres->isForRelations()) {
        for (unsigned i=0; i<10; i++) {
            std::cerr << "i";
            Gen.randomizeIdentity(Aset, scard);
            Gen.randomizeIdentity(Bset, scard);

            compare(Aset, Bset, f1, f2, fres);
        }
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

    for (unsigned i=1; i<=MAX_SET_CARD; i*=2)
    {
        std::cout << "Testing set operations on sets of size " << i
                  << " out of " << SG.potential() << "\n";

        test_on_forests(i, F1, F1, F1);
        test_on_forests(i, F1, F1, F2);
        test_on_forests(i, F1, F2, F1);
        test_on_forests(i, F1, F2, F2);
        test_on_forests(i, F2, F1, F1);
        test_on_forests(i, F2, F1, F2);
        test_on_forests(i, F2, F2, F1);
        test_on_forests(i, F2, F2, F2);
    }
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
        std::cout << "Testing set operations on relations of size " << i
                  << " out of " << RG.potential() << "\n";

        test_on_forests(i, R1, R1, R1);
        test_on_forests(i, R1, R1, R2);
        test_on_forests(i, R1, R1, R3);

        test_on_forests(i, R1, R2, R1);
        test_on_forests(i, R1, R2, R2);
        test_on_forests(i, R1, R2, R3);

        test_on_forests(i, R1, R3, R1);
        test_on_forests(i, R1, R3, R2);
        test_on_forests(i, R1, R3, R3);
    // ---
        test_on_forests(i, R2, R1, R1);
        test_on_forests(i, R2, R1, R2);
        test_on_forests(i, R2, R1, R3);

        test_on_forests(i, R2, R2, R1);
        test_on_forests(i, R2, R2, R2);
        test_on_forests(i, R2, R2, R3);

        test_on_forests(i, R2, R3, R1);
        test_on_forests(i, R2, R3, R2);
        test_on_forests(i, R2, R3, R3);
    // ---
        test_on_forests(i, R3, R1, R1);
        test_on_forests(i, R3, R1, R2);
        test_on_forests(i, R3, R1, R3);

        test_on_forests(i, R3, R2, R1);
        test_on_forests(i, R3, R2, R2);
        test_on_forests(i, R3, R2, R3);

        test_on_forests(i, R3, R3, R1);
        test_on_forests(i, R3, R3, R2);
        test_on_forests(i, R3, R3, R3);
    }
}

int main(int argc, const char** argv)
{
    using namespace std;

    //
    // First argument: seed
    //
    long seed = 0;
    if (argv[1]) {
        seed = atol(argv[1]);
    }
    vectorgen_base::setSeed(seed);

    try {
        MEDDLY::initialize();

        //
        // Test sets
        //

#ifdef TEST_SETS
        domain* SD = SG.makeDomain();
        test_sets(SD);
        domain::destroy(SD);
#endif

        //
        // Test relations
        //

#ifdef TEST_RELATIONS
        domain* RD = RG.makeDomain();
        test_rels(RD);
        domain::destroy(RD);
#endif

        MEDDLY::cleanup();

        std::cout << "\nDone testing set operations:\n";
        std::cout << "    cardinality\n";
        std::cout << "    complement\n";
        std::cout << "    difference\n";
        std::cout << "    intersection\n";
        std::cout << "    union\n";
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
