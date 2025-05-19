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

vectorgen SG(MEDDLY::SET, 5, 4);
vectorgen RG(MEDDLY::RELATION, 5, 2);

const unsigned MIN_SET_CARD = 1;
const unsigned MAX_SET_CARD = 256;
const unsigned MULT_SET_CARD = 2;

const unsigned MIN_REL_CARD = 1;
const unsigned MAX_REL_CARD = 256;
const unsigned MULT_REL_CARD = 4;

using namespace MEDDLY;

// #define DEBUG_OPS

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
    long ecard = -1;
    apply(CARDINALITY, e, ecard);

    long scard = set_cardinality(set);

    if (ecard == scard) return;

    ostream_output out(std::cout);

    out << "\nCardinality mismatch\n";
    out << "  Explicit set: " << scard << "\n";
    out << "DD cardinality: " << ecard << "\n";

    // showMinterms(std::cout, set, e.getForest());
    e.showGraph(out);

    throw "mismatch";
}

void checkIndexing(const dd_edge &the_set, const dd_edge &the_indexes)
{
    //
    // Iterate through the set, and make sure indexes match.
    //
    long count = 0;
    for (dd_edge::iterator s = the_set.begin(); s; ++s)
    {
        long index;
        the_indexes.evaluate(*s, index);
        if (count != index) {
            throw "check indexing mismatch: set element missing index";
        }
        ++count;
    }
    //
    // Iterate through the index, and make sure everything
    // in the index is also in the set.
    //
    for (dd_edge::iterator s = the_indexes.begin(); s; ++s)
    {
        bool ok;
        the_set.evaluate(*s, ok);
        if (!ok) {
            throw "check indexing mismatch: index element not in set";
        }
    }
}

void compare(vectorgen &Gen,
        const std::vector <bool> &Aset, const std::vector <bool> &Bset,
        forest* f1, forest* f2, forest* fres, forest* findex)
{
    const unsigned POTENTIAL = Gen.potential();

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

    Gen.explicit2edge(Aset, Add);
    Gen.explicit2edge(Bset, Bdd);
    Gen.explicit2edge(AiBset, AiBdd);
    Gen.explicit2edge(AuBset, AuBdd);
    Gen.explicit2edge(AmBset, AmBdd);
    Gen.explicit2edge(cAset, cABdd);
    Gen.explicit2edge(Aset, ccABdd);

    checkCardinality(Add, Aset);
    checkCardinality(Bdd, Bset);

    if (findex) {
        dd_edge Aindex(findex);
        apply(CONVERT_TO_INDEX_SET, Add, Aindex);
        checkIndexing(Add, Aindex);

        dd_edge Bindex(findex);
        apply(CONVERT_TO_INDEX_SET, Bdd, Bindex);
        checkIndexing(Bdd, Bindex);
    }

    dd_edge AiBsym(fres), AuBsym(fres), AmBsym(fres),
            cAsym(fres), ccAsym(fres);

#ifdef DEBUG_OPS
    ostream_output out(std::cout);
    std::cout << "==================================================================\n";
    std::cout << "Sets:\n";
    std::cout << "    A: ";
    Gen.showSet(std::cout, Aset);
    std::cout << "\n  = ";
    Gen.showMinterms(std::cout, Aset);
    std::cout << "\n    B: ";
    Gen.showSet(std::cout, Bset);
    std::cout << "\n  = ";
    Gen.showMinterms(std::cout, Bset);
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


void test_on_forests(unsigned scard, forest* f1, forest* f2, forest* fres,
        forest* findex)
{
    if (!f1) throw "null f1";
    if (!f2) throw "null f2";
    if (!fres) throw "null fres";

    vectorgen &Gen = f1->isForRelations() ? RG : SG;

    std::cerr << "    " << shortNameOf(f1->getReductionRule())
              << " " << shortNameOf(f2->getReductionRule())
              << " : " << shortNameOf(fres->getReductionRule()) << ' ';

    std::vector <bool> Aset(Gen.potential());
    std::vector <bool> Bset(Gen.potential());

    std::vector <bool> values(1);
    values[0] = true;

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        Gen.randomizeVector(Aset, scard, values);
        Gen.randomizeVector(Bset, scard, values);

        compare(Gen, Aset, Bset, f1, f2, fres, findex);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        Gen.randomizeFully(Aset, scard, values);
        Gen.randomizeFully(Bset, scard, values);

        compare(Gen, Aset, Bset, f1, f2, fres, findex);
    }
    if (fres->isForRelations()) {
        for (unsigned i=0; i<10; i++) {
            std::cerr << "i";
            Gen.randomizeIdentity(Aset, scard, values);
            Gen.randomizeIdentity(Bset, scard, values);

            compare(Gen, Aset, Bset, f1, f2, fres, findex);
        }
    }
    std::cerr << std::endl;
}

void test_sets(domain* D)
{
    policies p;
    p.useDefaults(SET);

    forest* Findex = forest::create(D, SET, range_type::INTEGER,
                        edge_labeling::INDEX_SET, p);

    forest* F1 = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* F2 = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    for (unsigned i=MIN_SET_CARD; i<=MAX_SET_CARD; i*=MULT_SET_CARD)
    {
        std::cout << "Testing set operations on sets of size " << i
                  << " out of " << SG.potential() << "\n";

        test_on_forests(i, F1, F1, F1, Findex);
        test_on_forests(i, F1, F1, F2, Findex);
        test_on_forests(i, F1, F2, F1, Findex);
        test_on_forests(i, F1, F2, F2, Findex);
        test_on_forests(i, F2, F1, F1, Findex);
        test_on_forests(i, F2, F1, F2, Findex);
        test_on_forests(i, F2, F2, F1, Findex);
        test_on_forests(i, F2, F2, F2, Findex);
    }
    std::cout << "\nDone testing set operations:\n";
    std::cout << "    cardinality\n";
    std::cout << "    set indexing\n";
    std::cout << "    complement\n";
    std::cout << "    difference\n";
    std::cout << "    intersection\n";
    std::cout << "    union\n";
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


    for (unsigned i=MIN_REL_CARD; i<=MAX_REL_CARD; i*=MULT_REL_CARD)
    {
        std::cout << "Testing set operations on relations of size " << i
                  << " out of " << RG.potential() << "\n";

        test_on_forests(i, R1, R1, R1, nullptr);
        test_on_forests(i, R1, R1, R2, nullptr);
        test_on_forests(i, R1, R1, R3, nullptr);

        test_on_forests(i, R1, R2, R1, nullptr);
        test_on_forests(i, R1, R2, R2, nullptr);
        test_on_forests(i, R1, R2, R3, nullptr);

        test_on_forests(i, R1, R3, R1, nullptr);
        test_on_forests(i, R1, R3, R2, nullptr);
        test_on_forests(i, R1, R3, R3, nullptr);
    // ---
        test_on_forests(i, R2, R1, R1, nullptr);
        test_on_forests(i, R2, R1, R2, nullptr);
        test_on_forests(i, R2, R1, R3, nullptr);

        test_on_forests(i, R2, R2, R1, nullptr);
        test_on_forests(i, R2, R2, R2, nullptr);
        test_on_forests(i, R2, R2, R3, nullptr);

        test_on_forests(i, R2, R3, R1, nullptr);
        test_on_forests(i, R2, R3, R2, nullptr);
        test_on_forests(i, R2, R3, R3, nullptr);
    // ---
        test_on_forests(i, R3, R1, R1, nullptr);
        test_on_forests(i, R3, R1, R2, nullptr);
        test_on_forests(i, R3, R1, R3, nullptr);

        test_on_forests(i, R3, R2, R1, nullptr);
        test_on_forests(i, R3, R2, R2, nullptr);
        test_on_forests(i, R3, R2, R3, nullptr);

        test_on_forests(i, R3, R3, R1, nullptr);
        test_on_forests(i, R3, R3, R2, nullptr);
        test_on_forests(i, R3, R3, R3, nullptr);
    }
    std::cout << "\nDone testing relation operations:\n";
    std::cout << "    cardinality\n";
    std::cout << "    complement\n";
    std::cout << "    difference\n";
    std::cout << "    intersection\n";
    std::cout << "    union\n";
}

void usage(const char* arg0)
{
    /* Strip leading directory, if any: */
    const char* name = arg0;
    for (const char* ptr=arg0; *ptr; ptr++) {
        if ('/' == *ptr) name = ptr+1;
    }
    std::cerr << "\nUsage: " << name << "options seed\n\n";
    std::cerr << "Options:\n";
    std::cerr << "    --set:    Test sets (default).\n";
    std::cerr << "    --rel:    Test relations\n";
    std::cerr << "\n";

    exit(1);
}


int main(int argc, const char** argv)
{
    /*
     * Command-line options
     */
    bool sets = true;
    long seed = 0;

    for (int i=1; i<argc; i++) {

        if (0==strcmp("--set", argv[i])) {
            sets = true;
            continue;
        }
        if (0==strcmp("--rel", argv[i])) {
            sets = false;
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
        if (sets) {
            domain* SD = SG.makeDomain();
            test_sets(SD);
            domain::destroy(SD);
        } else {
            domain* RD = RG.makeDomain();
            test_rels(RD);
            domain::destroy(RD);
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
