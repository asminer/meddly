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

vectorgen SG(MEDDLY::SET, 5, 4, 1);
vectorgen RG(MEDDLY::RELATION, 5, 2, 1);

const unsigned MIN_SET_CARD = 1;
const unsigned MAX_SET_CARD = 256;
const unsigned MULT_SET_CARD = 4;


const unsigned MIN_REL_CARD = 1;
const unsigned MAX_REL_CARD = 256;
const unsigned MULT_REL_CARD = 16;

using namespace MEDDLY;

// #define DEBUG_MXDOPS

#define TEST_SETS
#define TEST_RELATIONS


template <typename T>
void EQ(const std::vector <T> &A, const std::vector <T> &B,
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

template <typename T>
void NE(const std::vector <T> &A, const std::vector <T> &B,
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

template <typename T>
void GE(const std::vector <T> &A, const std::vector <T> &B,
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

template <typename T>
void GT(const std::vector <T> &A, const std::vector <T> &B,
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

template <typename T>
void LE(const std::vector <T> &A, const std::vector <T> &B,
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

template <typename T>
void LT(const std::vector <T> &A, const std::vector <T> &B,
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

inline const char* getReductionString(const dd_edge &e)
{
    const forest* f = e.getForest();
    return shortNameOf(f->getReductionRule());
}

void checkEqual(const char* what, const dd_edge &in1, const dd_edge &in2,
        const dd_edge &e1, const dd_edge &e2, const std::vector<bool> &set)
{
    if (e1 == e2) return;

    ostream_output out(std::cout);

    out << "\nMismatch on " << what << "\n";
    out << "Input A (" << getReductionString(in1) << "):\n";
    in1.showGraph(out);
    out << "Input B (" << getReductionString(in2) << "):\n";
    in2.showGraph(out);
    out << "Expected Output (" << getReductionString(e2) << "):\n";
    e2.showGraph(out);
    out << "Obtained Output (" << getReductionString(e1) << "):\n";
    e1.showGraph(out);

    throw "mismatch";
}


template <typename T>
void compare(vectorgen &Gen,
        const std::vector <T> &Aset, const std::vector <T> &Bset,
        forest* f1, forest* f2, forest* fres)
{
    const unsigned POTENTIAL = Gen.potential();

    std::vector <bool> AeqBset(POTENTIAL);
    std::vector <bool> AneBset(POTENTIAL);
    std::vector <bool> AgtBset(POTENTIAL);
    std::vector <bool> AgeBset(POTENTIAL);
    std::vector <bool> AltBset(POTENTIAL);
    std::vector <bool> AleBset(POTENTIAL);

    dd_edge Add(f1), Bdd(f2),
            AeqBdd(fres), AneBdd(fres),
            AgtBdd(fres), AgeBdd(fres),
            AltBdd(fres), AleBdd(fres);

    EQ(Aset, Bset, AeqBset);
    NE(Aset, Bset, AneBset);
    GE(Aset, Bset, AgeBset);
    GT(Aset, Bset, AgtBset);
    LE(Aset, Bset, AleBset);
    LT(Aset, Bset, AltBset);

    Gen.explicit2edge(Aset, Add);
    Gen.explicit2edge(Bset, Bdd);
    Gen.explicit2edge(AeqBset, AeqBdd);
    Gen.explicit2edge(AneBset, AneBdd);
    Gen.explicit2edge(AgtBset, AgtBdd);
    Gen.explicit2edge(AgeBset, AgeBdd);
    Gen.explicit2edge(AltBset, AltBdd);
    Gen.explicit2edge(AleBset, AleBdd);

    dd_edge AeqBsym(fres), AneBsym(fres),
            AgtBsym(fres), AgeBsym(fres),
            AltBsym(fres), AleBsym(fres);

    apply(EQUAL,        Add, Bdd, AeqBsym);
    checkEqual("equal", Add, Bdd, AeqBsym, AeqBdd, AeqBset);

    apply(NOT_EQUAL,        Add, Bdd, AneBsym);
    checkEqual("not_equal", Add, Bdd, AneBsym, AneBdd, AneBset);

    apply(GREATER_THAN,        Add, Bdd, AgtBsym);
    checkEqual("greater_than", Add, Bdd, AgtBsym, AgtBdd, AgtBset);

    apply(GREATER_THAN_EQUAL,        Add, Bdd, AgeBsym);
    checkEqual("greater_than_equal", Add, Bdd, AgeBsym, AgeBdd, AgeBset);

    apply(LESS_THAN,        Add, Bdd, AltBsym);
    checkEqual("less_than", Add, Bdd, AltBsym, AltBdd, AltBset);

    apply(LESS_THAN_EQUAL,        Add, Bdd, AleBsym);
    checkEqual("less_than_equal", Add, Bdd, AleBsym, AleBdd, AleBset);
}

template <typename TYPE>
void test_on_functions(unsigned scard, forest* f1, forest* f2, forest* fres)
{
    if (!f1) throw "null f1";
    if (!f2) throw "null f2";
    if (!fres) throw "null fres";

    vectorgen &Gen = f1->isForRelations() ? RG : SG;

    std::cerr << "    " << shortNameOf(f1->getReductionRule())
              << ' ' << shortNameOf(f2->getReductionRule())
              << " : " << shortNameOf(fres->getReductionRule()) << ' ';

    std::vector <TYPE> Aset(Gen.potential());
    std::vector <TYPE> Bset(Gen.potential());

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        Gen.randomizeVector(Aset, scard);
        Gen.randomizeVector(Bset, scard);

        compare(Gen, Aset, Bset, f1, f2, fres);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        Gen.randomizeFully(Aset, scard);
        Gen.randomizeFully(Bset, scard);

        compare(Gen, Aset, Bset, f1, f2, fres);
    }
    if (fres->isForRelations()) {
        for (unsigned i=0; i<10; i++) {
            std::cerr << "i";
            Gen.randomizeIdentity(Aset, scard);
            Gen.randomizeIdentity(Bset, scard);

            compare(Gen, Aset, Bset, f1, f2, fres);
        }
    }
    std::cerr << std::endl;
}

template <typename TYPE>
void test_sets(domain* D, range_type rt)
{
    policies p;
    p.useDefaults(SET);

    p.setFullyReduced();

    forest* in_fully = forest::create(D, SET, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    forest* out_fully = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* in_quasi = forest::create(D, SET, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    forest* out_quasi = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);


    for (unsigned i=MIN_SET_CARD; i<=MAX_SET_CARD; i*=MULT_SET_CARD) {
        std::cout << "Testing ==, !=, >, >=, <, <= on "
                  << nameOf(rt) << " vectors; " << i << " / "
                  << SG.potential() << " nonzeroes\n";

        test_on_functions<TYPE>(i, in_fully, in_fully, out_fully);
        test_on_functions<TYPE>(i, in_fully, in_fully, out_quasi);
        test_on_functions<TYPE>(i, in_fully, in_quasi, out_fully);
        test_on_functions<TYPE>(i, in_fully, in_quasi, out_quasi);
        test_on_functions<TYPE>(i, in_quasi, in_fully, out_fully);
        test_on_functions<TYPE>(i, in_quasi, in_fully, out_quasi);
        test_on_functions<TYPE>(i, in_quasi, in_quasi, out_fully);
        test_on_functions<TYPE>(i, in_quasi, in_quasi, out_quasi);
    }
}

template <typename TYPE>
void test_rels(domain* D, range_type rt)
{
    policies p;
    p.useDefaults(RELATION);

    p.setFullyReduced();

    forest* in_fully = forest::create(D, RELATION, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    forest* out_fully = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* in_quasi = forest::create(D, RELATION, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    forest* out_quasi = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setIdentityReduced();

    forest* in_ident = forest::create(D, RELATION, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    forest* out_ident = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);



    for (unsigned i=MIN_REL_CARD; i<=MAX_REL_CARD; i*=MULT_REL_CARD) {
        std::cout << "Testing ==, !=, >, >=, <, <= on "
                  << nameOf(rt) << " matrices; " << i << " / "
                  << SG.potential() << " nonzeroes\n";

        test_on_functions<TYPE>(i, in_fully, in_fully, out_fully);
        test_on_functions<TYPE>(i, in_fully, in_fully, out_quasi);
        test_on_functions<TYPE>(i, in_fully, in_fully, out_ident);

        test_on_functions<TYPE>(i, in_fully, in_quasi, out_fully);
        test_on_functions<TYPE>(i, in_fully, in_quasi, out_quasi);
        test_on_functions<TYPE>(i, in_fully, in_quasi, out_ident);

        test_on_functions<TYPE>(i, in_fully, in_ident, out_fully);
        test_on_functions<TYPE>(i, in_fully, in_ident, out_quasi);
        test_on_functions<TYPE>(i, in_fully, in_ident, out_ident);

    //

        test_on_functions<TYPE>(i, in_ident, in_fully, out_fully);
        test_on_functions<TYPE>(i, in_ident, in_fully, out_quasi);
        test_on_functions<TYPE>(i, in_ident, in_fully, out_ident);

        test_on_functions<TYPE>(i, in_ident, in_quasi, out_fully);
        test_on_functions<TYPE>(i, in_ident, in_quasi, out_quasi);
        test_on_functions<TYPE>(i, in_ident, in_quasi, out_ident);

        test_on_functions<TYPE>(i, in_ident, in_ident, out_fully);
        test_on_functions<TYPE>(i, in_ident, in_ident, out_quasi);
        test_on_functions<TYPE>(i, in_ident, in_ident, out_ident);

    //

        test_on_functions<TYPE>(i, in_quasi, in_fully, out_fully);
        test_on_functions<TYPE>(i, in_quasi, in_fully, out_quasi);
        test_on_functions<TYPE>(i, in_quasi, in_fully, out_ident);

        test_on_functions<TYPE>(i, in_quasi, in_quasi, out_fully);
        test_on_functions<TYPE>(i, in_quasi, in_quasi, out_quasi);
        test_on_functions<TYPE>(i, in_quasi, in_quasi, out_ident);

        test_on_functions<TYPE>(i, in_quasi, in_ident, out_fully);
        test_on_functions<TYPE>(i, in_quasi, in_ident, out_quasi);
        test_on_functions<TYPE>(i, in_quasi, in_ident, out_ident);
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
    vectorgen::setSeed(seed);


    try {
        MEDDLY::initialize();

        //
        // Test sets
        //

#ifdef TEST_SETS
        domain* SD = SG.makeDomain();
        test_sets<int>(SD, range_type::INTEGER);
        test_sets<float>(SD, range_type::REAL);
        domain::destroy(SD);
#endif

        //
        // Test relations
        //

#ifdef TEST_RELATIONS
        domain* RD = RG.makeDomain();
        test_rels<int>(RD, range_type::INTEGER);
        test_rels<float>(RD, range_type::REAL);
        domain::destroy(RD);
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
