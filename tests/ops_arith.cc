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
    Tests element-wise arithmetic operations.
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

const char* OPS = "max, min";

using namespace MEDDLY;

// #define DEBUG_MXDOPS

#define TEST_SETS
#define TEST_RELATIONS

template <typename T>
void Max(const std::vector <T> &A, const std::vector <T> &B,
        std::vector<T> &C)
{
    if (A.size() != B.size()) {
        throw "EQ size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "EQ size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = MAX(A[i], B[i]);
    }
}

template <typename T>
void Min(const std::vector <T> &A, const std::vector <T> &B,
        std::vector<T> &C)
{
    if (A.size() != B.size()) {
        throw "NE size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "NE size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = MIN(A[i], B[i]);
    }
}

template <typename T>
void Plus(const std::vector <T> &A, const std::vector <T> &B,
        std::vector<T> &C)
{
    if (A.size() != B.size()) {
        throw "GE size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "GE size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = A[i] + B[i];
    }
}

template <typename T>
void Minus(const std::vector <T> &A, const std::vector <T> &B,
        std::vector<T> &C)
{
    if (A.size() != B.size()) {
        throw "GT size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "GT size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = A[i] - B[i];
    }
}

template <typename T>
void Times(const std::vector <T> &A, const std::vector <T> &B,
        std::vector<T> &C)
{
    if (A.size() != B.size()) {
        throw "GE size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "GE size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = A[i] * B[i];
    }
}



// TBD others

inline const char* getReductionString(const dd_edge &e)
{
    const forest* f = e.getForest();
    return shortNameOf(f->getReductionRule());
}

template <typename T>
void checkEqual(const char* what, const dd_edge &in1, const dd_edge &in2,
        const dd_edge &e1, const dd_edge &e2, const std::vector<T> &set)
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

    std::vector <T> AmaxBset(POTENTIAL);
    std::vector <T> AminBset(POTENTIAL);
    std::vector <T> AplusBset(POTENTIAL);
    std::vector <T> AminusBset(POTENTIAL);
    std::vector <T> AtimesBset(POTENTIAL);

    dd_edge Add(f1), Bdd(f2),
            AmaxBdd(fres), AminBdd(fres),
            AplusBdd(fres), AminusBdd(fres),
            AtimesBdd(fres);

    Max(Aset, Bset, AmaxBset);
    Min(Aset, Bset, AminBset);
    Plus(Aset, Bset, AplusBset);
    Minus(Aset, Bset, AminusBset);
    Times(Aset, Bset, AtimesBset);

    // TBD: for division, we need to add 1 to Bset or something

    Gen.explicit2edge(Aset, Add);
    Gen.explicit2edge(Bset, Bdd);
    Gen.explicit2edge(AmaxBset, AmaxBdd);
    Gen.explicit2edge(AminBset, AminBdd);
    Gen.explicit2edge(AplusBset, AplusBdd);
    Gen.explicit2edge(AminusBset, AminusBdd);
    Gen.explicit2edge(AtimesBset, AtimesBdd);

    dd_edge AmaxBsym(fres), AminBsym(fres),
            AplusBsym(fres), AminusBsym(fres),
            AtimesBsym(fres);

    apply(MAXIMUM, Add, Bdd, AmaxBsym);
    checkEqual("max", Add, Bdd, AmaxBsym, AmaxBdd, AmaxBset);

    apply(MINIMUM, Add, Bdd, AminBsym);
    checkEqual("min", Add, Bdd, AminBsym, AminBdd, AminBset);

    //apply(PLUS, Add, Bdd, AplusBsym);
    //checkEqual("plus", Add, Bdd, AplusBsym, AplusBdd, AplusBset);

    //apply(MINUS, Add, Bdd, AminusBsym);
    //checkEqual("minus", Add, Bdd, AminusBsym, AminusBdd, AminusBset);

    //apply(MULTIPLY, Add, Bdd, AtimesBsym);
    //checkEqual("times", Add, Bdd, AtimesBsym, AtimesBdd, AtimesBset);
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
void test_sets(domain* D, edge_labeling el, range_type rt)
{
    policies p;
    p.useDefaults(SET);

    p.setFullyReduced();

    forest* fully = forest::create(D, SET, rt, el, p);

    p.setQuasiReduced();

    forest* quasi = forest::create(D, SET, rt, el, p);


    for (unsigned i=MIN_SET_CARD; i<=MAX_SET_CARD; i*=MULT_SET_CARD) {
        std::cout << "Testing " << OPS << " on "
                  << nameOf(el) << " " << nameOf(rt) << " vectors; "
                  << i << " / " << SG.potential() << " nonzeroes\n";

        test_on_functions<TYPE>(i, fully, fully, fully);
        test_on_functions<TYPE>(i, fully, fully, quasi);
        test_on_functions<TYPE>(i, fully, quasi, fully);
        test_on_functions<TYPE>(i, fully, quasi, quasi);
        test_on_functions<TYPE>(i, quasi, fully, fully);
        test_on_functions<TYPE>(i, quasi, fully, quasi);
        test_on_functions<TYPE>(i, quasi, quasi, fully);
        test_on_functions<TYPE>(i, quasi, quasi, quasi);
    }
}

template <typename TYPE>
void test_rels(domain* D, edge_labeling el, range_type rt)
{
    policies p;
    p.useDefaults(RELATION);

    p.setFullyReduced();

    forest* fully = forest::create(D, RELATION, rt, el, p);

    p.setQuasiReduced();

    forest* quasi = forest::create(D, RELATION, rt, el, p);

    p.setIdentityReduced();

    forest* ident = forest::create(D, RELATION, rt, el, p);


    for (unsigned i=MIN_REL_CARD; i<=MAX_REL_CARD; i*=MULT_REL_CARD) {
        std::cout << "Testing " << OPS << " on "
                  << nameOf(el) << " " << nameOf(rt) << " matrices; "
                  << i << " / " << RG.potential() << " nonzeroes\n";

        test_on_functions<TYPE>(i, fully, fully, fully);
        test_on_functions<TYPE>(i, fully, fully, quasi);
        test_on_functions<TYPE>(i, fully, fully, ident);

        test_on_functions<TYPE>(i, fully, quasi, fully);
        test_on_functions<TYPE>(i, fully, quasi, quasi);
        test_on_functions<TYPE>(i, fully, quasi, ident);

        test_on_functions<TYPE>(i, fully, ident, fully);
        test_on_functions<TYPE>(i, fully, ident, quasi);
        test_on_functions<TYPE>(i, fully, ident, ident);

    //

        test_on_functions<TYPE>(i, ident, fully, fully);
        test_on_functions<TYPE>(i, ident, fully, quasi);
        test_on_functions<TYPE>(i, ident, fully, ident);

        test_on_functions<TYPE>(i, ident, quasi, fully);
        test_on_functions<TYPE>(i, ident, quasi, quasi);
        test_on_functions<TYPE>(i, ident, quasi, ident);

        test_on_functions<TYPE>(i, ident, ident, fully);
        test_on_functions<TYPE>(i, ident, ident, quasi);
        test_on_functions<TYPE>(i, ident, ident, ident);

    //

        test_on_functions<TYPE>(i, quasi, fully, fully);
        test_on_functions<TYPE>(i, quasi, fully, quasi);
        test_on_functions<TYPE>(i, quasi, fully, ident);

        test_on_functions<TYPE>(i, quasi, quasi, fully);
        test_on_functions<TYPE>(i, quasi, quasi, quasi);
        test_on_functions<TYPE>(i, quasi, quasi, ident);

        test_on_functions<TYPE>(i, quasi, ident, fully);
        test_on_functions<TYPE>(i, quasi, ident, quasi);
        test_on_functions<TYPE>(i, quasi, ident, ident);
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
        test_sets<int>(SD, edge_labeling::MULTI_TERMINAL, range_type::INTEGER);
        test_sets<float>(SD, edge_labeling::MULTI_TERMINAL, range_type::REAL);
        test_sets<long>(SD, edge_labeling::EVPLUS, range_type::INTEGER);
        domain::destroy(SD);
#endif

        //
        // Test relations
        //

#ifdef TEST_RELATIONS
        domain* RD = RG.makeDomain();
        test_rels<int>(RD, edge_labeling::MULTI_TERMINAL, range_type::INTEGER);
        test_rels<float>(RD, edge_labeling::MULTI_TERMINAL, range_type::REAL);
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
