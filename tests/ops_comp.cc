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
const unsigned MULT_SET_CARD = 4;

const unsigned MIN_REL_CARD = 1;
const unsigned MAX_REL_CARD = 256;
const unsigned MULT_REL_CARD = 16;

const char* OPS = "==, !=, >, >=, <, <=";

using namespace MEDDLY;

// #define DEBUG_MXDOPS

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

    std::vector <TYPE> values(4);
    values[0] =  6;
    values[1] =  4;
    values[2] =  2;
    values[3] = -2;

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        Gen.randomizeVector(Aset, scard, values);
        Gen.randomizeVector(Bset, scard, values);

        compare(Gen, Aset, Bset, f1, f2, fres);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        Gen.randomizeFully(Aset, scard, values);
        Gen.randomizeFully(Bset, scard, values);

        compare(Gen, Aset, Bset, f1, f2, fres);
    }
    if (fres->isForRelations()) {
        for (unsigned i=0; i<10; i++) {
            std::cerr << "i";
            Gen.randomizeIdentity(Aset, scard, values);
            Gen.randomizeIdentity(Bset, scard, values);

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

    forest* in_fully = forest::create(D, SET, rt, el, p);

    forest* out_fully = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* in_quasi = forest::create(D, SET, rt, el, p);

    forest* out_quasi = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);


    for (unsigned i=MIN_SET_CARD; i<=MAX_SET_CARD; i*=MULT_SET_CARD) {
        std::cout << "Testing " << OPS << " on "
                  << nameOf(el) << " " << nameOf(rt) << " vectors; "
                  << i << " / " << RG.potential() << " nonzeroes\n";

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
void test_rels(domain* D, edge_labeling el, range_type rt)
{
    policies p;
    p.useDefaults(RELATION);

    p.setFullyReduced();

    forest* in_fully = forest::create(D, RELATION, rt, el, p);

    forest* out_fully = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* in_quasi = forest::create(D, RELATION, rt, el, p);

    forest* out_quasi = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setIdentityReduced();

    forest* in_ident = forest::create(D, RELATION, rt, el, p);

    forest* out_ident = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);



    for (unsigned i=MIN_REL_CARD; i<=MAX_REL_CARD; i*=MULT_REL_CARD) {
        std::cout << "Testing " << OPS << " on "
                  << nameOf(el) << " " << nameOf(rt) << " matrices; "
                  << i << " / " << RG.potential() << " nonzeroes\n";

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
    std::cerr << "    --int:    Test int range/edge values (default)\n";
    std::cerr << "    --long:   Test long range/edge values\n";
    std::cerr << "    --float:  Test float range/edge values\n";
    std::cerr << "    --double: Test double range/edge values\n";
    std::cerr << "\n";
    std::cerr << "    --MT:     Multi-terminal (default)\n";
    std::cerr << "    --EVp:    EV+\n";
    std::cerr << "    --EVs:    EV*\n";
    std::cerr << "\n";

    exit(1);
}

int main(int argc, const char** argv)
{
    /*
     * Command-line options
     */
    bool sets = true;
    edge_labeling EL = edge_labeling::MULTI_TERMINAL;
    range_type    RT = range_type::INTEGER;
    char type = 'I';
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

        if (0==strcmp("--MT", argv[i])) {
            EL = edge_labeling::MULTI_TERMINAL;
            continue;
        }
        if (0==strcmp("--EVp", argv[i])) {
            EL = edge_labeling::EVPLUS;
            continue;
        }
        if (0==strcmp("--EVs", argv[i])) {
            EL = edge_labeling::EVTIMES;
            continue;
        }

        if (0==strcmp("--int", argv[i])) {
            RT = range_type::INTEGER;
            type = 'I';
            continue;
        }
        if (0==strcmp("--long", argv[i])) {
            RT = range_type::INTEGER;
            type = 'L';
            continue;
        }
        if (0==strcmp("--float", argv[i])) {
            RT = range_type::REAL;
            type = 'F';
            continue;
        }
        if (0==strcmp("--double", argv[i])) {
            RT = range_type::REAL;
            type = 'D';
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

        domain* D = sets ? SG.makeDomain() : RG.makeDomain();
        if (sets) {
            switch (type) {
                case 'I':   test_sets<int>(D, EL, RT);
                            break;

                case 'L':   test_sets<long>(D, EL, RT);
                            break;

                case 'F':   test_sets<float>(D, EL, RT);
                            break;

                case 'D':   test_sets<double>(D, EL, RT);
                            break;

                default:    throw "Unknown type";
            }
        } else {
            switch (type) {
                case 'I':   test_rels<int>(D, EL, RT);
                            break;

                case 'L':   test_rels<long>(D, EL, RT);
                            break;

                case 'F':   test_rels<float>(D, EL, RT);
                            break;

                case 'D':   test_rels<double>(D, EL, RT);
                            break;

                default:    throw "Unknown type";
            }
        }
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
