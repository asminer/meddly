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

vectorgen SG(MEDDLY::SET, 5, 4);
vectorgen RG(MEDDLY::RELATION, 5, 2);

const unsigned MIN_SET_CARD = 1;
const unsigned MAX_SET_CARD = 256;
const unsigned MULT_SET_CARD = 4;


const unsigned MIN_REL_CARD = 1;
const unsigned MAX_REL_CARD = 256;
const unsigned MULT_REL_CARD = 16;

const char* OPS = "max, min";

using namespace MEDDLY;

// #define DEBUG_MXDOPS

template <typename T>
void Max(const std::vector <MEDDLY::rangeval> &A,
        const std::vector <MEDDLY::rangeval> &B,
        std::vector<MEDDLY::rangeval> &C)
{
    if (A.size() != B.size()) {
        throw "Max size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "Max size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = MAX(T(A[i]), T(B[i]));
    }
}

template <typename T>
void Min(const std::vector <MEDDLY::rangeval> &A,
        const std::vector <MEDDLY::rangeval> &B,
        std::vector<MEDDLY::rangeval> &C)
{
    if (A.size() != B.size()) {
        throw "Min size mismatch A,B";
    }
    if (A.size() != C.size()) {
        throw "Min size mismatch A,C";
    }
    for (unsigned i=0; i<C.size(); i++) {
        C[i] = MIN(T(A[i]), T(B[i]));
    }
}


inline const char* getReductionString(const dd_edge &e)
{
    const forest* f = e.getForest();
    return shortNameOf(f->getReductionRule());
}

void checkEqual(const char* what, const dd_edge &in1, const dd_edge &in2,
        const dd_edge &e1, const dd_edge &e2)
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
        const std::vector <MEDDLY::rangeval> &Aset,
        const std::vector <MEDDLY::rangeval> &Bset,
        forest* f1, forest* f2, forest* fres)
{
    const unsigned POTENTIAL = Gen.potential();

    std::vector <MEDDLY::rangeval> AmaxBset(POTENTIAL);
    std::vector <MEDDLY::rangeval> AminBset(POTENTIAL);

    dd_edge Add(f1), Bdd(f2),
            AmaxBdd(fres), AminBdd(fres);

    Max<T>(Aset, Bset, AmaxBset);
    Min<T>(Aset, Bset, AminBset);

    Gen.explicit2edgeMax(Aset, Add, T(0));
    Gen.explicit2edgeMax(Bset, Bdd, T(0));
    Gen.explicit2edgeMax(AmaxBset, AmaxBdd, T(0));
    Gen.explicit2edgeMax(AminBset, AminBdd, T(0));

    dd_edge AmaxBsym(fres), AminBsym(fres);

    apply(MAXIMUM, Add, Bdd, AmaxBsym);
    checkEqual("max", Add, Bdd, AmaxBsym, AmaxBdd);

    apply(MINIMUM, Add, Bdd, AminBsym);
    checkEqual("min", Add, Bdd, AminBsym, AminBdd);
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

    std::vector <MEDDLY::rangeval> Aset(Gen.potential());
    std::vector <MEDDLY::rangeval> Bset(Gen.potential());

    std::vector <MEDDLY::rangeval> values(5);
    values[0] =  TYPE(0);
    values[1] =  TYPE(6);
    values[2] =  TYPE(4);
    values[3] =  TYPE(2);
    values[4] =  TYPE(-2);

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        Gen.randomizeVector(Aset, scard, values);
        Gen.randomizeVector(Bset, scard, values);

        compare<TYPE>(Gen, Aset, Bset, f1, f2, fres);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        Gen.randomizeFully(Aset, scard, values);
        Gen.randomizeFully(Bset, scard, values);

        compare<TYPE>(Gen, Aset, Bset, f1, f2, fres);
    }
    if (fres->isForRelations()) {
        for (unsigned i=0; i<10; i++) {
            std::cerr << "i";
            Gen.randomizeIdentity(Aset, scard, values);
            Gen.randomizeIdentity(Bset, scard, values);

            compare<TYPE>(Gen, Aset, Bset, f1, f2, fres);
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
