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
    Tests element-wise user-defined unary operations.
*/


#include <cstdlib>
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

using namespace MEDDLY;

// #define DEBUG_MXDOPS

void MyAbs(const rangeval &x, rangeval &y)
{
    if (x.isPlusInfinity()) {
        y = x;
        return;
    }
    if (x.isInteger()) {
        const long L = x;
        y = (L<0) ? -L : L;
    } else {
        const double D = x;
        y = (D<0) ? -D : D;
    }
}

user_unary_factory MyAbsFact("MyAbs", MyAbs);

void MyNeg(const rangeval &x, rangeval &y)
{
    if (x.isPlusInfinity()) {
        y = x;
        return;
    }
    if (x.isInteger()) {
        y = -long(x);
    } else {
        y = -double(x);
    }
}

user_unary_factory MyNegFact("MyNeg", MyNeg);

void MyEven(const rangeval &x, rangeval &y)
{
    if (x.isPlusInfinity()) {
        y = false;
        return;
    }
    if (x.isInteger()) {
        y = ( long(x) % 2 == 0 );
    } else {
        double hx = double(x) / 2.0;
        y = ( double(long(hx)) == hx );
    }
}

user_unary_factory MyEvenFact("MyEven", MyEven);



void compare(vectorgen &Gen, const std::vector <MEDDLY::rangeval> &A,
        rangeval deflt, user_unary_factory &UUF, forest* fin, forest* fout)
{
    ostream_output out(std::cout);

    const unsigned POTENTIAL = Gen.potential();
    std::vector <rangeval> FA(POTENTIAL);
    if (FA.size() != A.size()) {
        throw "size mismatch";
    }
    dd_edge Add(fin), FAdd(fout), FAsym(fout);

    Gen.explicit2edgeMin(A, Add, deflt);

    // Build output vector explicitly
    user_defined_unary f = UUF.getFunc();
    for (unsigned i=0; i<POTENTIAL; i++) {
        f(A[i], FA[i]);
    }

    // Output vector -> MDD
    rangeval fdeflt;
    f(deflt, fdeflt);
    Gen.explicit2edgeMin(FA, FAdd, fdeflt);

    // Build output vector "symbolically"
    apply(UUF, Add, FAsym);

    // compare outputs
    if (FAsym == FAdd) return;

    out << "\nMismatch on " << UUF.getName() << "\n";

    out << "Input (" << shortNameOf(fin->getReductionRule()) << "):\n";
    Add.showGraph(out);
    out << "Expected Output (" << shortNameOf(fout->getReductionRule()) << "):\n";
    FAdd.showGraph(out);
    out << "Obtained Output (" << shortNameOf(fout->getReductionRule()) << "):\n";
    FAsym.showGraph(out);

    out << "Input vector [";
    for (unsigned i=0; i<A.size(); i++) {
        if (i) out << ", ";
        A[i].write(out);
    }
    out << "\n";


    throw "mismatch";


}

template <typename TYPE>
void test_on_functions(unsigned scard, forest* fin, forest* fout, forest* fbool)
{
    if (!fin) throw "null fin";
    if (!fout) throw "null fout";

    vectorgen &Gen = fin->isForRelations() ? RG : SG;

    std::cerr << "    " << shortNameOf(fin->getReductionRule())
              << " : " << shortNameOf(fout->getReductionRule()) << ' ';

    std::vector <MEDDLY::rangeval> Aset(Gen.potential());

    std::vector <MEDDLY::rangeval> values(7);
    if (fin->isEVPlus()) {
        values[0] = MEDDLY::rangeval(MEDDLY::range_special::PLUS_INFINITY,
                            MEDDLY::range_type::INTEGER);
    } else if (fin->isEVTimes()) {
        values[0] =  TYPE(1);
    } else {
        values[0] =  TYPE(0);
    }
    values[1] =  TYPE(4);
    values[2] =  TYPE(3);
    values[3] =  TYPE(2);
    values[4] =  TYPE(-2);
    values[5] =  TYPE(-3);
    values[6] =  TYPE(-4);

    for (unsigned i=0; i<10; i++) {
        std::cerr << '.';
        Gen.randomizeVector(Aset, scard, values);

        compare(Gen, Aset, values[0], MyAbsFact, fin, fout);
        compare(Gen, Aset, values[0], MyNegFact, fin, fout);
        compare(Gen, Aset, values[0], MyEvenFact, fin, fbool);
    }
    for (unsigned i=0; i<10; i++) {
        std::cerr << "x";
        Gen.randomizeFully(Aset, scard, values);

        compare(Gen, Aset, values[0], MyAbsFact, fin, fout);
        compare(Gen, Aset, values[0], MyNegFact, fin, fout);
        compare(Gen, Aset, values[0], MyEvenFact, fin, fbool);
    }
    if (fout->isForRelations()) {
        for (unsigned i=0; i<10; i++) {
            std::cerr << "i";
            Gen.randomizeIdentity(Aset, scard, values);

            compare(Gen, Aset, values[0], MyAbsFact, fin, fout);
            compare(Gen, Aset, values[0], MyNegFact, fin, fout);
            compare(Gen, Aset, values[0], MyEvenFact, fin, fbool);
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
    forest* fbool = forest::create(D, SET, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* quasi = forest::create(D, SET, rt, el, p);
    forest* qbool = forest::create(D, SET, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);


    for (unsigned i=MIN_SET_CARD; i<=MAX_SET_CARD; i*=MULT_SET_CARD) {
        std::cout << "Testing on "
                  << nameOf(el) << " " << nameOf(rt) << " vectors; "
                  << i << " / " << SG.potential() << " nonzeroes\n";

        test_on_functions<TYPE>(i, fully, fully, fbool);
        test_on_functions<TYPE>(i, fully, quasi, qbool);
        test_on_functions<TYPE>(i, quasi, fully, fbool);
        test_on_functions<TYPE>(i, quasi, quasi, qbool);
    }
}

template <typename TYPE>
void test_rels(domain* D, edge_labeling el, range_type rt)
{
    policies p;
    p.useDefaults(RELATION);

    p.setFullyReduced();

    forest* fully = forest::create(D, RELATION, rt, el, p);
    forest* fbool = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* quasi = forest::create(D, RELATION, rt, el, p);
    forest* qbool = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);

    p.setIdentityReduced();

    forest* ident = forest::create(D, RELATION, rt, el, p);
    forest* ibool = forest::create(D, RELATION, range_type::BOOLEAN,
                        edge_labeling::MULTI_TERMINAL, p);


    for (unsigned i=MIN_REL_CARD; i<=MAX_REL_CARD; i*=MULT_REL_CARD) {
        std::cout << "Testing on "
                  << nameOf(el) << " " << nameOf(rt) << " matrices; "
                  << i << " / " << RG.potential() << " nonzeroes\n";

        test_on_functions<TYPE>(i, fully, fully, fbool);
        test_on_functions<TYPE>(i, fully, quasi, qbool);
        test_on_functions<TYPE>(i, fully, ident, ibool);

    //

        test_on_functions<TYPE>(i, ident, fully, fbool);
        test_on_functions<TYPE>(i, ident, quasi, qbool);
        test_on_functions<TYPE>(i, ident, ident, ibool);

    //

        test_on_functions<TYPE>(i, quasi, fully, fbool);
        test_on_functions<TYPE>(i, quasi, quasi, qbool);
        test_on_functions<TYPE>(i, quasi, ident, ibool);
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
