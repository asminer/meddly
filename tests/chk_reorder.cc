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

#include "../src/meddly.h"
#include <iostream>
#include <iomanip>

#define TEST_SETS
#define TEST_RELS

const int DOMSIZE = 4;       // DO NOT change
const int VARS = 10;

const int NICE_ORDER[] = { 0, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4 };

/*
 *
 * RNG stuff
 *
 */

long seed;

double Random()
{
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

int Equilikely(int a, int b)
{
    return (a + (int) ((b - a + 1) * Random()));
}

unsigned Equilikely(unsigned a, unsigned b)
{
    return (a + (unsigned) ((b - a + 1) * Random()));
}

void shuffle(int* A, unsigned n)
{
    for (unsigned i=1; i<n; i++) {
        unsigned j = Equilikely(i, n);
        if (i != j) {
            MEDDLY::SWAP(A[i], A[j]);
        }
    }
}

void randomizeMinterm(MEDDLY::minterm &m, MEDDLY::range_type r)
{
    const int vals[9] = { -1, 0, 0, 1, 1, 2, 2, 3, 3 };
    const int unvals[52] = {
        -1, -1, -1, -1,     // 4 (x,x) pairs
        -1, -1, -1, -1,     // 4 (x,i) pairs
        -1, -1, -1, -1,     // 4 (x, normal) pairs
        0,  1,  2,  3,      // 4 (normal, x) pairs
        0,  1,  2,  3,      // 4 (normal, i) pairs
        0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3,  // 16 normal pairs
        0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3   // 16 normal pairs
    };

    const int prvals[52] = {
        -1, -1, -1, -1,     // 4 (x,x) pairs
        -2, -2, -2, -2,     // 4 (x,i) pairs
        0,  1,  2,  3,      // 4 (x, normal) pairs
        -1, -1, -1, -1,     // 4 (normal, x) pairs
        -2, -2, -2, -2,     // 4 (normal, i) pairs
        0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  // 16 normal pairs
        0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3   // 16 normal pairs
    };


    if (m.isForRelations()) {
        for (unsigned i=1; i<=m.getNumVars(); i++) {
            int index = Equilikely(0, 51);
            m.setVars(i, unvals[index], prvals[index]);
        }
    } else {
        for (unsigned i=1; i<=m.getNumVars(); i++) {
            int index = Equilikely(0, 8);
            m.setVar(i, vals[index]);
        }
    }

    static int val = 0;

    val = 1 + (val % 5);

    if (MEDDLY::range_type::INTEGER == r) {
        m.setValue(val);
    }
    if (MEDDLY::range_type::REAL == r) {
        m.setValue( float(val) );
    }
}

void reorderMinterm(const MEDDLY::minterm &in, MEDDLY::minterm &out,
        const int* reord)
{
    if (in.isForRelations()) {
        for (unsigned i=in.getNumVars(); i; i--) {
            out.setVars(i, in.from(reord[i]), in.to(reord[i]));
        }
    } else {
        for (unsigned i=in.getNumVars(); i; i--) {
            out.setVar(i, in.from(reord[i]));
        }
    }
    out.setValue( in.getValue() );
}

using namespace MEDDLY;

inline const char* name(range_type rt)
{
    switch (rt) {
        case range_type::BOOLEAN:
            return "boolean";

        case range_type::INTEGER:
            return "integer";

        case range_type::REAL:
            return "real";

        default:
            throw "unknown range type";
    }
}

inline const char* name(edge_labeling el)
{
    switch (el) {
        case edge_labeling::MULTI_TERMINAL:
            return "MT";

        case edge_labeling::EVPLUS:
            return "EV+";

        case edge_labeling::EVTIMES:
            return "EV*";

        default:
            throw "unknown edge labeling";
    }
}

/*
 *
 * Tests for set-type forests
 *
 */

void compare(domain* D, range_type rt, edge_labeling el, const policies &p,
        minterm_coll &mtc, const int* order)
{

    forest* F = forest::create(D, mtc.isForRelations(), rt, el, p);
    rangeval zero;
    F->getValueForEdge(F->getTransparentEdge(), F->getTransparentNode(), zero);

    dd_edge e1(F);
    mtc.buildFunctionMax(zero, e1);
    F->reorderVariables(order);

    // Compare against mtc, reordered
    minterm_coll mt2(mtc.size(), F);
    for (unsigned i=0; i<mtc.size(); i++) {
        reorderMinterm(mtc.at(i), mt2.unused(), order);
        mt2.pushUnused();
    }

    dd_edge e2(F);
    mt2.buildFunctionMax(zero, e2);

    if (e1 == e2) {
        forest::destroy(F);
        return;
    }

    ostream_output out(std::cout);

    out << "After reorder:\n";
    e1.showGraph(out);

    out << "According to us:\n";
    e2.showGraph(out);

    throw "mismatch";
}

void test_sets(domain *D, range_type r, edge_labeling el)
{
    //
    // Policies
    //
    policies pf;
    pf.useDefaults(SET);
    pf.setPessimistic();
    pf.setFullyReduced();

    policies pq;
    pq.useDefaults(SET);
    pq.setPessimistic();
    pq.setQuasiReduced();

    //
    // Random variable order
    //
    int randord[VARS+1];
    for (unsigned i=VARS; i; --i)
    {
        randord[i] = i;
    }
    randord[0] = 0;
    shuffle(randord, VARS);

    //
    // Minterms
    //
    minterm_coll mtlist(128, D, SET);

    std::cout << "Checking reordering for " << name(el) << " "
              << name(r) << " vectors\n";

    for (unsigned s=1; s <= mtlist.maxsize(); s *= 2)
    {
        std::cout << std::setw(5) << s << ": ";
        std::cout.flush();
        while (mtlist.size() < s) {
            randomizeMinterm(mtlist.unused(), r);
            mtlist.pushUnused();
        }

        compare(D, r, el, pq, mtlist, NICE_ORDER);
        std::cout << "qn ";
        std::cout.flush();
        compare(D, r, el, pq, mtlist, randord);
        std::cout << "qr ";
        std::cout.flush();

        compare(D, r, el, pf, mtlist, NICE_ORDER);
        std::cout << "fn ";
        std::cout.flush();
        compare(D, r, el, pf, mtlist, randord);
        std::cout << "fr ";
        std::cout.flush();

        std::cout << std::endl;
    }
}


/*
 *
 * Tests for relation-type forests
 *
 */

void test_rels(domain *D, range_type r, edge_labeling el)
{
    //
    // Policies
    //
    policies pf;
    pf.useDefaults(RELATION);
    pf.setPessimistic();
    pf.setFullyReduced();

    policies pq;
    pq.useDefaults(RELATION);
    pq.setPessimistic();
    pq.setQuasiReduced();

    policies pi;
    pi.useDefaults(RELATION);
    pi.setPessimistic();
    pi.setIdentityReduced();

    //
    // Random variable order
    //
    int randord[VARS+1];
    for (unsigned i=VARS; i; --i)
    {
        randord[i] = i;
    }
    randord[0] = 0;
    shuffle(randord, VARS);

    //
    // Minterms
    //
    minterm_coll mtlist(128, D, RELATION);

    std::cout << "Checking reordering for " << name(el) << " "
              << name(r) << " matrices\n";

    for (unsigned s=1; s <= mtlist.maxsize(); s *= 2)
    {
        std::cout << std::setw(5) << s << ": ";
        std::cout.flush();
        while (mtlist.size() < s) {
            randomizeMinterm(mtlist.unused(), r);
            mtlist.pushUnused();
        }

        compare(D, r, el, pq, mtlist, NICE_ORDER);
        std::cout << "qn ";
        std::cout.flush();
        compare(D, r, el, pq, mtlist, randord);
        std::cout << "qr ";
        std::cout.flush();

        compare(D, r, el, pf, mtlist, NICE_ORDER);
        std::cout << "fn ";
        std::cout.flush();
        compare(D, r, el, pf, mtlist, randord);
        std::cout << "fr ";
        std::cout.flush();

        compare(D, r, el, pi, mtlist, NICE_ORDER);
        std::cout << "in ";
        std::cout.flush();
        compare(D, r, el, pi, mtlist, randord);
        std::cout << "ir ";
        std::cout.flush();

        std::cout << std::endl;
    }
}

/*
 *
 * Main
 *
 */

int main(int argc, const char** argv)
{
    using namespace std;

    //
    // First argument: seed
    //
    if (argv[1]) {
        seed = atol(argv[1]);
    } else {
        seed = time(NULL);
    }
    if (seed < 0) seed *= -1;
    if (0==seed)  seed = 12345;
    cout << "Using rng seed " << seed << "\n";

    try {
        MEDDLY::initialize();
        //
        // Build domain - once
        //
        int bs[VARS];
        for (unsigned i=0; i<VARS; i++) {
            bs[i] = DOMSIZE;
        }
        domain* D = domain::createBottomUp(bs, VARS);

#ifdef TEST_SETS
        test_sets(D, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
        // test_sets(D, range_type::INTEGER, edge_labeling::MULTI_TERMINAL);
        // test_sets(D, range_type::INTEGER, edge_labeling::EVPLUS);
        // test_sets(D, range_type::REAL, edge_labeling::MULTI_TERMINAL);
#endif
#ifdef TEST_RELS
        test_rels(D, range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
        // test_rels(D, range_type::INTEGER, edge_labeling::MULTI_TERMINAL);
        // test_rels(D, range_type::INTEGER, edge_labeling::EVPLUS);
        // test_rels(D, range_type::REAL, edge_labeling::MULTI_TERMINAL);
        // test_rels(D, range_type::REAL, edge_labeling::EVTIMES);
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
