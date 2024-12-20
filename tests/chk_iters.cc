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
#include <time.h>
#include <iostream>

#define TEST_SETS
#define TEST_RELS
// #define SHOW_MINTERMS

const int DOMSIZE = 4;       // DO NOT change
const int SETVARS = 10;
const int RELVARS = 6;

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

using namespace MEDDLY;

void randomizeMinterm(minterm &m, range_type rt)
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

    static int count=0;
    count = (1+count%4);
    switch (rt) {
        case range_type::INTEGER:
                m.setTerm(count);
                break;

        case range_type::REAL:
                m.setTerm(count * 1.1f);
                break;

        default:
                m.setTerm(true);
    }
}


/*
 *
 * Tests for set-type forests
 *
 */

void test_sets(char reduction, range_type rt, edge_labeling el)
{
    //
    // Build domain - once
    //
    int bs[SETVARS];
    for (unsigned i=0; i<SETVARS; i++) {
        bs[i] = DOMSIZE;
    }
    domain* D = domain::createBottomUp(bs, SETVARS);

    //
    // Build the forest
    //
    policies p;
    p.useDefaults(SET);

    if ('f' == reduction) {
        p.setFullyReduced();
    } else {
        p.setQuasiReduced();
    }

    forest* F = forest::create(D, SET, rt, el, p);

    ostream_output out(std::cout);
    out << "Checking iterators over ";
    out << nameOf(F->getReductionRule()) << " "
        << nameOf(el) << "MDDs of " << nameOf(rt) << "s\n";
    out.flush();


    //
    // Set up minterm collection
    // and evaluation minterm
    //
#ifdef SHOW_MINTERMS
    minterm_coll mtcoll(4, F);
#else
    minterm_coll mtcoll(16, F);
    minterm_coll mtctwo(64, F);
#endif
    minterm eval(F);

    //
    // Build the collection
    //

    while (mtcoll.size() < mtcoll.maxsize()) {
        randomizeMinterm(mtcoll.unused(), rt);
        mtcoll.pushUnused();
    }

    //
    // Set up ddedges
    //
    dd_edge E(F);
    mtcoll.buildFunction(E);

#ifdef SHOW_MINTERMS
    out << "\nMinterms:\n";
    mtcoll.show(out);
    out << "Iterating over ";
    minterm mask(F);
    mask.setAllVars(DONT_CARE);
    mask.setVar(1, 2);
    mask.show(out);
    out << "\n";
    unsigned long count = 0;
    for (dd_edge::iterator i = E.begin(&mask);
            i != E.end();
            i++)
    {
        out.put(count, 5);
        out << ": ";
        (*i).show(out);
        out << "\n";
        out.flush();
        ++count;
    }

    out << "Done, " << count << " minterms\n";
#else
    dd_edge E2(F), tmp(F);
    // F->createEdge(false, E2);

    out << "Iterating...\n";
    for (dd_edge::iterator i = E.begin(); i; ++i)
    {
        mtctwo.unused().setFrom(*i);
        mtctwo.pushUnused();
        if (mtctwo.isFull()) {
            mtctwo.buildFunction(tmp);
            E2 += tmp;
            mtctwo.clear();
        }
    }
    mtctwo.buildFunction(tmp);
    E2 += tmp;

    if (E != E2) {
        throw "Function mismatch";
    }
    out << "    matches\n";
#endif

    domain::destroy(D);
}



void test_rels(char reduction, range_type rt, edge_labeling el)
{
    //
    // Build domain - once
    //
    int bs[RELVARS];
    for (unsigned i=0; i<RELVARS; i++) {
        bs[i] = DOMSIZE;
    }
    domain* D = domain::createBottomUp(bs, RELVARS);

    //
    // Build the forest
    //
    policies p;
    p.useDefaults(RELATION);

    if ('f' == reduction) {
        p.setFullyReduced();
    } else if ('i' == reduction) {
        p.setIdentityReduced();
    } else {
        p.setQuasiReduced();
    }

    forest* F = forest::create(D, RELATION, rt, el, p);

    ostream_output out(std::cout);
    out << "Checking iterators over ";
    out << nameOf(F->getReductionRule()) << " "
        << nameOf(el) << "MxDs of " << nameOf(rt) << "s\n";
    out.flush();


    //
    // Set up minterm collection
    // and evaluation minterm
    //
#ifdef SHOW_MINTERMS
    minterm_coll mtcoll(4, F);
#else
    minterm_coll mtcoll(16, F);
    minterm_coll mtctwo(64, F);
#endif
    minterm eval(F);

    //
    // Build the collection
    //

    while (mtcoll.size() < mtcoll.maxsize()) {
        randomizeMinterm(mtcoll.unused(), rt);
        mtcoll.pushUnused();
    }

    //
    // Set up ddedges
    //
    dd_edge E(F);
    mtcoll.buildFunction(E);

#ifdef SHOW_MINTERMS
    out << "\nMinterms:\n";
    mtcoll.show(out);
    out << "Iterating over ";
    minterm mask(F);
    mask.setAllVars(DONT_CARE, DONT_CARE);
    // mask.setVars(1, 2, 1);
    mask.show(out);
    out << "\n";
    unsigned long count = 0;
    for (dd_edge::iterator i = E.begin(&mask);
            i != E.end();
            i++)
    {
        out.put(count, 5);
        out << ": ";
        (*i).show(out);
        out << "\n";
        out.flush();
        ++count;
    }

    out << "Done, " << count << " minterms\n";
#else
    dd_edge E2(F), tmp(F);
    // F->createEdge(false, E2);

    out << "Iterating...\n";
    for (dd_edge::iterator i = E.begin();
            i != E.end();
            i++)
    {
        mtctwo.unused().setFrom(*i);
        mtctwo.pushUnused();
        if (mtctwo.isFull()) {
            mtctwo.buildFunction(tmp);
            E2 += tmp;
            mtctwo.clear();
        }
    }
    mtctwo.buildFunction(tmp);
    E2 += tmp;

    if (E != E2) {
        throw "Function mismatch";
    }
    out << "    matches\n";
#endif

    domain::destroy(D);
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
#ifdef TEST_SETS
        test_sets('q', range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
        test_sets('f', range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
        test_sets('q', range_type::INTEGER, edge_labeling::MULTI_TERMINAL);
        test_sets('f', range_type::INTEGER, edge_labeling::MULTI_TERMINAL);
        test_sets('q', range_type::REAL, edge_labeling::MULTI_TERMINAL);
        test_sets('f', range_type::REAL, edge_labeling::MULTI_TERMINAL);
#endif
#ifdef TEST_RELS
        test_rels('q', range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
        test_rels('f', range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
        test_rels('i', range_type::BOOLEAN, edge_labeling::MULTI_TERMINAL);
        test_rels('q', range_type::INTEGER, edge_labeling::MULTI_TERMINAL);
        test_rels('f', range_type::INTEGER, edge_labeling::MULTI_TERMINAL);
        test_rels('i', range_type::INTEGER, edge_labeling::MULTI_TERMINAL);
        test_rels('q', range_type::REAL, edge_labeling::MULTI_TERMINAL);
        test_rels('f', range_type::REAL, edge_labeling::MULTI_TERMINAL);
        test_rels('i', range_type::REAL, edge_labeling::MULTI_TERMINAL);
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
