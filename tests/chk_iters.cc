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

#include "randomize.h"

//
// Minterm generators :)
//

vectorgen SG(MEDDLY::SET, 12, 4);
vectorgen RG(MEDDLY::RELATION, 10, 4);

#define TEST_SETS
#define TEST_RELS
// #define SHOW_MINTERMS

using namespace MEDDLY;

/*
 *
 * Tests for set-type forests
 *
 */

void test_sets(domain* D, char reduction, range_type rt, edge_labeling el,
        const std::vector <rangeval> &vals)
{
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
    rangeval zero;
    F->getValueForEdge(F->getTransparentEdge(), F->getTransparentNode(), zero);

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
        SG.randomizeMinterm(mtcoll.unused(), rt);
        mtcoll.unused().setValue( vals[mtcoll.size() % vals.size()] );
        mtcoll.pushUnused();
    }

    //
    // Set up ddedges
    //
    dd_edge E(F);
    mtcoll.buildFunctionMax(zero, E);

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
            mtctwo.buildFunctionMax(zero, tmp);
            E2 += tmp;
            mtctwo.clear();
        }
    }
    mtctwo.buildFunctionMax(zero, tmp);
    E2 += tmp;

    if (E != E2) {
        throw "Function mismatch";
    }
    out << "    matches\n";
#endif
}



void test_rels(domain* D, char reduction, range_type rt, edge_labeling el,
        const std::vector <rangeval> &vals)
{
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
    rangeval zero;
    F->getValueForEdge(F->getTransparentEdge(), F->getTransparentNode(), zero);

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
        RG.randomizeMinterm(mtcoll.unused(), rt);
        mtcoll.unused().setValue( vals[mtcoll.size() % vals.size()] );
        mtcoll.pushUnused();
    }

    //
    // Set up ddedges
    //
    dd_edge E(F);
    mtcoll.buildFunctionMax(zero, E);

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
            mtctwo.buildFunctionMax(zero, tmp);
            E2 += tmp;
            mtctwo.clear();
        }
    }
    mtctwo.buildFunctionMax(zero, tmp);
    E2 += tmp;

    if (E != E2) {
        throw "Function mismatch";
    }
    out << "    matches\n";
#endif
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
    long seed = 0;
    if (argv[1]) {
        seed = atol(argv[1]);
    }
    vectorgen::setSeed(seed);

    //
    // Initialize function values
    //
    std::vector<rangeval> boolvals(1);
    boolvals[0] = true;

    std::vector<rangeval> intvals(6);
    intvals[0] =  6;
    intvals[1] =  4;
    intvals[2] =  2;
    intvals[3] = -2;
    intvals[4] = -4;
    intvals[5] = -6;

    std::vector<rangeval> realvals(6);
    realvals[0] =  6.0;
    realvals[1] =  4.0;
    realvals[2] =  2.0;
    realvals[3] = -2.0;
    realvals[4] = -4.0;
    realvals[5] = -6.0;

    try {
        MEDDLY::initialize();
#ifdef TEST_SETS
        domain* SD = SG.makeDomain();
        test_sets(SD, 'q', range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL, boolvals);
        test_sets(SD, 'f', range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL, boolvals);
        test_sets(SD, 'q', range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL, intvals);
        test_sets(SD, 'f', range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL, intvals);
        test_sets(SD, 'q', range_type::REAL,
                edge_labeling::MULTI_TERMINAL, realvals);
        test_sets(SD, 'f', range_type::REAL,
                edge_labeling::MULTI_TERMINAL, realvals);
        domain::destroy(SD);
#endif
#ifdef TEST_RELS
        domain* RD = RG.makeDomain();
        test_rels(RD, 'q', range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL, boolvals);
        test_rels(RD, 'f', range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL, boolvals);
        test_rels(RD, 'i', range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL, boolvals);
        test_rels(RD, 'q', range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL, intvals);
        test_rels(RD, 'f', range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL, intvals);
        test_rels(RD, 'i', range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL, intvals);
        test_rels(RD, 'q', range_type::REAL,
                edge_labeling::MULTI_TERMINAL, realvals);
        test_rels(RD, 'f', range_type::REAL,
                edge_labeling::MULTI_TERMINAL, realvals);
        test_rels(RD, 'i', range_type::REAL,
                edge_labeling::MULTI_TERMINAL, realvals);
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
