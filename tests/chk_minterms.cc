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

#include "randomize.h"

//
// Minterm generators :)
//

vectorgen SG(MEDDLY::SET, 8, 4, 5);
vectorgen RG(MEDDLY::RELATION, 5, 3, 5);

const unsigned NUM_MINTERMS = 9;

using namespace MEDDLY;

/*
 *
 * Testing for sets
 *
 */

bool matches_set(const minterm &m, const minterm& vars)
{
    for (unsigned i=1; i<=vars.getNumVars(); i++) {
        if (m.from(i) == DONT_CARE) continue;
        if (m.from(i) != vars.from(i)) return false;
    }
    return true;
}


void check_equal_set(const dd_edge &E, const minterm &m, const rangeval &deflt)
{
    //
    // Make sure E encodes the function: if vars match m, then
    // the value for m; otherwise the default value.
    //
    // Do this by evaluating the function for all variable assignments.
    //
    minterm vars(E.getForest());
    rangeval Eval;
    do {
        E.evaluate(vars, Eval);

        const rangeval &Fval = matches_set(m, vars) ? m.getValue() : deflt;

        if (Fval != Eval) {
            MEDDLY::ostream_output out(std::cout);
            out << "\nMismatch!";
            out << "\n    MDD says   ";
            Eval.write(out);
            out << "\n    Should be ";
            Fval.write(out);
            out << "\n    Var assignments: ";
            vars.show(out);
            out << "\n    Generated using: ";
            m.show(out);
            out << "deflt ";
            deflt.write(out);
            out << "\n    Func. MDD encoding:\n";
            E.showGraph(out);
            out.flush();
            throw "mismatch";
        }

    } while (SG.nextMinterm(vars));
}

void test_sets(forest* F, rangeval deflt, rangeval fval)
{
    //
    // Display what we're doing
    //
    ostream_output out(std::cout);

    out << "Checking " << shortNameOf(F->getRangeType())
        << " " << shortNameOf(F->getReductionRule())
        << " " << nameOf(F->getEdgeLabeling()) << " MDD with ";
    deflt.write(out);
    out << "off, ";
    fval.write(out);
    out << "on\n    ";
    out.flush();

    //
    // Set up minterm, edge
    //
    minterm m(F);
    dd_edge E(F);

    //
    // Check all don't cares
    //

    out << "x ";
    out.flush();
    m.setAllVars(DONT_CARE);
    m.setValue(fval);
#ifdef SHOW_MINTERMS
    out << "\n  minterm: ";
    m.show(out);
    out << "\n";
    out.flush();
#endif
    m.buildFunction(deflt, E);
    check_equal_set(E, m, deflt);

    //
    // Check a few random minterms
    //
    for (unsigned i=0; i<NUM_MINTERMS; i++)
    {
        out << "r ";
        out.flush();
        SG.randomizeMinterm(m, range_type::BOOLEAN);
        m.setValue(fval);
#ifdef SHOW_MINTERMS
        out << "\n  minterm: ";
        m.show(out);
        out << "\n";
        out.flush();
#endif
        m.buildFunction(deflt, E);
        check_equal_set(E, m, deflt);
    }
    out << "\n";
}

void test_sets_with_policies(domain *D, policies p)
{
    //
    // Boolean
    //
    forest* Fbool = forest::create(D, SET, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL, p);

    test_sets(Fbool, false, true);
    test_sets(Fbool, true, false);

    forest::destroy(Fbool);

    //
    // Integer
    //
    forest* Fint  = forest::create(D, SET, range_type::INTEGER,
            edge_labeling::MULTI_TERMINAL, p);

    test_sets(Fint,  0L, 2L);
    test_sets(Fint, -1L, 1L);
    test_sets(Fint, -2L, 0L);

    //
    // Real
    //
    forest* Freal = forest::create(D, SET, range_type::REAL,
            edge_labeling::MULTI_TERMINAL, p);

    test_sets(Freal, 0.0, 2.5);
    test_sets(Freal, -1.25, 1.25);
    test_sets(Freal, -2.5, 0.0);

    //
    // EV+ Integer
    //
    forest* Fevp  = forest::create(D, SET, range_type::INTEGER,
            edge_labeling::EVPLUS, p);

    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);

    test_sets(Fevp,  0L, 2L);
    test_sets(Fevp, -1L, 1L);
    test_sets(Fevp, -2L, 0L);
    test_sets(Fevp,  infty, 0L);
    test_sets(Fevp,  infty, 8L);
    test_sets(Fevp,  0L, infty);
    test_sets(Fevp,  8L, infty);
}

/*
 *
 * Testing for relations
 *
 */

bool matches_rel(const minterm &m, const minterm& vars)
{
    for (unsigned i=1; i<=vars.getNumVars(); i++) {
        if ((m.from(i) != DONT_CARE) && (m.from(i) != vars.from(i)))  {
            return false;
        }
        if (m.to(i) == DONT_CARE) continue;
        if (m.to(i) == DONT_CHANGE) {
            if (vars.to(i) == vars.from(i)) continue;
        }
        if (m.to(i) != vars.to(i)) return false;
    }
    return true;
}

void check_equal_rel(const dd_edge &E, const minterm &m, const rangeval &deflt)
{
    //
    // Make sure E encodes the function: if vars match m, then
    // the value for m; otherwise the default value.
    //
    // Do this by evaluating the function for all variable assignments.
    //
    minterm vars(E.getForest());
    rangeval Eval;
    do {
        E.evaluate(vars, Eval);

        const rangeval &Fval = matches_rel(m, vars) ? m.getValue() : deflt;

        if (Fval != Eval) {
            MEDDLY::ostream_output out(std::cout);
            out << "\nMismatch!";
            out << "\n    MxD says   ";
            Eval.write(out);
            out << "\n    Should be ";
            Fval.write(out);
            out << "\n    Var assignments: ";
            vars.show(out);
            out << "\n    Generated using: ";
            m.show(out);
            out << "deflt ";
            deflt.write(out);
            out << "\n    Func. MxD encoding:\n";
            E.showGraph(out);
            out.flush();
            throw "mismatch";
        }

    } while (RG.nextMinterm(vars));
}

void test_rels(forest* F, rangeval deflt, rangeval fval)
{
    //
    // Display what we're doing
    //
    ostream_output out(std::cout);

    out << "Checking " << shortNameOf(F->getRangeType())
        << " " << shortNameOf(F->getReductionRule())
        << " " << nameOf(F->getEdgeLabeling()) << " MxD with ";
    deflt.write(out);
    out << "off, ";
    fval.write(out);
    out << "on\n    ";
    out.flush();

    //
    // Set up minterm, edge
    //
    minterm m(F);
    dd_edge E(F);

    //
    // Check all don't cares
    //

    out << "x ";
    out.flush();
    m.setAllVars(DONT_CARE, DONT_CARE);
    m.setValue(fval);
#ifdef SHOW_MINTERMS
    out << "\n  minterm: ";
    m.show(out);
    out << "\n";
    out.flush();
#endif
    m.buildFunction(deflt, E);
    check_equal_rel(E, m, deflt);

    //
    // Check all 'identity'
    //

    out << "i ";
    out.flush();
    m.setAllVars(DONT_CARE, DONT_CHANGE);
    m.setValue(fval);
#ifdef SHOW_MINTERMS
    out << "\n  minterm: ";
    m.show(out);
    out << "\n";
    out.flush();
#endif
    m.buildFunction(deflt, E);
    check_equal_rel(E, m, deflt);

    //
    // Check a few random minterms
    //
    for (unsigned i=0; i<NUM_MINTERMS; i++)
    {
        out << "r ";
        out.flush();
        RG.randomizeMinterm(m, range_type::BOOLEAN);
        m.setValue(fval);
#ifdef SHOW_MINTERMS
        out << "\n  minterm: ";
        m.show(out);
        out << "\n";
        out.flush();
#endif
        m.buildFunction(deflt, E);
        check_equal_rel(E, m, deflt);
    }
    out << "\n";
}

void test_rels_with_policies(domain *D, policies p)
{
    //
    // Boolean
    //
    forest* Fbool = forest::create(D, RELATION, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL, p);

    test_rels(Fbool, false, true);
    test_rels(Fbool, true, false);

    forest::destroy(Fbool);

    //
    // Integer
    //
    forest* Fint  = forest::create(D, RELATION, range_type::INTEGER,
            edge_labeling::MULTI_TERMINAL, p);

    test_rels(Fint,  0L, 2L);
    test_rels(Fint, -1L, 1L);
    test_rels(Fint, -2L, 0L);

    //
    // EV+ Integer
    //
    forest* Fevp  = forest::create(D, RELATION, range_type::INTEGER,
            edge_labeling::EVPLUS, p);

    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);

    test_rels(Fevp,  0L, 2L);
    test_rels(Fevp, -1L, 1L);
    test_rels(Fevp, -2L, 0L);
    test_rels(Fevp,  infty, 0L);
    test_rels(Fevp,  infty, 8L);
    test_rels(Fevp,  0L, infty);
    test_rels(Fevp,  8L, infty);

    //
    // Real
    //
    forest* Freal = forest::create(D, RELATION, range_type::REAL,
            edge_labeling::MULTI_TERMINAL, p);

    test_rels(Freal, 0.0, 2.5);
    test_rels(Freal, -1.25, 1.25);
    test_rels(Freal, -2.5, 0.0);

    //
    // EV* real
    //
    forest* Fevt = forest::create(D, RELATION, range_type::REAL,
            edge_labeling::EVTIMES, p);

    test_rels(Fevt, 0.0, 2.5);
    test_rels(Fevt, 0.0, 5.0);
    test_rels(Fevt, -1.25, 1.25);
    test_rels(Fevt, -2.5, 0.0);
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

    try {
        MEDDLY::initialize();
#ifdef TEST_SETS
        domain* SD = SG.makeDomain();
        policies Sp;
        Sp.useDefaults(SET);

        // Test fully-reduced set forests
        Sp.setFullyReduced();
        test_sets_with_policies(SD, Sp);

        // Test quasi-reduced set forests
        Sp.setQuasiReduced();
        test_sets_with_policies(SD, Sp);

        domain::destroy(SD);
#endif
#ifdef TEST_RELS
        domain* RD = RG.makeDomain();
        policies Rp;
        Rp.useDefaults(RELATION);

        // Test fully-reduced set forests
        Rp.setFullyReduced();
        test_rels_with_policies(RD, Rp);

        // Test identity-reduced set forests
        Rp.setIdentityReduced();
        test_rels_with_policies(RD, Rp);

        // Test quasi-reduced set forests
        Rp.setQuasiReduced();
        test_rels_with_policies(RD, Rp);

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

