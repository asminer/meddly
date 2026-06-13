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

#define TEST_SETS
#define TEST_RELS
// #define SHOW_MINTERMS

#include "randomize.h"

//
// Minterm generators :)
//

vectorgen SG(MEDDLY::SET, 8, 4);
vectorgen RG(MEDDLY::RELATION, 5, 3);

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
    vars.setAllVars(0);
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
    vars.setAllVars(0, 0);
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
    std::cerr << "    --quasi:  Test quasi reduced\n";
    std::cerr << "    --fully:  Test fully reduced (default).\n";
    std::cerr << "    --ident:  Test identity reduced (relations only)\n";
    std::cerr << "\n";

    exit(1);
}

void setReductionLetter(policies &p, char letter)
{
    switch (letter)
    {
        case 'Q':
            p.setQuasiReduced();
            return;

        case 'F':
            p.setFullyReduced();
            return;

        case 'I':
            p.setIdentityReduced();
            return;

        default:
            throw "Unknown reduction";
    }
}

/*
 *
 * Main
 *
 */

int main(int argc, const char** argv)
{
    bool sets = true;
    char reduction = 'F';
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
        if (0==strcmp("--quasi", argv[i])) {
            reduction = 'Q';
            continue;
        }
        if (0==strcmp("--fully", argv[i])) {
            reduction = 'F';
            continue;
        }
        if (0==strcmp("--ident", argv[i])) {
            reduction = 'I';
            continue;
        }

        if ((argv[i][0] < '0') || (argv[i][0] > '9')) {
            usage(argv[0]);
        }

        seed = atol(argv[i]);
    }

    if (('I' == reduction) && sets) {
        std::cerr << "Cannot use identity with sets.\n";
        usage(argv[0]);
    }

    //
    // Set seed
    //
    vectorgen::setSeed(seed);

    try {
        MEDDLY::initialize();

        if (sets) {
            domain* D = SG.makeDomain();
            policies p;
            p.useDefaults(SET);
            setReductionLetter(p, reduction);
            test_sets_with_policies(D, p);
            domain::destroy(D);
        } else {
            domain* D = RG.makeDomain();
            policies p;
            p.useDefaults(RELATION);
            setReductionLetter(p, reduction);
            test_rels_with_policies(D, p);
            domain::destroy(D);
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

