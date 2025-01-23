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

vectorgen SG(MEDDLY::SET, 8, 4);
vectorgen RG(MEDDLY::RELATION, 5, 3);

const unsigned SET_COLL_SIZE = 64;
const unsigned REL_COLL_SIZE = 32;

using namespace MEDDLY;

/*
 *
 * Comparisons
 *
 */

bool operator>=(const rangeval &v1, const rangeval &v2)
{
    if (v1.getType() != v2.getType()) {
        throw "type mismatch on our >=";
    }
    switch (v1.getType()) {
        case range_type::BOOLEAN:
            return bool(v1) >= bool(v2);

        case range_type::INTEGER:
            if (v1.isPlusInfinity()) return true;
            if (v2.isPlusInfinity()) return false;
            return long(v1) >= long(v2);

        case range_type::REAL:
            return double(v1) >= double(v2);

        default:
            throw "unhandled type in our >=";
    }
}

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

template <bool MAX>
void eval_set(const rangeval &deflt, const minterm_coll &MC,
        const minterm &vars, rangeval &ans)
{
    ans = deflt;
    for (unsigned i=0; i<MC.size(); i++) {
        if (matches_set(MC.at(i), vars)) {
            const rangeval& mcv = MC.at(i).getValue();
            if (MAX) {
                if (mcv >= ans) {
                    ans = mcv;
                }
            } else {
                if (ans >= mcv) {
                    ans = mcv;
                }
            }
        } // matches_set
    } // for i
}


template <bool MAX>
void check_equal_set(const dd_edge &E,
        const rangeval &deflt, const minterm_coll &MC)
{
    //
    // Make sure E encodes the max of minterms in MC.
    //
    // Do this by evaluating the function for all variable assignments.
    //
    minterm vars(E.getForest());
    rangeval Eval, Fval;
    vars.setAllVars(0);
    do {
        E.evaluate(vars, Eval);

        eval_set<MAX>(deflt, MC, vars, Fval);

        if (Fval != Eval) {
            MEDDLY::ostream_output out(std::cout);
            out << "\nMismatch!";
            out << "\n    MDD says   ";
            Eval.write(out);
            out << "\n    Collection (using "
                << (MAX ? "max" : "min") << ") says ";
            Fval.write(out);
            out << "\n    Var assignments: ";
            vars.show(out);
            out << "\n    Collection:\n";
            MC.show(out, "    ");
            out << "deflt ";
            deflt.write(out);
            out << "\n    Func. MDD encoding:\n";
            E.showGraph(out);
            out.flush();
            throw "mismatch";
        }

    } while (SG.nextMinterm(vars));
}

template <bool MAX>
void test_sets(forest* F, rangeval deflt, std::vector <rangeval> fvals)
{
    //
    // Display what we're doing
    //
    ostream_output out(std::cout);

    out << "Checking " << shortNameOf(F->getRangeType())
        << " " << shortNameOf(F->getReductionRule())
        << " " << nameOf(F->getEdgeLabeling()) << " MDD with default ";
    deflt.write(out);
    out << "\n";
    out.flush();

    //
    // Set up collection, edge
    //
    minterm_coll MC(SET_COLL_SIZE, F);
    dd_edge E(F);

    //
    // Check all don't cares
    //

    out << "    x";
    out.flush();

    MC.unused().setAllVars(DONT_CARE);
    MC.unused().setValue(fvals[0]);
    MC.pushUnused();

#ifdef SHOW_MINTERMS
    out << "\n  minterms:\n";
    MC.show(out, "  ");
    out << "\n";
    out.flush();
#endif
    if (MAX) {
        MC.buildFunctionMax(deflt, E);
    } else {
        MC.buildFunctionMin(deflt, E);
    }
    out << ".";
    out.flush();
    check_equal_set<MAX>(E, deflt, MC);

    //
    // Check collections of increasing size
    //
    MC.clear();
    for (unsigned i=1; i<=SET_COLL_SIZE; i*=2)
    {
        out << i;
        out.flush();

        // more random minterms
        while (MC.size() < i) {
            SG.randomizeMinterm(MC.unused(), range_type::BOOLEAN);
            MC.unused().setValue( fvals[ MC.size() % fvals.size() ] );
            MC.pushUnused();
        }
#ifdef SHOW_MINTERMS
        out << "\n  minterms:\n";
        MC.show(out, "  ");
        out << "\n";
        out.flush();
#endif
        if (MAX) {
            MC.buildFunctionMax(deflt, E);
        } else {
            MC.buildFunctionMin(deflt, E);
        }
        out << ".";
        out.flush();
        check_equal_set<MAX>(E, deflt, MC);
    }
    out << "\n";
}

void test_sets_with_policies(domain *D, policies p)
{
    //
    // Function values to use
    //
    std::vector <rangeval> allT;
    allT.push_back(true);

    std::vector <rangeval> allF;
    allF.push_back(false);

    std::vector <rangeval> intvals;
    intvals.push_back(2);
    intvals.push_back(3);
    intvals.push_back(5);

    std::vector <rangeval> realvals;
    realvals.push_back(2.5);
    realvals.push_back(3.5);
    realvals.push_back(5.5);

    //
    // Boolean
    //
    forest* Fbool = forest::create(D, SET, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL, p);

    test_sets<true >(Fbool, false, allT);
    test_sets<false>(Fbool, true,  allF);

    forest::destroy(Fbool);

    //
    // Integer
    //
    forest* Fint  = forest::create(D, SET, range_type::INTEGER,
            edge_labeling::MULTI_TERMINAL, p);

    test_sets<true >(Fint,  0L, intvals);
    test_sets<true >(Fint,  1L, intvals);
    test_sets<false>(Fint,  7L, intvals);

    forest::destroy(Fint);

    //
    // Real
    //
    forest* Freal = forest::create(D, SET, range_type::REAL,
            edge_labeling::MULTI_TERMINAL, p);

    test_sets<true >(Freal, 0.0, realvals);
    test_sets<true >(Freal, -1.25, realvals);
    test_sets<false>(Freal, 8.75, realvals);

    //
    // EV+ Integer
    //

    forest* Fevp  = forest::create(D, SET, range_type::INTEGER,
            edge_labeling::EVPLUS, p);

    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);

    test_sets<true >(Fevp,  0L, intvals);
    test_sets<true >(Fevp,  1L, intvals);
    test_sets<false>(Fevp,  7L, intvals);
    test_sets<false>(Fevp,  infty, intvals);

    forest::destroy(Fevp);
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

template <bool MAX>
void eval_rel(const rangeval &deflt, const minterm_coll &MC,
        const minterm &vars, rangeval &ans)
{
    ans = deflt;
    for (unsigned i=0; i<MC.size(); i++) {
        if (matches_rel(MC.at(i), vars)) {
            const rangeval& mcv = MC.at(i).getValue();
            if (MAX) {
                if (mcv >= ans) {
                    ans = mcv;
                }
            } else {
                if (ans >= mcv) {
                    ans = mcv;
                }
            }
        } // matches_rel
    } // for i
}

template <bool MAX>
void check_equal_rel(const dd_edge &E,
        const rangeval &deflt, const minterm_coll &MC)
{
    //
    // Make sure E encodes the max of minterms in MC.
    //
    // Do this by evaluating the function for all variable assignments.
    //
    minterm vars(E.getForest());
    rangeval Eval, Fval;
    vars.setAllVars(0, 0);
    do {
        E.evaluate(vars, Eval);

        eval_rel<MAX>(deflt, MC, vars, Fval);

        if (Fval != Eval) {
            MEDDLY::ostream_output out(std::cout);
            out << "\nMismatch!";
            out << "\n    MxD says   ";
            Eval.write(out);
            out << "\n    Collection (using "
                << (MAX ? "max" : "min") << ") says ";
            Fval.write(out);
            out << "\n    Var assignments: ";
            vars.show(out);
            out << "\n    Collection:\n";
            MC.show(out, "    ");
            out << "deflt ";
            deflt.write(out);
            out << "\n    Func. MxD encoding:\n";
            E.showGraph(out);
            out.flush();
            throw "mismatch";
        }

    } while (RG.nextMinterm(vars));
}

template <bool MAX>
void test_rels(forest* F, rangeval deflt, std::vector <rangeval> fvals)
{
    //
    // Display what we're doing
    //
    ostream_output out(std::cout);

    out << "Checking " << shortNameOf(F->getRangeType())
        << " " << shortNameOf(F->getReductionRule())
        << " " << nameOf(F->getEdgeLabeling()) << " MxD with default ";
    deflt.write(out);
    out << "\n";
    out.flush();

    //
    // Set up collection, edge
    //
    minterm_coll MC(REL_COLL_SIZE, F);
    dd_edge E(F);

    //
    // Check all don't cares
    //
    out << "    x";
    out.flush();

    MC.unused().setAllVars(DONT_CARE, DONT_CARE);
    MC.unused().setValue(fvals[0]);
    MC.pushUnused();

#ifdef SHOW_MINTERMS
    out << "\n  minterms:\n";
    MC.show(out, "  ");
    out << "\n";
    out.flush();
#endif
    if (MAX) {
        MC.buildFunctionMax(deflt, E);
    } else {
        MC.buildFunctionMin(deflt, E);
    }
    out << ".";
    out.flush();
    check_equal_rel<MAX>(E, deflt, MC);

    //
    // Check all don't change
    //
    out << "i";
    out.flush();

    MC.clear();
    MC.unused().setAllVars(DONT_CARE, DONT_CHANGE);
    MC.unused().setValue(fvals[0]);
    MC.pushUnused();

#ifdef SHOW_MINTERMS
    out << "\n  minterms:\n";
    MC.show(out, "  ");
    out << "\n";
    out.flush();
#endif
    if (MAX) {
        MC.buildFunctionMax(deflt, E);
    } else {
        MC.buildFunctionMin(deflt, E);
    }
    out << ".";
    out.flush();
    check_equal_rel<MAX>(E, deflt, MC);

    //
    // Check collections of increasing size
    //
    MC.clear();
    for (unsigned i=1; i<=REL_COLL_SIZE; i*=2)
    {
        out << i;
        out.flush();

        // more random minterms
        while (MC.size() < i) {
            RG.randomizeMinterm(MC.unused(), range_type::BOOLEAN);
            MC.unused().setValue( fvals[ MC.size() % fvals.size() ] );
            MC.pushUnused();
        }
#ifdef SHOW_MINTERMS
        out << "\n  minterms:\n";
        MC.show(out, "  ");
        out << "\n";
        out.flush();
#endif
        if (MAX) {
            MC.buildFunctionMax(deflt, E);
        } else {
            MC.buildFunctionMin(deflt, E);
        }
        out << ".";
        out.flush();
        check_equal_rel<MAX>(E, deflt, MC);
    }
    out << "\n";
}

void test_rels_with_policies(domain *D, policies p)
{
    //
    // Function values to use
    //
    std::vector <rangeval> allT;
    allT.push_back(true);

    std::vector <rangeval> allF;
    allF.push_back(false);

    std::vector <rangeval> intvals;
    intvals.push_back(2);
    intvals.push_back(3);
    intvals.push_back(5);

    std::vector <rangeval> realvals;
    realvals.push_back(2.5);
    realvals.push_back(3.5);
    realvals.push_back(5.5);

    //
    // Boolean
    //
    forest* Fbool = forest::create(D, RELATION, range_type::BOOLEAN,
            edge_labeling::MULTI_TERMINAL, p);

    test_rels<true >(Fbool, false, allT);
    test_rels<false>(Fbool, true,  allF);

    forest::destroy(Fbool);

    //
    // Integer
    //
    forest* Fint  = forest::create(D, RELATION, range_type::INTEGER,
            edge_labeling::MULTI_TERMINAL, p);

    test_rels<true >(Fint,  0L, intvals);
    test_rels<true >(Fint,  1L, intvals);
    test_rels<false>(Fint,  7L, intvals);

    forest::destroy(Fint);

    //
    // Real
    //
    forest* Freal = forest::create(D, RELATION, range_type::REAL,
            edge_labeling::MULTI_TERMINAL, p);

    test_rels<true >(Freal, 0.0, realvals);
    test_rels<true >(Freal, -1.25, realvals);
    test_rels<false>(Freal, 8.75, realvals);

    //
    // EV+ Integer
    //

    forest* Fevp  = forest::create(D, RELATION, range_type::INTEGER,
            edge_labeling::EVPLUS, p);

    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);

    test_rels<true >(Fevp,  0L, intvals);
    test_rels<true >(Fevp,  1L, intvals);
    test_rels<false>(Fevp,  7L, intvals);
    test_rels<false>(Fevp,  infty, intvals);

    forest::destroy(Fevp);
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

