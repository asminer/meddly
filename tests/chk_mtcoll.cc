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

// #define TEST_SETS
#define TEST_RELS
// #define SHOW_MINTERMS

#include "randomize.h"

//
// Minterm generators :)
//

vectorgen SG(MEDDLY::SET, 8, 4, 5);
vectorgen RG(MEDDLY::RELATION, 5, 3, 5);

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

    /*
    forest* Fevp  = forest::create(D, SET, range_type::INTEGER,
            edge_labeling::EVPLUS, p);

    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);

    test_sets<true >(Fevp,  0L, intvals);
    test_sets<true >(Fevp,  1L, intvals);
    test_sets<false>(Fevp,  7L, intvals);
    test_sets<false>(Fevp,  infty, intvals);

    forest::destroy(Fevp);
    */
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

    // test_rels<true >(Fbool, false, allT);
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

    /*
    forest* Fevp  = forest::create(D, RELATION, range_type::INTEGER,
            edge_labeling::EVPLUS, p);

    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);

    test_rels<true >(Fevp,  0L, intvals);
    test_rels<true >(Fevp,  1L, intvals);
    test_rels<false>(Fevp,  7L, intvals);
    test_rels<false>(Fevp,  infty, intvals);

    forest::destroy(Fevp);
    */
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


//
// OLD BELOW
//

#if 0

/*
 *  common for sets and relations
 */

template <typename RTYPE>
inline void show(MEDDLY::output &out, RTYPE v)
{
    out << v;
}

template <>
inline void show(MEDDLY::output &out, bool v)
{
    out << (v ? "true" : "false");
}

template <typename RTYPE>
inline bool different(RTYPE a, RTYPE b)
{
    return a != b;
}

template <>
inline bool different(double a, double b)
{
    if (a < b) {
        return (b-a) > 1e-5;
    } else {
        return (a-b) > 1e-5;
    }
}


template <typename RTYPE>
void mismatch(const minterm &eval,
        char mddtype, const MEDDLY::dd_edge &E, RTYPE val_mdd,
        const minterm_coll &mtcoll, RTYPE mtcoll_answer)
{
    const char* MDD = eval.isForRelations() ? "MxD" : "MDD";
    MEDDLY::ostream_output out(std::cout);
    out << "\nMismatch on ";
    eval.show(out);
    out << "\n  " << mddtype << MDD << ": ";
    show(out, val_mdd);
    out << "\nmtlist: ";
    show(out, mtcoll_answer);
    out << "\n\n";
    out << "Minterm list:\n";
    mtcoll.show(out);

    out << mddtype << MDD << ":\n";
    E.showGraph(out);
    out.flush();

    throw "mismatch";
}


/*
 *
 * Manipulate minterms for sets
 *
 */




/*
 *
 * Tests for set-type forests
 *
 */

template <typename RTYPE>
void test_sets(domain* D, range_type rt)
{
    //
    // Build the various types of boolean forests
    //
    policies p;
    p.useDefaults(SET);

    p.setFullyReduced();

    forest* Ff = forest::create(D, SET, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* Fq = forest::create(D, SET, rt,
                    edge_labeling::MULTI_TERMINAL, p);


    //
    // Set up minterm collection
    // and evaluation minterm
    //
    minterm_coll mtcoll(64, D, SET);
    minterm eval(D, SET);

    //
    // For various collections of minterms,
    //  (1) build sets
    //  (2) explicitly verify sets against minterms

    ostream_output out(std::cout);

    out << "Checking vectors of " << nameOf(rt)
        << " built from minterm collections:\n";
    out.flush();

    char     testtype[] = { 'f', 'r', 'r', 'r', 'r', 'r', 'r', 'r', 0 };
    unsigned mtsizes[] =  {   0,   1,   2,   4,   8,  16,  32,  64    };

    for (unsigned i=0; testtype[i]; ++i) {
        //
        // Build the collection
        //
        if ('f' == testtype[i]) {
            mtcoll.clear();
            mtcoll.unused().setAllVars(DONT_CARE);
            RTYPE val;
            SG.vno2val(3, val);
            mtcoll.unused().setValue(val);
            mtcoll.pushUnused();
            out << " fully: ";
            out.flush();
        } else {
            if (1 == mtsizes[i]) mtcoll.clear();

            while (mtcoll.size() < mtsizes[i]) {
                SG.randomizeMinterm(mtcoll.unused(), rt);
                mtcoll.pushUnused();
            }

            out << "    ";
            out.put((unsigned long) mtsizes[i], 2);
            out << ": ";
            out.flush();
        }

#ifdef SHOW_MINTERMS
        out << "\nMinterms:\n";
        mtcoll.show(out);
#endif

        //
        // Set up ddedges
        //
        dd_edge Ef(Ff);
        dd_edge Eq(Fq);

        mtcoll.buildFunctionMax(RTYPE(0), Eq);
        out << "q ";
        out.flush();

        mtcoll.buildFunctionMax(RTYPE(0), Ef);
        out << "f ";
        out.flush();

        //
        // Brute force: compare functions
        //
        eval.setAllVars(0);
        do {
            RTYPE in_mtcoll, qval, fval;
            evaluate_set(eval, mtcoll, in_mtcoll);
            Eq.evaluate(eval, qval);
            Ef.evaluate(eval, fval);

            if (different(qval, in_mtcoll)) {
                mismatch(eval, 'Q', Eq, qval, mtcoll, in_mtcoll);
            }
            if (different(fval, in_mtcoll)) {
                mismatch(eval, 'F', Ef, fval, mtcoll, in_mtcoll);
            }
        } while (SG.nextMinterm(eval));
        out << "=\n";
        out.flush();

    } // for mtsize
}

/*
 *
 * Manipulate minterms for relations
 *
 */

bool matches_rel(const minterm &m, const minterm& va)
{
    for (unsigned i=1; i<=va.getNumVars(); i++) {
        if ((m.from(i) != DONT_CARE) && (m.from(i) != va.from(i)))  {
            return false;
        }
        if (m.to(i) == DONT_CARE) continue;
        if (m.to(i) == DONT_CHANGE) {
            if (va.to(i) == va.from(i)) continue;
        }
        if (m.to(i) != va.to(i)) return false;
    }
    return true;
}

template <typename RTYPE>
void evaluate_rel(const minterm &m, const minterm_coll &MC, RTYPE &ans)
{
    ans = 0;
    for (unsigned i=0; i<MC.size(); i++) {
        if (matches_rel(MC.at(i), m)) {
            ans = MAX(ans, RTYPE(MC.at(i).getValue()) );
        }
    }
}

/*
 *
 * Tests for relation-type forests
 *
 */



template <typename RTYPE>
void test_rels(domain* D, range_type rt)
{
    //
    // Build the various types of boolean forests
    //
    policies p;
    p.useDefaults(RELATION);

    p.setFullyReduced();

    forest* Ff = forest::create(D, RELATION, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setIdentityReduced();

    forest* Fi = forest::create(D, RELATION, rt,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* Fq = forest::create(D, RELATION, rt,
                    edge_labeling::MULTI_TERMINAL, p);


    //
    // Set up minterm collection
    // and evaluation minterm
    //
    minterm_coll mtcoll(32, D, RELATION);
    minterm eval(D, RELATION);

    //
    // For various collections of minterms,
    //  (1) build sets
    //  (2) explicitly verify sets against minterms

    ostream_output out(std::cout);

    out << "Checking matrices of " << nameOf(rt)
        << " built from minterm collections:\n";
    out.flush();

    char     testtype[] = { 'f', 'i', 'r', 'r', 'r', 'r', 'r', 'r',  0 };
    unsigned mtsizes[] =  {   0,   0,   1,   2,   4,   8,  16,  32     };

    for (unsigned i=0; testtype[i]; ++i) {
        //
        // Build the collection
        //
        if ('f' == testtype[i]) {
            mtcoll.clear();
            mtcoll.unused().setAllVars(DONT_CARE, DONT_CARE);
            RTYPE val;
            RG.vno2val(3, val);
            mtcoll.unused().setValue(val);
            mtcoll.pushUnused();
            out << " fully: ";
            out.flush();
        } else if ('i' == testtype[i]) {
            mtcoll.clear();
            mtcoll.unused().setAllVars(DONT_CARE, DONT_CHANGE);
            RTYPE val;
            RG.vno2val(5, val);
            mtcoll.unused().setValue(val);
            mtcoll.pushUnused();
            out << " ident: ";
            out.flush();

        } else {
            if (1 == mtsizes[i]) mtcoll.clear();

            while (mtcoll.size() < mtsizes[i]) {
                RG.randomizeMinterm(mtcoll.unused(), rt);
                mtcoll.pushUnused();
            }

            out << "    ";
            out.put((unsigned long) mtsizes[i], 2);
            out << ": ";
            out.flush();
        }

#ifdef SHOW_MINTERMS
        out << "\nMinterms:\n";
        mtcoll.show(out);
#endif

        //
        // Set up ddedges
        //
        dd_edge Ef(Ff);
        dd_edge Ei(Fi);
        dd_edge Eq(Fq);

        mtcoll.buildFunctionMax(RTYPE(0), Eq);
        out << "q ";
        out.flush();

        mtcoll.buildFunctionMax(RTYPE(0), Ef);
        out << "f ";
        out.flush();

        mtcoll.buildFunctionMax(RTYPE(0), Ei);
        out << "i ";
        out.flush();

        //
        // Brute force: compare functions
        //
        eval.setAllVars(0, 0);
        do {
            RTYPE in_mtcoll, qval, fval, ival;
            evaluate_rel(eval, mtcoll, in_mtcoll);
            Eq.evaluate(eval, qval);
            Ef.evaluate(eval, fval);
            Ei.evaluate(eval, ival);

            if (different(qval, in_mtcoll)) {
                mismatch(eval, 'Q', Eq, qval, mtcoll, in_mtcoll);
            }
            if (different(fval, in_mtcoll)) {
                mismatch(eval, 'F', Ef, fval, mtcoll, in_mtcoll);
            }
            if (different(ival, in_mtcoll)) {
                mismatch(eval, 'I', Ei, ival, mtcoll, in_mtcoll);
            }
        } while (RG.nextMinterm(eval));
        out << "=\n";
        out.flush();

    } // for mtsize
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
        test_sets<bool>(SD, range_type::BOOLEAN);
        test_sets<long>(SD, range_type::INTEGER);
        test_sets<double>(SD, range_type::REAL);
        domain::destroy(SD);
#endif
#ifdef TEST_RELS
        domain* RD = RG.makeDomain();
        test_rels<bool>(RD, range_type::BOOLEAN);
        test_rels<long>(RD, range_type::INTEGER);
        test_rels<double>(RD, range_type::REAL);
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
#endif
