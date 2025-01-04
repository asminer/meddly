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
#define SHOW_MINTERMS

#include "randomize.h"

//
// Minterm generators :)
//

vectorgen SG(MEDDLY::SET, 8, 4, 5);
vectorgen RG(MEDDLY::RELATION, 5, 3, 5);

const unsigned NUM_MINTERMS = 9;

using namespace MEDDLY;

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
            out << "\n    DD says   ";
            Eval.write(out);
            out << "\n    Should be ";
            Fval.write(out);
            out << "\n    Var assignments: ";
            vars.show(out);
            out << "\n    Generated using: ";
            m.show(out);
            out << "deflt ";
            deflt.write(out);
            out << "\n    Func. DD encoding:\n";
            E.showGraph(out);
            out.flush();
            throw "mismatch";
        }

    } while (SG.nextMinterm(vars));
}

void test_sets(forest* F, rangeval deflt)
{
    ostream_output out(std::cout);

    out << "Checking " << nameOf(F->getRangeType())
        << " " << shortNameOf(F->getReductionRule())
        << " " << nameOf(F->getEdgeLabeling()) << " mdd, default ";
    deflt.write(out);
    out << "\n    ";
    out.flush();

    //
    // Determine non-default value
    //
    rangeval fval;
    switch (deflt.getType()) {
        case range_type::BOOLEAN:
            fval = ! bool(deflt);
            break;

        case range_type::INTEGER:
            fval = 2 + long(deflt);
            break;

        case range_type::REAL:
            fval = 2.5 + double(deflt);
            break;

        default:
            throw "unknown type";
    }

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
    m.buildFunction(E); // TBD - default!
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
        m.buildFunction(E); // TBD - default!
        check_equal_set(E, m, deflt);
    }
    out << "\n";
}

/*
 *
 * Testing for sets
 *
 */

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

        //
        // Build fully-reduced set forests
        //
        Sp.setFullyReduced();

        forest* S_fr_mdd_b = forest::create(SD, SET, range_type::BOOLEAN,
                                edge_labeling::MULTI_TERMINAL, Sp);

        forest* S_fr_mdd_i = forest::create(SD, SET, range_type::INTEGER,
                                edge_labeling::MULTI_TERMINAL, Sp);

        forest* S_fr_mdd_r = forest::create(SD, SET, range_type::REAL,
                                edge_labeling::MULTI_TERMINAL, Sp);

        //
        // Build quasi-reduced set forests
        //

        Sp.setQuasiReduced();

        forest* S_qr_mdd_b = forest::create(SD, SET, range_type::BOOLEAN,
                                edge_labeling::MULTI_TERMINAL, Sp);

        forest* S_qr_mdd_i = forest::create(SD, SET, range_type::INTEGER,
                                edge_labeling::MULTI_TERMINAL, Sp);

        forest* S_qr_mdd_r = forest::create(SD, SET, range_type::REAL,
                                edge_labeling::MULTI_TERMINAL, Sp);

        //
        // Run tests
        //

        test_sets(S_fr_mdd_b, false);
        test_sets(S_qr_mdd_b, false);

        test_sets(S_fr_mdd_b, true);
        test_sets(S_qr_mdd_b, true);

        // TBD - integer and real

        domain::destroy(SD);
#endif
#ifdef TEST_RELS
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
// OLD BELOW HERE
//

#if 0

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

template <typename RTYPE>
void mismatch(const minterm &eval,
        char mddtype, const MEDDLY::dd_edge &E, RTYPE val_mdd,
        const minterm &mt, RTYPE mtcoll_answer)
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
    out << "Minterm:\n  ";
    mt.show(out);

    out << "\n" << mddtype << MDD << ":\n";
    E.showGraph(out);
    out.flush();

    throw "mismatch";
}


/*
 *
 * Manipulate minterms for sets
 *
 */

bool matches_set(const minterm &m, const minterm& va)
{
    for (unsigned i=1; i<=va.getNumVars(); i++) {
        if (m.from(i) == DONT_CARE) continue;
        if (m.from(i) != va.from(i)) return false;
    }
    return true;
}

template <typename RTYPE>
void evaluate_set(const minterm &m, const minterm_coll &MC, RTYPE &ans)
{
    ans = 0;
    for (unsigned i=0; i<MC.size(); i++) {
        if (matches_set(MC.at(i), m)) {
            ans = MAX(ans, RTYPE(MC.at(i).getValue()) );
        }
    }
}




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

        mtcoll.buildFunction(Eq);
        out << "q ";
        out.flush();

        mtcoll.buildFunction(Ef);
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

    out << "Checking vectors of " << nameOf(rt)
        << " built from single minterms:\n";
    out.flush();
    for (unsigned i=0; i<5; ++i) {
        //
        // Set up ddedges
        //
        dd_edge Ef(Ff);
        dd_edge Eq(Fq);

        mtcoll.at(i).buildFunction(Eq);
        out << "        q ";
        out.flush();

        mtcoll.at(i).buildFunction(Ef);
        out << "f ";
        out.flush();

        //
        // Brute force: compare functions
        //
        eval.setAllVars(0);
        do {
            RTYPE in_mtcoll, qval, fval;
            if (matches_set(mtcoll.at(i), eval)) {
                in_mtcoll = mtcoll.at(i).getValue();
            } else {
                in_mtcoll = 0;
            }
            Eq.evaluate(eval, qval);
            Ef.evaluate(eval, fval);

            if (different(qval, in_mtcoll)) {
                mismatch(eval, 'Q', Eq, qval, mtcoll.at(i), in_mtcoll);
            }
            if (different(fval, in_mtcoll)) {
                mismatch(eval, 'F', Ef, fval, mtcoll.at(i), in_mtcoll);
            }
        } while (SG.nextMinterm(eval));
        out << "=\n";
        out.flush();
    }
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

        mtcoll.buildFunction(Eq);
        out << "q ";
        out.flush();

        mtcoll.buildFunction(Ef);
        out << "f ";
        out.flush();

        mtcoll.buildFunction(Ei);
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

    out << "Checking matrices of " << nameOf(rt)
        << " built from single minterms:\n";
    out.flush();
    for (unsigned i=0; i<5; ++i) {
        //
        // Set up ddedges
        //
        dd_edge Ef(Ff);
        dd_edge Eq(Fq);
        dd_edge Ei(Fi);

        mtcoll.at(i).buildFunction(Eq);
        out << "        q ";
        out.flush();

        mtcoll.at(i).buildFunction(Ef);
        out << "f ";
        out.flush();

        mtcoll.at(i).buildFunction(Ei);
        out << "i ";
        out.flush();

        //
        // Brute force: compare functions
        //
        eval.setAllVars(0, 0);
        do {
            RTYPE in_mtcoll, qval, fval, ival;
            if (matches_rel(mtcoll.at(i), eval)) {
                in_mtcoll = mtcoll.at(i).getValue();
            } else {
                in_mtcoll = 0;
            }
            Eq.evaluate(eval, qval);
            Ef.evaluate(eval, fval);
            Ei.evaluate(eval, ival);

            if (different(qval, in_mtcoll)) {
                mismatch(eval, 'Q', Eq, qval, mtcoll.at(i), in_mtcoll);
            }
            if (different(fval, in_mtcoll)) {
                mismatch(eval, 'F', Ef, fval, mtcoll.at(i), in_mtcoll);
            }
            if (different(ival, in_mtcoll)) {
                mismatch(eval, 'I', Ei, ival, mtcoll.at(i), in_mtcoll);
            }
        } while (RG.nextMinterm(eval));
        out << "=\n";
        out.flush();
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
