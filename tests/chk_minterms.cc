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
const int SETVARS = 8;
const int RELVARS = 4;

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

void nextTermValue(minterm &m, range_type rt)
{
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
}



/*
 *
 * Manipulate minterms for sets
 *
 */

bool nextSetMinterm(minterm& m)
{
    for (unsigned i=1; i<=m.getNumVars(); i++) {
        if (++m.from(i) < DOMSIZE) {
            return true;
        }
        m.from(i) = 0;
    }
    return false;
}

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
            RTYPE tmp;
            MC.at(i).getTerm().getValue(tmp);
            ans = MAX(ans, tmp);
        }
    }
}




/*
 *
 * Tests for set-type forests
 *
 */

template <typename RTYPE>
void test_sets(range_type rt)
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
            nextTermValue(mtcoll.unused(), rt);
            mtcoll.pushUnused();
            out << " fully: ";
            out.flush();
        } else {
            if (1 == mtsizes[i]) mtcoll.clear();

            while (mtcoll.size() < mtsizes[i]) {
                randomizeMinterm(mtcoll.unused(), rt);
                nextTermValue(mtcoll.unused(), rt);
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
        } while (nextSetMinterm(eval));
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
                mtcoll.at(i).getTerm().getValue(in_mtcoll);
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
        } while (nextSetMinterm(eval));
        out << "=\n";
        out.flush();
    }

    domain::destroy(D);
}

/*
 *
 * Manipulate minterms for relations
 *
 */

bool nextRelMinterm(minterm& m)
{
    for (unsigned i=1; i<=m.getNumVars(); i++) {
        if (++m.from(i) < DOMSIZE) {
            return true;
        }
        m.from(i) = 0;
        if (++m.to(i) < DOMSIZE) {
            return true;
        }
        m.to(i) = 0;
    }
    return false;
}

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
            RTYPE tmp;
            MC.at(i).getTerm().getValue(tmp);
            ans = MAX(ans, tmp);
        }
    }
}

/*
 *
 * Tests for relation-type forests
 *
 */



template <typename RTYPE>
void test_rels(range_type rt)
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
            nextTermValue(mtcoll.unused(), rt);
            mtcoll.pushUnused();
            out << " fully: ";
            out.flush();
        } else if ('i' == testtype[i]) {
            mtcoll.clear();
            mtcoll.unused().setAllVars(DONT_CARE, DONT_CHANGE);
            nextTermValue(mtcoll.unused(), rt);
            mtcoll.pushUnused();
            out << " ident: ";
            out.flush();

        } else {
            if (1 == mtsizes[i]) mtcoll.clear();

            while (mtcoll.size() < mtsizes[i]) {
                randomizeMinterm(mtcoll.unused(), rt);
                nextTermValue(mtcoll.unused(), rt);
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
        } while (nextRelMinterm(eval));
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
                mtcoll.at(i).getTerm().getValue(in_mtcoll);
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
        } while (nextRelMinterm(eval));
        out << "=\n";
        out.flush();
    }

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
        test_sets<bool>(range_type::BOOLEAN);
        test_sets<long>(range_type::INTEGER);
        test_sets<double>(range_type::REAL);
#endif
#ifdef TEST_RELS
        test_rels<bool>(range_type::BOOLEAN);
        test_rels<long>(range_type::INTEGER);
        test_rels<double>(range_type::REAL);
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
