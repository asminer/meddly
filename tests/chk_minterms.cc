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

// #define SHOW_MINTERMS

const unsigned MAXTERMS = 32;

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

/*
 *  common for sets and relations
 */

bool evaluate(const minterm &m, const minterm_coll &MC)
{
    for (unsigned i=0; i<MC.size(); i++) {
        if (MC.at(i)->contains(m)) return true;
    }
    return false;
}

void check_equal(char mddtype, const MEDDLY::dd_edge &E,
        minterm &eval, minterm_coll &mtcoll)
{
    const char* MDD = eval.isForRelations() ? "MxD" : "MDD";
    bool val_mdd, val_mts;
    val_mts = evaluate(eval, mtcoll);
    E.evaluate(eval, val_mdd);

    if (val_mdd != val_mts) {
        MEDDLY::ostream_output out(std::cout);
        out << "\nMismatch on ";
        eval.show(out);
        out << "\n  " << mddtype << MDD << ": "
            << (val_mdd ? "true" : "false") << "\n";
        out << "mtlist: " << (val_mts ? "true" : "false") << "\n";
        out << "\n";
        out << "Minterm list:\n";
        mtcoll.show(out, nullptr, "\n");

        out << mddtype << MDD << ":\n";
        E.showGraph(out);
        out.flush();

        throw "mismatch";
    }
}


/*
 *
 * Manipulate minterms for sets
 *
 */

void randomSetMinterm(minterm* m)
{
    const int vals[9] = { -1, 0, 0, 1, 1, 2, 2, 3, 3 };

    for (unsigned i=1; i<=m->getNumVars(); i++) {
        int index = Equilikely(0, 8);
        m->setVar(i, vals[index]);
    }
}

void zeroSetMinterm(minterm& m)
{
    for (unsigned i=1; i<=m.getNumVars(); i++) {
        m.setVar(i, 0);
    }
}

bool nextSetMinterm(minterm& m)
{
    for (unsigned i=1; i<=m.getNumVars(); i++) {
        int v = m.getVar(i);
        v++;
        if (v < DOMSIZE) {
            m.setVar(i, v);
            return true;
        }
        m.setVar(i, 0);
    }
    return false;
}


/*
 *
 * Tests for set-type forests
 *
 */

void test_sets()
{
    using namespace MEDDLY;

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

    forest* Ff = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* Fq = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);


    //
    // Set up minterm collection
    // and evaluation minterm
    //
    minterm_coll mtcoll(D, SET);
    minterm eval(D, SET, FULL_ONLY);

    //
    // For various collections of minterms,
    //  (1) build sets
    //  (2) explicitly verify sets against minterms

    ostream_output out(std::cout);
    out << "Checking sets built from minterms:\n";
    out.flush();
    for (unsigned mtsize=1; mtsize<=MAXTERMS; mtsize*=2)
    {
        //
        // Add minterms to the collection
        //
        while (mtcoll.size() < mtsize) {
            minterm* m = mtcoll.makeTemp();
            randomSetMinterm(m);
            mtcoll.addToCollection(m);
        }

#ifdef SHOW_MINTERMS
        mtcoll.sort();
        mtcoll.show(out, nullptr, "\n");
#endif
        out << "    ";
        out.put((unsigned long) mtsize, 2);
        out << ": ";
        out.flush();


        //
        // Set up ddedges
        //
        dd_edge Ef(Ff);
        dd_edge Eq(Fq);
        Ff->createEdge(false, Ef);
        Fq->createEdge(false, Eq);

        //
        // Build unions
        //
        apply(UNION, Eq, mtcoll, Eq);
        out << "q ";
        out.flush();

        apply(UNION, Ef, mtcoll, Ef);
        out << "f ";
        out.flush();

        //
        // Brute force: compare functions
        //
        zeroSetMinterm(eval);
        do {
            check_equal('Q', Eq, eval, mtcoll);
            check_equal('F', Ef, eval, mtcoll);
        } while (nextSetMinterm(eval));
        out << "=\n";
        out.flush();

    } // for mtsize

    domain::destroy(D);
}

/*
 *
 * Manipulate minterms for relations
 *
 */

void randomRelMinterm(MEDDLY::minterm* m)
{
    const int unvals[52] =
    { -1, -1, -1, -1,   // 4 (x,x) pairs
      -1, -1, -1, -1,   // 4 (x,i) pairs
      -1, -1, -1, -1,   // 4 (x, normal) pairs
       0,  1,  2,  3,   // 4 (normal, x) pairs
       0,  1,  2,  3,   // 4 (normal, i) pairs
       0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3,   // 16 normal pairs
       0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3};  // 16 normal pairs

    const int prvals[52] =
    { -1, -1, -1, -1,   // 4 (x,x) pairs
      -2, -2, -2, -2,   // 4 (x,i) pairs
       0,  1,  2,  3,   // 4 (x, normal) pairs
      -1, -1, -1, -1,   // 4 (normal, x) pairs
      -2, -2, -2, -2,   // 4 (normal, i) pairs
       0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,   // 16 normal pairs
       0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3};  // 16 normal pairs

    for (unsigned i=1; i<=m->getNumVars(); i++) {
        int index = Equilikely(0, 51);

        m->setVars(i, unvals[index], prvals[index]);
    }
}

void zeroRelMinterm(minterm& m)
{
    for (unsigned i=1; i<=m.getNumVars(); i++) {
        m.setVars(i, 0, 0);
    }
}

bool nextRelMinterm(minterm& m)
{
    for (unsigned i=1; i<=m.getNumVars(); i++) {
        int v = m.getFrom(i);
        v++;
        if (v < DOMSIZE) {
            m.setFrom(i, v);
            return true;
        }
        m.setFrom(i, 0);
        v = m.getTo(i);
        v++;
        if (v < DOMSIZE) {
            m.setTo(i, v);
            return true;
        }
        m.setTo(i, 0);
    }
    return false;
}



/*
 *
 * Tests for relation-type forests
 *
 */



void test_rels()
{
    using namespace MEDDLY;

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
    p.useDefaults(SET);

    forest* Ff = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    p.setQuasiReduced();

    forest* Fq = forest::create(D, RELATION, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);


    //
    // Set up minterm collection
    // and evaluation minterm
    //
    minterm_coll mtcoll(D, RELATION);
    minterm eval(D, RELATION, FULL_ONLY);

    //
    // For various collections of minterms,
    //  (1) build sets
    //  (2) explicitly verify sets against minterms

    ostream_output out(std::cout);
    out << "Checking relations built from minterms:\n";
    out.flush();

        //
        // TBD
        //
    for (unsigned mtsize=1; mtsize<=MAXTERMS; ++mtsize)
    {
        minterm* m = mtcoll.makeTemp();
        randomRelMinterm(m);
        mtcoll.addToCollection(m);
    }

    mtcoll.sort();
    mtcoll.show(out, nullptr, "\n");

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
        if (seed < 0) seed *= -1;
        if (0==seed)  seed = 12345; // probably never happen
    }
    cout << "Using rng seed " << seed << "\n";

    try {
        MEDDLY::initialize();
        test_sets();
        test_rels();
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
