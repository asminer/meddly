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

#define DEBUG_MINTERM_COLL

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

bool evaluate(const minterm &m, const minterm_coll &MC)
{
    for (unsigned i=0; i<MC.size(); i++) {
        if (MC.at(i)->contains(m)) return true;
    }
    return false;
}


void zeroMinterm(minterm& m)
{
    for (unsigned i=1; i<=m.getNumVars(); i++) {
        m.setVar(i, 0);
    }
}

bool nextMinterm(minterm& m)
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

void check_set_eq(char mddtype, const MEDDLY::dd_edge &E,
        minterm &eval, minterm_coll &mtcoll)
{
    bool val_mdd, val_mts;
    val_mts = evaluate(eval, mtcoll);
    E.evaluate(eval, val_mdd);

    if (val_mdd != val_mts) {
        MEDDLY::ostream_output out(std::cout);
        out << "\nMismatch on ";
        eval.show(out);
        out << "\n  " << mddtype << "MDD: " << (val_mdd ? "true" : "false") << "\n";
        out << "mtlist: " << (val_mts ? "true" : "false") << "\n";

        out << "\n";
        out << "Minterm list:\n";
        mtcoll.show(out, nullptr, "\n");

        out << mddtype << "MDD:\n";
        E.showGraph(out);
        out.flush();

        throw "mismatch";
    }
}

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

    std::cout << "Checking sets built from minterms:\n";
    std::cout.flush();
    for (unsigned mtsize=1; mtsize<=MAXTERMS; mtsize*=2)
    {
        std::cout << "    ";
        if (mtsize<10) std::cout << ' ';
        std::cout << mtsize << ": ";
        std::cout.flush();

        //
        // Add minterms to the collection
        //
        while (mtcoll.size() < mtsize) {
            minterm* m = mtcoll.makeTemp();
            randomSetMinterm(m);
            mtcoll.addToCollection(m);
        }

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
        std::cout << "q ";
        std::cout.flush();

        apply(UNION, Ef, mtcoll, Ef);
        std::cout << "f ";
        std::cout.flush();

        //
        // Brute force: compare functions
        //
        zeroMinterm(eval);
        do {
            check_set_eq('Q', Eq, eval, mtcoll);
            check_set_eq('F', Ef, eval, mtcoll);
        } while (nextMinterm(eval));
        std::cout << "=" << std::endl;

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
    //
    // Separated
    //
    //
    //

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

void randomRelMinterm(int* un, int* pr, unsigned vars)
{
    //
    // Separated
    //
    //
    //

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

    for (unsigned i=1; i<=vars; i++) {
        int index = Equilikely(0, 51);
        un[i] = unvals[index];
        pr[i] = prvals[index];
#ifdef RESTRICT_DONT_CHANGE
        if (pr[i] == -2) {
            un[i] = -1;
        }
#endif
    }
}

void showMinterm(const int* un, const int* pr, unsigned vars)
{
    using namespace std;

    cout << "[ bot";
    for (unsigned i=1; i<=vars; i++) {
        cout << ", ";
        if (un[i] < 0)  cout << 'x';
        else            cout << un[i];
        cout << "->";
        if (pr[i] == -2)    cout << 'i';
        else if (pr[i] < 0) cout << 'x';
        else                cout << pr[i];
    }
    cout << "]";
}

bool mintermMatches(const int* mtun, const int* mtpr,
        const int* valun, const int* valpr, unsigned vars)
{
    for (unsigned i=1; i<=vars; i++) {
        if ((mtun[i] != -1) && (mtun[i] != valun[i])) return false;
        if (mtpr[i] == -1) continue;
        if (mtpr[i] == -2) {
            if (valpr[i] == valun[i]) continue;
        }
        if (mtpr[i] != valpr[i]) return false;
    }
    return true;
}

bool evaluate(const int* vun, const int* vpr, unsigned vars,
        int** mtun, int** mtpr, unsigned nmt)
{
    for (unsigned i=0; i<nmt; i++) {
        if (mintermMatches(mtun[i], mtpr[i], vun, vpr, vars)) return true;
    }
    return false;
}

void zeroMinterm(int* vun, int* vpr, unsigned vars)
{
    for (unsigned i=1; i<=vars; i++) {
        vun[i] = 0;
        vpr[i] = 0;
    }
}

bool nextMinterm(int* vun, int* vpr, unsigned vars)
{
    for (unsigned i=1; i<=vars; i++) {
        vpr[i]++;
        if (vpr[i] < DOMSIZE) return true;
        vpr[i] = 0;
        vun[i]++;
        if (vun[i] < DOMSIZE) return true;
        vun[i] = 0;
    }
    return false;
}

/*
 *
 * Tests for relation-type forests
 *
 */

void check_rel_eq(char mxdtype, const MEDDLY::dd_edge &E,
        const int* uneval, const int* preval,
        int** unlist, int** prlist, const unsigned mtsize)
{
    bool val_mxd;
    bool val_mts = evaluate(uneval, preval, RELVARS, unlist, prlist, mtsize);
    const MEDDLY::forest* F = E.getForest();
    F->evaluate(E, uneval, preval, val_mxd);

    if (val_mxd != val_mts) {
        std::cout << "\nMismatch on ";
        showMinterm(uneval, preval, RELVARS);
        std::cout << "\n  " << mxdtype << "MxD: " << (val_mxd ? "true" : "false") << "\n";
        std::cout << "mtlist: " << (val_mts ? "true" : "false") << "\n";

        std::cout << "\n";
        std::cout << "Minterm list:\n";

        for (unsigned n=0; n<mtsize; n++) {
            std::cout << "    ";
            showMinterm(unlist[n], prlist[n], RELVARS);
            std::cout << "\n";
        }

        std::cout << mxdtype << "Mxd:\n";
        MEDDLY::ostream_output out(std::cout);
        E.showGraph(out);

        throw "mismatch";
    }

}


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
    // Build random minterms
    //
    int* unlist[MAXTERMS];
    int* prlist[MAXTERMS];
    for (unsigned i=0; i<MAXTERMS; i++) {
        unlist[i] = new int[1+RELVARS];
        prlist[i] = new int[1+RELVARS];
        randomRelMinterm(unlist[i], prlist[i], RELVARS);
    }

    int uneval[1+RELVARS];
    int preval[1+RELVARS];

#ifdef DEBUG_MINTERM_COLL
    minterm_coll mtcoll(D, RELATION);
    for (unsigned i=0; i<MAXTERMS; i++) {
        minterm* m = mtcoll.makeTemp();
        randomRelMinterm(m);
        mtcoll.addToCollection(m);
    }
    ostream_output out(std::cout);
    out << "Built relation minterm collection:\n";
    mtcoll.show(out, nullptr, "\n");
    out.flush();
    out << "Sorting...\n";
    mtcoll.sort();
    out << "Sorted collection:\n";
    mtcoll.show(out, nullptr, "\n");
    out.flush();
    return;
#endif


    //
    // For various collections of minterms,
    //  (1) build relations
    //  (2) explicitly verify sets against minterms

    std::cout << "Checking relations built from minterms:\n";
    std::cout.flush();
    for (unsigned mtsize=1; mtsize<=MAXTERMS; mtsize*=2)
    {
        std::cout << "    ";
        if (mtsize<10) std::cout << ' ';
        std::cout << mtsize << ": ";
        std::cout.flush();
        dd_edge Ef(Ff);
        dd_edge Eq(Fq);

        Fq->createEdge(unlist, prlist, mtsize, Eq);
        std::cout << "q ";
        std::cout.flush();
        Ff->createEdge(unlist, prlist, mtsize, Ef);
        std::cout << "f ";
        std::cout.flush();

        zeroMinterm(uneval, preval, RELVARS);
        do {
            check_rel_eq('Q', Eq, uneval, preval, unlist, prlist, mtsize);
            check_rel_eq('F', Ef, uneval, preval, unlist, prlist, mtsize);
        } while (nextMinterm(uneval, preval, RELVARS));
        std::cout << "=" << std::endl;
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
        if (seed < 0) seed *= -1;
        if (0==seed)  seed = 12345; // probably never happen
    }
    cout << "Using rng seed " << seed << "\n";

    try {
        MEDDLY::initialize();
        test_sets();
    //    test_rels();
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
