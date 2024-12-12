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
// #define TEST_RELS
#define SHOW_MINTERMS

const int DOMSIZE = 4;       // DO NOT change
const int SETVARS = 10;
const int RELVARS = 5;

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

void randomizeMinterm(MEDDLY::minterm &m)
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
 * Tests for set-type forests
 *
 */

void test_sets(char reduction)
{
    using namespace MEDDLY;

    ostream_output out(std::cout);

    out << "Checking iterators over ";
    out << ( ('f' == reduction) ? "fully" : "quasi" );
    out << "-reduced sets:\n";
    out.flush();


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

    if ('f' == reduction) {
        p.setFullyReduced();
    } else {
        p.setQuasiReduced();
    }

    forest* F = forest::create(D, SET, range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL, p);

    //
    // Set up minterm collection
    // and evaluation minterm
    //
    minterm_coll mtcoll(16, D, SET);
    minterm eval(D, SET);

    //
    // Build the collection
    //

    while (mtcoll.size() < mtcoll.maxsize()) {
        randomizeMinterm(mtcoll.unused());
        mtcoll.pushUnused();
    }

#ifdef SHOW_MINTERMS
    out << "\nMinterms:\n";
    mtcoll.show(out);
#endif

    //
    // Set up ddedges
    //
    dd_edge E(F);
    mtcoll.buildFunction(E);

    out << "Iterating...\n";

    unsigned long count = 0;
    for (dd_edge::iterator i = E.begin(); i != E.end(); i++) {
        out.put(count, 5);
        out << ": ";
        (*i).show(out);
        out << "\n";
        out.flush();
        ++count;
    }

    domain::destroy(D);
}



void test_rels()
{
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
        test_sets('q');
        test_sets('f');
#endif
#ifdef TEST_RELS
        test_rels();
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
