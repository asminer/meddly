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

const unsigned MAXTERMS = 32;
// const unsigned MAXTERMS = 4;

const int DOMSIZE = 4;       // DO NOT change
const int SETVARS = 10;
const int RELVARS = 6;

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


/*
 *
 * Manipulate minterms for sets
 *
 */

void randomSetMinterm(int* mt, unsigned vars)
{
    const int vals[9] = { -1, 0, 0, 1, 1, 2, 2, 3, 3 };

    for (unsigned i=1; i<=vars; i++) {
        int index = Equilikely(0, 8);
        mt[i] = vals[index];
    }
}

void showMinterm(const int* mt, unsigned vars)
{
    using namespace std;

    cout << "[ bot";
    for (unsigned i=1; i<=vars; i++) {
        cout << ", ";
        if (mt[i] < 0)  cout << 'x';
        else            cout << mt[i];
    }
    cout << "]";
}

bool mintermMatches(const int* mt, const int* vals, unsigned vars)
{
    for (unsigned i=1; i<=vars; i++) {
        if (mt[i] == -1) continue;  // don't care
        if (mt[i] != vals[i]) return false;
    }
    return true;
}

bool evaluate(const int* vals, unsigned vars, int** mt, unsigned nmt)
{
    for (unsigned i=0; i<nmt; i++) {
        if (mintermMatches(mt[i], vals, vars)) return true;
    }
    return false;
}

void zeroMinterm(int* mt, unsigned vars)
{
    for (unsigned i=1; i<=vars; i++) {
        mt[i] = 0;
    }
}

bool nextMinterm(int* mt, unsigned vars)
{
    for (unsigned i=1; i<=vars; i++) {
        mt[i]++;
        if (mt[i] < DOMSIZE) return true;
        mt[i] = 0;
    }
    return false;
}

/*
 *
 * Manipulate minterms for relations
 *
 */

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

    int unterm[1+RELVARS];
    int prterm[1+RELVARS];

    for (unsigned n=0; n<MAXTERMS; n++) {
        randomRelMinterm(unterm, prterm, RELVARS);
        cout << "    ";
        showMinterm(unterm, prterm, RELVARS);
        cout << "\n";
    }


    /*

    int* mtlist[MAXTERMS];
    for (unsigned i=0; i<MAXTERMS; i++) {
        mtlist[i] = new int[1+SETVARS];
    }

    cout << "Minterms:\n";
    for (unsigned n=0; n<MAXTERMS; n++) {
        randomSetMinterm(mtlist[n], SETVARS);
        cout << "    ";
        showMinterm(mtlist[n], SETVARS);
        cout << "\n";
    }

    cout << "Set:\n";
    int vars[1+SETVARS];
    zeroMinterm(vars, SETVARS);
    unsigned card = 0;
    do {
        if (evaluate(vars, SETVARS, mtlist, MAXTERMS)) {
            ++card;
            cout << "    ";
            showMinterm(vars, SETVARS);
            cout << "\n";

        }
    } while (nextMinterm(vars, SETVARS));
    cout << "Cardinality " << card << "\n";

    */

    //

    return 0;
}
