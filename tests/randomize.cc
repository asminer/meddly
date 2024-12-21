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

#include "randomize.h"
#include <time.h>
#include <iostream>

long vectorgen_base::seed;

vectorgen_base::vectorgen_base(unsigned v, unsigned rd)
    : VARS(v), RELDOM(rd)
{
    if (rd < 2) {
        throw "RELDOM too small (min 2)";
    }
    if (rd > 15) {
        throw "RELDOM too large";
    }
    if (v < 2) {
        throw "VARS too small (min 2)";
    }
    if (v > 15) {
        throw "VARS too large";
    }

    //
    // Compute POTENTIAL = (rd*rd) ^ v
    //

    const unsigned rdrd = rd*rd;
    POTENTIAL = 1;
    while (v--) {
        POTENTIAL *= rdrd;
    }

    if (POTENTIAL > 99999999) {
        throw "POTENTIAL too large";
    }
}

void vectorgen_base::setSeed(long s, bool print)
{
    if (s<1) {
        s = time(NULL);
        // on the off chance that time(NULL) returns 0 or negative:
        if (s < 0) s *= -1;
        if (0 == s) s = 12345;
    }
    seed = s;
    if (print) {
        std::cout << "Using rng seed " << seed << "\n";
    }
}

double vectorgen_base::random()
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

void vectorgen_base::index2minterm(unsigned x, MEDDLY::minterm &m)
{
    if (m.getNumVars() != vars()) {
        throw "minterm mismatch in call to index2minterm";
    }
    if (m.isForSets()) {
        for (unsigned j=1; j<=m.getNumVars(); j++) {
            m.setVar(j, x % setdom());
            x /= setdom();
        }
    } else {
        for (unsigned j=1; j<=m.getNumVars(); j++) {
            const int from = x % reldom();
            x /= reldom();
            const int to = x % reldom();
            x /= reldom();
            m.setVars(j, from, to);
        }
    }
}

void vectorgen_base::buildIdentityFromIndex(unsigned ndx,
        unsigned k1, unsigned k2, std::vector <unsigned> &ilist) const
{
}

void vectorgen_base::buildFullyFromIndex(unsigned ndx,
        unsigned k1, unsigned k2, std::vector <unsigned> &ilist) const
{
    if (k1 < k2) {
        MEDDLY::SWAP(k1, k2);
    }
    if (k1 == k2) {
        k2 = 0;
    }

    //
    // Conceptually, ndx is equivalent to a VARS digit number
    // in base SETDOM. We want to get the values of those numbers
    // with digits k1 and k2 going over all possible digits.
    //
    // To do this, we determine the value when digits k1 and k2 are 0
    // (variable mindex), and multipliers mult1 and mult2,
    // such that the numbers we want are all those of the form
    //      mindex + mult1 * digit1 + mult2 * digit2
    //
    unsigned mindex = ndx;
    unsigned mult   = 1;
    unsigned mult1  = 0;
    unsigned mult2  = 0;
    unsigned x = ndx;
    for (unsigned j=1; j<=vars(); j++) {
        if (k1 == j) {
            mult1 = mult;
            unsigned d = x % setdom();
            mindex -= d * mult;
        }
        if (k2 == j) {
            mult2 = mult;
            unsigned d = x % setdom();
            mindex -= d * mult;
        }
        x /= setdom();
        mult *= setdom();
    }


    //
    // Ready to build the list. Clear it and then
    // loop over pairs
    //
    ilist.clear();
    for (unsigned v1=0; v1<setdom(); v1++) {
        for (unsigned v2=0; v2<setdom(); v2++) {
            ilist.push_back( mindex + mult1*v1 + mult2*v2 );
            if (0==mult2) break;
        }
        if (0==mult1) break;
    }
}

