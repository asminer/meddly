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

vectorgen_base::vectorgen_base(bool sr, unsigned v, unsigned d, unsigned r)
    : is_for_relations(sr), VARS(v), DOM(d), RANGE(r)
{
    if (d < 2) {
        throw "DOM too small (min 2)";
    }
    if (d > 15) {
        throw "DOM too large";
    }
    if (v < 2) {
        throw "VARS too small (min 2)";
    }
    if (v > 15) {
        throw "VARS too large";
    }

    //
    // Compute POTENTIAL = rd ^ v for sets, rd ^ 2v for relations
    //

    POTENTIAL = 1;
    if (is_for_relations) {
        v *= 2;
    }
    while (v--) {
        POTENTIAL *= DOM;
    }

    if (POTENTIAL > 99999999) {
        throw "POTENTIAL too large";
    }

    current_terminal = 0;
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

MEDDLY::domain* vectorgen_base::makeDomain()
{
    int bs[vars()];
    for (unsigned i=0; i<vars(); i++) {
        bs[i] = dom();
    }
    return MEDDLY::domain::createBottomUp(bs, vars());
}

void vectorgen_base::randomizeMinterm(MEDDLY::minterm &m, MEDDLY::range_type rt)
{
    const int X = MEDDLY::DONT_CARE;
    const int I = MEDDLY::DONT_CHANGE;

    if (m.isForRelations()) {
        for (unsigned i=vars(); i; --i) {
            int val = Equilikely_I(-5, 2*dom()-1);
            switch (val) {
                case -5:    // (x,x) pair
                            m.setVars(i, X, I);
                            continue;

                case -4:    // (x,i) pair
                            m.setVars(i, X, I);
                            continue;

                case -3:    // (x,normal) pair
                            m.setVars(i, X, Equilikely_I(0, dom()-1));
                            continue;

                case -2:    // (normal,x) pair
                            m.setVars(i, Equilikely_I(0, dom()-1), X);
                            continue;

                case -1:    // (normal,i) pair
                            m.setVars(i, Equilikely_I(0, dom()-1), I);
                            continue;

                default:    // (normal, normal) pair
                            m.setVars(i, val/2, Equilikely_I(0, dom()-1));

            }
        }
    } else {
        ///
        /// Roll a die with -1 on one side,
        /// each possible value from 0..setdom()-1
        /// on two sides.
        ///
        for (unsigned i=vars(); i; --i) {
            int val = Equilikely_I(-1, 2*dom()-1);
            if (val>=0) {
                val /= 2;
            } else {
                val = X;
            }
            m.setVar(i, val);
        }
    }

    //
    // Randomize the terminal
    //
    current_terminal = 1 + (current_terminal % RANGE);

    int iterm;
    float fterm;
    switch (rt) {
        case MEDDLY::range_type::INTEGER:
            vno2val(current_terminal, iterm);
            m.setTerm(iterm);
            break;

        case MEDDLY::range_type::REAL:
            vno2val(current_terminal, fterm);
            m.setTerm(fterm);
            break;

        default:
            m.setTerm(true);
    }
}

void vectorgen_base::index2minterm(unsigned x, MEDDLY::minterm &m)
{
    if (m.getNumVars() != vars()) {
        throw "minterm mismatch in call to index2minterm";
    }
    if (m.isForSets()) {
        for (unsigned j=1; j<=m.getNumVars(); j++) {
            m.setVar(j, x % dom());
            x /= dom();
        }
    } else {
        for (unsigned j=1; j<=m.getNumVars(); j++) {
            const int from = x % dom();
            x /= dom();
            const int to = x % dom();
            x /= dom();
            m.setVars(j, from, to);
        }
    }
}

void vectorgen_base::buildIdentityFromIndex(unsigned ndx,
        unsigned k1, unsigned k2, std::vector <unsigned> &ilist) const
{
    if (isForSets()) {
        throw "buildIdentityFromIndex called on set";
    }
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
    const unsigned dom_squared = dom() * dom();
    for (unsigned j=1; j<=vars(); j++) {
        if (k1 == j) {
            mult1 = mult;
            unsigned d = x % dom_squared;
            mindex -= d * mult;
        }
        if (k2 == j) {
            mult2 = mult;
            unsigned d = x % dom_squared;
            mindex -= d * mult;
        }
        x /= dom_squared;
        mult *= dom_squared;
    }


    //
    // Ready to build the list. Clear it and then
    // loop over pairs (v1, v2) with
    //      v1 = (r1, r1) = r1 * dom() + r1,
    // and similarly for v2.
    //
    ilist.clear();
    for (unsigned r1=0; r1<dom(); r1++) {
        const unsigned v1 = r1 * dom() + r1;
        for (unsigned r2=0; r2<dom(); r2++) {
            const unsigned v2 = r2 * dom() + r2;
            ilist.push_back( mindex + mult1*v1 + mult2*v2 );
            if (0==mult2) break;
        }
        if (0==mult1) break;
    }
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
    const unsigned pairdom = isForSets() ? dom() : dom() * dom();
    for (unsigned j=1; j<=vars(); j++) {
        if (k1 == j) {
            mult1 = mult;
            unsigned d = x % pairdom;
            mindex -= d * mult;
        }
        if (k2 == j) {
            mult2 = mult;
            unsigned d = x % pairdom;
            mindex -= d * mult;
        }
        x /= pairdom;
        mult *= pairdom;
    }


    //
    // Ready to build the list. Clear it and then
    // loop over pairs (v1, v2) with v1 in {0, .., setdom()-1}
    // and v2 in {0, ..., setdom()-1}
    //
    ilist.clear();
    for (unsigned v1=0; v1<pairdom; v1++) {
        for (unsigned v2=0; v2<pairdom; v2++) {
            ilist.push_back( mindex + mult1*v1 + mult2*v2 );
            if (0==mult2) break;
        }
        if (0==mult1) break;
    }
}

