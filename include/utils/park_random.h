/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#ifndef PARK_RANDOM_H
#define PARK_RANDOM_H

long seed = -1;

inline void SeedRNG(long s)
{
    seed = s;
}

inline double Random()
{
    static const long MODULUS = 2147483647L;
    static const long MULTIPLIER = 48271L;
    static const long Q = MODULUS / MULTIPLIER;
    static const long R = MODULUS % MULTIPLIER;

    long t = MULTIPLIER * (seed % Q) - R * (seed / Q);
    if (t > 0) {
        seed = t;
    } else {
        seed = t + MODULUS;
    }
    return ((double) seed / MODULUS);
}

inline unsigned Equilikely(unsigned a, unsigned b)
{
    return (a + (unsigned) ((b - a + 1) * Random()));
}

#endif

