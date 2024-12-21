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

#ifndef RANDOMIZE_H
#define RANDOMIZE_H

#include "../src/meddly.h"
#include <vector>

class vectorgen_base {
    public:
        ///
        /// Constructor.
        ///     @param  v   number of domain variables.
        ///     @param  rd  variable domain for relations;
        ///                 square this for sets.
        ///
        /// The explicit representation will be a vector of
        /// dimension POTENTIAL = (rd*rd) ^ v
        /// so these need to be very small integers.
        ///
        vectorgen_base(unsigned v=5, unsigned rd=2);

        static void setSeed(long s, bool print=true);

        static double random();

        static inline unsigned Equilikely(unsigned a, unsigned b)
        {
            return (a + (unsigned) ((b - a + 1) * random()));
        }

        inline unsigned potential() const {
            return POTENTIAL;
        }
        inline unsigned vars() const {
            return VARS;
        }
        inline unsigned reldom() const {
            return RELDOM;
        }
        inline unsigned setdom() const {
            return RELDOM * RELDOM;
        }


        ///
        /// Converts an index into the corresponding minterm.
        ///
        ///     @param  x   Index to convert.
        ///                 Only x % POTENTIAL is significant;
        ///                 any excess is ignored.
        ///
        ///     @param  m   Minterm to fill.
        ///                 Should have VARS variables.
        ///
        void index2minterm(unsigned x, MEDDLY::minterm &m);

        ///
        /// Converts ndx into an array of digits, and then loops
        /// digits corresponding to k1 and k1', and k2 and k2',
        /// such that digit k1 = digit k1', and digit k2 = digit k2'.
        /// The resulting digits are converted back to indexes,
        /// and added to a list.
        ///
        ///     @param  ndx     Index to start from.
        ///     @param  k1      Some variable
        ///     @param  k2      Another variable; can equal k1
        ///     @param  ilist   List of indexes, will be cleared.
        ///                     On output will contain ndx and
        ///                     any other indexes obtained.
        ///
        void buildIdentityFromIndex(unsigned ndx, unsigned k1, unsigned k2,
                std::vector <unsigned> &ilist) const;

        ///
        /// Converts ndx into an array of digits, and then loops
        /// digits k1 and k2 through all possible values. The resulting
        /// digits are converted back to indexes, and added to a list.
        ///
        ///     @param  ndx     Index to start from.
        ///     @param  k1      Some variable
        ///     @param  k2      Another variable; can equal k1
        ///     @param  ilist   List of indexes, will be cleared.
        ///                     On output will contain ndx and
        ///                     any other indexes obtained.
        ///
        void buildFullyFromIndex(unsigned ndx, unsigned k1, unsigned k2,
                std::vector <unsigned> &ilist) const;

    protected:
        static inline void vno2val(unsigned v, bool &val)
        {
            val = true;
        }
        static inline void vno2val(unsigned v, int &val)
        {
            val = int(v);
        }
        static inline void vno2val(unsigned v, float &val)
        {
            val = v / 2.0;
        }

    private:
        const unsigned VARS;
        const unsigned RELDOM;
        unsigned POTENTIAL;

        static long seed;
};


template <typename TYPE>
class vectorgen : public vectorgen_base {
    public:
        vectorgen(unsigned v=5, unsigned rd=2) : vectorgen_base(v, rd)
        {}

        ///
        /// Randomly generate an explicit vector.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-zero elements
        ///     @param  nv      Number of distinct values to use
        ///                     for the non-zero elements.
        ///                     Will be evenly distributed.
        ///
        void randomizeVector(std::vector <TYPE> &elems, unsigned card,
                unsigned nv) const
        {
        }

        ///
        /// Randomly generate an explicit vector, with fully patterns.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-zero fully patterns.
        ///                     Because patterns might overlap,
        ///                     and the size of each pattern is more
        ///                     than one, this will not be the
        ///                     cardinality of the final set/vector.
        ///     @param  nv      Number of distinct values to use
        ///                     for the non-zero elements.
        ///                     Will be evenly distributed.
        ///
        void randomizeFully(std::vector <TYPE> &elems, unsigned card,
                unsigned nv) const
        {
        }

        ///
        /// Randomly generate an explicit vector, with identity patterns.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-zero identity patterns.
        ///                     Because patterns might overlap,
        ///                     and the size of each pattern is more
        ///                     than one, this will not be the
        ///                     cardinality of the final set/vector.
        ///     @param  nv      Number of distinct values to use
        ///                     for the non-zero elements.
        ///                     Will be evenly distributed.
        ///
        void randomizeIdentity(std::vector <TYPE> &elems, unsigned card,
                unsigned nv) const
        {
        }


};

#endif
