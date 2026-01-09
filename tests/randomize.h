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

class vectorgen {
    public:
        ///
        /// Constructor.
        ///     @param  sr      Set (false) or relation (true).
        ///     @param  vars    number of domain variables.
        ///     @param  dom     domain; size for each variable.
        ///
        /// The explicit representation will be a vector of
        /// dimension POTENTIAL = (rd*rd) ^ v
        /// so these need to be very small integers.
        ///
        vectorgen(bool sr, unsigned vars, unsigned dom);

        static void setSeed(long s, bool print=true);

        inline static long getSeed()
        {
            return seed;
        }

        static double random();

        static inline unsigned Equilikely_U(unsigned a, unsigned b)
        {
            return (a + (unsigned) ((b - a + 1) * random()));
        }

        static inline int Equilikely_I(int a, int b)
        {
            return (a + (int) ((b - a + 1) * random()));
        }

        inline unsigned potential() const {
            return POTENTIAL;
        }
        inline unsigned vars() const {
            return VARS;
        }
        inline unsigned dom() const {
            return DOM;
        }
        inline bool isForRelations() const {
            return is_for_relations;
        }
        inline bool isForSets() const {
            return !is_for_relations;
        }

        ///
        /// Build a domain based on the variables and domain
        ///
        MEDDLY::domain* makeDomain();

        ///
        /// Increment a minterm.
        ///
        inline bool nextMinterm(MEDDLY::minterm &m) const
        {
            if (m.isForSets()) {
                for (unsigned i=1; i<=m.getNumVars(); ++i)
                {
                    if (++m.from(i) < dom()) {
                        return true;
                    }
                    m.from(i) = 0;
                }
                return false;
            } else {
                for (unsigned i=1; i<=m.getNumVars(); ++i)
                {
                    if (++m.to(i) < dom()) {
                        return true;
                    }
                    m.to(i) = 0;
                    if (++m.from(i) < dom()) {
                        return true;
                    }
                    m.from(i) = 0;
                }
                return false;
            }
        }

        ///
        /// Generate a random minterm.
        /// There may be don't care or don't change values.
        ///
        void randomizeMinterm(MEDDLY::minterm &m, MEDDLY::range_type rt);


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
        void index2minterm(unsigned x, MEDDLY::minterm &m) const;

        ///
        /// Converts a pair of indexes (from, to) to a relation minterm.
        ///
        ///     @param  fr  Index of source state.
        ///                 Only fr % POTENTIAL is significant;
        ///                 any excess is ignored.
        ///
        ///     @param  to  Index of destination state.
        ///                 Only to % POTENTIAL is significant;
        ///                 any excess is ignored.
        ///
        ///     @param  m   Minterm to fill.
        ///                 Should have VARS variables and be for relations.
        ///
        void indexes2minterm(unsigned fr, unsigned to, MEDDLY::minterm &m) const;


        ///
        /// Calculates the index of a minterm.
        /// Essentially, the reverse of index2minterm().
        ///
        ///     @param  m   Minterm in question.
        ///                 Should have VARS variables.
        ///
        unsigned minterm2index(const MEDDLY::minterm &m) const;

        ///
        /// Calculates a pair of indexes (from, to) from a minterm.
        /// Essentially, the reverse of indexes2minterm().
        ///
        ///     @param  m   Minterm in question.
        ///                 Should have VARS variables.
        ///
        ///     @param  fr  On return: Index of source state.
        ///
        ///     @param  to  On return: Index of destination state.
        ///
        void minterm2indexes(const MEDDLY::minterm &m, unsigned &fr, unsigned &to) const;

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

        ///
        /// Randomly generate an explicit vector.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-default elements
        ///     @param  values  Values to sample from. Element 0 is the
        ///                     default value; all others are for non-default.
        ///
        void randomizeVector(std::vector <MEDDLY::rangeval> &elems,
                unsigned card, const std::vector <MEDDLY::rangeval> &values);

        ///
        /// Randomly generate an explicit vector, with fully patterns.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-default fully patterns.
        ///                     Because patterns might overlap,
        ///                     and the size of each pattern is more
        ///                     than one, this will not be the
        ///                     cardinality of the final set/vector.
        ///     @param  values  Values to sample from. Element 0 is the
        ///                     default value; all others are for non-default.
        ///
        void randomizeFully(std::vector <MEDDLY::rangeval> &elems,
                unsigned card, const std::vector <MEDDLY::rangeval> &values);

        ///
        /// Randomly generate an explicit vector, with identity patterns.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-default identity patterns.
        ///                     Because patterns might overlap,
        ///                     and the size of each pattern is more
        ///                     than one, this will not be the
        ///                     cardinality of the final set/vector.
        ///     @param  values  Values to sample from. Element 0 is the
        ///                     default value; all others are for non-default.
        ///
        void randomizeIdentity(std::vector <MEDDLY::rangeval> &elems,
                unsigned card, const std::vector <MEDDLY::rangeval> &values);

        //
        // Build a dd_edge from an explicit vector, using MAXIMUM
        // (use deflt as the smallest element in x).
        //
        void explicit2edgeMax(const std::vector<MEDDLY::rangeval> &x,
            MEDDLY::dd_edge &s, MEDDLY::rangeval deflt) const;

        //
        // Build a dd_edge from an explicit vector, using MINIMUM
        // (use deflt as the largest element in x).
        //
        void explicit2edgeMin(const std::vector<MEDDLY::rangeval> &x,
            MEDDLY::dd_edge &s, MEDDLY::rangeval deflt) const;

        //
        // Display an explicit vector, as a set
        //
        void showSet(std::ostream &out,
                const std::vector <MEDDLY::rangeval> &elems) const;

        //
        // Display minterms corresponding to an explicit vector
        //
        void showMinterms(std::ostream &out, const MEDDLY::domain* D,
                const std::vector <MEDDLY::rangeval> &elems) const;

    private:
        const bool is_for_relations;
        const unsigned VARS;
        const unsigned DOM;
        unsigned POTENTIAL;
        std::vector <unsigned> ilist;

        static long seed;

        unsigned current_terminal;
};

#endif
