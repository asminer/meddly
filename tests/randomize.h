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

    public:
        // template member functions

        ///
        /// Randomly generate an explicit vector.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-zero elements
        ///     @param  values  Values to sample from
        ///
        template <typename TYPE>
        void randomizeVector(std::vector <TYPE> &elems, unsigned card,
                const std::vector <TYPE> &values)
        {
            if (elems.size() != potential()) {
                throw "vector size mismatch in randomizeVector";
            }
            unsigned v=0;
            for (unsigned i=0; i<elems.size(); i++) {
                if (i<card) {
                    elems[i] = values[v];
                    v = (v+1) % values.size();
                } else {
                    elems[i] = 0;
                }
            }
            //
            // Shuffle the array
            //
            for (unsigned i=0; i<elems.size()-1; i++) {
                unsigned j = Equilikely_U(i, elems.size()-1);
                if (elems[i] != elems[j]) {
                    TYPE t = elems[i];
                    elems[i] = elems[j];
                    elems[j] = t;
                }
            }
        }

        ///
        /// Randomly generate an explicit vector, with fully patterns.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-zero fully patterns.
        ///                     Because patterns might overlap,
        ///                     and the size of each pattern is more
        ///                     than one, this will not be the
        ///                     cardinality of the final set/vector.
        ///     @param  values  Values to sample from
        ///
        template <typename TYPE>
        void randomizeFully(std::vector <TYPE> &elems, unsigned card,
                const std::vector <TYPE> &values)
        {
            if (elems.size() != potential()) {
                throw "vector size mismatch in randomizeFully";
            }
            for (unsigned i=0; i<elems.size(); i++) {
                elems[i] = 0;
            }
            unsigned v=0;
            for (unsigned i=0; i<card; i++) {
                unsigned x = Equilikely_U(0, elems.size()-1);
                unsigned k1 = Equilikely_U(1, vars());
                unsigned k2 = Equilikely_U(1, vars());

                ilist.clear();
                buildFullyFromIndex(x, k1, k2, ilist);
                for (unsigned z=0; z<ilist.size(); z++) {
                    elems[ ilist[z] ] = values[v];
                }
                v = (v+1) % values.size();
            }
        }

        ///
        /// Randomly generate an explicit vector, with identity patterns.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-zero identity patterns.
        ///                     Because patterns might overlap,
        ///                     and the size of each pattern is more
        ///                     than one, this will not be the
        ///                     cardinality of the final set/vector.
        ///     @param  values  Values to sample from
        ///
        template <typename TYPE>
        void randomizeIdentity(std::vector <TYPE> &elems, unsigned card,
                const std::vector <TYPE> &values)
        {
            if (elems.size() != potential()) {
                throw "vector size mismatch in randomizeIdentity";
            }
            for (unsigned i=0; i<elems.size(); i++) {
                elems[i] = 0;
            }
            unsigned v=0;
            for (unsigned i=0; i<card; i++) {
                unsigned x = Equilikely_U(0, elems.size()-1);
                unsigned k1 = Equilikely_U(1, vars());
                unsigned k2 = Equilikely_U(1, vars());

                ilist.clear();
                buildIdentityFromIndex(x, k1, k2, ilist);
                for (unsigned z=0; z<ilist.size(); z++) {
                    elems[ ilist[z] ] = values[v];
                }

                v = (1+v) % values.size();
            }
        }

        //
        // Build a dd_edge from an explicit vector
        //
        template <typename TYPE>
        void explicit2edge(const std::vector<TYPE> &x, MEDDLY::dd_edge &s)
            const
        {
            MEDDLY::forest* F = s.getForest();
            if (!F) throw "null forest in explicit2edge";
            // Determine x 'cardinality'
            unsigned card = 0;
            for (unsigned i=0; i<x.size(); i++) {
                if (x[i]) ++card;
            }

            if (0==card) {
                TYPE zero = 0;
                F->createConstant(zero, s);
                return;
            }

            MEDDLY::minterm_coll mtlist(card, F);
            for (unsigned i=0; i<x.size(); i++) {
                if (!x[i]) continue;
                index2minterm(i, mtlist.unused());
                TYPE val = x[i];
                mtlist.unused().setValue(val);
                mtlist.pushUnused();
            }
            mtlist.buildFunctionMax(TYPE(0), s);
        }

        //
        // Display an explicit vector
        //
        template <typename TYPE>
        void showSet(std::ostream &out, const std::vector <TYPE> &elems)
            const
        {
            out << "{ ";
            bool printed = false;
            for (unsigned i=0; i<elems.size(); i++) {
                if (!elems[i]) continue;
                if (printed) out << ", ";
                out << i;
                showElem(out, elems[i]);
                printed = true;
            }
            out << " }";
        }

        //
        // Display minterms corresponding to an explicit vector
        //
        template <typename TYPE>
        void showMinterms(std::ostream &out, const std::vector <TYPE> &elems)
            const
        {
            if (!_D) throw "null domain, showMinterms";
            MEDDLY::minterm mt(_D, isForRelations());
            MEDDLY::ostream_output mout(out);
            mout << "{ ";
            bool printed = false;
            for (unsigned i=0; i<elems.size(); i++) {
                if (!elems[i]) continue;

                if (printed) mout << ",\n      ";
                printed = true;

                index2minterm(i, mt);
                mt.show(mout);
            }
            mout << " }";

        }

    protected:
        template <typename TYPE>
        static void showElem(std::ostream &out, TYPE elem);

    private:
        const bool is_for_relations;
        const unsigned VARS;
        const unsigned DOM;
        unsigned POTENTIAL;
        std::vector <unsigned> ilist;
        const MEDDLY::domain* _D;

        static long seed;

        unsigned current_terminal;
};

template <>
inline void vectorgen::showElem(std::ostream &out, bool elem)
{
}

template <typename TYPE>
inline void vectorgen::showElem(std::ostream &out, TYPE elem)
{
    out << ":" << elem;
}

#endif
