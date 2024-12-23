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
        ///     @param  sr      Set (false) or relation (true).
        ///     @param  vars    number of domain variables.
        ///     @param  dom     domain; size for each variable.
        ///     @param  range   number of distinct range values.
        ///
        /// The explicit representation will be a vector of
        /// dimension POTENTIAL = (rd*rd) ^ v
        /// so these need to be very small integers.
        ///
        vectorgen_base(bool sr, unsigned vars, unsigned dom, unsigned range);

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
        inline unsigned range() const {
            return RANGE;
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

        static inline void vno2val(unsigned v, bool &val)
        {
            val = true;
        }
        static inline void vno2val(unsigned v, int &val)
        {
            val = int(v);
        }
        static inline void vno2val(unsigned v, long &val)
        {
            val = long(v);
        }
        static inline void vno2val(unsigned v, float &val)
        {
            val = v / 2.0;
        }
        static inline void vno2val(unsigned v, double &val)
        {
            val = v / 2.0;
        }

    protected:
        std::vector <unsigned> ilist;
        const MEDDLY::domain* _D;

    private:
        const bool is_for_relations;
        const unsigned VARS;
        const unsigned DOM;
        const unsigned RANGE;
        unsigned POTENTIAL;

        static long seed;

        unsigned current_terminal;
};


template <typename TYPE>
class vectorgen : public vectorgen_base {
    public:
        vectorgen(bool sr, unsigned vars, unsigned dom, unsigned range)
            : vectorgen_base(sr, vars, dom, range) { }

        ///
        /// Randomly generate an explicit vector.
        ///     @param  elems   Vector of elements; will be cleared.
        ///     @param  card    Number of non-zero elements
        ///
        void randomizeVector(std::vector <TYPE> &elems, unsigned card)
        {
            if (elems.size() != potential()) {
                throw "vector size mismatch in randomizeVector";
            }
            unsigned v=0;
            for (unsigned i=0; i<elems.size(); i++) {
                if (i<card) {
                    v = 1+ (v % range());
                    TYPE val;
                    vno2val(v, val);
                    elems[i] = val;
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
        ///
        void randomizeFully(std::vector <TYPE> &elems, unsigned card)
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
                v = 1+ (v % range());
                TYPE val;
                vno2val(v, val);

                ilist.clear();
                buildFullyFromIndex(x, k1, k2, ilist);
                for (unsigned z=0; z<ilist.size(); z++) {
                    elems[ ilist[z] ] = val;
                }
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
        ///
        void randomizeIdentity(std::vector <TYPE> &elems, unsigned card)
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
                v = 1+ (v % range());
                TYPE val;
                vno2val(v, val);

                ilist.clear();
                buildIdentityFromIndex(x, k1, k2, ilist);
                for (unsigned z=0; z<ilist.size(); z++) {
                    elems[ ilist[z] ] = val;
                }
            }
        }

        //
        // Build a dd_edge from an explicit vector
        //
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
                F->createEdge(zero, s);
                return;
            }

            MEDDLY::minterm_coll mtlist(card, F);
            for (unsigned i=0; i<x.size(); i++) {
                if (!x[i]) continue;
                index2minterm(i, mtlist.unused());
                TYPE val = x[i];
                mtlist.unused().setTerm(val);
                mtlist.pushUnused();
            }
            mtlist.buildFunction(s);
        }

        //
        // Display an explicit vector
        //
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
        static void showElem(std::ostream &out, TYPE elem);
};

template <>
inline void vectorgen<bool>::showElem(std::ostream &out, bool elem)
{
}

template <typename TYPE>
inline void vectorgen<TYPE>::showElem(std::ostream &out, TYPE elem)
{
    out << ":" << elem;
}

#endif
