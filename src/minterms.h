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

#ifndef MEDDLY_MINTERMS_H
#define MEDDLY_MINTERMS_H

//
// Idea: define minterm objects, for both full and sparse storage,
// and for set vs relation.
//
// For now: defines that should eventually be part of minterm objects
//

//
// A minterm should be a std::vector of ints, plus a rangeval.
//
// Should we define collections of minterms here?
//

#include "defines.h"
#include "minterms.h"
#include "policies.h"
#include "terminal.h"
#include "dd_edge.h"

namespace MEDDLY {

    /** Special value for minterms: don't care what this variable does.
      I.e., do the same thing for all possible assignments for a variable.
      */
    const int DONT_CARE  = -1;

    /** Special value for primed minterms: don't change this variable.
      Forces the primed value to equal the unprimed value for a variable.
      Undefined for unprimed variables.
      */
    const int DONT_CHANGE = -2;

    class domain;
    class output;

    class minterm;
    class minterm_coll;

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                         minterm  class                         *
// *                                                                *
// ******************************************************************


/**
    Minterm object.
    We assign values to each variable, or we allow variables
    to be unset using a value of DONT_CARE, or in a relation
    we can restruct variables to be unchanged using a value
    of DONT_CARE.

    Unlike a proper minterm, however, we can set the
    function value to something other than true for the
    given set of variable assignments.
*/
class MEDDLY::minterm {
    public:
        minterm(const domain* D, set_or_rel sr);
        ~minterm();

        inline bool isForSets()             const   { return !for_relations; }
        inline bool isForRelations()        const   { return for_relations; }

        inline unsigned getNumVars()        const   { return num_vars; }

        inline const domain* getDomain()    const   { return _D; }

        /// Set the terminal value
        inline void setTerm(const terminal &t)      { termval = t; }

        /// Get the terminal value
        inline const terminal& getTerm()    const   { return termval; }

        // set access

        inline int getVar(unsigned i) const {
            MEDDLY_DCASSERT(isForSets());
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            MEDDLY_DCASSERT(_from);
            return _from[i];
        }
        inline void setVar(unsigned i, int from) {
            MEDDLY_DCASSERT(isForSets());
            MEDDLY_DCASSERT(from >= 0 || from == DONT_CARE);
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            MEDDLY_DCASSERT(_from);
            _from[i] = from;
        }

        // relation access

        inline void getVars(unsigned i, int &from, int &to) const {
            MEDDLY_DCASSERT(isForRelations());
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            MEDDLY_DCASSERT(_from);
            MEDDLY_DCASSERT(_to);
            from = _from[i];
            to   = _to[i];
        }
        inline void setVars(unsigned i, int from, int to) {
            MEDDLY_DCASSERT(isForRelations());
            MEDDLY_DCASSERT(from >= 0 || from == DONT_CARE);
            MEDDLY_DCASSERT(to >= 0 || to == DONT_CARE || to == DONT_CHANGE);
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            MEDDLY_DCASSERT(_from);
            MEDDLY_DCASSERT(_to);
            if (DONT_CHANGE == to) {
                if (from >= 0) {
                    to = from;
                }
            }
            _from[i] = from;
            _to[i] = to;
        }

        // both access

        inline int from(unsigned i) const {
            MEDDLY_DCASSERT(_from);
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            return _from[i];
        }
        inline int& from(unsigned i) {
            MEDDLY_DCASSERT(_from);
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            return _from[i];
        }
        inline int to(unsigned i) const {
            MEDDLY_DCASSERT(_to);
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            return _to[i];
        }
        inline int& to(unsigned i) {
            MEDDLY_DCASSERT(_to);
            CHECK_RANGE(__FILE__, __LINE__, 1u, i, num_vars+1);
            return _to[i];
        }

        inline int var(int i) const {
            return (i>0) ? from(i) : to(-i);
        }

        //
        // For convenience and debugging
        void show(output &s) const;

    private:
        const domain* _D;

        int* _from;
        int* _to;

        terminal termval;

        unsigned num_vars;
        bool for_relations;
};

// ******************************************************************
// *                                                                *
// *                       minterm_coll class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::minterm_coll {
    public:
        /// Build a collection for a given domain.
        ///     @param  maxsize Maximum size of the collection
        ///     @param  D       The domain
        ///     @param  sr      Set or relation
        minterm_coll(unsigned maxsize,
                const domain* D, set_or_rel sr);

        ~minterm_coll();

        inline bool isForRelations()        const   { return for_relations; }
        inline unsigned getNumVars()        const   { return num_vars; }
        inline const domain* getDomain()    const   { return _D; }

        /// Current collection size
        inline unsigned size() const { return first_unused; }

        /// Collection max size
        inline unsigned maxsize() const { return max_coll_size; }

        /// Clear the collection
        inline void clear() { first_unused = 0; }

        /// Get an element
        inline minterm& at(unsigned i) {
            CHECK_RANGE(__FILE__, __LINE__, 0u, i, 1+first_unused);
            MEDDLY_DCASSERT(_mtlist);
            MEDDLY_DCASSERT(_mtlist[i]);
            return *(_mtlist[i]);
        }

        /// Get minterm i
        inline const minterm& at(unsigned i) const {
            CHECK_RANGE(__FILE__, __LINE__, 0u, i, 1+first_unused);
            MEDDLY_DCASSERT(_mtlist);
            MEDDLY_DCASSERT(_mtlist[i]);
            return *(_mtlist[i]);
        }

        /// Get x_k (unprimed) from minterm i
        inline int unprimed(unsigned i, int k) const {
            CHECK_RANGE(__FILE__, __LINE__, 1, k, int(num_vars)+1);
            MEDDLY_DCASSERT(_mtlist);
            MEDDLY_DCASSERT(_mtlist[i]);
            return _mtlist[i]->from(k);
        }

        /// Get x'_k (primed) from minterm i
        inline int primed(unsigned i, int k) const {
            CHECK_RANGE(__FILE__, __LINE__, 1, k, int(num_vars)+1);
            MEDDLY_DCASSERT(_mtlist);
            MEDDLY_DCASSERT(_mtlist[i]);
            return _mtlist[i]->to(k);
        }

        /// Get the first unused element
        inline minterm& unused() {
            return at(first_unused);
        }

        /// Add the unused element to the collection.
        /// Returns true on success.
        inline bool pushUnused() {
            if (first_unused < maxsize()) {
                ++first_unused;
                return true;
            } else {
                return false;
            }
        }

        /// Swap minterms i and j
        inline void swap(unsigned i, unsigned j)
        {
            CHECK_RANGE(__FILE__, __LINE__, 0u, i, first_unused);
            CHECK_RANGE(__FILE__, __LINE__, 0u, j, first_unused);
            MEDDLY_DCASSERT(_mtlist);
            MEDDLY_DCASSERT(i != j);
            SWAP(_mtlist[i], _mtlist[j]);
        }


        /** Create an edge from this collection.
            The default value is 0.
            If the collection contains overlapping minterms
            for a particular variable assignment,
            the maximum value is taken.

                @param  e   On input: should be attached to the
                            forest we want to create the function in.
                            On output: the function.
        */
        void buildFunction(dd_edge &e);


        //
        // For convenience and debugging
        void show(output &s, const char* pre, const char* post) const;

        //
        // Used by fbuilder.
        // TBD: move there?
        //


#ifdef ALLOW_MINTERM_OPS
        /** Partial sort, as used by explicit minterm operations.
            Based on level L, and the range of minterms [low, hi),
            collect all elements together with variable L equal to item #low.
            on output, lend is such that [low, lend) has equal elements
            where low < lend <= hi

            As a special case, if L>0 and the collection is for relations,
            if the first item has value "DONT_CARE", then we ALSO look at
            the primed value but only to distinguish "DONT_CHANGE" and
            not equal to "DONT_CHANGE".

                @param  L       Variable number; use negative for primed
                @param  low     First minterm to check
                @param  hi      One past the last minterm to check
                @param  lend    On output, one past the last minterm
                                in the group equal to low.

                @return The first value
        */
        inline int collect_first(int L, unsigned low, unsigned hi,
                unsigned& lend)
        {
            if (hi-low < 2) {
                // Only one element, no checking or reordering required.
                lend = hi;
                return at(low).var(L);
            }
            return _collect_first(L, low, hi, lend);
        }
#endif


    private:


#ifdef ALLOW_MINTERM_OPS
    private:
        int _collect_first(int L, unsigned low, unsigned hi, unsigned& lend);
#endif


        // OLD STUFF
    private:

        /**
            Sort minterms, as needed for explicit minterm operations.
            Based on level L, and the range of minterms [low, hi),
            sort the minterms and store the interval points in
            vector index such that the interval [ index[i], index[i+1] )
            gives the range of indexes with equal values, after the sort.
            If variable L has the same value for the whole interval
            [low, high), then index will have size 2 with
            index[0] = low and index[1] = high.
            Actually, no matter what, we will have index[0] = low
            and index[size()-1] = high.

                @param  L       Variable number
                @param  low     Low end of range
                @param  high    High end of range
                @param  index   On output: sub-intervals.
                @param  value   On output: value for each sub-interval.
        */
        /*
        void sortOnVariable(int L, unsigned low, unsigned high,
                std::vector<unsigned> &index, std::vector<int> &value);
                */


    private:
        minterm** _mtlist;
        const unsigned max_coll_size;

        unsigned first_unused;

        const domain* _D;

        unsigned num_vars;
        bool for_relations;
};

#endif
