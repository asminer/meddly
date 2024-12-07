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
        minterm(const forest* F);
        ~minterm();

        inline bool isForSets()             const   { return !for_relations; }
        inline bool isForRelations()        const   { return for_relations; }

        inline unsigned getNumVars()        const   { return num_vars; }

        inline const domain* getDomain()    const   { return _D; }

        /// Set the terminal value
        inline void setTerm(const terminal &t)      { termval = t; }

        /// Get the terminal value
        inline const terminal& getTerm()    const   { return termval; }

        /// Get the memory usage, in bytes
        inline long memoryUsage()           const
        {
            const long me = sizeof(minterm);
            const long arr = sizeof(int) * (1+num_vars);
            return _to  ? ( me + 2* arr) : ( me + arr);
        }

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
        // Set from an array.
        //
        inline void setAll(const int* v, const terminal &t)
        {
            MEDDLY_DCASSERT(isForSets());
            MEDDLY_DCASSERT(v);
            MEDDLY_DCASSERT(_from);
            for (unsigned i=num_vars; i; --i) {
                _from[i] = v[i];
            }
            termval = t;
        }

        //
        // Set from unprimed and primed arrays.
        //
        inline void setAll(const int* un, const int* pr, const terminal &t)
        {
            MEDDLY_DCASSERT(isForRelations());
            MEDDLY_DCASSERT(un);
            MEDDLY_DCASSERT(pr);
            MEDDLY_DCASSERT(_from);
            MEDDLY_DCASSERT(_to);
            for (unsigned i=num_vars; i; --i) {
                _from[i] = un[i];
                _to[i] = pr[i];
            }
            termval = t;
        }

        //
        // For convenience and debugging
        void show(output &s) const;

        /** Create an edge from this minterm.
            The default value is 0.

                @param  e   On input: should be attached to the
                            forest we want to create the function in.
                            On output: the function.
        */
        void buildFunction(dd_edge &e) const;

    private:
        /// PRE:    _D, termval, and for_relations are set.
        ///
        /// POST:   _D, termval, and for_relations are unchanged;
        ///         num_vars, _from, and _to are set.
        ///
        void initVectors();

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

        /// Build a collection for a given forest.
        /// The collection will work for any forest with the same
        /// domain, and with the same dimension (i.e., sets
        /// vs. relations).
        ///     @param  maxsize Maximum size of the collection
        ///     @param  F       The forest
        minterm_coll(unsigned maxsize, const forest* F);


        ~minterm_coll();

        inline bool isForRelations()        const   { return for_relations; }
        inline unsigned getNumVars()        const   { return num_vars; }
        inline const domain* getDomain()    const   { return _D; }

        /// Memory used
        inline long memoryUsage()           const
        {
            long bytes = sizeof(minterm_coll);
            if (max_coll_size) {
                MEDDLY_DCASSERT(_mtlist);
                MEDDLY_DCASSERT(_mtlist[0]);
                bytes += max_coll_size * _mtlist[0]->memoryUsage();
            }
            return bytes;
        }

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
        inline void pushUnused() {
            MEDDLY_DCASSERT(first_unused < maxsize());
            ++first_unused;
        }

        /// Expand the collection if necessary.
        /// Call this right after pushUnused(), if you want to
        /// allow expansions.
        inline void expandIfNecessary() {
            if (first_unused < maxsize()) return;
            _expand();
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
            The default value is 0 except for EVPLUS which is infinity.

                @param  e   On input: should be attached to the
                            forest we want to create the function in.
                            On output: the function.

                @param  minimize    If there are overlapping minterms in the
                                    collection, then we will take the minimum
                                    value if this paramter is true; otherwise
                                    we take the maximum value.
        */
        void buildFunction(dd_edge &e, bool minimize=false);


        //
        // For convenience and debugging: show all minterms.
        void show(output &s, const char* pre=nullptr, const char* post="\n")
            const;

    protected:
        void _expand();

    private:
        minterm** _mtlist;
        unsigned max_coll_size;

        unsigned first_unused;

        const domain* _D;

        unsigned num_vars;
        bool for_relations;
};

#endif
