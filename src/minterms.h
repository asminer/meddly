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

#include <vector>

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

    Storage-wise, a minterm may be stored in sparse or full format.
    In sparse format, any unspecified variable is assumed
    to be DONT_CARE (any possible value).
    In both formats, for relations we must set the unprimed
    and primed variables together.

    Unlike a proper minterm, however, we can set the
    function value to something other than true for the
    given set of variable assignments.
*/
class MEDDLY::minterm {
    public:
        minterm(const domain* D, set_or_rel sr, node_storage_flags fs);

        inline bool isSparse()              const   { return sparse; }
        inline bool isFull()                const   { return !sparse; }

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
#ifdef DEVELOPMENT_CODE
            return _from.at(i);
#else
            return _from[i];
#endif
        }
        inline void setVar(unsigned i, int from) {
            MEDDLY_DCASSERT(isForSets());
            MEDDLY_DCASSERT(from >= 0 || from == DONT_CARE);
#ifdef DEVELOPMENT_CODE
            _from.at(i) = from;
#else
            _from[i] = from;
#endif
        }

        // relation access

        inline void getVars(unsigned i, int &from, int &to) const {
            MEDDLY_DCASSERT(isForRelations());
#ifdef DEVELOPMENT_CODE
            from = _from.at(i);
            to   = _to.at(i);
#else
            from = _from[i];
            to   = _to[i];
#endif
        }
        inline void setVars(unsigned i, int from, int to) {
            MEDDLY_DCASSERT(isForRelations());
            MEDDLY_DCASSERT(from >= 0 || from == DONT_CARE);
            MEDDLY_DCASSERT(to >= 0 || to == DONT_CARE || to == DONT_CHANGE);
            if (DONT_CHANGE == to) {
                if (from >= 0) {
                    to = from;
                }
            }
#ifdef DEVELOPMENT_CODE
            _from.at(i) = from;
            _to.at(i) = to;
#else
            _from[i] = from;
            _to[i] = to;
#endif
        }

        // both access

        inline int var(int i) const {
            MEDDLY_DCASSERT(i>0 || isForRelations());
#ifdef DEVELOPMENT_CODE
            return (i>0) ? _from.at(i) : _to.at(-i);
#else
            return (i>0) ? _from[i] : _to[-i];
#endif
        }

        inline int from(unsigned i) const {
#ifdef DEVELOPMENT_CODE
            return _from.at(i);
#else
            return _from[i];
#endif
        }
        inline int& from(unsigned i) {
#ifdef DEVELOPMENT_CODE
            return _from.at(i);
#else
            return _from[i];
#endif
        }
        inline int to(unsigned i) const {
            MEDDLY_DCASSERT(isForRelations());
#ifdef DEVELOPMENT_CODE
            return _to.at(i);
#else
            return _to[i];
#endif
        }
        inline int& to(unsigned i) {
            MEDDLY_DCASSERT(isForRelations());
#ifdef DEVELOPMENT_CODE
            return _to.at(i);
#else
            return _to[i];
#endif
        }

        // sparse access

        inline void append(unsigned ndx, int from) {
            if (DONT_CARE == from) return;
            MEDDLY_DCASSERT(from >= 0);
            MEDDLY_DCASSERT(isForSets());
            MEDDLY_DCASSERT(isSparse());
            MEDDLY_DCASSERT( _index.empty() || _index.back() < ndx );
            _index.push_back(ndx);
            _from.push_back(from);

        }

        inline void append(unsigned ndx, int from, int to) {
            if ((DONT_CARE == from) && (DONT_CARE == to)) return;
            MEDDLY_DCASSERT(from >= 0 || from == DONT_CARE);
            MEDDLY_DCASSERT(to >= 0 || to == DONT_CARE || to == DONT_CHANGE);
            MEDDLY_DCASSERT(isForSets());
            MEDDLY_DCASSERT(isSparse());
            MEDDLY_DCASSERT( _index.empty() || _index.back() < ndx );
            if (DONT_CHANGE == to) {
                if (from >= 0) {
                    to = from;
                }
            }
            _index.push_back(ndx);
            _from.push_back(from);
            _to.push_back(to);
        }

        inline unsigned size() const {
            MEDDLY_DCASSERT(isSparse());
            return _index.size();
        }

        inline int whichVar(unsigned nz) const {
            MEDDLY_DCASSERT(isSparse());
#ifdef DEVELOPMENT_CODE
            return _index.at(nz);
#else
            return _index[nz];
#endif
        }

        //
        // Check that we are sparse, and have the same variables
        // set as the given list L.
        //
        bool hasSparseVars(const std::vector <unsigned> &vars) const;

        //
        // For convenience and debugging
        void show(output &s) const;

    private:
        const domain* _D;

        std::vector<int> _from;
        std::vector<int> _to;
        std::vector<unsigned> _index;

        terminal termval;

        unsigned num_vars;
        bool for_relations;
        bool sparse;

    private:
        friend class MEDDLY::minterm_coll;
        minterm* next;
};

// ******************************************************************
// *                                                                *
// *                       minterm_coll class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::minterm_coll {
    public:
        /// Build a collection for a given domain.
        /// Each minterm is stored in full.
        ///     @param  D   The domain
        ///     @param  sr  Set or relation
        minterm_coll(const domain* D, set_or_rel sr);

        /// Build a collection of sparse minterms,
        /// for a given domain and with the same set
        /// of variables used.
        ///     @param  D       The domain
        ///     @param  sr      Set or relation
        ///     @param  varlist Array list of variables used,
        ///                     in decreasing order.
        minterm_coll(const domain* D, set_or_rel sr,
                const std::vector<unsigned> &varlist);

        inline bool isForRelations()        const   { return for_relations; }
        inline unsigned getNumVars()        const   { return num_vars; }
        inline const domain* getDomain()    const   { return _D; }

        /// Make a minterm for someone to build
        inline minterm* makeTemp()
        {
            if (freelist) {
                minterm* p = freelist;
                freelist = freelist->next;
                return p;
            }
            return new minterm(_D, for_relations,
                    sparse ? SPARSE_ONLY : FULL_ONLY);
        }

        /// Add a minterm to the collection
        inline void addToCollection(minterm* m)
        {
            if (m) {
                MEDDLY_DCASSERT(_check(m));
                _mtlist.push_back(m);
                sorted = false;
            }
        }

        /// Recycle a minterm
        inline void recycle(minterm* m)
        {
            if (m) {
                MEDDLY_DCASSERT(m->getDomain() == _D);
                m->next = freelist;
                freelist = m;
            }
        }


        /// Sort the collection into a good order
        /// for building a MDD / MXD
        inline void sort()
        {
            if (sorted) return;
            if (_mtlist.size()) {
                _sort(num_vars, false, 0, _mtlist.size());
            }
            sorted = true;
        }

        /// Clear the collection
        void clear();


        /// Current collection size
        inline unsigned size() const { return _mtlist.size(); }

        /// Get an element
        inline const minterm* at(unsigned i) const {
#ifdef DEVELOPMENT_CODE
            return _mtlist.at(i);
#else
            return _mtlist[i];
#endif
        }

        //
        // For convenience and debugging
        void show(output &s, const char* pre, const char* post) const;

    private:
        /// Make sure we can add the given minterm to this collection.
        bool _check(minterm*) const;

        /// recursive radix sort helper.
        ///     @param  k       Variable to sort on
        ///     @param  primed  Are we on the primed variable?
        ///     @param  low     Smallest index
        ///     @param  high    One past largest index
        ///
        /// We sort on elements [low, low+1, ..., high-1]
        ///
        void _sort(unsigned k, bool primed, unsigned low, unsigned high);

        ///
        /// Helper for _sort
        ///
        ///     For all items between lo (included) and hi (not included),
        ///     move any item with value val for variable x_k (or x'_k,
        ///     if pr is true) to the end of the array.
        ///
        ///     Returns the index of the first element equal to val,
        ///     after the move.
        ///     If this equals hi, then there were no elements.
        ///
        ///     Remembers the largest value not equal to val
        ///     that was seen.
        ///
        unsigned moveValuesToEnd(int val, int &max, unsigned k, bool pr,
                unsigned lo, unsigned hi);

    private:
        std::vector<minterm*> _mtlist;

        std::vector<unsigned> _varlist;

        const domain* _D;
        minterm* freelist;

        unsigned num_vars;
        bool for_relations;
        bool sparse;
        bool sorted;
};

#endif
