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

#include "minterms.h"
#include "policies.h"
#include "terminal.h"

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


    class forest;
    class minterm;

};  // namespace MEDDLY


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
        minterm(const forest* parent, node_storage_flags fs);

        inline bool isSparse()              const   { return sparse; }
        inline bool isFull()                const   { return !sparse; }

        inline bool isForSets()             const   { return !for_relations; }
        inline bool isForRelations()        const   { return for_relations; }

        inline unsigned numVars()           const   { return num_vars; }

        /// Return true iff there are no don't cares/don't change
        inline bool areAllAssigned()        const   { return all_assigned; }

        /// Set the terminal value
        inline void setTerm(const terminal &t)      { termval = t; }

        /// Get the terminal value
        inline const terminal& getTerm()    const   { return termval; }

        // access for sets

        inline void getVar(unsigned i, int &from) const {
            MEDDLY_DCASSERT(isForSets());
#ifdef DEVELOPMENT_CODE
            from = _from.at(i);
#else
            from = _from[i];
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
            if (from<0) {
                all_assigned = false;
            }
        }

        // access for relations

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
#ifdef DEVELOPMENT_CODE
            _from.at(i) = from;
            _to.at(i) = to;
#else
            _from[i] = from;
            _to[i] = to;
#endif
            if (from<0 || to<0) {
                all_assigned = false;
            }
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
        // check if we are equal to m,
        // except we may contain more don't cares
        //
        bool contains(const minterm& m) const;

    private:
        inline bool matches(unsigned i, const minterm& m, unsigned j) const
        {
            MEDDLY_DCASSERT(isForRelations() == m.isForRelations());

            if ((DONT_CARE != _from.at(i)) && (_from.at(i) == m._from.at(j)))
                return false;

            if (!isForRelations())  return true;

            if (DONT_CARE == _to.at(i)) return true;
            if (_to.at(i) == m._to.at(j)) return true;

            return (DONT_CHANGE == _to.at(i))
                && (m._to.at(j) == m._from.at(j));
        }

        inline bool isDontCare(unsigned i) const
        {
            if (DONT_CARE != _from.at(i)) return false;
            if (!isForRelations()) return true;
            return DONT_CARE == _to.at(i);
        }

    private:
        std::vector<int> _from;
        std::vector<int> _to;
        std::vector<unsigned> _index;

        terminal termval;

        unsigned num_vars;
        bool for_relations;
        bool sparse;
        bool all_assigned;
};

#endif
