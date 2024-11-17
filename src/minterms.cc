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

#include "minterms.h"

#include "io.h"
#include "forest.h"

// ******************************************************************
// *                                                                *
// *                        minterm  methods                        *
// *                                                                *
// ******************************************************************


MEDDLY::minterm::minterm(const forest* parent, node_storage_flags fs)
    : termval(true)
{
    _parent = parent;
    if (parent) {
        num_vars = parent->getNumVariables();
        for_relations = parent->isForRelations();
    } else {
        num_vars = 0;
        for_relations = false;
    }
    sparse = (fs == SPARSE_ONLY);
    if (!sparse) {
        _from.resize(1+num_vars, DONT_CARE);
        if (for_relations) {
            _to.resize(1+num_vars, DONT_CARE);
        }
    }
}

bool MEDDLY::minterm::contains(const minterm& m) const
{
    if ( (m.isForRelations() != isForRelations()) ||
         (m.numVars() != numVars()) )
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    if ( m.getParent() != getParent() )
    {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    if (isFull()) {
        //
        // We use full storage
        //
        if (m.isFull()) {
            //
            // m is also full
            //
            for (unsigned i=1; i<=num_vars; i++) {
                if (!matches(i, m, i)) return false;
            }
            return true;
        } else {
            //
            // m is sparse; missing entries are don't care.
            // first check that all entries in m match
            //
            unsigned z=0;
            for (unsigned i=1; i<=num_vars; i++) {
                if (z<m.size() && m.whichVar(z)==i) {
                    if (!matches(i, m, z)) return false;
                    z++;
                } else {
                    if (!isDontCare(i)) return false;
                }
            }
            return true;
        }
    } else {
        //
        // We use sparse storage;
        //

        if (m.isFull()) {
            //
            // m is full. we can skip the missing don't care entries
            //
            for (unsigned z=0; z<size(); z++) {
                const unsigned i = _index.at(z);
                if (!matches(z, m, i)) return false;
            }
            return true;
        } else {
            //
            // m is also sparse.
            // Make sure all of our entries match m
            //
            unsigned mz=0;
            for (unsigned z=0; z<size(); z++) {
                while (mz<m.size() && m.whichVar(mz) < whichVar(z)) {
                    // skip all extra entries in m
                    ++mz;
                }
                if (mz >= m.size()) return false;
                if (m.whichVar(mz) > whichVar(z)) return false;
                MEDDLY_DCASSERT(m.whichVar(mz) == whichVar(z));

                if (!matches(z, m, mz)) return false;
            }
            return true;
        }
    }
}

bool MEDDLY::minterm::hasSparseVars(const std::vector <unsigned> &vars) const
{
    if (isFull()) return false;
    if (_index.size() != vars.size()) return false;
    for (unsigned i=0; i<vars.size(); i++) {
        if (_index[i] != vars[i]) return false;
    }
    return true;
}

void MEDDLY::minterm::show(output &s) const
{
    if (isFull()) {
        //
        // Full storage
        //
        s.put("[ top");
        for (unsigned i=num_vars; i; i--) {
            s.put(", ");
            if (_from[i] < 0)   s.put('x');
            else                s.put(_from[i]);
            if (isForSets()) continue;
            s.put("->");
            if (DONT_CARE == _to[i])        s.put('x');
            else if (DONT_CHANGE == _to[i]) s.put('i');
            else                            s.put(_to[i]);
        }
        s.put(", bot]");
    } else {
        //
        // Sparse storage
        //
        s.put('(');
        for (unsigned z=0; z<size(); z++) {
            if (z) s.put(", ");
            s.put(whichVar(z));
            s.put(':');
            if (_from[z] < 0)   s.put('x');
            else                s.put(_from[z]);
            if (isForSets()) continue;
            s.put("->");
            if (DONT_CARE == _to[z])        s.put('x');
            else if (DONT_CHANGE == _to[z]) s.put('i');
            else                            s.put(_to[z]);
        }
        s.put(')');
    }
}

// ******************************************************************
// *                                                                *
// *                      minterm_coll methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::minterm_coll::minterm_coll(const forest* parent)
{
    _parent = parent;
    if (parent) {
        num_vars = parent->getNumVariables();
        for_relations = parent->isForRelations();
    } else {
        num_vars = 0;
        for_relations = false;
    }
    sparse = false;
    freelist = nullptr;
    sorted = false;
}

MEDDLY::minterm_coll::minterm_coll(const forest* parent,
        const std::vector<unsigned> &varlist)
{
    _parent = parent;
    if (parent) {
        num_vars = parent->getNumVariables();
        for_relations = parent->isForRelations();
    } else {
        num_vars = 0;
        for_relations = false;
    }
    sparse = true;
    freelist = nullptr;
    sorted = false;

    _varlist.resize(varlist.size());
    for (unsigned i=0; i<varlist.size(); i++) {
        _varlist[i] = varlist[i];
    }
}

void MEDDLY::minterm_coll::clear()
{
    while (_mtlist.size()) {
        recycle( _mtlist.back() );
        _mtlist.back() = nullptr;
        _mtlist.pop_back();
    }
    sorted = false;
}

void MEDDLY::minterm_coll::show(output &s,
        const char* pre, const char* post) const
{
    for (unsigned i=0; i<_mtlist.size(); i++) {
        s.put(pre);
        _mtlist[i]->show(s);
        s.put(post);
    }
}

bool MEDDLY::minterm_coll::_check(minterm* m) const
{
    MEDDLY_DCASSERT(m->getParent() == _parent);
    if (sparse) {
        MEDDLY_DCASSERT(!m->isFull());
        MEDDLY_DCASSERT(m->hasSparseVars(_varlist));
    } else {
        MEDDLY_DCASSERT(m->isFull());
    }
    return true;
}

void MEDDLY::minterm_coll::_sort(unsigned k, unsigned low, unsigned high)
{
    if (0==k) return;
    if (low+1 == high) return;
    MEDDLY_DCASSERT(low < high);
}

