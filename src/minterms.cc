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

#include "forest.h"

MEDDLY::minterm::minterm(const forest* parent, node_storage_flags fs)
    : termval(true)
{
    if (parent) {
        num_vars = parent->getNumVariables();
        for_relations = parent->isForRelations();
    } else {
        num_vars = 0;
        for_relations = false;
    }
    sparse = (fs == SPARSE_ONLY);
    all_assigned = !sparse;
}

bool MEDDLY::minterm::contains(const minterm& m) const
{
    if ( (m.isForRelations() != isForRelations()) ||
         (m.numVars() != numVars()) )
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
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

