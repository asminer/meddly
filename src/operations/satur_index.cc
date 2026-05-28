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

#include "satur_index.h"
#include "../variable.h"
#include "../domain.h"
#include "../forest.h"
#include "../rel_node.h"

#define DEBUG_BASE

// **********************************************************************
// *                                                                    *
// *                     sat_index_explorer methods                     *
// *                                                                    *
// **********************************************************************

MEDDLY::sat_index_explorer::sat_index_explorer(forest* F, int lvl, bool fwd)
    : For(F), level(lvl), forwd(fwd)
{
    MEDDLY_DCASSERT(For);
    var = For->getDomain()->getVar(level);
    RN = nullptr;
    U = unpacked_node::New(For, SPARSE_ONLY);
}

MEDDLY::sat_index_explorer::~sat_index_explorer()
{
    if (RN) {
        For->doneRelNode(RN);
        RN = nullptr;
    }
    unpacked_node::Recycle(U);
}

void MEDDLY::sat_index_explorer::restart(node_handle n)
{
    //
    // Clear old
    //
    for (unsigned i=0; i<diagonals.size(); i++) {
        diagonals[i] = 0;
    }
    for (unsigned i=0; i<rows.size(); i++) {
        if (!rows[i].explored) continue;
        rows[i].elements.clear();
        rows[i].explored = false;
    }
    if (RN) {
        For->doneRelNode(RN);
    }

    //
    // Build new
    //
    RN = For->buildRelNode(n);

    // update size
    expandRows(var->getBound(!forwd));

    //
    // For forward exploration, build rows on the fly as needed.
    //
    if (forwd) return;

    //
    // For backward exploration, build everything at the beginning.
    // We explore the relation node elements [i,j]
    // and store the transpose [j,i] internally.
    //
    row_element elem;
    for (unsigned i=0; i<var->getBound(false); i++) {
        if (!RN->outgoing(i, *U)) continue;
        elem.index = i;

        for (unsigned z=0; z<U->getSize(); z++) {
            const unsigned j = U->index(z);
            if (j==i) {
                // Set diagonal
                diagonals[i] = U->down(z);
            } else {
                // append (i, down) to row j
                elem.down = U->down(z);
                rows[j].elements.push_back(elem);
            }
        } // for z
    } // for row i
    for (unsigned i=0; i<rows.size(); i++) {
        rows[i].explored = true;
    }

    finishAllRows();
}

void MEDDLY::sat_index_explorer::clear()
{
    // Do nothing
}

void MEDDLY::sat_index_explorer::exploreRow(unsigned i)
{
    MEDDLY_DCASSERT(forwd);

    unsigned max_index = i;
    if (RN->outgoing(i, *U)) {
        for (unsigned z=0; z<U->getSize(); z++) {
            const unsigned j = U->index(z);
            if (j==i) {
                // Set diagonal
                diagonals[i] = U->down(z);
            } else {
                // append (j, down) to row i
                row_element elem;
                elem.index = j;
                elem.down = U->down(z);
                rows[i].elements.push_back(elem);
                max_index = MAX(max_index, j);
            }
        } // for z
    }
    rows[i].explored = true;
    if (max_index > rows.size()) {
        expandRows(1+max_index);
    }

    finishRow(i);
}

void MEDDLY::sat_index_explorer::finishRow(unsigned)
{
    // Do nothing
}

void MEDDLY::sat_index_explorer::finishAllRows()
{
    // Do nothing
}


void MEDDLY::sat_index_explorer::expandRows(unsigned newsz)
{
    if (newsz <= rows.size()) return;

    const unsigned oldsize = rows.size();

    rows.resize(newsz);
    diagonals.resize(newsz);

    finishExpandRows(oldsize, newsz);
}

void MEDDLY::sat_index_explorer::finishExpandRows(unsigned, unsigned)
{
    // Do nothing
}

void MEDDLY::sat_index_explorer::finishUpdate(unsigned)
{
    // Do nothing
}


// **********************************************************************
// *                                                                    *
// *                             Front  end                             *
// *                                                                    *
// **********************************************************************

MEDDLY::sat_index_explorer* MEDDLY::makeSatIndexExplorer(char which,
        forest* F, int level, bool forwd)
{
    // TBD
    return nullptr;
}

