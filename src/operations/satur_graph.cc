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

#include "satur_graph.h"

#include "../io.h"
#include "../forest.h"

// #define TRACE

// ******************************************************************
// *                                                                *
// *                      satur_graph  methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::satur_graph::satur_graph()
{
    For = nullptr;
    var = nullptr;
    RN = nullptr;
    U = nullptr;
    node = 0;
}

MEDDLY::satur_graph::~satur_graph()
{
    if (RN) {
        MEDDLY_DCASSERT(For);
        For->doneRelNode(RN);
    }
    unpacked_node::Recycle(U);
}

void MEDDLY::satur_graph::attach(forest* F, int lvl, bool fwd)
{
    MEDDLY_DCASSERT(!For);
    For = F;
    level = lvl;
    forwd = fwd;

    var = For->getDomain()->getVar(level);
    U = unpacked_node::New(For, SPARSE_ONLY);
}

void MEDDLY::satur_graph::show(output &s) const
{
    if (!For) return;
    s.put("graph at level ");
    s.put(level);
    s.put(", node ");
    s.put(node);
    s.put("\n");
    for (unsigned i=0; i<rowptr.size(); i++) {
        s.put("Row ");
        s.put(i);
        s.put(":");
        if (isRowUnexplored(i)) {
            s.put(" unexplored\n");
            continue;
        }
        s.put("\n");
        if (diagonals[i]) {
            s.put("    diag: ");
            s.put(diagonals[i]);
            s.put("\n");
        }
        for (int z=rowptr[i]; elements[z].index >= 0; ++z) {
            MEDDLY_DCASSERT(z>=0);
            s.put("    ");
            s.put(elements[z].index);
            s.put(": ");
            s.put(elements[z].down);
            s.put("\n");
        }
    }
}

void MEDDLY::satur_graph::_restart(node_handle n)
{
    node = n;

#ifdef TRACE
    std::cout << "saturating level " << level
              << " (mxd " << n << ")\n";
#endif
    //
    // Clear old
    //
    for (unsigned i=0; i<diagonals.size(); i++) {
        diagonals[i] = 0;
    }
    elements.clear();
    for (unsigned i=0; i<rowptr.size(); i++) {
        rowptr[i] = -1;
    }
    if (RN) {
        For->doneRelNode(RN);
        RN = nullptr;
    }

    if (0==n) return;

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
    // For backward exploration, build the rows as usual
    // and then transpose the matrix
    //
    for (unsigned i=0; i<var->getBound(false); i++) {
        exploreRow(i);
    }
    transpose();
}

void MEDDLY::satur_graph::exploreRow(unsigned i)
{
    unsigned max_index = i;
    rowptr[i] = elements.size();

    if (RN->outgoing(i, *U)) {
        for (unsigned z=0; z<U->getSize(); z++) {
            const unsigned j = U->index(z);
            if (j==i) {
                // Set diagonal
                diagonals[i] = U->down(z);
            } else {
                // append (j, down) to row i
                elements.push_back(sparse_element(j, U->down(z)));
                max_index = MAX(max_index, j);
            }
        } // for z
    }
    // null terminate the list of elements
    elements.push_back(sparse_element(-1, 0));
    expandRows(1+max_index);
}

void MEDDLY::satur_graph::transpose()
{
    // TBD
}
