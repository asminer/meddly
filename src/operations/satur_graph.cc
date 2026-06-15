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
// #define DEBUG_TRANSPOSE

#ifdef DEBUG_TRANSPOSE

inline void showVector(MEDDLY::output &s, const char* name,
        const std::vector <int> &A)
{
    s.put(name);
    s.put('[');
    for (unsigned i=0; i<A.size(); i++) {
        if (i) s.put(", ");
        s.put(A[i]);
    }
    s.put(']');
    s.put('\n');
}

inline void showVector(MEDDLY::output &s, const char* name,
        const std::vector <MEDDLY::satur_graph::sparse_element> &A)
{
    s.put(name);
    s.put('[');
    for (unsigned i=0; i<A.size(); i++) {
        if (i) s.put(", ");
        s.put('(');
        s.put(A[i].index);
        s.put(':');
        s.put(A[i].down);
        s.put(')');
    }
    s.put(']');
    s.put('\n');
}

#endif

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
    // Backward exploration: build the entire transpose, now
    //
    if (!forwd) buildTranspose();

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

void MEDDLY::satur_graph::buildTranspose()
{
#ifdef DEBUG_TRANSPOSE
    ostream_output out(std::cout);
    out.put("Transposing relation:\n");
    RN->show(out);
#endif

    //
    // Go through and count number of elements per column
    //
    MEDDLY_DCASSERT(var->getBound(false) >= 0);
    MEDDLY_DCASSERT(var->getBound(true) >= 0);
    const unsigned num_rows = var->getBound(false);
    const unsigned num_cols = var->getBound(true);
    std::vector <int> col_counts(num_cols, 0);
    for (unsigned i=0; i<num_rows; i++) {
        if (RN->outgoing(i, *U)) {
            for (unsigned z=0; z<U->getSize(); z++) {
                const unsigned j = U->index(z);
                if (i != j) {
                    // count column j unless this is a diagonal entry
                    ++col_counts[j];
                }
            }
        }
    } // for i
#ifdef DEBUG_TRANSPOSE
    showVector(out, "col_counts: ", col_counts);
#endif

    //
    // Build the rowptr array
    //
    int acc = 0;
    for (unsigned i=0; i<rowptr.size(); i++) {
        rowptr[i] = acc;
        acc += col_counts[i] + 1;   // we null-terminate every row
    }

#ifdef DEBUG_TRANSPOSE
    showVector(out, "rowptr: ", rowptr);
    out.put("total size of entry array: ");
    out.put(long(acc));
    out.put('\n');
#endif

    //
    // Fill elements with 0
    //
    elements.resize(acc, sparse_element(-1, 0));

    //
    // Go through and add elements in transpose order.
    // This will update the rowptr array, so we will correct it after.
    //
    for (unsigned i=0; i<num_rows; i++) {
        if (RN->outgoing(i, *U)) {
            for (unsigned z=0; z<U->getSize(); z++) {
                const unsigned j = U->index(z);
                if (i != j) {
                    elements[ rowptr[j]++ ] = sparse_element(i, U->down(z));
                } else {
                    diagonals[i] = U->down(z);
                }
            }
        }
    } // for i

    //
    // Rebuild the rowptr array
    //
    acc = 0;
    for (unsigned i=0; i<rowptr.size(); i++) {
        rowptr[i] = acc;
        acc += col_counts[i] + 1;   // we null-terminate every row
    }

#ifdef DEBUG_TRANSPOSE
    out.put("Done building\n");
    show(out);
#endif
}

