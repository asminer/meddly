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

#include "node_marker.h"


MEDDLY::node_marker::node_marker(bool permanent, node_headers &H,
        const node_storage *nm, expert_forest* F)
    : marked(permanent, &H), nodeHead(H)
{
    nodeMan = nm;
    For = F;

    MEDDLY_DCASSERT(nodeMan);
    MEDDLY_DCASSERT(For);
}

MEDDLY::node_marker::~node_marker()
{
    // No other data structures to destroy
}

//
// Public methods
//

size_t MEDDLY::node_marker::countEdges() const
{
    size_t ec = 0;
    unpacked_node* M = unpacked_node::New();
    node_handle i=0;
    while( (i=marked.firstOne(i+1)) < marked.getSize() )
    {
        For->unpackNode(M, i, FULL_ONLY);
        ec += M->getSize();
    }
    unpacked_node::recycle(M);
    return ec;
}

size_t MEDDLY::node_marker::countNonzeroEdges() const
{
    size_t ec = 0;
    unpacked_node* M = unpacked_node::New();
    node_handle i=0;
    while( (i=marked.firstOne(i+1)) < marked.getSize() )
    {
        For->unpackNode(M, i, SPARSE_ONLY);
        ec += M->getNNZs();
    }
    unpacked_node::recycle(M);
    return ec;
}

void MEDDLY::node_marker::showByLevels(output &s) const
{
    unpacked_node* M = unpacked_node::New();

    const unsigned lwid = digits(For->getNumVariables());
    const unsigned nwid = digits(getSize());

    for (int k=For->getNumVariables(); k; k = forest::downLevel(k)) {

        bool level_printed = false;

        node_handle i=0;
        while( (i=marked.firstOne(i+1)) < marked.getSize() )
        {
            if (For->getNodeLevel(i) != k) continue;
            //
            // Node i is at level k.
            //

            if (!level_printed) {
                //
                // Show level name
                //
                s << "Level: ";
                s.put(k, lwid);
                s << " Var: ";
                const variable* v = For->getDomain()->getVar(
                    For->getVarByLevel(ABS(k))
                );
                char primed = (k>0) ? ' ' : '\'';
                if (v->getName()) {
                    s << v->getName() << primed << '\n';
                } else {
                    s << For->getVarByLevel(ABS(k)) << primed << '\n';
                }
                level_printed = true;
            }

            //
            // Show the node
            //
            s << "    node:";
            s.put(i, nwid);
            s.put(' ');
            For->unpackNode(M, i, FULL_OR_SPARSE);
            M->show(s, true);
            s.put('\n');

        } // for i

    } // for k
    unpacked_node::recycle(M);
}

//
// Private methods
//

void MEDDLY::node_marker::_mark(node_handle p)
{
    CHECK_RANGE(__FILE__, __LINE__, 1, p, (node_handle) marked.getSize());
    MEDDLY_DCASSERT(nodeMan);

    marked.set(p, true);
    nodeMan->markDownPointers( *this, nodeHead.getNodeAddress(p) );
}


