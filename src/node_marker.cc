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
        const node_storage *nm)
    : marked(permanent, &H), nodeHead(H)
{
    nodeMan = nm;
}

MEDDLY::node_marker::~node_marker()
{
    // No other data structures to destroy
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
