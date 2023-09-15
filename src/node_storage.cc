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

#include "defines.h"
#include "node_storage.h"
#include "io.h"
#include "operators.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                   node_storage_style methods                   *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_storage_style::node_storage_style(const char* n)
{
    name = n;
}

MEDDLY::node_storage_style::~node_storage_style()
{
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_storage methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_storage::node_storage(const char* n, expert_forest* f)
{
    style_name = n;
    parent = f;
}

MEDDLY::node_storage::~node_storage()
{
    // nothing, derived classes must handle everything
}

void MEDDLY::node_storage::dumpInternal(output &s, unsigned flags) const
{
    dumpInternalInfo(s);
    s << "Data array by record:\n";
    for (node_address a=firstNodeAddress(); a > 0; ) {
        s.flush();
        a = dumpInternalNode(s, a, flags);
    } // for a
    dumpInternalTail(s);
    s.flush();
}

