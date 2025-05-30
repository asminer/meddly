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

#include "relforest.h"

// ******************************************************************
// *                                                                *
// *                       relforest  methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::relforest::relforest(domain *_D)
{
    D = _D;
}

MEDDLY::relforest::~relforest()
{
}

unsigned MEDDLY::relforest::levelOf(node_handle ID) const
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

bool MEDDLY::relforest::outgoing(node_handle ID, unsigned i, unpacked &u)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

bool MEDDLY::relforest::incoming(node_handle ID, unsigned i, unpacked &u)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

