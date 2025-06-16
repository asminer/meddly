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
// *                       relforest  members                       *
// *                                                                *
// ******************************************************************

std::vector <MEDDLY::relforest*> MEDDLY::relforest::all_forests;

// ******************************************************************
// *                                                                *
// *                    relforest::node  methods                    *
// *                                                                *
// ******************************************************************

bool MEDDLY::relforest::node::outgoing(unsigned i, slice &u)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

bool MEDDLY::relforest::node::incoming(unsigned i, slice &u)
{
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

MEDDLY::relforest::node::node(const relforest* _parent) : parent(_parent)
{
#ifdef DEVELOPMENT_CODE
    my_ID = ~0;
    my_level = ~0;
#endif
}

MEDDLY::relforest::node::~node()
{
}


// ******************************************************************
// *                                                                *
// *                       relforest  methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::relforest::relforest(domain *_D)
{
    D = _D;

    registerForest(this);
}

MEDDLY::relforest::~relforest()
{
    unregisterForest(this);
}

