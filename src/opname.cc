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
#include "opname.h"
#include "error.h"
#include "dd_edge.h"
#include "forest.h"

#include "initializer.h"

#include "oper.h"
#include "oper_unary.h"
#include "oper_binary.h"

// These are needed for builtin opname initialization

// ******************************************************************
// *                                                                *
// *                         opname methods                         *
// *                                                                *
// ******************************************************************

MEDDLY::opname::opname(const char* n)
{
    name = n;
}

MEDDLY::opname::~opname()
{
}

// ******************************************************************
// *                                                                *
// *                      unary_opname methods                      *
// *                                                                *
// ******************************************************************

/*
MEDDLY::unary_opname::unary_opname(const char* n)
    : opname(n), cache(n)
{
}

MEDDLY::unary_opname::~unary_opname()
{
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::getOperation(const dd_edge &ar, const dd_edge &res)
{
    unary_operation* match = cache.find(ar.getForest(), res.getForest());
    if (match) return match;

    return cache.add(
            buildOperation(cache, ar.getForest(), res.getForest())
    );
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::getOperation(const dd_edge &ar, opnd_type res)
{
    unary_operation* match = cache.find(ar.getForest(), res);
    if (match) return match;

    return cache.add(
            buildOperation(cache, ar.getForest(), res)
    );
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(unary_list &, forest *, forest *)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(unary_list &, forest *, opnd_type)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}
*/

// ******************************************************************
// *                                                                *
// *                     binary_opname  methods                     *
// *                                                                *
// ******************************************************************

/*
MEDDLY::binary_opname::binary_opname(const char* n)
    : opname(n), cache(n)
{
}

MEDDLY::binary_opname::~binary_opname()
{
}

MEDDLY::binary_operation*
MEDDLY::binary_opname::getOperation(const dd_edge &a1, const dd_edge &a2,
        const dd_edge &res)
{
    binary_operation* match = cache.find(
            a1.getForest(), a2.getForest(), res.getForest()
    );
    if (match) return match;

    return cache.add(
            buildOperation(cache, a1.getForest(), a2.getForest(), res.getForest())
    );
}
*/

// ******************************************************************
// *                                                                *
// *                   specialized_opname methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::specialized_opname::arguments::arguments()
{
    setAutoDestroy(true);
}

MEDDLY::specialized_opname::arguments::~arguments()
{
}


MEDDLY::specialized_opname::specialized_opname(const char* n) : opname(n)
{
}

MEDDLY::specialized_opname::~specialized_opname()
{
}


