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

#include "oper_minterm.h"
#include "error.h"
#include "forest.h"
#include "minterms.h"

// ******************************************************************
// *                                                                *
// *                   minterm_operation  methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::minterm_operation::minterm_operation(forest* arg1,
        minterm_coll &a2, forest* res) : arg2(a2)
{
    arg1F = arg1;
    resF  = res;
}

MEDDLY::minterm_operation::~minterm_operation()
{
}

void MEDDLY::minterm_operation::compute(const dd_edge &ar1, dd_edge &res)
{
    if (0==arg2.size()) {
        res = ar1;
        return;
    }

    node_handle resp;
    compute(resF->getMaxLevelIndex(), ~0,
            ar1.getEdgeValue(), ar1.getNode(),
            0, arg2.size(),
            res.setEdgeValue(), resp);
    res.set(resp);
}

