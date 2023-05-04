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

#include "oper_binary.h"
#include "error.h"
#include "forest.h"

// ******************************************************************
// *                    binary_operation methods                    *
// ******************************************************************

MEDDLY::binary_operation::binary_operation(binary_opname* op,
    unsigned et_slots, expert_forest* arg1, expert_forest* arg2,
    expert_forest* res) : operation(op, et_slots)
{
    arg1F = arg1;
    arg2F = arg2;
    resF = res;

    registerInForest(arg1F);
    registerInForest(arg2F);
    registerInForest(resF);

    can_commute = false;
}

MEDDLY::binary_operation::~binary_operation()
{
    unregisterInForest(arg1F);
    unregisterInForest(arg2F);
    unregisterInForest(resF);
}

bool MEDDLY::binary_operation::matches(const dd_edge &arg1,
        const dd_edge &arg2, const dd_edge &res) const
{
    return (arg1.getForest() == arg1F
        && arg2.getForest() == arg2F
        && res.getForest() == resF);
}


void MEDDLY::binary_operation::compute(const dd_edge &ar1,
    const dd_edge &ar2, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    computeDDEdge(ar1, ar2, res, true);
}

void MEDDLY::binary_operation::computeTemp(const dd_edge &ar1,
    const dd_edge &ar2, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    computeDDEdge(ar1, ar2, res, false);
}

bool MEDDLY::binary_operation::checkForestCompatibility() const
{
    auto o1 = arg1F->variableOrder();
    auto o2 = arg2F->variableOrder();
    auto o3 = resF->variableOrder();
    return o1->is_compatible_with(*o2) && o1->is_compatible_with(*o3);
}


