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

MEDDLY::binary_operation::binary_operation(binary_list& owner,
    unsigned et_slots, forest* arg1, forest* arg2,
    forest* res) : operation(owner.getName(), et_slots), parent(owner)
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

    parent.removeOperation(this);
}

/*
bool MEDDLY::binary_operation::matches(const dd_edge &arg1,
        const dd_edge &arg2, const dd_edge &res) const
{
    return
        arg1.isAttachedTo(arg1F) &&
        arg2.isAttachedTo(arg2F) &&
        res.isAttachedTo(resF);
}
*/

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

// ******************************************************************
// *                      binary_list  methods                      *
// ******************************************************************

MEDDLY::binary_list::binary_list(const char* n)
{
    name = n;
    front = nullptr;
}

MEDDLY::binary_operation*
MEDDLY::binary_list::mtfBinary(const forest* arg1F,
        const forest* arg2F, const forest* resF)
{
    binary_operation* prev = front;
    binary_operation* curr = front->next;
    while (curr) {
        if ((curr->arg1F == arg1F) && (curr->arg2F == arg2F) && (curr->resF == resF))
        {
            // Move to front
            prev->next = curr->next;
            curr->next = front;
            front = curr;
            return curr;
        }
        prev = curr;
        curr = curr->next;
    }
    return nullptr;
}


void MEDDLY::binary_list::searchRemove(binary_operation* bop)
{
    if (!front) return;
    binary_operation* prev = front;
    binary_operation* curr = front->next;
    while (curr) {
        if (curr == bop) {
            prev->next = curr->next;
            return;
        }
        prev = curr;
        curr = curr->next;
    }
}


