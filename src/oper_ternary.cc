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

#include "oper_ternary.h"
#include "error.h"
#include "forest.h"

// ******************************************************************
// *                    ternary_operation methods                    *
// ******************************************************************

MEDDLY::ternary_operation::ternary_operation(ternary_list& owner,
    unsigned et_slots, forest* arg1, forest* arg2, forest* arg3,
    forest* res) : operation(owner.getName(), et_slots), parent(owner)
{
    arg1F = arg1;
    arg2F = arg2;
    arg3F = arg3;
    resF = res;

    registerInForest(arg1F);
    registerInForest(arg2F);
    registerInForest(arg3F);
    registerInForest(resF);
}

MEDDLY::ternary_operation::~ternary_operation()
{
    unregisterInForest(arg1F);
    unregisterInForest(arg2F);
    unregisterInForest(arg3F);
    unregisterInForest(resF);

    parent.remove(this);
}

void MEDDLY::ternary_operation::compute(const dd_edge &ar1,
        const dd_edge &ar2, const dd_edge &ar3, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    computeDDEdge(ar1, ar2, ar3, res, true);
}

#ifdef ALLOW_DEPRECATED_0_17_6

void MEDDLY::ternary_operation::computeTemp(const dd_edge &ar1,
        const dd_edge &ar2, const dd_edge &ar3, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    computeDDEdge(ar1, ar2, ar3, res, false);
}

#endif

// ******************************************************************
// *                      ternary_list methods                      *
// ******************************************************************

MEDDLY::ternary_list::ternary_list(const char* n)
{
    reset(n);
}

void MEDDLY::ternary_list::reset(const char* n)
{
    front = nullptr;
    name = n;
}

MEDDLY::ternary_operation*
MEDDLY::ternary_list::mtfTernary(const forest* arg1F,
        const forest* arg2F, const forest* arg3F, const forest* resF)
{
    ternary_operation* prev = front;
    ternary_operation* curr = front->next;
    while (curr) {
        if ((curr->arg1F == arg1F) &&
            (curr->arg2F == arg2F) &&
            (curr->arg3F == arg3F) &&
            (curr->resF == resF))
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


void MEDDLY::ternary_list::searchRemove(ternary_operation* top)
{
    if (!front) return;
    ternary_operation* prev = front;
    ternary_operation* curr = front->next;
    while (curr) {
        if (curr == top) {
            prev->next = curr->next;
            return;
        }
        prev = curr;
        curr = curr->next;
    }
}



