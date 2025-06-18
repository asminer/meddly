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

#include "oper_binrel.h"
#include "error.h"
#include "forest.h"

#define DEBUG_THROWS

// ******************************************************************
// *                    binrel_operation methods                    *
// ******************************************************************

MEDDLY::binrel_operation::binrel_operation(forest* arg1, relforest* arg2,
        forest* res): operation()
{
    parent = nullptr;

    arg1F = arg1;
    arg2F = arg2;
    resF = res;

    registerInForest(arg1F);
    registerInForest(arg2F);
    registerInForest(resF);
}


MEDDLY::binrel_operation::~binrel_operation()
{
    unregisterInForest(arg1F);
    unregisterInForest(arg2F);
    unregisterInForest(resF);

    if (parent) parent->remove(this);
}


void MEDDLY::binrel_operation::compute(const dd_edge &ar1,
        const relforest::edge &ar2, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    node_handle resp;
    compute(resF->getMaxLevelIndex(), ~0,
            ar1.getEdgeValue(), ar1.getNode(),
            ar2.value, ar2.root,
            res.setEdgeValue(), resp);
    res.set(resp);
#ifdef DEVELOPMENT_CODE
    resF->validateIncounts(true, __FILE__, __LINE__, getName());
#endif
}

// ******************************************************************
// *                      binrel_list  methods                      *
// ******************************************************************

MEDDLY::binrel_list::binrel_list(const char* n)
{
    reset(n);
}

void MEDDLY::binrel_list::reset(const char* n)
{
    front = nullptr;
    name = n;
}

MEDDLY::binrel_operation*
MEDDLY::binrel_list::mtfBinary(const forest* arg1F,
        const relforest* arg2F, const forest* resF)
{
    binrel_operation* prev = front;
    binrel_operation* curr = front->next;
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


void MEDDLY::binrel_list::searchRemove(binrel_operation* bop)
{
    if (!front) return;
    binrel_operation* prev = front;
    binrel_operation* curr = front->next;
    while (curr) {
        if (curr == bop) {
            prev->next = curr->next;
            return;
        }
        prev = curr;
        curr = curr->next;
    }
}


