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

#include "oper_unary.h"
#include "error.h"
#include "forest.h"
#include "dd_edge.h"

// ******************************************************************
// *                    unary_operation  methods                    *
// ******************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

MEDDLY::unary_operation::unary_operation(unary_list& owner,
    unsigned et_slots, forest* arg, forest* res)
    : operation(owner.getName(), et_slots)
{
    parent = &owner;

    argF = arg;
    resultType = opnd_type::FOREST;
    resF = res;

    registerInForest(argF);
    registerInForest(resF);
    new_style = false;
}

MEDDLY::unary_operation::unary_operation(unary_list& owner,
    unsigned et_slots, forest* arg, opnd_type res)
    : operation(owner.getName(), et_slots)
{
    parent = &owner;
    argF = arg;
    resultType = res;
    resF = nullptr;

    registerInForest(argF);
    new_style = false;
}

#endif

MEDDLY::unary_operation::unary_operation(forest* arg, forest* res)
    : operation()
{
    parent = nullptr;

    argF = arg;
    resultType = opnd_type::FOREST;
    resF = res;

    registerInForest(argF);
    registerInForest(resF);

#ifdef ALLOW_DEPRECATED_0_17_6
    new_style = true;
#endif
}

MEDDLY::unary_operation::unary_operation(forest* arg, opnd_type res)
    : operation()
{
    parent = nullptr;

    argF = arg;
    resultType = res;
    resF = nullptr;

    registerInForest(argF);

#ifdef ALLOW_DEPRECATED_0_17_6
    new_style = true;
#endif
}

MEDDLY::unary_operation::~unary_operation()
{
    unregisterInForest(argF);
    unregisterInForest(resF);
    if (parent) parent->remove(this);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
#ifdef ALLOW_DEPRECATED_0_17_6
    if (new_style) {
        node_handle resp;
        compute(arg.getEdgeValue(), arg.getNode(),
                resF->getMaxLevelIndex(),
                res.setEdgeValue(), resp);
        res.set(resp);
    } else {
        computeDDEdge(arg, res, true);
    }
#else
    node_handle resp;
    compute(arg.getEdgeValue(), arg.getNode(),
            resF->getMaxLevelIndex(),
            res.setEdgeValue(), resp);
    res.set(resp);
#endif
}

#ifdef ALLOW_DEPRECATED_0_17_6
void MEDDLY::unary_operation::computeTemp(const dd_edge &arg, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }

    if (new_style) {
    } else {
        computeDDEdge(arg, res, false);
        node_handle resp;
        compute(arg.getEdgeValue(), arg.getNode(),
                0,
                res.setEdgeValue(), resp);
        res.set(resp);
    }
}

void MEDDLY::unary_operation::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}
#endif

// new compute methods

void MEDDLY::unary_operation::compute(const edge_value &av, node_handle ap,
                int L, edge_value &cv, node_handle &cp)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::unary_operation::compute(int L, const edge_value &av,
                node_handle ap, oper_item &res)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


// ******************************************************************
// *                       unary_list methods                       *
// ******************************************************************

MEDDLY::unary_list::unary_list(const char* n)
{
    reset(n);
}

void MEDDLY::unary_list::reset(const char* n)
{
    front = nullptr;
    name = n;
}

MEDDLY::unary_operation*
MEDDLY::unary_list::mtfUnary(const forest* argF, const forest* resF)
{
    unary_operation* prev = front;
    unary_operation* curr = front->next;
    while (curr) {
        if ((curr->argF == argF) && (curr->resF == resF)) {
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

MEDDLY::unary_operation*
MEDDLY::unary_list::mtfUnary(const forest* argF, opnd_type resType)
{
    unary_operation* prev = front;
    unary_operation* curr = front->next;
    while (curr) {
        if ((curr->argF == argF) && (curr->resultType == resType)) {
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


void MEDDLY::unary_list::searchRemove(unary_operation* uop)
{
    if (!front) return;
    unary_operation* prev = front;
    unary_operation* curr = front->next;
    while (curr) {
        if (curr == uop) {
            prev->next = curr->next;
            return;
        }
        prev = curr;
        curr = curr->next;
    }
}

