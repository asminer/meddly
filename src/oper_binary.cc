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
#include "forest_levels.h"

#define DEBUG_THROWS

// ******************************************************************
// *                    binary_operation methods                    *
// ******************************************************************

#ifdef ALLOW_OLD_BINARY_0_17_6

MEDDLY::binary_operation::binary_operation(binary_list& owner,
    unsigned et_slots, forest* arg1, forest* arg2,
    forest* res) : operation(owner.getName(), et_slots)
{
    factory = nullptr;
    parent = &owner;

    arg1F = arg1;
    arg2F = arg2;
    resF = res;

    registerInForest(arg1F);
    registerInForest(arg2F);
    registerInForest(resF);
    can_commute = false;
    new_style = false;
}

#endif

MEDDLY::binary_operation::binary_operation(forest* arg1, forest* arg2,
        forest* res): operation()
{
    factory = nullptr;
#ifdef ALLOW_OLD_BINARY_0_17_6
    parent = nullptr;
#endif

    arg1F = arg1;
    arg2F = arg2;
    resF = res;

    registerInForest(arg1F);
    registerInForest(arg2F);
    registerInForest(resF);

#ifdef ALLOW_OLD_BINARY_0_17_6
    can_commute = false;
    new_style = true;
#endif
}


MEDDLY::binary_operation::~binary_operation()
{
    unregisterInForest(arg1F);
    unregisterInForest(arg2F);
    unregisterInForest(resF);

    if (factory) factory->remove(this);
#ifdef ALLOW_OLD_BINARY_0_17_6
    if (parent) parent->remove(this);
#endif
}


void MEDDLY::binary_operation::compute(const dd_edge &ar1,
        const dd_edge &ar2, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
#ifdef ALLOW_OLD_BINARY_0_17_6
    if (new_style) {
        node_handle resp;
        compute(resF->getMaxLevelIndex(), ~0,
                ar1.getEdgeValue(), ar1.getNode(),
                ar2.getEdgeValue(), ar2.getNode(),
                res.setEdgeValue(), resp);
        res.set(resp);
    } else {
        computeDDEdge(ar1, ar2, res, true);
   }
#else
    node_handle resp;
    compute(resF->getMaxLevelIndex(), ~0,
            ar1.getEdgeValue(), ar1.getNode(),
            ar2.getEdgeValue(), ar2.getNode(),
            res.setEdgeValue(), resp);
    res.set(resp);
#endif
#ifdef DEVELOPMENT_CODE
    resF->validateIncounts(true, __FILE__, __LINE__, getName());
#endif
}

#ifdef ALLOW_OLD_BINARY_0_17_6

void MEDDLY::binary_operation::computeTemp(const dd_edge &ar1,
        const dd_edge &ar2, dd_edge &res)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    if (new_style) {
        node_handle resp;
        // THIS vvvvv is an ugly hack
        int toplevel = arg1F->isForRelations()
            ?  MXD_levels::topUnprimed(ar1.getLevel(), ar2.getLevel())
            :  MDD_levels::topLevel(ar1.getLevel(), ar2.getLevel());
        compute(toplevel, ~0,
                ar1.getEdgeValue(), ar1.getNode(),
                ar2.getEdgeValue(), ar2.getNode(),
                res.setEdgeValue(), resp);
        res.set(resp);
    } else {
        computeDDEdge(ar1, ar2, res, false);
    }
}


void MEDDLY::binary_operation::computeDDEdge(const dd_edge &ar1,
        const dd_edge &ar2, dd_edge &res, bool userFlag)
{
#ifdef DEBUG_THROWS
    std::cerr << "default computeDDEdge for binary operation " << getName() << "\n";
    FAIL(__FILE__, __LINE__);
#endif
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::binary_operation::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
#ifdef DEBUG_THROWS
    std::cerr << "default compute for binary operation " << getName() << "\n";
    FAIL(__FILE__, __LINE__);
#endif
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

#endif

// ******************************************************************
// *                     binary_factory methods                     *
// ******************************************************************

void MEDDLY::binary_factory::cleanup()
{
    _cleanup();
}

MEDDLY::binary_operation*
MEDDLY::binary_factory::build_new(forest* a, forest* b, forest* c)
{
    return nullptr;
}

void MEDDLY::binary_factory::_setup(const char* file, const char* name,
        const char* doc)
{
    _file = file;
    _name = name;
    _doc = doc;
    front = nullptr;
}

void MEDDLY::binary_factory::_cleanup()
{
    MEDDLY_DCASSERT(nullptr == front);
}

bool
MEDDLY::binary_factory::mtfBinary(const forest* arg1F,
        const forest* arg2F, const forest* resF)
{
    binary_operation* prev = front;
    binary_operation* curr = front->next;
    while (curr) {
        if ((curr->arg1F == arg1F) && (curr->arg2F == arg2F)
                && (curr->resF == resF))
        {
            // Move to front
            prev->next = curr->next;
            curr->next = front;
            front = curr;
            return true;
        }
        prev = curr;
        curr = curr->next;
    }
    return false;
}


void MEDDLY::binary_factory::searchRemove(binary_operation* bop)
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



// ******************************************************************
// *                      binary_list  methods                      *
// ******************************************************************

#ifdef ALLOW_OLD_BINARY_0_17_6

MEDDLY::binary_list::binary_list(const char* n)
{
    reset(n);
}

void MEDDLY::binary_list::reset(const char* n)
{
    front = nullptr;
    name = n;
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

#endif
