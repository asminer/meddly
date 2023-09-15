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
#include "operations/mpz_object.h"  // for mpz wrapper
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
// *                         Helper methods                         *
// ******************************************************************

template <class RT>
static MEDDLY::unary_operation*
FindUnary(MEDDLY::operation* list, const MEDDLY::dd_edge &arg, const RT &res)
{
    MEDDLY::unary_operation* match = nullptr;
    MEDDLY::operation* prev = nullptr;
    MEDDLY::operation* curr;
    for (curr = list; curr; curr = curr->getNext()) {
        match = static_cast <MEDDLY::unary_operation*> (curr);
        if (match->matches(arg, res)) {
            // Move to front, unless it's already there
            if (prev) {
                prev->setNext(curr->getNext());
                curr->setNext(list);
                list = curr;
            }
            return match;
        }
    }
    return nullptr;
}

// ******************************************************************
// *                                                                *
// *                     gmp  wrapper functions                     *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP

MEDDLY::ct_object& MEDDLY::get_mpz_wrapper()
{
    static MEDDLY::mpz_object foo;
    return foo;
}

void MEDDLY::unwrap(const ct_object &x, mpz_t &value)
{
    using namespace MEDDLY;
    const mpz_object &mx = static_cast <const mpz_object &> (x);
    mx.copyInto(value);
}

#endif


// ******************************************************************
// *                                                                *
// *                         opname methods                         *
// *                                                                *
// ******************************************************************

MEDDLY::opname::opname(const char* n)
{
    name = n;
//     index = next_index;
    // next_index++;

    cache = nullptr;
}

MEDDLY::opname::~opname()
{
    // library must be closing
}

void
MEDDLY::opname::removeOperationFromCache(operation* op)
{
    if (!op) return;
    operation* prev = nullptr;
    operation* curr;
    for (curr = cache; curr; curr = curr->getNext()) {
        if (curr == op) break;
        prev = curr;
    }
    if (!curr) return;

    // Remove curr from list
    if (prev) {
        prev->setNext(curr->getNext());
    } else {
        cache = curr->getNext();
    }
    curr->setNext(nullptr);
}


// ******************************************************************
// *                                                                *
// *                      unary_opname methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname::unary_opname(const char* n) : opname(n)
{
}

MEDDLY::unary_opname::~unary_opname()
{
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::getOperation(const dd_edge &ar, const dd_edge &res)
{
    unary_operation* match = FindUnary(cache, ar, res);
    if (match) {
        cache = match;
        return match;
    }

    //
    // Not found; build a new one and add it to the front
    //
    match = buildOperation(
                dynamic_cast<expert_forest*>(ar.getForest()),
                dynamic_cast<expert_forest*>(res.getForest())
            );
    if (match) {
        match->setNext(cache);
        cache = match;
    }
    return match;
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::getOperation(const dd_edge &ar, opnd_type res)
{
    unary_operation* match = FindUnary(cache, ar, res);
    if (match) {
        cache = match;
        return match;
    }

    //
    // Not found; build a new one and add it to the front
    //
    match = buildOperation(
                dynamic_cast<expert_forest*>(ar.getForest()),
                res
            );
    if (match) {
        match->setNext(cache);
        cache = match;
    }
    return match;
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(expert_forest *, expert_forest *)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(expert_forest *, opnd_type)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


// ******************************************************************
// *                                                                *
// *                     binary_opname  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname::binary_opname(const char* n) : opname(n)
{
}

MEDDLY::binary_opname::~binary_opname()
{
}

MEDDLY::binary_operation*
MEDDLY::binary_opname::getOperation(const dd_edge &a1, const dd_edge &a2,
        const dd_edge &res)
{
    binary_operation* match = nullptr;
    operation* prev = nullptr;
    operation* curr;
    for (curr = cache; curr; curr = curr->getNext()) {
        match = static_cast <binary_operation*> (curr);
        if (match->matches(a1, a2, res)) {
            // Move to front, unless it's already there
            if (prev) {
                prev->setNext(curr->getNext());
                curr->setNext(cache);
                cache = curr;
            }
            return match;
        }
    }

    //
    // Not found; build a new one and add it to the front
    //
    match = buildOperation(
                static_cast<expert_forest*>(a1.getForest()),
                static_cast<expert_forest*>(a2.getForest()),
                static_cast<expert_forest*>(res.getForest())
            );
    match->setNext(cache);
    cache = match;
    return match;
}

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


// ******************************************************************
// *                                                                *
// *                        apply  functions                        *
// *                                                                *
// ******************************************************************

void MEDDLY::apply(unary_handle code, const dd_edge &a, dd_edge &c)
{
    if (nullptr == code) {
        throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }
    unary_opname* opn = code();
    if (nullptr == opn) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
    unary_operation* op = opn->getOperation(a, c);
    op->computeTemp(a, c);
}

void MEDDLY::apply(unary_handle code, const dd_edge &a, long &c)
{
    if (nullptr == code) {
        throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }
    unary_opname* opn = code();
    if (nullptr == opn) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
    unary_operation* op = opn->getOperation(a, opnd_type::INTEGER);
    op->compute(a, c);
}

void MEDDLY::apply(unary_handle code, const dd_edge &a, double &c)
{
    if (nullptr == code) {
        throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }
    unary_opname* opn = code();
    if (nullptr == opn) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
    unary_operation* op = opn->getOperation(a, opnd_type::REAL);
    op->compute(a, c);
}

void MEDDLY::apply(unary_handle code, const dd_edge &a,
    opnd_type cr, ct_object &c)
{
    if (nullptr == code) {
        throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }
    unary_opname* opn = code();
    if (nullptr == opn) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
    unary_operation* op = opn->getOperation(a, cr);
    op->compute(a, c);
}

void MEDDLY::apply(binary_handle code, const dd_edge &a,
    const dd_edge &b, dd_edge &c)
{
    if (nullptr == code) {
        throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }
    binary_opname* opn = code();
    if (nullptr == opn) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
    binary_operation* op = opn->getOperation(a, b, c);
    op->computeTemp(a, b, c);
}


