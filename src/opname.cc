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

#include "initializer.h"

#include "oper.h"
#include "oper_unary.h"
#include "oper_binary.h"

// These are needed for builtin opname initialization

//
// Unary operations
//
#include "operations/copy.h"
#include "operations/cardinality.h"
#include "operations/complement.h"
#include "operations/maxmin_range.h"
#include "operations/mdd2index.h"
#include "operations/cycle.h"
#include "operations/select.h"


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
// *                       opname_init  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::opname_init : public initializer_list {
    public:
        opname_init(initializer_list* p);
    protected:
        virtual void setup();
        virtual void cleanup();
    private:
        inline void mydelete(unary_opname* &p) {
            delete p;
            p = nullptr;
        }
};

MEDDLY::opname_init::opname_init(initializer_list* p)
    : initializer_list(p)
{
    // opname::next_index = 0;

    opname::_COPY       = nullptr;
    opname::_CARD       = nullptr;
    opname::_COMPL      = nullptr;
    opname::_MAXRANGE   = nullptr;
    opname::_MINRANGE   = nullptr;
    opname::_MDD2INDEX  = nullptr;
    opname::_CYCLE      = nullptr;
    opname::_SELECT     = nullptr;
}

void MEDDLY::opname_init::setup()
{
    //
    // Unary ops
    //
    opname::_COPY       = initializeCopy();
    opname::_CARD       = initializeCardinality();
    opname::_COMPL      = initializeComplement();
    opname::_MAXRANGE   = initializeMaxRange();
    opname::_MINRANGE   = initializeMaxRange();
    opname::_MDD2INDEX  = initializeMDD2INDEX();
    opname::_CYCLE      = initializeCycle();
    opname::_SELECT     = initializeSelect();

}

void MEDDLY::opname_init::cleanup()
{
    mydelete(opname::_COPY);
    mydelete(opname::_CARD);
    mydelete(opname::_COMPL);
    mydelete(opname::_MAXRANGE);
    mydelete(opname::_MINRANGE);
    mydelete(opname::_MDD2INDEX);
    mydelete(opname::_CYCLE);
    mydelete(opname::_SELECT);

    // opname::next_index = 0;
}

// ******************************************************************
// *                                                                *
// *                         opname methods                         *
// *                                                                *
// ******************************************************************

// size_t MEDDLY::opname::next_index;

MEDDLY::unary_opname* MEDDLY::opname::_COPY;
MEDDLY::unary_opname* MEDDLY::opname::_CARD;
MEDDLY::unary_opname* MEDDLY::opname::_COMPL;
MEDDLY::unary_opname* MEDDLY::opname::_MAXRANGE;
MEDDLY::unary_opname* MEDDLY::opname::_MINRANGE;
MEDDLY::unary_opname* MEDDLY::opname::_MDD2INDEX;
MEDDLY::unary_opname* MEDDLY::opname::_CYCLE;
MEDDLY::unary_opname* MEDDLY::opname::_SELECT;

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

MEDDLY::initializer_list*
MEDDLY::opname::makeInitializer(initializer_list* prev)
{
    return new opname_init(prev);
}

void
MEDDLY::opname::removeFromCache(operation* op)
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
    match = buildOperation(ar, res);
    match->setNext(cache);
    cache = match;
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
    match = buildOperation(ar, res);
    match->setNext(cache);
    cache = match;
    return match;
}

void MEDDLY::unary_opname::removeOperationFromCache(unary_operation* op)
{
    removeFromCache(op);
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(const dd_edge &, const dd_edge &) const
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(const dd_edge &, opnd_type) const
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
    match = buildOperation(a1, a2, res);
    match->setNext(cache);
    cache = match;
    return match;
}

void MEDDLY::binary_opname::removeOperationFromCache(binary_operation* op)
{
    removeFromCache(op);
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
// *                  front end:  unary operations                  *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::COPY()
{
    return opname::COPY();
}

MEDDLY::unary_opname* MEDDLY::CARDINALITY()
{
    return opname::CARDINALITY();
}

MEDDLY::unary_opname* MEDDLY::COMPLEMENT()
{
    return opname::COMPLEMENT();
}

MEDDLY::unary_opname* MEDDLY::MAX_RANGE()
{
    return opname::MAX_RANGE();
}

MEDDLY::unary_opname* MEDDLY::MIN_RANGE()
{
    return opname::MIN_RANGE();
}

MEDDLY::unary_opname* MEDDLY::CONVERT_TO_INDEX_SET()
{
    return opname::CONVERT_TO_INDEX_SET();
}

MEDDLY::unary_opname* MEDDLY::CYCLE()
{
    return opname::CYCLE();
}

MEDDLY::unary_opname* MEDDLY::SELECT()
{
    return opname::SELECT();
}

// ******************************************************************
// *                                                                *
// *                        apply  functions                        *
// *                                                                *
// ******************************************************************

void MEDDLY::apply(unary_opname* (*code)(), const dd_edge &a, dd_edge &c)
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

void MEDDLY::apply(unary_opname* (*code)(), const dd_edge &a, long &c)
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

void MEDDLY::apply(unary_opname* (*code)(), const dd_edge &a, double &c)
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

void MEDDLY::apply(unary_opname* (*code)(), const dd_edge &a,
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

void MEDDLY::apply(binary_opname* (*code)(), const dd_edge &a,
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


