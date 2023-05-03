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
//
// Binary operations
//
#include "operations/union.h"
#include "operations/intersection.h"
#include "operations/difference.h"
#include "operations/cross.h"
#include "operations/maxmin.h"
#include "operations/plus.h"
#include "operations/minus.h"
#include "operations/multiply.h"
#include "operations/divide.h"
#include "operations/modulo.h"
#include "operations/comp_eq.h"
#include "operations/comp_ne.h"
#include "operations/comp_lt.h"
#include "operations/comp_le.h"
#include "operations/comp_gt.h"
#include "operations/comp_ge.h"
#include "operations/prepostplus.h"
#include "operations/prepostimage.h"
#include "operations/reach_bfs.h"
#include "operations/reach_dfs.h"
#include "operations/vect_matr.h"
#include "operations/mm_mult.h"

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
        inline void mydelete(binary_opname* &p) {
            delete p;
            p = nullptr;
        }
};

MEDDLY::opname_init::opname_init(initializer_list* p)
    : initializer_list(p)
{
    //
    // Unary ops
    //
    opname::_COPY       = nullptr;
    opname::_CARD       = nullptr;
    opname::_COMPL      = nullptr;
    opname::_MAXRANGE   = nullptr;
    opname::_MINRANGE   = nullptr;
    opname::_MDD2INDEX  = nullptr;
    opname::_CYCLE      = nullptr;
    opname::_SELECT     = nullptr;
    //
    // Binary ops
    //
    opname::_UNION          = nullptr;
    opname::_INTERSECT      = nullptr;
    opname::_DIFFERENCE     = nullptr;
    opname::_CROSS          = nullptr;
    opname::_MIN            = nullptr;
    opname::_MAX            = nullptr;
    opname::_PLUS           = nullptr;
    opname::_MINUS          = nullptr;
    opname::_MULTIPLY       = nullptr;
    opname::_DIVIDE         = nullptr;
    opname::_MODULO         = nullptr;
    opname::_EQ             = nullptr;
    opname::_NE             = nullptr;
    opname::_LT             = nullptr;
    opname::_LE             = nullptr;
    opname::_GT             = nullptr;
    opname::_GE             = nullptr;
    opname::_PRE_PLUS       = nullptr;
    opname::_POST_PLUS      = nullptr;
    opname::_PRE_IMAGE      = nullptr;
    opname::_POST_IMAGE     = nullptr;
    opname::_TC_POST_IMAGE  = nullptr;
    opname::_FORWARD_DFS    = nullptr;
    opname::_FORWARD_BFS    = nullptr;
    opname::_BACKWARD_DFS   = nullptr;
    opname::_BACKWARD_BFS   = nullptr;
    opname::_VM_MULTIPLY    = nullptr;
    opname::_MV_MULTIPLY    = nullptr;
    opname::_MM_MULTIPLY    = nullptr;
}

void MEDDLY::opname_init::setup()
{
    //
    // Unary ops
    //
    opname::_COPY       =   initializeCopy()        ;
    opname::_CARD       =   initializeCardinality() ;
    opname::_COMPL      =   initializeComplement()  ;
    opname::_MAXRANGE   =   initializeMaxRange()    ;
    opname::_MINRANGE   =   initializeMaxRange()    ;
    opname::_MDD2INDEX  =   initializeMDD2INDEX()   ;
    opname::_CYCLE      =   initializeCycle()       ;
    opname::_SELECT     =   initializeSelect()      ;
    //
    // Binary ops
    //
    opname::_UNION          =   initializeUnion()           ;
    opname::_INTERSECT      =   initializeIntersection()    ;
    opname::_DIFFERENCE     =   initializeDifference()      ;
    opname::_CROSS          =   initializeCross()           ;
    opname::_MAX            =   initializeMaximum()         ;
    opname::_MIN            =   initializeMinimum()         ;
    opname::_PLUS           =   initializePlus()            ;
    opname::_MINUS          =   initializeMinus()           ;
    opname::_MULTIPLY       =   initializeMultiply()        ;
    opname::_DIVIDE         =   initializeDivide()          ;
    opname::_MODULO         =   initializeModulo()          ;
    opname::_EQ             =   initializeEQ()              ;
    opname::_NE             =   initializeNE()              ;
    opname::_LT             =   initializeLT()              ;
    opname::_LE             =   initializeLE()              ;
    opname::_GT             =   initializeGT()              ;
    opname::_GE             =   initializeGE()              ;
    opname::_PRE_PLUS       =   initializePrePlus()         ;
    opname::_POST_PLUS      =   initializePostPlus()        ;
    opname::_PRE_IMAGE      =   initializePreImage()        ;
    opname::_POST_IMAGE     =   initializePostImage()       ;
    opname::_TC_POST_IMAGE  =   initializeTCPostImage()     ;
    opname::_FORWARD_DFS    =   initializeForwardDFS()      ;
    opname::_FORWARD_BFS    =   initializeForwardBFS()      ;
    opname::_BACKWARD_DFS   =   initializeBackwardDFS()     ;
    opname::_BACKWARD_BFS   =   initializeBackwardBFS()     ;
    opname::_VM_MULTIPLY    =   initializeVMmult()          ;
    opname::_MV_MULTIPLY    =   initializeMVmult()          ;
    opname::_MM_MULTIPLY    =   initializeMMMultiply()      ;
}

void MEDDLY::opname_init::cleanup()
{
    //
    // Unary ops
    //
    mydelete(opname::_COPY);
    mydelete(opname::_CARD);
    mydelete(opname::_COMPL);
    mydelete(opname::_MAXRANGE);
    mydelete(opname::_MINRANGE);
    mydelete(opname::_MDD2INDEX);
    mydelete(opname::_CYCLE);
    mydelete(opname::_SELECT);
    //
    // Binary ops
    //
    mydelete(opname::_UNION);
    mydelete(opname::_INTERSECT);
    mydelete(opname::_DIFFERENCE);
    mydelete(opname::_CROSS);
    mydelete(opname::_MIN);
    mydelete(opname::_MAX);
    mydelete(opname::_PLUS);
    mydelete(opname::_MINUS);
    mydelete(opname::_MULTIPLY);
    mydelete(opname::_DIVIDE);
    mydelete(opname::_MODULO);
    mydelete(opname::_EQ);
    mydelete(opname::_NE);
    mydelete(opname::_LT);
    mydelete(opname::_LE);
    mydelete(opname::_GT);
    mydelete(opname::_GE);
    mydelete(opname::_PRE_PLUS);
    mydelete(opname::_POST_PLUS);
    mydelete(opname::_PRE_IMAGE);
    mydelete(opname::_POST_IMAGE);
    mydelete(opname::_TC_POST_IMAGE);
    mydelete(opname::_FORWARD_DFS);
    mydelete(opname::_FORWARD_BFS);
    mydelete(opname::_BACKWARD_DFS);
    mydelete(opname::_BACKWARD_BFS);
    mydelete(opname::_VM_MULTIPLY);
    mydelete(opname::_MV_MULTIPLY);
    mydelete(opname::_MM_MULTIPLY);
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
MEDDLY::unary_opname::buildOperation(const dd_edge &, const dd_edge &)
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(const dd_edge &, opnd_type)
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
// *                  front end: binary operations                  *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::UNION()
{
    return opname::UNION();
}

MEDDLY::binary_opname* MEDDLY::INTERSECTION()
{
    return opname::INTERSECTION();
}

MEDDLY::binary_opname* MEDDLY::DIFFERENCE()
{
    return opname::DIFFERENCE();
}

MEDDLY::binary_opname* MEDDLY::CROSS()
{
    return opname::CROSS();
}

MEDDLY::binary_opname* MEDDLY::MINIMUM()
{
    return opname::MINIMUM();
}

MEDDLY::binary_opname* MEDDLY::MAXIMUM()
{
    return opname::MAXIMUM();
}

MEDDLY::binary_opname* MEDDLY::PLUS()
{
    return opname::PLUS();
}

MEDDLY::binary_opname* MEDDLY::MINUS()
{
    return opname::MINUS();
}

MEDDLY::binary_opname* MEDDLY::MULTIPLY()
{
    return opname::MULTIPLY();
}

MEDDLY::binary_opname* MEDDLY::DIVIDE()
{
    return opname::DIVIDE();
}

MEDDLY::binary_opname* MEDDLY::MODULO()
{
    return opname::MODULO();
}

MEDDLY::binary_opname* MEDDLY::EQUAL()
{
    return opname::EQUAL();
}

MEDDLY::binary_opname* MEDDLY::NOT_EQUAL()
{
    return opname::NOT_EQUAL();
}

MEDDLY::binary_opname* MEDDLY::LESS_THAN()
{
    return opname::LESS_THAN();
}

MEDDLY::binary_opname* MEDDLY::LESS_THAN_EQUAL()
{
    return opname::LESS_THAN_EQUAL();
}

MEDDLY::binary_opname* MEDDLY::GREATER_THAN()
{
    return opname::GREATER_THAN();
}

MEDDLY::binary_opname* MEDDLY::GREATER_THAN_EQUAL()
{
    return opname::GREATER_THAN_EQUAL();
}

MEDDLY::binary_opname* MEDDLY::PRE_PLUS()
{
    return opname::PRE_PLUS();
}

MEDDLY::binary_opname* MEDDLY::POST_PLUS()
{
    return opname::POST_PLUS();
}

MEDDLY::binary_opname* MEDDLY::PRE_IMAGE()
{
    return opname::PRE_IMAGE();
}

MEDDLY::binary_opname* MEDDLY::POST_IMAGE()
{
    return opname::POST_IMAGE();
}

MEDDLY::binary_opname* MEDDLY::TC_POST_IMAGE()
{
    return opname::TC_POST_IMAGE();
}

MEDDLY::binary_opname* MEDDLY::REACHABLE_STATES_DFS()
{
    return opname::REACHABLE_STATES_DFS();
}

MEDDLY::binary_opname* MEDDLY::REACHABLE_STATES_BFS()
{
    return opname::REACHABLE_STATES_BFS();
}

MEDDLY::binary_opname* MEDDLY::REVERSE_REACHABLE_DFS()
{
    return opname::REVERSE_REACHABLE_DFS();
}

MEDDLY::binary_opname* MEDDLY::REVERSE_REACHABLE_BFS()
{
    return opname::REVERSE_REACHABLE_BFS();
}

MEDDLY::binary_opname* MEDDLY::VM_MULTIPLY()
{
    return opname::VM_MULTIPLY();
}

MEDDLY::binary_opname* MEDDLY::MV_MULTIPLY()
{
    return opname::MV_MULTIPLY();
}

MEDDLY::binary_opname* MEDDLY::MM_MULTIPLY()
{
    return opname::MM_MULTIPLY();
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


