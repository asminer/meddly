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

#include "ops_builtin.h"

#include "opname.h"
#include "opname_numer.h"
#include "opname_satur.h"

#include "initializer.h"

namespace MEDDLY {
    class builtin_init;
};

//
// For initializing unary operations
//
#include "operations/copy.h"
#include "operations/cardinality.h"
#include "operations/complement.h"
#include "operations/maxmin_range.h"
#include "operations/mdd2index.h"
#include "operations/cycle.h"
#include "operations/select.h"
//
// For initializing binary operations
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
//
// For initializing numerical operations
//
#include "operations/vect_matr.h"



// ******************************************************************
// *                                                                *
// *                       builtin_init class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::builtin_init : public initializer_list {
    public:
        builtin_init(initializer_list* p);
    protected:
        virtual void setup();
        virtual void cleanup();
    private:
        template <class OPN>
        inline void mydelete(OPN* &p) {
            delete p;
            p = nullptr;
        }
    public:
        static unary_opname* _CARD;
        static unary_opname* _COMPL;
        static unary_opname* _MDD2INDEX;
        static unary_opname* _COPY;
        static unary_opname* _CYCLE;
        static unary_opname* _MAXRANGE;
        static unary_opname* _MINRANGE;
        static unary_opname* _SELECT;
    public:
        static binary_opname* _UNION;
        static binary_opname* _INTERSECT;
        static binary_opname* _DIFFERENCE;
        static binary_opname* _CROSS;
        static binary_opname* _MIN;
        static binary_opname* _MAX;
        static binary_opname* _PLUS;
        static binary_opname* _MINUS;
        static binary_opname* _MULTIPLY;
        static binary_opname* _DIVIDE;
        static binary_opname* _MODULO;
        static binary_opname* _EQ;
        static binary_opname* _NE;
        static binary_opname* _LT;
        static binary_opname* _LE;
        static binary_opname* _GT;
        static binary_opname* _GE;
        static binary_opname* _PRE_PLUS;
        static binary_opname* _POST_PLUS;
        static binary_opname* _PRE_IMAGE;
        static binary_opname* _POST_IMAGE;
        static binary_opname* _TC_POST_IMAGE;
        static binary_opname* _FORWARD_DFS;
        static binary_opname* _FORWARD_BFS;
        static binary_opname* _BACKWARD_DFS;
        static binary_opname* _BACKWARD_BFS;
        static binary_opname* _VM_MULTIPLY;
        static binary_opname* _MV_MULTIPLY;
        static binary_opname* _MM_MULTIPLY;
    public:
        static numerical_opname* _EXPLVECT_MATR_MULT;
        static numerical_opname* _MATR_EXPLVECT_MULT;
};

// ******************************************************************

MEDDLY::builtin_init::builtin_init(initializer_list* p)
    : initializer_list(p)
{
    //
    // Unary ops
    //
    _COPY       = nullptr;
    _CARD       = nullptr;
    _COMPL      = nullptr;
    _MAXRANGE   = nullptr;
    _MINRANGE   = nullptr;
    _MDD2INDEX  = nullptr;
    _CYCLE      = nullptr;
    _SELECT     = nullptr;
    //
    // Binary ops
    //
    _UNION          = nullptr;
    _INTERSECT      = nullptr;
    _DIFFERENCE     = nullptr;
    _CROSS          = nullptr;
    _MIN            = nullptr;
    _MAX            = nullptr;
    _PLUS           = nullptr;
    _MINUS          = nullptr;
    _MULTIPLY       = nullptr;
    _DIVIDE         = nullptr;
    _MODULO         = nullptr;
    _EQ             = nullptr;
    _NE             = nullptr;
    _LT             = nullptr;
    _LE             = nullptr;
    _GT             = nullptr;
    _GE             = nullptr;
    _PRE_PLUS       = nullptr;
    _POST_PLUS      = nullptr;
    _PRE_IMAGE      = nullptr;
    _POST_IMAGE     = nullptr;
    _TC_POST_IMAGE  = nullptr;
    _FORWARD_DFS    = nullptr;
    _FORWARD_BFS    = nullptr;
    _BACKWARD_DFS   = nullptr;
    _BACKWARD_BFS   = nullptr;
    _VM_MULTIPLY    = nullptr;
    _MV_MULTIPLY    = nullptr;
    _MM_MULTIPLY    = nullptr;
    //
    // Numerical ops
    //
    _EXPLVECT_MATR_MULT = nullptr;
    _MATR_EXPLVECT_MULT = nullptr;
}

void MEDDLY::builtin_init::setup()
{
    //
    // Unary ops
    //
    _COPY       =   initializeCopy()        ;
    _CARD       =   initializeCardinality() ;
    _COMPL      =   initializeComplement()  ;
    _MAXRANGE   =   initializeMaxRange()    ;
    _MINRANGE   =   initializeMaxRange()    ;
    _MDD2INDEX  =   initializeMDD2INDEX()   ;
    _CYCLE      =   initializeCycle()       ;
    _SELECT     =   initializeSelect()      ;
    //
    // Binary ops
    //
    _UNION          =   initializeUnion()           ;
    _INTERSECT      =   initializeIntersection()    ;
    _DIFFERENCE     =   initializeDifference()      ;
    _CROSS          =   initializeCross()           ;
    _MAX            =   initializeMaximum()         ;
    _MIN            =   initializeMinimum()         ;
    _PLUS           =   initializePlus()            ;
    _MINUS          =   initializeMinus()           ;
    _MULTIPLY       =   initializeMultiply()        ;
    _DIVIDE         =   initializeDivide()          ;
    _MODULO         =   initializeModulo()          ;
    _EQ             =   initializeEQ()              ;
    _NE             =   initializeNE()              ;
    _LT             =   initializeLT()              ;
    _LE             =   initializeLE()              ;
    _GT             =   initializeGT()              ;
    _GE             =   initializeGE()              ;
    _PRE_PLUS       =   initializePrePlus()         ;
    _POST_PLUS      =   initializePostPlus()        ;
    _PRE_IMAGE      =   initializePreImage()        ;
    _POST_IMAGE     =   initializePostImage()       ;
    _TC_POST_IMAGE  =   initializeTCPostImage()     ;
    _FORWARD_DFS    =   initializeForwardDFS()      ;
    _FORWARD_BFS    =   initializeForwardBFS()      ;
    _BACKWARD_DFS   =   initializeBackwardDFS()     ;
    _BACKWARD_BFS   =   initializeBackwardBFS()     ;
    _VM_MULTIPLY    =   initializeVMmult()          ;
    _MV_MULTIPLY    =   initializeMVmult()          ;
    _MM_MULTIPLY    =   initializeMMMultiply()      ;
    //
    // Numerical ops
    //
    _EXPLVECT_MATR_MULT = initExplVectorMatrixMult()  ;
    _MATR_EXPLVECT_MULT = initMatrixExplVectorMult()  ;
}

void MEDDLY::builtin_init::cleanup()
{
    //
    // Unary ops
    //
    mydelete(_COPY);
    mydelete(_CARD);
    mydelete(_COMPL);
    mydelete(_MAXRANGE);
    mydelete(_MINRANGE);
    mydelete(_MDD2INDEX);
    mydelete(_CYCLE);
    mydelete(_SELECT);
    //
    // Binary ops
    //
    mydelete(_UNION);
    mydelete(_INTERSECT);
    mydelete(_DIFFERENCE);
    mydelete(_CROSS);
    mydelete(_MIN);
    mydelete(_MAX);
    mydelete(_PLUS);
    mydelete(_MINUS);
    mydelete(_MULTIPLY);
    mydelete(_DIVIDE);
    mydelete(_MODULO);
    mydelete(_EQ);
    mydelete(_NE);
    mydelete(_LT);
    mydelete(_LE);
    mydelete(_GT);
    mydelete(_GE);
    mydelete(_PRE_PLUS);
    mydelete(_POST_PLUS);
    mydelete(_PRE_IMAGE);
    mydelete(_POST_IMAGE);
    mydelete(_TC_POST_IMAGE);
    mydelete(_FORWARD_DFS);
    mydelete(_FORWARD_BFS);
    mydelete(_BACKWARD_DFS);
    mydelete(_BACKWARD_BFS);
    mydelete(_VM_MULTIPLY);
    mydelete(_MV_MULTIPLY);
    mydelete(_MM_MULTIPLY);
    //
    // Numerical ops
    //
    mydelete(_EXPLVECT_MATR_MULT);
    mydelete(_MATR_EXPLVECT_MULT);
}

// ******************************************************************
// *                                                                *
// *                  front end:  unary operations                  *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::CARDINALITY() {
    return builtin_init::_CARD;
}
MEDDLY::unary_opname* MEDDLY::COMPLEMENT() {
    return builtin_init::_COMPL;
}
MEDDLY::unary_opname* MEDDLY::CONVERT_TO_INDEX_SET() {
    return builtin_init::_MDD2INDEX;
}
MEDDLY::unary_opname* MEDDLY::COPY() {
    return builtin_init::_COPY;
}
MEDDLY::unary_opname* MEDDLY::CYCLE() {
    return builtin_init::_CYCLE;
}
MEDDLY::unary_opname* MEDDLY::MAX_RANGE() {
    return builtin_init::_MAXRANGE;
}
MEDDLY::unary_opname* MEDDLY::MIN_RANGE() {
    return builtin_init::_MINRANGE;
}
MEDDLY::unary_opname* MEDDLY::SELECT() {
    return builtin_init::_SELECT;
}

// ******************************************************************
// *                                                                *
// *                  front end: binary operations                  *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::UNION() {
    return builtin_init::_UNION;
}
MEDDLY::binary_opname* MEDDLY::INTERSECTION() {
    return builtin_init::_INTERSECT;
}
MEDDLY::binary_opname* MEDDLY::DIFFERENCE() {
    return builtin_init::_DIFFERENCE;
}
MEDDLY::binary_opname* MEDDLY::CROSS() {
    return builtin_init::_CROSS;
}
MEDDLY::binary_opname* MEDDLY::MINIMUM() {
    return builtin_init::_MIN;
}
MEDDLY::binary_opname* MEDDLY::MAXIMUM() {
    return builtin_init::_MAX;
}
MEDDLY::binary_opname* MEDDLY::PLUS() {
    return builtin_init::_PLUS;
}
MEDDLY::binary_opname* MEDDLY::MINUS() {
    return builtin_init::_MINUS;
}
MEDDLY::binary_opname* MEDDLY::MULTIPLY() {
    return builtin_init::_MULTIPLY;
}
MEDDLY::binary_opname* MEDDLY::DIVIDE() {
    return builtin_init::_DIVIDE;
}
MEDDLY::binary_opname* MEDDLY::MODULO() {
    return builtin_init::_MODULO;
}
MEDDLY::binary_opname* MEDDLY::EQUAL() {
    return builtin_init::_EQ;
}
MEDDLY::binary_opname* MEDDLY::NOT_EQUAL() {
    return builtin_init::_NE;
}
MEDDLY::binary_opname* MEDDLY::LESS_THAN() {
    return builtin_init::_LT;
}
MEDDLY::binary_opname* MEDDLY::LESS_THAN_EQUAL() {
    return builtin_init::_LE;
}
MEDDLY::binary_opname* MEDDLY::GREATER_THAN() {
    return builtin_init::_GT;
}
MEDDLY::binary_opname* MEDDLY::GREATER_THAN_EQUAL() {
    return builtin_init::_GE;
}
MEDDLY::binary_opname* MEDDLY::PRE_PLUS() {
    return builtin_init::_PRE_PLUS;
}
MEDDLY::binary_opname* MEDDLY::POST_PLUS() {
    return builtin_init::_POST_PLUS;
}
MEDDLY::binary_opname* MEDDLY::PRE_IMAGE() {
    return builtin_init::_PRE_IMAGE;
}
MEDDLY::binary_opname* MEDDLY::POST_IMAGE() {
    return builtin_init::_POST_IMAGE;
}
MEDDLY::binary_opname* MEDDLY::TC_POST_IMAGE() {
    return builtin_init::_TC_POST_IMAGE;
}
MEDDLY::binary_opname* MEDDLY::REACHABLE_STATES_DFS() {
    return builtin_init::_FORWARD_DFS;
}
MEDDLY::binary_opname* MEDDLY::REACHABLE_STATES_BFS() {
    return builtin_init::_FORWARD_BFS;
}
MEDDLY::binary_opname* MEDDLY::REVERSE_REACHABLE_DFS() {
    return builtin_init::_BACKWARD_DFS;
}
MEDDLY::binary_opname* MEDDLY::REVERSE_REACHABLE_BFS() {
    return builtin_init::_BACKWARD_BFS;
}
MEDDLY::binary_opname* MEDDLY::VM_MULTIPLY() {
    return builtin_init::_VM_MULTIPLY;
}
MEDDLY::binary_opname* MEDDLY::MV_MULTIPLY() {
    return builtin_init::_MV_MULTIPLY;
}
MEDDLY::binary_opname* MEDDLY::MM_MULTIPLY() {
    return builtin_init::_MM_MULTIPLY;
}

// ******************************************************************
// *                                                                *
// *                front end:  numerical operations                *
// *                                                                *
// ******************************************************************

MEDDLY::numerical_opname* MEDDLY::EXPLVECT_MATR_MULT() {
    return builtin_init::_EXPLVECT_MATR_MULT;
}
MEDDLY::numerical_opname* MEDDLY::MATR_EXPLVECT_MULT() {
    return builtin_init::_MATR_EXPLVECT_MULT;
}

// ******************************************************************
// *                                                                *
// *                   front end:  initialization                   *
// *                                                                *
// ******************************************************************

MEDDLY::initializer_list*
MEDDLY::makeBuiltinInitializer(initializer_list* prev)
{
    return new builtin_init(prev);
}

