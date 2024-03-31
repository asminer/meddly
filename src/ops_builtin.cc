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
#include "operations/vect_matr.h"
#include "operations/mm_mult.h"
//
// For initializing numerical operations
//
#include "operations/vect_matr.h"
//
// For initializing saturation-line operations
//
#include "operations/reach_bfs.h"
#include "operations/reach_dfs.h"
#include "operations/sat_pregen.h"
#include "operations/sat_otf.h"
#include "operations/sat_impl.h"
#include "operations/sat_hyb.h"
#include "operations/constrained.h"
#include "operations/transitive_closure.h"



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
//
// Working to remove everything below here
//
    private:
        template <class OPN>
        inline void mydelete(OPN* &p) {
            delete p;
            p = nullptr;
        }
    public:
        // static satpregen_opname* _SATURATION_FORWARD;
        // static satpregen_opname* _SATURATION_BACKWARD;
        static satotf_opname* _SATURATION_OTF_FORWARD;
        static satimpl_opname* _SATURATION_IMPL_FORWARD;
        static sathyb_opname* _SATURATION_HYB_FORWARD;

        static constrained_opname* _CONSTRAINED_BACKWARD_BFS;
        static constrained_opname* _CONSTRAINED_FORWARD_DFS;
        static constrained_opname* _CONSTRAINED_BACKWARD_DFS;
        static constrained_opname* _TRANSITIVE_CLOSURE_DFS;
};

// ******************************************************************

// MEDDLY::satpregen_opname* MEDDLY::builtin_init::_SATURATION_FORWARD;
// MEDDLY::satpregen_opname* MEDDLY::builtin_init::_SATURATION_BACKWARD;
MEDDLY::satotf_opname* MEDDLY::builtin_init::_SATURATION_OTF_FORWARD;
MEDDLY::satimpl_opname* MEDDLY::builtin_init::_SATURATION_IMPL_FORWARD;
MEDDLY::sathyb_opname* MEDDLY::builtin_init::_SATURATION_HYB_FORWARD;

MEDDLY::constrained_opname* MEDDLY::builtin_init::_CONSTRAINED_BACKWARD_BFS;
MEDDLY::constrained_opname* MEDDLY::builtin_init::_CONSTRAINED_FORWARD_DFS;
MEDDLY::constrained_opname* MEDDLY::builtin_init::_CONSTRAINED_BACKWARD_DFS;
MEDDLY::constrained_opname* MEDDLY::builtin_init::_TRANSITIVE_CLOSURE_DFS;

// ******************************************************************

MEDDLY::builtin_init::builtin_init(initializer_list* p)
    : initializer_list(p)
{
    //
    // Saturation-like ops
    //
    // _SATURATION_FORWARD         = nullptr;
    // _SATURATION_BACKWARD        = nullptr;
    _SATURATION_OTF_FORWARD     = nullptr;
    _SATURATION_IMPL_FORWARD    = nullptr;
    _SATURATION_HYB_FORWARD     = nullptr;

    _CONSTRAINED_BACKWARD_BFS   = nullptr;
    _CONSTRAINED_FORWARD_DFS    = nullptr;
    _CONSTRAINED_BACKWARD_DFS   = nullptr;
    _TRANSITIVE_CLOSURE_DFS     = nullptr;
}

void MEDDLY::builtin_init::setup()
{
    //
    // Unary ops
    //
    COPY_init();
    CARDINALITY_init();
    COMPLEMENT_init();
    CONVERT_TO_INDEX_SET_init();
    CYCLE_init();
    MAX_RANGE_init();
    MIN_RANGE_init();
    SELECT_init();

    //
    // Binary ops
    //
    UNION_init();
    INTERSECTION_init();
    DIFFERENCE_init();
    CROSS_init();
    MAXIMUM_init();
    MINIMUM_init();
    PLUS_init();
    MINUS_init();
    MULTIPLY_init();
    DIVIDE_init();
    MODULO_init();

    EQUAL_init();
    NOT_EQUAL_init();
    LESS_THAN_init();
    LESS_THAN_EQUAL_init();
    GREATER_THAN_init();
    GREATER_THAN_EQUAL_init();

    PRE_PLUS_init();
    POST_PLUS_init();
    PRE_IMAGE_init();
    POST_IMAGE_init();
    TC_POST_IMAGE_init();
    REACHABLE_STATES_BFS_init();
    REACHABLE_STATES_DFS_init();
    REVERSE_REACHABLE_BFS_init();
    REVERSE_REACHABLE_DFS_init();

    VM_MULTIPLY_init();
    MV_MULTIPLY_init();
    MM_MULTIPLY_init();

    //
    // Saturation-like ops
    //
    // _SATURATION_FORWARD         =   initSaturationForward()         ;
    // _SATURATION_BACKWARD        =   initSaturationBackward()        ;
    _SATURATION_OTF_FORWARD     =   initOtfSaturationForward()      ;
    _SATURATION_IMPL_FORWARD    =   initImplSaturationForward()     ;
    _SATURATION_HYB_FORWARD     =   initHybSaturationForward()      ;

    _CONSTRAINED_BACKWARD_BFS   =   initConstrainedBFSBackward()    ;
    _CONSTRAINED_FORWARD_DFS    =   initConstrainedDFSForward()     ;
    _CONSTRAINED_BACKWARD_DFS   =   initConstrainedDFSBackward()    ;
    _TRANSITIVE_CLOSURE_DFS     =   initTransitiveClosureDFS()      ;
}

void MEDDLY::builtin_init::cleanup()
{
    //
    // Unary ops
    //
    COPY_done();
    CARDINALITY_done();
    COMPLEMENT_done();
    CONVERT_TO_INDEX_SET_done();
    CYCLE_done();
    MAX_RANGE_done();
    MIN_RANGE_done();
    SELECT_done();

    //
    // Binary ops
    //
    UNION_done();
    INTERSECTION_done();
    DIFFERENCE_done();
    CROSS_done();
    MAXIMUM_done();
    MINIMUM_done();
    PLUS_done();
    MINUS_done();
    MULTIPLY_done();
    DIVIDE_done();
    MODULO_done();

    EQUAL_done();
    NOT_EQUAL_done();
    LESS_THAN_done();
    LESS_THAN_EQUAL_done();
    GREATER_THAN_done();
    GREATER_THAN_EQUAL_done();

    PRE_PLUS_done();
    POST_PLUS_done();
    PRE_IMAGE_done();
    POST_IMAGE_done();
    TC_POST_IMAGE_done();
    REACHABLE_STATES_BFS_done();
    REACHABLE_STATES_DFS_done();
    REVERSE_REACHABLE_BFS_done();
    REVERSE_REACHABLE_DFS_done();

    VM_MULTIPLY_done();
    MV_MULTIPLY_done();
    MM_MULTIPLY_done();

    //
    // Saturation-like ops
    //
    // mydelete(_SATURATION_FORWARD);
    // mydelete(_SATURATION_BACKWARD);
    mydelete(_SATURATION_OTF_FORWARD);
    mydelete(_SATURATION_IMPL_FORWARD);
    mydelete(_SATURATION_HYB_FORWARD);

    mydelete(_CONSTRAINED_BACKWARD_BFS);
    mydelete(_CONSTRAINED_FORWARD_DFS);
    mydelete(_CONSTRAINED_BACKWARD_DFS);
    mydelete(_TRANSITIVE_CLOSURE_DFS);
}


// ******************************************************************
// *                                                                *
// *             front end:  saturation-like operations             *
// *                                                                *
// ******************************************************************

/*
MEDDLY::satpregen_opname* MEDDLY::SATURATION_FORWARD() {
    return builtin_init::_SATURATION_FORWARD;
}
MEDDLY::satpregen_opname* MEDDLY::SATURATION_BACKWARD() {
    return builtin_init::_SATURATION_BACKWARD;
}
*/

MEDDLY::satotf_opname* MEDDLY::SATURATION_OTF_FORWARD() {
    return builtin_init::_SATURATION_OTF_FORWARD;
}
MEDDLY::satimpl_opname* MEDDLY::SATURATION_IMPL_FORWARD() {
    return builtin_init::_SATURATION_IMPL_FORWARD;
}
MEDDLY::sathyb_opname* MEDDLY::SATURATION_HYB_FORWARD() {
    return builtin_init::_SATURATION_HYB_FORWARD;
}
MEDDLY::constrained_opname* MEDDLY::CONSTRAINED_BACKWARD_BFS() {
    return builtin_init::_CONSTRAINED_BACKWARD_BFS;
}
MEDDLY::constrained_opname* MEDDLY::CONSTRAINED_FORWARD_DFS() {
    return builtin_init::_CONSTRAINED_FORWARD_DFS;
}
MEDDLY::constrained_opname* MEDDLY::CONSTRAINED_BACKWARD_DFS() {
    return builtin_init::_CONSTRAINED_BACKWARD_DFS;
}
MEDDLY::constrained_opname* MEDDLY::TRANSITIVE_CLOSURE_DFS() {
    return builtin_init::_TRANSITIVE_CLOSURE_DFS;
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

