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

#ifdef ALLOW_DEPRECATED_0_17_8
#include "operations/select.h"
#endif

//
// For initializing binary operations
//
#include "operations/union.h"
#include "operations/intersection.h"
#include "operations/difference.h"
#include "operations/compare.h"

#include "operations/arith_max.h"
#include "operations/arith_min.h"
#include "operations/arith_plus.h"
#include "operations/arith_minus.h"
#include "operations/arith_mult.h"
#include "operations/arith_div.h"
#include "operations/arith_mod.h"

#include "operations/cross.h"
// #include "operations/maxmin.h"
// #include "operations/plus.h"
// #include "operations/minus.h"
// #include "operations/multiply.h"
// #include "operations/divide.h"
// #include "operations/modulo.h"
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
// *                      builtin_init methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::builtin_init::builtin_init(initializer_list* p)
    : initializer_list(p)
{
    //
    // Add all unary factories
    //
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
#ifdef ALLOW_DEPRECATED_0_17_8
    SELECT_init();
#endif

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

    CONSTRAINED_BACKWARD_BFS_init();
    CONSTRAINED_FORWARD_DFS_init();
    CONSTRAINED_BACKWARD_DFS_init();
    TRANSITIVE_CLOSURE_DFS_init();
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
#ifdef ALLOW_DEPRECATED_0_17_8
    SELECT_done();
#endif

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
    // Ternary ops
    //

    CONSTRAINED_BACKWARD_BFS_done();
    CONSTRAINED_FORWARD_DFS_done();
    CONSTRAINED_BACKWARD_DFS_done();
    TRANSITIVE_CLOSURE_DFS_done();
}

