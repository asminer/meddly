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

#include "meddly/ops_builtin.h"
#include "meddly/oper_unary.h"
#include "meddly/oper_binary.h"

#include "meddly/initializer.h"

namespace MEDDLY {
    class builtin_init;
};

//
// For initializing unary operations
//
#include "meddly/operations/copy.h"
#include "meddly/operations/cardinality.h"
#include "meddly/operations/complement.h"
#include "meddly/operations/maxmin_range.h"
#include "meddly/operations/mdd2index.h"
#include "meddly/operations/cycle.h"

//
// For initializing binary operations
//
#include "meddly/operations/union.h"
#include "meddly/operations/intersection.h"
#include "meddly/operations/difference.h"
#include "meddly/operations/compare.h"

#include "meddly/operations/arith_max.h"
#include "meddly/operations/arith_min.h"
#include "meddly/operations/arith_distmin.h"
#include "meddly/operations/arith_plus.h"
#include "meddly/operations/arith_minus.h"
#include "meddly/operations/arith_mult.h"
#include "meddly/operations/arith_div.h"
#include "meddly/operations/arith_mod.h"

#include "meddly/operations/cross.h"
#include "meddly/operations/prepostplus.h"
#include "meddly/operations/prepost_sets.h"
#include "meddly/operations/prepost_rels.h"
#include "meddly/operations/reach_trad.h"
#include "meddly/operations/vect_matr.h"
#include "meddly/operations/mm_mult.h"
//
// For initializing numerical operations
//
#include "meddly/operations/vect_matr.h"
//
// For initializing saturation-line operations
//
#include "meddly/operations/reach_dfs.h"
#include "meddly/operations/sat_pregen.h"
#include "meddly/operations/sat_otf.h"
#include "meddly/operations/sat_impl.h"
#include "meddly/operations/sat_hyb.h"
#include "meddly/operations/constrained.h"
#include "meddly/operations/transitive_closure.h"



// ******************************************************************
// *                                                                *
// *                      builtin_init methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::builtin_init::builtin_init(initializer_list* p)
    : initializer_list(p)
{
    //
    // Add all unary factories here
    //
    all_unary.push_back( &BOGUS_UNARY()             );
    all_unary.push_back( &CARDINALITY()             );
    all_unary.push_back( &COMPLEMENT()              );
    all_unary.push_back( &CONVERT_TO_INDEX_SET()    );
    all_unary.push_back( &COPY()                    );
    all_unary.push_back( &CYCLE()                   );
    all_unary.push_back( &DIST_INC()                );
    all_unary.push_back( &MAX_RANGE()               );
    all_unary.push_back( &MIN_RANGE()               );

    //
    // Add all binary factories here
    //
    all_binary.push_back( &BOGUS_BINARY()           );
    all_binary.push_back( &CROSS()                  );
    all_binary.push_back( &DIFFERENCE()             );
    all_binary.push_back( &DIST_MIN()               );
    all_binary.push_back( &DIVIDE()                 );
    all_binary.push_back( &EQUAL()                  );
    all_binary.push_back( &GREATER_THAN()           );
    all_binary.push_back( &GREATER_THAN_EQUAL()     );
    all_binary.push_back( &INTERSECTION()           );
    all_binary.push_back( &LESS_THAN()              );
    all_binary.push_back( &LESS_THAN_EQUAL()        );
    all_binary.push_back( &MAXIMUM()                );
    all_binary.push_back( &MINIMUM()                );
    all_binary.push_back( &MINUS()                  );
    all_binary.push_back( &MODULO()                 );
    all_binary.push_back( &MULTIPLY()               );
    all_binary.push_back( &MV_MULTIPLY()            );
    all_binary.push_back( &NOT_EQUAL()              );
    all_binary.push_back( &PLUS()                   );
    all_binary.push_back( &POST_IMAGE()             );
    all_binary.push_back( &PRE_IMAGE()              );

    all_binary.push_back( &REACHABLE_SATUR(true, 1)     );
    all_binary.push_back( &REACHABLE_SATUR(false, 1)    );
    all_binary.push_back( &REACHABLE_SATUR(true, 2)     );
    all_binary.push_back( &REACHABLE_SATUR(false, 2)    );
    all_binary.push_back( &REACHABLE_TRAD_FS(true)      );
    all_binary.push_back( &REACHABLE_TRAD_FS(false)     );
    all_binary.push_back( &REACHABLE_TRAD_NOFS(true)    );
    all_binary.push_back( &REACHABLE_TRAD_NOFS(false)   );

    all_binary.push_back( &UNION()                  );
    all_binary.push_back( &VM_MULTIPLY()            );
}

void MEDDLY::builtin_init::setup()
{
    //
    // Unary ops
    //
    for (unsigned i=0; i<all_unary.size(); i++) {
        if (all_unary[i]) {
            all_unary[i]->setup();
        }
    }
    //
    // Binary ops
    //
    for (unsigned i=0; i<all_binary.size(); i++) {
        if (all_binary[i]) {
            all_binary[i]->setup();
        }
    }

    // OLD BELOW HERE

    //
    // Binary ops
    //

    PRE_PLUS_init();
    POST_PLUS_init();
    TC_POST_IMAGE_init();

#ifdef ALLOW_DEPRECATED_0_18_1
    REACHABLE_STATES_DFS_init();
    REVERSE_REACHABLE_DFS_init();
#endif

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
    for (unsigned i=0; i<all_unary.size(); i++) {
        if (all_unary[i]) {
            all_unary[i]->cleanup();
        }
    }
    //
    // Binary ops
    //
    for (unsigned i=0; i<all_binary.size(); i++) {
        if (all_binary[i]) {
            all_binary[i]->cleanup();
        }
    }

    // OLD BELOW HERE

    //
    // Binary ops
    //
    PRE_PLUS_done();
    POST_PLUS_done();
    TC_POST_IMAGE_done();

#ifdef ALLOW_DEPRECATED_0_18_1
    REACHABLE_STATES_DFS_done();
    REVERSE_REACHABLE_DFS_done();
#endif

    MM_MULTIPLY_done();

    //
    // Ternary ops
    //

    CONSTRAINED_BACKWARD_BFS_done();
    CONSTRAINED_FORWARD_DFS_done();
    CONSTRAINED_BACKWARD_DFS_done();
    TRANSITIVE_CLOSURE_DFS_done();
}

