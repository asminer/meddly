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

#include "../defines.h"
#include "../forest.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"
#include "reach_sat1.h"

/*
namespace MEDDLY {

};
*/

// ******************************************************************
// *                                                                *
// *                  reachset_sat1_factory  class                  *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <bool FWD>
    class reachset_sat1_factory : public binary_factory {
        public:
            virtual void setup();

            virtual binary_operation*
                build_new(forest* a, forest* b, forest* c);
    };
};

// ******************************************************************

template <bool FWD>
void MEDDLY::reachset_sat1_factory <FWD>::setup()
{
    if (FWD) {
        _setup(__FILE__, "REACHABLE_SAT1(true)", "Build forward reachability set using saturation (citation? TBD). The first argument is the set of initial states, and the second argument is the transition relation.");
    } else {
        _setup(__FILE__, "REACHABLE_SAT1(false)", "Build backward reachability set using saturation (citation? TBD). The first argument is the set of initial states, and the second argument is the transition relation.");
    }
}

template <bool FWD>
MEDDLY::binary_operation*
MEDDLY::reachset_sat1_factory <FWD>::build_new(forest* a, forest* b, forest* c)
{
    return nullptr;
    /*
    binary_operation* imageOp = nullptr;
    binary_operation* unionOp = nullptr;

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {

        switch (c->getRangeType()) {
            case range_type::BOOLEAN:
                imageOp = FWD ? MEDDLY::build(POST_IMAGE, a, b, c)
                              : MEDDLY::build(PRE_IMAGE,  a, b, c);

                unionOp = MEDDLY::build(UNION, c, c, c);
                return new reachset_no_frontier(imageOp, unionOp);

            case range_type::INTEGER:
                imageOp = FWD ? MEDDLY::build(POST_IMAGE, a, b, c)
                              : MEDDLY::build(PRE_IMAGE,  a, b, c);

                unionOp = MEDDLY::build(DIST_MIN, c, c, c);
                return new reachset_no_frontier(imageOp, unionOp);

            default:
                return nullptr;
        }
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS) {
        if (c->getRangeType() == range_type::INTEGER) {
            imageOp = FWD ? MEDDLY::build(POST_IMAGE, a, b, c)
                          : MEDDLY::build(PRE_IMAGE,  a, b, c);

            unionOp = MEDDLY::build(MINIMUM, c, c, c);
            return new reachset_no_frontier(imageOp, unionOp);
        }
    }

    return nullptr;
    */
}




// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::REACHABLE_SAT1(bool fwd)
{
    static reachset_sat1_factory<true>  forwd;
    static reachset_sat1_factory<false> bckwd;

    if (fwd) return forwd;
    else     return bckwd;
}

