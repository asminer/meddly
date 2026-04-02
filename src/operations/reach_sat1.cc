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

// ************************************************************************
// ************************************************************************

/*
    Template class for saturation (version 1) operations,
    when the relation is MT-based.

    Template parameters:

        EOP: one of the EdgeOp classes in forests_edgerules.h.

        FORWD: if true, do forward reachability, otherwise
                        do backward reachability.

        ATYPE: arithmetic type class, must provide the following methods.

            /// Get the operation name, for display purposes
            static const char* name(bool forwd);

            /// Return true if the given edge is unreachable
            static bool isUnreachable(const edge_value &cv, node_handle c);

            /// Set an edge to be unreachable
            static void setUnreachable(edge_value &cv, node_handle &c);

            /// Set all edges of an unpacked node to unreachable
            static void setAllUnreachable(unpacked_node *u);

            /// Get the accumulate operation for result nodes.
            static binary_operation* accumulateOp(const forest* resF);

            /// Apply the operation when b is a terminal node.
            static void apply(
                const forest* fa, const edge_value &av, node_handle a,
                const forest* fb, node_handle b,
                const forest* fc, edge_value &cv, node_handle &c
            );

*/
namespace MEDDLY {

    template <class EOP, bool FORWD, class ATYPE>
    class saturation1_set_mtrel : public binary_operation {
        public:
            saturation1_set_mtrel(forest* arg1, forest* arg2,
                    forest* res);

            virtual ~saturation1_set_mtrel();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    const edge_value &bv, node_handle bp,
                    edge_value &cv, node_handle &cp);

        protected:
            void _compute(int L, const edge_value &av, node_handle A,
                    node_handle B, edge_value &cv, node_handle &C);

        private:
            inline const edge_value &edgeval(unpacked_node *U, unsigned i) const
            {
                if (EOP::hasEdgeValues()) {
                    return U->edgeval(i);
                } else {
                    return nothing;
                }
            }

            //
            // Correctly do C[i] = C[i] + <v, p>
            //
            inline void addToCi(int nextL, unpacked_node *C, unsigned i,
                    const edge_value &v, node_handle p)
            {
                // This case should be caught already
                MEDDLY_DCASSERT( !ATYPE::isUnreachable(v, p) );
                if (ATYPE::isUnreachable(edgeval(C, i), C->down(i)))
                {
                    C->setFull(i, v, p);
                    return;
                }
                edge_value  newdv;
                node_handle newdp;
                accumulateOp->compute(nextL, ~0,
                    edgeval(C, i), C->down(i), v, p, newdv, newdp
                );
                resF->unlinkNode(p);
                resF->unlinkNode(C->down(i));
                C->setFull(i, newdv, newdp);
            }

        private:
            ct_entry_type* ct;
            binary_operation* accumulateOp;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            edge_value nothing;
            bool forced_by_levels;

    }; // class
}; // namespace MEDDLY


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

