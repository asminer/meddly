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
#include "reach_sat1.h"

#include "../forest.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"
#include "../ct_vector.h"
#include "../forest_levels.h"
#include "../forest_edgerules.h"

#include "prepost_common.h"

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


            static inline const char* opName() {
                return FORWD ? "fwd-sat1" : "bck-sat1";
            }

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

// ************************************************************************

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::saturation1_set_mtrel<EOP, FORWD, ATYPE>
    ::saturation1_set_mtrel(forest* arg1, forest* arg2, forest* res)
    : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
    checkLabelings(__FILE__, __LINE__,
        res->getEdgeLabeling(),
        edge_labeling::MULTI_TERMINAL,
        res->getEdgeLabeling()
    );

    if (arg1F->getRangeType() != resF->getRangeType()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    //
    // Addition operation for the vector-matrix multiply.
    //  union for pre/post image, addition otherwise.
    //
    accumulateOp = ATYPE::accumulateOp(res);
    MEDDLY_DCASSERT(accumulateOp);

    //
    // Do we need to recurse by levels and store level info in the CT?
    // YES, if the set and relation are both fully-reduced.
    // (If the set is quasi reduced, we will recurse by levels anyway.)
    // (If the relation is identity-reduced, we can skip levels.)
    //
    forced_by_levels = arg1->isFullyReduced() && arg2->isFullyReduced();

    // TBD

    //
    // Build compute table key and result types.
    //
    ct = new ct_entry_type(opName());
    if (forced_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    if (EOP::hasEdgeValues()) {
        ct->setResult(EOP::edgeValueTypeLetter(), res);
    } else {
        ct->setResult(res);
    }
    ct->doneBuilding();
}

template <class EOP, bool FORWD, class ATYPE>
MEDDLY::saturation1_set_mtrel<EOP, FORWD, ATYPE>::~saturation1_set_mtrel()
{
    ct->markForDestroy();
}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::saturation1_set_mtrel<EOP, FORWD, ATYPE>
    ::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << opName() << " #" << top_count << " begin\n";
#endif

    // TBD: split the relation here

    MEDDLY_DCASSERT(bv.isVoid());
    _compute(L, av, ap, bp, cv, cp);

#ifdef TRACE
    out << opName() << " #" << top_count << " end\n";
#endif
}

template <class EOP, bool FORWD, class ATYPE>
void MEDDLY::saturation1_set_mtrel<EOP, FORWD, ATYPE>::_compute(int L,
        const edge_value &av, node_handle A,
        node_handle B, edge_value &cv, node_handle &C)
{
    // TBD
    MEDDLY_DCASSERT(false);
}

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
    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {

        switch (c->getRangeType()) {
            case range_type::BOOLEAN:
                return new saturation1_set_mtrel<EdgeOp_none, FWD,
                            mt_prepost>(a, b, c);

            case range_type::INTEGER:
                if (c->isFullyReduced())  {
                    return new saturation1_set_mtrel<EdgeOp_none, FWD,
                            mt_distance>(a, b, c);
                }

            default:
                return nullptr;
        }
    }

    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {

        switch (a->getEdgeType()) {
            case edge_type::INT:
                return new saturation1_set_mtrel<EdgeOp_plus<int>, FWD,
                            ev_prepost<int> > (a, b, c);

            case edge_type::LONG:
                return new saturation1_set_mtrel<EdgeOp_plus<long>, FWD,
                            ev_prepost<long> > (a, b, c);

            default:
                return nullptr;
        };
    }
    return nullptr;
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

