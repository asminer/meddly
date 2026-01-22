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
#include "../ct_vector.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../forest_levels.h"
#include "../forest_edgerules.h"
#include "user_unary.h"

#ifdef TRACE_USER
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                      user_unary_op  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EdgeOp>
    class user_unary_op : public unary_operation {
        public:
            user_unary_op(forest* arg, forest* res,
                    user_defined_unary F);

            virtual ~user_unary_op();

            virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

        protected:
            void _compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

        private:
            inline const edge_value &edgeval(unpacked_node *U, unsigned i) const
            {
                if (EdgeOp::hasEdgeValues()) {
                    return U->edgeval(i);
                } else {
                    return nothing;
                }
            }

        private:
            user_defined_unary F;
            ct_entry_type* ct;
#ifdef TRACE_USER
            ostream_output out;
            unsigned top_count;
#endif
            edge_value nothing;
            bool preserves_sparse;
    };
};

// ******************************************************************
// *                     user_unary_op  methods                     *
// ******************************************************************

template <class EdgeOp>
MEDDLY::user_unary_op<EdgeOp>::user_unary_op(
        forest* arg, forest* res, user_defined_unary _F)
        : unary_operation(arg, res), F(_F)
#ifdef TRACE_USER
            , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);

    ct = new ct_entry_type("user_unary");
    if (arg->isMultiTerminal()) {
        ct->setFixed(arg);
    } else {
        ct->setFixed(arg->getEdgeType(), arg);
    }
    if (res->isMultiTerminal()) {
        ct->setResult(res);
    } else {
        ct->setResult(res->getEdgeType(), res);
    }
    ct->doneBuilding();

    //
    // Does F(0) == 0?
    // In the generic sense.
    //
    edge_value  zev;
    node_handle zdp;
    rangeval zero;
    argF->getTransparentEdge(zev, zdp);
    argF->getValueForEdge(zev, zdp, zero);

    rangeval Fzero;
    F(zero, Fzero);

    resF->getEdgeForValue(Fzero, zev, zdp);
    preserves_sparse = resF->isTransparentEdge(zev, zdp);
}

// ******************************************************************

template <class EdgeOp>
MEDDLY::user_unary_op<EdgeOp>::~user_unary_op()
{
    ct->markForDestroy();
}

// ******************************************************************

template <class EdgeOp>
void MEDDLY::user_unary_op<EdgeOp>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
#ifdef TRACE_USER
    ++top_count;
    out.indentation(0);
    out << "Starting top-level user_unary_op::compute #" << top_count << "\n";
#endif
    _compute(L, in, av, ap, cv, cp);
#ifdef TRACE_USER
    out.indentation(0);
    out << "Finishing top-level user_unary_op::compute #" << top_count << "\n";
#endif
}

// ******************************************************************

template <class EdgeOp>
void MEDDLY::user_unary_op<EdgeOp>::_compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
    //
    // Terminal case
    //
    if (argF->isTerminalNode(ap)) {
        rangeval A, C;
        argF->getValueForEdge(av, ap, A);
        F(A, C);
        resF->getEdgeForValue(C, cv, cp);
        return;
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(ap);
#ifdef TRACE_USER
    out << "user_unary_op::sparse_compute(";
    argF->showEdge(out, av, ap);
    out << ")\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    if (argF->isMultiTerminal()) {
        key[0].setN(ap);
    } else {
        key[0].set(av);
        key[1].setN(ap);
    }
    if (ct->findCT(key, res)) {
        //
        // compute table 'hit'
        //
        if (resF->isMultiTerminal()) {
            cv.set();
            cp = resF->linkNode(res[0].getN());
        } else {
            res[0].get(cv);
            cp = resF->linkNode(res[1].getN());
        }
#ifdef TRACE_USER
        out << "  CT hit ";
        resF->showEdge(out, cv, cp);
        out << "\n";
        out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
        //
        // done: compute table 'hit'
        //
    } else {
        //
        // compute table 'miss'; do computation.
        // Initialize unpacked nodes
        //
        unpacked_node* Au = nullptr;
        unpacked_node* Cu = nullptr;
        if (preserves_sparse) {
            Au = unpacked_node::newFromNode(argF, ap, SPARSE_ONLY);
            Cu = unpacked_node::newWritable(resF, Alevel, SPARSE_ONLY);
        } else {
            Au = unpacked_node::newFromNode(argF, ap, FULL_ONLY);
            Cu = unpacked_node::newWritable(resF, Alevel, FULL_ONLY);
            MEDDLY_DCASSERT(Au->getSize() == Cu->getSize());
        }
#ifdef TRACE_USER
        out << "A: ";
        Au->show(out, true);
        out.indent_more();
        out.put('\n');
#endif
        const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

        if (preserves_sparse) {
            //
            // Nodes are sparse; visit children
            //
            unsigned cz = 0;
            for (unsigned z=0; z<Au->getSize(); z++) {
                const unsigned i = Au->index(z);
                edge_value v;
                node_handle d;
                _compute(Cnextlevel, i,
                        EdgeOp::applyOp(av, edgeval(Au, z)), Au->down(z),
                        v, d);

                if (resF->isTransparentEdge(v, d)) continue;
                if (resF->isMultiTerminal()) {
                    MEDDLY_DCASSERT(v.isVoid());
                    Cu->setSparse(cz, i, d);
                } else {
                    Cu->setSparse(cz, i, v, d);
                }
                ++cz;
            } // for z
            Cu->shrink(cz);
        } else {
            //
            // Nodes are full; visit children
            //
            for (unsigned i=0; i<Au->getSize(); i++) {
                edge_value v;
                node_handle d;
                _compute(Cnextlevel, i,
                        EdgeOp::applyOp(av, edgeval(Au, i)), Au->down(i),
                        v, d);

                if (resF->isTransparentEdge(v, d)) continue;
                if (resF->isMultiTerminal()) {
                    MEDDLY_DCASSERT(v.isVoid());
                    Cu->setFull(i, d);
                } else {
                    Cu->setFull(i, v, d);
                }
            } // for i
        }


        //
        // Recurse over child (sparse) edges, pushing values down
        //

#ifdef TRACE_USER
        out.indent_less();
        out.put('\n');
        out << "user_unary_op::sparse_compute(";
        argF->showEdge(out, av, ap);
        out << ") done\nA: ";
        Au->show(out, true);
        out << "\nC: ";
        Cu->show(out, true);
        out << "\n";
#endif

        //
        // Reduce
        //
        resF->createReducedNode(Cu, cv, cp);
#ifdef TRACE_USER
        out << "reduced to ";
        resF->showEdge(out, cv, cp);
        out << ": ";
        resF->showNode(out, cp, SHOW_DETAILS);
        out << "\n";
#endif

        //
        // Add to CT
        //
        if (resF->isMultiTerminal()) {
            MEDDLY_DCASSERT(cv.isVoid());
            res[0].setN(cp);
        } else {
            res[0].set(cv);
            res[1].setN(cp);
        }
        ct->addCT(key, res);

        //
        // Cleanup
        //
        unpacked_node::Recycle(Au);
        // Cu is recycled when we reduce it :)

        //
        // done: compute table 'miss'
        //
    }
    //
    // Adjust result for singletons, added identities/redundants.
    // Do this for both CT hits and misses.
    //
    const int Clevel = resF->getNodeLevel(cp);
    if (Clevel == L) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        cp = resF->redirectSingleton(in, cp);
        return;
    }

    //
    // Add nodes but only from Alevel;
    // if Clevel is below Alevel it means nodes were eliminated
    // in the result forest.
    //
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, Alevel, L, in);
    } else {
        cp = resF->makeRedundantsTo(cp, Alevel, L);
    }
}

// ******************************************************************
// *                                                                *
// *                   user_unary_factory methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::user_unary_factory::user_unary_factory(const char* n,
        user_defined_unary _F)
{
    name = n;
    F = _F;
    setup();
}

void MEDDLY::user_unary_factory::setup()
{
    _setup(__FILE__, "user-defined", "User-defined unary operation");
}

MEDDLY::unary_operation*
MEDDLY::user_unary_factory::build_new(forest* arg, forest* res)
{
    if (arg->isMultiTerminal()) {
        return new user_unary_op<EdgeOp_none>(arg, res, F);
    }

    if (arg->isEVTimes()) {
        return new user_unary_op< EdgeOp_times<float> >(arg, res, F);
    }

    return new user_unary_op< EdgeOp_plus<long> >(arg, res, F);
}

