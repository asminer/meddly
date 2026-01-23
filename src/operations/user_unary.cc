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

// #define TRACE
// #define TRACE_FI

#ifdef TRACE
#include "../operators.h"
#endif
#ifdef TRACE_FI
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

            inline void F_identity(edge_value &ev, node_handle &p,
                    int K, int L, unsigned in)
            {
                if (0==L) return;
                if (K==L) return;
                MEDDLY_DCASSERT(MXD_levels::topLevel(L, K) == L);
                if (preserves_sparse) {
                    p = resF->makeIdentitiesTo(p, K, L, in);
                } else {
                    _F_identity(ev, p, K, L, in);
                }
            }

#ifdef TRACE_FI
            inline void showCount(output &out, const char* w, node_handle p)
                const
            {
                out << "    " << w << " " << p;
                if (p>0) {
                    out << " incount " << resF->getNodeInCount(p);
                }
                out << '\n';
            }
#endif

            /*
                Build a translated identity pattern:
                    [ p f f ... f ]
                    [ f p f ... f ]
                    [ f f p ... f ]
                    [ :         : ]
                    [ f f f ... p ]
                where f is Fzero.

                    @param  ev  Edge value on diagonal
                    @param  p   Down pointer on diagonal
                    @param  K   Start the pattern above this level
                    @param  L   Stop the pattern at this level
                    @param  in  Incoming edge index; needed only if L is primed

                    On output, <ev, p> encodes the desired function.
            */
            void _F_identity(edge_value &ev, node_handle& p,
                    int K, int L, unsigned in);

        private:
            user_defined_unary F;
            ct_entry_type* ct;
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            edge_value nothing;
            rangeval Fzero;
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
#ifdef TRACE
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
#ifdef TRACE
    ++top_count;
    out.indentation(0);
    out << "Starting top-level user_unary_op::compute #" << top_count << "\n";
#endif
    _compute(L, in, av, ap, cv, cp);
#ifdef TRACE
    out.indentation(0);
    out << "Finishing top-level user_unary_op::compute #" << top_count << "\n";
#endif
}

// ******************************************************************

template <class EdgeOp>
void MEDDLY::user_unary_op<EdgeOp>::_compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
    const int Alevel = argF->getNodeLevel(ap);

    //
    // Terminal case
    //
    if (argF->isTerminalNode(ap)) {
        rangeval A, C;
        argF->getValueForEdge(av, ap, A);
        F(A, C);
        resF->getEdgeForValue(C, cv, cp);
        if (argF->isIdentityReduced()) {
            F_identity(cv, cp, Alevel, L, in);
        } else {
            cp = resF->makeRedundantsTo(cp, Alevel, L);
        }
        return;
    }

#ifdef TRACE
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
#ifdef TRACE
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
#ifdef TRACE
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

#ifdef TRACE
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
#ifdef TRACE
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
        F_identity(cv, cp, Alevel, L, in);
    } else {
        cp = resF->makeRedundantsTo(cp, Alevel, L);
    }
}

// ******************************************************************

template <class EdgeOp>
void MEDDLY::user_unary_op<EdgeOp>::_F_identity(edge_value &ev, node_handle &p,
    int K, int L, unsigned in)
{
#ifdef TRACE_FI
    ostream_output trout(std::cout);
    trout << "F_identity(";
    resF->showEdge(trout, ev, p);
    trout << ", " << K << ", "
          << L << ", " << in << ")\n";
#endif
    MEDDLY_DCASSERT(L!=0);
    unpacked_node* Uun;
    unpacked_node* Upr;

    if (K<0) {
        //
        // Add a redundant layer as needed;
        // this also handles the only case when we need to check
        // if p is a singleton node :)
        //
#ifdef TRACE_FI
        trout << "  add redundant layer\n";
        const node_handle oldp = p;
        showCount(trout, "old p", oldp);
#endif
        p = resF->makeRedundantsTo(p, K, MXD_levels::upLevel(K));
#ifdef TRACE_FI
        trout << "    new p " << p << "\n";
        showCount(trout, "old p", oldp);
#endif
        K = MXD_levels::upLevel(K);
        if (L == K) return;
    }

    MEDDLY_DCASSERT(K>=0);
    const int Lstop = (L<0) ? MXD_levels::downLevel(L) : L;

    edge_value  Fzero_ev;
    node_handle Fzero_p;

    resF->getEdgeForValue(Fzero, Fzero_ev, Fzero_p);
    Fzero_p = resF->makeRedundantsTo(Fzero_p, 0, K);

#ifdef TRACE_FI
    trout << "    F(0) is ";
    resF->showEdge(trout, Fzero_ev, Fzero_p);
    trout << ": ";
    resF->showNode(trout, Fzero_p, SHOW_DETAILS);
    trout << "\n      p " << p << " ";
    resF->showNode(trout, p, SHOW_DETAILS);
    trout << "\n";
#endif


    //
    // Proceed in unprimed, primed pairs
    //
    for (K++; K<=Lstop; K++) {
        Uun = unpacked_node::newWritable(resF, K, FULL_ONLY);

        for (unsigned i=0; i<Uun->getSize(); i++) {
            Upr = unpacked_node::newWritable(resF, -K, FULL_ONLY);
            for (unsigned j=0; j<Uun->getSize(); j++) {
                if (j==i) {
                    Upr->setFull(j, ev, resF->linkNode(p));
                } else {
                    Upr->setFull(j, Fzero_ev, resF->linkNode(Fzero_p));
                }
            }
            node_handle hp;
            edge_value  hv;
            resF->createReducedNode(Upr, hv, hp);
            Uun->setFull(i, hv, hp);
        }
#ifdef TRACE_FI
        trout << "After level " << K << "\n";
#endif

        resF->unlinkNode(p);
#ifdef TRACE_FI
        showCount(trout, "p", p);
        resF->showNode(trout, p, SHOW_DETAILS);
        trout << "\n";
        trout << "\n      p := reduce ";
        Uun->show(trout, true);
        trout << "\n";
#endif
        resF->createReducedNode(Uun, ev, p);


#ifdef TRACE_FI
        trout << "    F(0) node " << Fzero_p << " ";
        resF->showNode(trout, Fzero_p, SHOW_DETAILS);
        trout << "\n";
#endif
        Fzero_p = resF->makeRedundantsTo(Fzero_p, K-1, K);
#ifdef TRACE_FI
        trout << "Before next\n";
        trout << "    F(0) node " << Fzero_p << " ";
        resF->showNode(trout, Fzero_p, SHOW_DETAILS);
        trout << "\n";
        showCount(trout, "p", p);
        resF->showNode(trout, p, SHOW_DETAILS);
        trout << "\n";
#endif
    } // for k

    //
    // Add top primed node, if L is negative
    //
    if (L<0) {
        MEDDLY_DCASSERT(-K == L);
        Upr = unpacked_node::newWritable(resF, -K, FULL_ONLY);
        for (unsigned j=0; j<Upr->getSize(); j++) {
            if (j == in) {
                Upr->setFull(j, ev, p);
            } else {
                Upr->setFull(j, Fzero_ev, resF->linkNode(Fzero_p));
            }
        }
        resF->createReducedNode(Upr, ev, p);
    }
    resF->unlinkNode(Fzero_p);
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

