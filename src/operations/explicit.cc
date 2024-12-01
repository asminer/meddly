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

#ifdef ALLOW_MINTERM_OPS

#include "../defines.h"
#include "explicit.h"
#include "../ops_builtin.h"
#include "../oper_binary.h"
#include "../oper_minterm.h"
#include "../minterms.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

/*
    To use the generic template classes below,
    class OP must provide the following static methods:


        bool stop_recursion_on(const edge_value &av, node_handle ap);

            If true, we can stop the recursion and return <av, ap>


        void finalize(const edge_value &av, node_handle ap,
                const minterm_coll &mc, unsigned low, unsigned high,
                edge_value &cv, node_handle &cp);

            Determine the terminal edge for the given input edge
            <av, ap> and minterms in mc with indexes in [low, high).
            Store the result in <cv, cp>.


        void get_null_edge(edge_value &av, node_handle &ap);

            Get the edge e such that e op x is x for all x.
            For op=union this is false. For op=intersection this
            is true.  For op=minimize this is infinity. For
            op=maximize this is negative infinity.


        bool is_null_edge(const edge_value &av, node_handle ap);

            Return true if (av, ap) is the edge e such that
            e op x is x for all x.
 */

// ******************************************************************
// *                                                                *
// *           template class for "set" minterm operators           *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class OP, class EdgeOp>
    class minterm_op : public minterm_operation {
        public:
            minterm_op(binary_builtin Union,
                    forest* arg1, minterm_coll &arg2, forest* res);

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

        private:
            // Set version of compute
            void setCompute(int L, const edge_value &av, node_handle ap,
                    unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

            // Relation version of compute
            void relCompute(int L, const edge_value &av, node_handle ap,
                    unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

            /// Build a chain of nodes for a single "set" minterm.
            ///     @param L    Last level to build a node.
            ///     @param m    Minterm to build from
            ///     @param cv   Input: bottom value
            ///                 Output: top value
            ///     @param cp   Input: bottom node
            ///                 Output: top node
            ///
            void setPathToBottom(int L, const minterm &m,
                    edge_value &cv, node_handle &cp)
            {
#ifdef TRACE
                out << "    setPathToBottom level " << L << "\n";
#endif
                MEDDLY_DCASSERT(m.isForSets());
                for (int k=1; k<=L; k++) {
                    if (DONT_CARE == m.from(k)) {
                        cp = resF->makeRedundantsTo(cp, k-1, k);
                        continue;
                    }
                    // make a singleton node
                    unpacked_node* nb = unpacked_node::newSparse(resF, k, 1);
                    MEDDLY_DCASSERT(m.from(k) >= 0);
                    nb->setSparse(0, m.from(k), cv, cp);
                    resF->createReducedNode(nb, cv, cp);
                } // for k
            }

            /// Build a chain of nodes for a single "relation" minterm.
            /// This is for identity-reduced forests:
            ///     we don't need to build identity patterns (yay)
            ///     but we need to check if edges to singletons
            ///     are legal (boo).
            ///
            ///     @param L    Last level to build a node.
            ///     @param m    Minterm to build from
            ///     @param cv   Input: bottom value
            ///                 Output: top value
            ///     @param cp   Input: bottom node
            ///                 Output: top node
            ///
            void identityPathToBottom(int L, const minterm &m,
                    edge_value &cv, node_handle &cp)
            {
#ifdef TRACE
                out << "    identityPathToBottom level " << L << "\n";
#endif
                MEDDLY_DCASSERT(m.isForRelations());
                MEDDLY_DCASSERT(resF->isIdentityReduced());
                MEDDLY_DCASSERT(L>0);
                for (int k=1; k<=L; k++) {
                    //
                    // Check for identity pattern at levels (k, k')
                    //
                    if (DONT_CHANGE == m.to(k)) {
                        MEDDLY_DCASSERT(DONT_CARE == m.from(k));
                        continue;
                    }

                    //
                    // Build node at primed level, unless skipped?
                    //
                    if (DONT_CARE == m.to(k)) {
                        if (DONT_CARE == m.from(k)) {
                            cp = resF->makeRedundantsTo(cp, k-1, k);
                            continue;
                        }
                        cp = resF->makeRedundantsTo(cp, k-1, -k);
                    } else {
                        if (m.from(k) != m.to(k)) {
                            //
                            // The singleton node at the primed
                            // level can be pointed at, so
                            // build the node.
                            //
                            unpacked_node* nb
                                = unpacked_node::newSparse(resF, k, 1);
                            nb->setSparse(0, m.to(k), cv, cp);
                            resF->createReducedNode(nb, cv, cp);
                        }
                    }

                    //
                    // Build node at unprimed level, unless skipped?
                    //
                    if (DONT_CARE != m.from(k)) {
                        unpacked_node* nb
                            = unpacked_node::newSparse(resF, k, 1);
                        nb->setSparse(0, m.from(k), cv, cp);
                        resF->createReducedNode(nb, cv, cp);
                    }
                } // for k
            }

            /// Build a chain of nodes for a single "relation" minterm.
            /// This is for non-identity-reduced forests:
            ///     we need to build identity patterns (boo)
            ///     but no need to check incoming edges to singletons (yay).
            ///
            ///     @param L    Last level to build a node.
            ///     @param m    Minterm to build from
            ///     @param cv   Input: bottom value
            ///                 Output: top value
            ///     @param cp   Input: bottom node
            ///                 Output: top node
            ///
            void relPathToBottom(int L, const minterm &m,
                    edge_value &cv, node_handle &cp)
            {
#ifdef TRACE
                out << "    relPathToBottom level " << L << "\n";
#endif
                MEDDLY_DCASSERT(m.isForRelations());
                MEDDLY_DCASSERT(!resF->isIdentityReduced());
                MEDDLY_DCASSERT(L>0);
                for (int k=1; k<=L; k++) {
                    //
                    // Check for identity pattern at levels (k, k')
                    //
                    if (DONT_CHANGE == m.to(k)) {
                        MEDDLY_DCASSERT(DONT_CARE == m.from(k));
                        cp = resF->makeIdentitiesTo(cp, k-1, k, ~0);
                        continue;
                    }

                    //
                    // Build node at primed level
                    //
                    if (DONT_CARE == m.to(k)) {
                        if (DONT_CARE == m.from(k)) {
                            cp = resF->makeRedundantsTo(cp, k-1, k);
                            continue;
                        }
                        cp = resF->makeRedundantsTo(cp, k-1, -k);
                    } else {
                        unpacked_node* nb
                            = unpacked_node::newSparse(resF, k, 1);
                        nb->setSparse(0, m.to(k), cv, cp);
                        resF->createReducedNode(nb, cv, cp);
                    }

                    //
                    // Build node at unprimed level
                    //
                    if (DONT_CARE == m.from(k)) {
                        cp = resF->makeRedundantsTo(cp, -k, k);
                    } else {
                        unpacked_node* nb
                            = unpacked_node::newSparse(resF, k, 1);
                        nb->setSparse(0, m.from(k), cv, cp);
                        resF->createReducedNode(nb, cv, cp);
                    }
                } // for k
            }

        private:
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            binary_operation* union_op;
    };
};

// ******************************************************************

template <class OP, class EdgeOp>
MEDDLY::minterm_op<OP,EdgeOp>::minterm_op(binary_builtin Union,
        forest* arg1, minterm_coll &arg2, forest* res)
        : minterm_operation(arg1, arg2, res)
#ifdef TRACE
            , out(std::cout), top_count(0)
#endif
{
    if (arg1 != res) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    if (arg1->getDomain() != arg2.getDomain()) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    if (arg1->isForRelations() != arg2.isForRelations())
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    union_op = Union(arg1, arg1, arg1);
    MEDDLY_DCASSERT(union_op);
}

// ******************************************************************

template <class OP, class EdgeOp>
void MEDDLY::minterm_op<OP,EdgeOp>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        unsigned low, unsigned high,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << "generic_minterm_op #" << top_count << " begin\n";
    arg2.show(out, nullptr, "\n");

    out << "compute L=" << L << ", in=" << in << ", ap=" << ap
        << ", [" << low << ", " << high << ")\n";
#endif

    if (arg2.isForRelations()) {
        relCompute(L, av, ap, low, high, cv, cp);
    } else {
        setCompute(L, av, ap, low, high, cv, cp);
    }
    return;

    //
    // OLD HERE
    //

#if 0

    if (OP::stop_recursion_on(av, ap)) {
        cv = av;
        cp = ap;
        return;
    }

    if (0==L) {
        OP::finalize(av, ap, arg2, low, high, cv, cp);
        return;
    }

    unpacked_node* Cu;

    //
    // Special, fast case: edge (av, ap) is the identity for
    // this operation, and we're down to only one element
    // in the collection. Use a loop to build from the bottom up.
    //
    MEDDLY_DCASSERT(high > low);
    if ( (high - low == 1) && OP::is_null_edge(av, ap) )
    {
#ifdef TRACE
        out << "    one midterm loop from bottom to level " << L << "\n";
#endif
        OP::finalize(av, ap, arg2, low, high, cv, cp);
        setPathToBottom(L, arg2.at(low), cv, cp);
        return;
    }

    // TBD HERE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // NO CT :)

    //
    // Initialize "don't care" result,
    // and our "don't change" result (relations only).
    //
    edge_value dnc_val, ident_val;
    node_handle dnc_node, ident_node;
    bool has_dont_care = false;
    bool has_identity = false;


    //
    // Allocate unpacked result node,
    // and copy from ap
    //
    const int Lnext = resF->isForRelations()
                        ? MXD_levels::downLevel(L)
                        : MDD_levels::downLevel(L);

    const int Alevel = arg1F->getNodeLevel(ap);
    if (Alevel != L) {
        if ((L<0) && (resF->isIdentityReduced()) && ap) {
            Cu = unpacked_node::newIdentity(arg1F, L, in, ap, FULL_ONLY);
        } else {
            Cu = unpacked_node::newRedundant(arg1F, L, ap, FULL_ONLY);
        }
    } else {
        Cu = arg1F->newUnpacked(ap, FULL_ONLY);
    }
    Cu->allowWrites(resF);

    //
    // Apply av to edges in Cu,
    // and increase reference counts
    //
    for (unsigned i=0; i<Cu->getSize(); i++)
    {
        Cu->setEdgeval(i, EdgeOp::accumulate(av, Cu->edgeval(i)) );
        resF->linkNode(Cu->down(i));
    }

    unsigned mid;
    do {
        //
        // Break off next interval
        //
        int val = arg2.collect_first(L, low, high, mid);
#ifdef TRACE
        out << "    processing " << val << " on [" << low << ", " << mid << ")\n";
        out.indent_more();
        out.put('\n');
#endif

        if (DONT_CARE == val) {
            //
            // We have to deal with either DONT_CARE,
            // or a DONT_CARE, DONT_CHANGE pair group.
            // First, determine which one.
            //
            if (resF->isForRelations() && L>0
                    && (DONT_CHANGE == arg2.at(low).var(-L)))
            {
                //
                // Build the "don't change" portion
                // Note this is down TWO levels, at the next unprimed level
                //
                MEDDLY_DCASSERT(!has_identity);
                OP::get_null_edge(ident_val, ident_node);
                ident_node = resF->makeRedundantsTo(ident_node, 0, L-1);
                compute(L-1, ~0, ident_val, ident_node, low, mid,
                            ident_val, ident_node);
                ident_node = resF->makeIdentitiesTo(ident_node, L-1, L, ~0);
                has_identity = true;
            } else {
                //
                // Build the "don't care" portion
                //
                MEDDLY_DCASSERT(!has_dont_care);
                OP::get_null_edge(dnc_val, dnc_node);
                dnc_node = resF->makeRedundantsTo(dnc_node, 0, Lnext);
                compute(Lnext, ~0, dnc_val, dnc_node, low, mid,
                            dnc_val, dnc_node);
                dnc_node = resF->makeRedundantsTo(dnc_node, Lnext, L);
                has_dont_care = true;
            }

        } else {
            //
            // Ordinary value
            // Just recurse
            //
            edge_value tv;
            node_handle tp;
            const unsigned i = unsigned(val);
            compute(Lnext, i, Cu->edgeval(i), Cu->down(i), low, mid, tv, tp);
            resF->unlinkNode(Cu->down(i));
            Cu->setFull(i, tv, tp);
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "    done val " << val << " on [" << low << ", " << mid << ")\n";
#endif
        low = mid;
    } while (mid < high);

    resF->createReducedNode(Cu, cv, cp);
    if (has_dont_care) {
        union_op->compute(L, in, dnc_val, dnc_node, cv, cp, cv, cp);
    }
    if (has_identity) {
        union_op->compute(L, in, ident_val, ident_node, cv, cp, cv, cp);
    }
#endif
}

// ******************************************************************

template <class OP, class EdgeOp>
void MEDDLY::minterm_op<OP,EdgeOp>::setCompute(int L,
        const edge_value &av, node_handle ap,
        unsigned low, unsigned high,
        edge_value &cv, node_handle &cp)
{
    if (OP::stop_recursion_on(av, ap)) {
        cv = av;
        cp = ap;
        return;
    }

    if (0==L) {
        OP::finalize(av, ap, arg2, low, high, cv, cp);
        return;
    }

    MEDDLY_DCASSERT(high > low);

    //
    // Special case for only one minterm left
    //
    if ( (high - low == 1) && OP::is_null_edge(av, ap) )
    {
        OP::finalize(av, ap, arg2, low, high, cv, cp);
        setPathToBottom(L, arg2.at(low), cv, cp);
        return;
    }

    //
    // Determine minterm value range
    //
    int minV, maxV;
    arg2.getMinMax(L, low, high, minV, maxV);

    //
    // Special case: all values are DONT_CARE
    //
    if (DONT_CARE == maxV) {
        edge_value zv;
        node_handle zp;
        OP::get_null_edge(zv, zp);
        setCompute(L-1, zv, zp, low, high, cv, cp);
        cp = resF->makeRedundantsTo(cp, L-1, L);
        if (!OP::is_null_edge(av, ap)) {
            union_op->compute(L, ~0, av, ap, cv, cp, cv, cp);
        }
        return;
    }

    //
    // Start unpacked result node, from <av, ap>
    //
    unpacked_node* Cu;
    const int Alevel = arg1F->getNodeLevel(ap);
    if (Alevel != L || OP::is_null_edge(av, ap)) {
        Cu = unpacked_node::newRedundant(arg1F, L, ap, FULL_ONLY);
    } else {
        Cu = arg1F->newUnpacked(ap, FULL_ONLY);
    }
    Cu->allowWrites(resF);

    //
    // Apply av to edges in Cu,
    // and increase reference counts
    //
    for (unsigned i=0; i<Cu->getSize(); i++)
    {
        Cu->setEdgeval(i, EdgeOp::accumulate(av, Cu->edgeval(i)) );
        resF->linkNode(Cu->down(i));
    }

    //
    // Vars for recursion
    //
    edge_value dnc_val, cz_val;
    node_handle dnc_node, cz_node;
    bool has_dont_care = false;
    unsigned mid = low;

    //
    // Recursion loop
    //
    while (minV < maxV) {
        const int currV = minV;
        arg2.moveValueToFront(L, minV, low, high, mid);

        if (DONT_CARE == currV) {
            //
            // Build a redundant node at level L
            // but recurse to figure out where it points to.
            //
            edge_value zv;
            node_handle zp;
            OP::get_null_edge(zv, zp);
            setCompute(L-1, zv, zp, low, high, dnc_val, dnc_node);
            dnc_node = resF->makeRedundantsTo(dnc_node, L-1, L);
            has_dont_care = true;

            // Make sure dnc_node isn't recycled yet
            Cu->setTempRoot(dnc_node);

        } else {
            //
            // Recurse and add result to unpacked node Cu
            //
            setCompute(L-1, Cu->edgeval(currV), Cu->down(currV),
                    low, mid, cz_val, cz_node);
            resF->unlinkNode(Cu->down(currV));
            Cu->setFull(currV, cz_val, cz_node);
        }

        low = mid;
    }
    MEDDLY_DCASSERT(maxV>=0);
    MEDDLY_DCASSERT(minV == maxV);
    //
    // Last value
    //
    setCompute(L-1, Cu->edgeval(minV), Cu->down(minV),
            mid, high, cz_val, cz_node);
    resF->unlinkNode(Cu->down(minV));
    Cu->setFull(minV, cz_val, cz_node);

    //
    // Reduce result node
    //
    resF->createReducedNode(Cu, cv, cp);
    if (has_dont_care) {
        union_op->compute(L, ~0, dnc_val, dnc_node, cv, cp, cv, cp);
    }
}

// ******************************************************************

template <class OP, class EdgeOp>
void MEDDLY::minterm_op<OP,EdgeOp>::relCompute(int L,
        const edge_value &av, node_handle ap,
        unsigned low, unsigned high,
        edge_value &cv, node_handle &cp)
{
}

// ******************************************************************
// *                                                                *
// *                        operator classes                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    //
    // union operation for multi-terminal forests
    //
    struct mt_union_templ {
        static inline bool stop_recursion_on(const edge_value &av,
                node_handle ap)
        {
            MEDDLY_DCASSERT(av.isVoid());

            if (!forest::isTerminalNode(ap)) return false;

            return ap != 0;
        }

        static inline void finalize(const edge_value &av, node_handle ap,
                        const minterm_coll &mc, unsigned low, unsigned high,
                        edge_value &cv, node_handle &cp)
        {
            MEDDLY_DCASSERT(av.isVoid());
            cv.set();
            for (unsigned i=low; i<high; i++) {
                terminal t = mc.at(i).getTerm();
                cp = t.getHandle();
                if (cp) return;
            }
        }

        static inline void get_null_edge(edge_value &av, node_handle &ap)
        {
            av.set();
            terminal t(false);
            ap = t.getHandle();
        }

        static inline bool is_null_edge(const edge_value &av, node_handle ap)
        {
            MEDDLY_DCASSERT(av.isVoid());
            terminal t(false);
            return (ap == t.getHandle());
        }

    };
};


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::minterm_operation*
MEDDLY::UNION(forest* a, minterm_coll &b, forest* c)
{
    if (!a || !c) {
        return nullptr;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return new minterm_op<mt_union_templ,EdgeOp_none> (UNION, a, b, c);
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

#endif // allow_minterm_ops
