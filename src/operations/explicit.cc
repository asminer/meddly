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
 */

// ******************************************************************
// *                                                                *
// *           template class for "set" minterm operators           *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class OP, class EdgeOp, bool RELS>
    class set_minterm_op : public minterm_operation {
        public:
            set_minterm_op(binary_builtin Union,
                    forest* arg1, minterm_coll &arg2, forest* res);

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

        private:
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
            binary_operation* union_op;
    };
};

// ******************************************************************

template <class OP, class EdgeOp, bool RELS>
MEDDLY::set_minterm_op<OP,EdgeOp,RELS>::set_minterm_op(binary_builtin Union,
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

    if ( (arg1->isForRelations() != RELS) ||
         (arg2.isForRelations() != RELS) )
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    union_op = Union(arg1, arg1, arg1);
    MEDDLY_DCASSERT(union_op);
}

// ******************************************************************

template <class OP, class EdgeOp, bool RELS>
void MEDDLY::set_minterm_op<OP,EdgeOp,RELS>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        unsigned low, unsigned high,
        edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    if (L == resF->getNumVariables()) {
        out.indentation(0);
        ++top_count;
        out << "generic_minterm_op #" << top_count << " begin\n";
        arg2.show(out, nullptr, "\n");
    }

    out << "compute L=" << L << ", in=" << in << ", ap=" << ap
        << ", [" << low << ", " << high << ")\n";
#endif

    if (OP::stop_recursion_on(av, ap)) {
        cv = av;
        cp = ap;
        return;
    }

    if (0==L) {
        OP::finalize(av, ap, arg2, low, high, cv, cp);
        return;
    }

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
    const int Lnext =
        RELS ? MXD_levels::downLevel(L)
             : MDD_levels::downLevel(L);
    const int Alevel = arg1F->getNodeLevel(ap);
    unpacked_node* Cu;
    if (Alevel != L) {
        if ((L<0) && (resF->isIdentityReduced())) {
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

    //
    // Break interval [low, high) into subintervals
    // based on value of variable x_L
    //
    unsigned dl = low;
    while (dl < high) {
        //
        // find end of subinterval
        //
        const int
            lowval = (L>0) ? arg2.at(dl)->getFrom(L) : arg2.at(dl)->getTo(-L);
        unsigned dh;
        for (dh=dl+1; dh<high; dh++) {
            int hv = (L>0) ? arg2.at(dh)->getFrom(L) : arg2.at(dh)->getTo(-L);
            if (hv != lowval) break;
        }
        //
        // Current interval is [dl, dh)
        // and all elements have value: lowval.
        //
        // Ready to recurse
        //
        edge_value tv;
        node_handle tp;

        if (lowval == DONT_CARE) {
            //
            // If we're at an unprimed level and we're a relation,
            // split the DONT_CARE region into
            // [dl, identstart), where the primed level is not DONT_CHANGE
            // [identstart, dh), where the primed level is DONT_CHANGE
            //
            unsigned identstart = dh;
            if (RELS && L>0) {
                for (identstart=dl; identstart<dh; identstart++) {
                    if ( DONT_CHANGE == arg2.at(identstart)->getTo(L) ) break;
                }

                //
                // Build the "don't change" portion
                // Note this is down TWO levels, at the next unprimed level
                //
                if (identstart < dh) {
                    MEDDLY_DCASSERT(!has_identity);
#ifdef TRACE
                    out << "    [" << identstart << ", " << dh << ") identity\n";
                    out.indent_more();
                    out.put('\n');
#endif
                    OP::get_null_edge(ident_val, ident_node);
                    ident_node = resF->makeRedundantsTo(ident_node, 0, L-1);
                    compute(L-1, ~0, ident_val, ident_node, identstart, dh,
                            ident_val, ident_node);
                    ident_node = resF->makeIdentitiesTo(ident_node, L-1, L, ~0);
                    has_identity = true;
#ifdef TRACE
                    out.indent_less();
                    out.put('\n');
                    out << "    [" << identstart << ", " << dh << ") done\n";
#endif
                }

            }

            //
            // Build the "don't care" portion
            //
            if (dl < identstart) {
                MEDDLY_DCASSERT(!has_dont_care);
#ifdef TRACE
                out << "    [" << dl << ", " << identstart << ") don't care\n";
                out.indent_more();
                out.put('\n');
#endif
                OP::get_null_edge(dnc_val, dnc_node);
                dnc_node = resF->makeRedundantsTo(dnc_node, 0, Lnext);
                compute(Lnext, ~0, dnc_val, dnc_node, dl, identstart,
                        dnc_val, dnc_node);
                dnc_node = resF->makeRedundantsTo(dnc_node, Lnext, L);
                has_dont_care = true;
#ifdef TRACE
                out.indent_less();
                out.put('\n');
                out << "    [" << dl << ", " << identstart << ") done\n";
#endif
            }

        } else {
#ifdef TRACE
            out << "    [" << dl << ", " << dh << ") value " << lowval << "\n";
            out.indent_more();
            out.put('\n');
#endif

            const unsigned i = unsigned(lowval);
            compute(Lnext, i, Cu->edgeval(i), Cu->down(i), dl, dh, tv, tp);
            resF->unlinkNode(Cu->down(i));
            Cu->setFull(i, tv, tp);

#ifdef TRACE
            out.indent_less();
            out.put('\n');
            out << "    [" << dl << ", " << dh << ") value " << lowval
                << " done\n";
#endif
        }

        //
        // Advance interval
        //
        dl = dh;

    } // while more subintervals

    resF->createReducedNode(Cu, cv, cp);
    if (has_dont_care) {
        union_op->compute(L, in, dnc_val, dnc_node, cv, cp, cv, cp);
    }
    if (RELS && has_identity) {
        union_op->compute(L, in, ident_val, ident_node, cv, cp, cv, cp);
    }
}

// ******************************************************************
// *                                                                *
// *        template class  for "relation" minterm operators        *
// *                                                                *
// ******************************************************************

// TBD

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
                terminal t = mc.at(i)->getTerm();
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
        if (c->isForRelations()) {
            return new set_minterm_op<mt_union_templ, EdgeOp_none, true>
                (UNION, a, b, c);
        } else {
            return new set_minterm_op<mt_union_templ, EdgeOp_none, false>
                (UNION, a, b, c);
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

