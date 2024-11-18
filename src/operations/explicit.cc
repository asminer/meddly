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
#include "../oper_minterm.h"
#include "../minterms.h"

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif


// ******************************************************************
// *                                                                *
// *          generic template class for minterm operators          *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class OP, class EdgeOp>
    class generic_minterm_op : public minterm_operation {
        /*
            OP must provide the following static methods:


                bool stop_recursion_on(const edge_value &av, node_handle ap);

                    If true, we can stop the recursion and return <av, ap>


                void finalize(const edge_value &av, node_handle ap,
                        const minterm_coll &mc, unsigned low, unsigned high,
                        edge_value &cv, node_handle &cp);

                    Determine the terminal edge for the given input edge
                    <av, ap> and minterms in mc with indexes in [low, high).
                    Store the result in <cv, cp>.
         */
        public:
            generic_minterm_op(forest* arg1, minterm_coll &arg2,
                    forest* res);

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

        private:
#ifdef TRACE
            ostream_output out;
            unsigned top_count;
#endif
    };
};

// ******************************************************************

template <class OP, class EdgeOp>
MEDDLY::generic_minterm_op<OP,EdgeOp>::generic_minterm_op(forest* arg1,
        minterm_coll &arg2, forest* res)
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

    if (arg1->isForRelations() != arg2.isForRelations()) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
}

// ******************************************************************

template <class OP, class EdgeOp>
void MEDDLY::generic_minterm_op<OP,EdgeOp>::compute(int L, unsigned in,
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
    // Allocate unpacked result node,
    // and copy from ap
    //
    const int Lnext = resF->isForRelations()
        ? MXD_levels::downLevel(L)
        : MDD_levels::downLevel(L)
    ;
    const int Alevel = arg1F->getNodeLevel(ap);
    unpacked_node* Cu;
    if (Alevel != L) {
        if (arg1F->isIdentityReduced() && L<0) {
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
        for (dh=low+1; dh<high; dh++) {
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

#ifdef TRACE
        out << "    [" << dl << ", " << dh << ") value " << lowval << "\n";
        out.indent_more();
        out.put('\n');
#endif

        if (lowval == DONT_CARE) {
            //
            // Special case: DONT_CARE
            //
            for (unsigned i=0; i<Cu->getSize(); i++) {
                compute(Lnext, i, Cu->edgeval(i), Cu->down(i), dl, dh, tv, tp);
                resF->unlinkNode(Cu->down(i));
                Cu->setFull(i, tv, tp);
            } // for i
        } else {
            //
            // Normal recursion, except for DONT_CHANGE
            // we set the index to the incoming index
            //
            if (lowval<0) {
                MEDDLY_DCASSERT(L<0);
                MEDDLY_DCASSERT(in != ~0);
            }
            const unsigned i = (lowval<0) ? in : unsigned(lowval);
            compute(Lnext, i, Cu->edgeval(i), Cu->down(i), dl, dh, tv, tp);
            resF->unlinkNode(Cu->down(i));
            Cu->setFull(i, tv, tp);
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "    [" << dl << ", " << dh << ") value " << lowval << " done\n";
#endif

        //
        // Advance interval
        //
        dl = dh;

    } // while more subintervals

    resF->createReducedNode(Cu, cv, cp);
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
                terminal t = mc.at(i)->getTerm();
                cp = t.getHandle();
                if (cp) return;
            }
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
        return new generic_minterm_op<mt_union_templ, EdgeOp_none> (a, b, c);
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

