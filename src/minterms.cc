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

#include "minterms.h"

#include "io.h"
#include "domain.h"
#include "forest.h"
#include "oper_binary.h"
#include "ops_builtin.h"

// #define DEBUG_COLLECT_FIRST
// #include "operators.h"

// ******************************************************************
// *                                                                *
// *                         fbuilder class                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    /*
        struct OP must provide the following static methods:

        void finalize(const minterm_coll &mc, unsigned low, unsigned high,
                edge_value &cv, node_handle &cp);

            Determine the terminal edge for the minterms in mc with
            indexes in [low, high). Store the result in <cv, cp>.

    */
    template <class OP, class EdgeOp>
    class fbuilder {
        public:
            fbuilder(forest* f, minterm_coll &mtlist, binary_builtin Union)
                : mtc(mtlist)
            {
                F = f;
                union_op = Union(f, f, f);
                MEDDLY_DCASSERT(union_op);
            }

            void compute(int L, unsigned in,
                    unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

        private:
            forest* F;
            minterm_coll &mtc;
            binary_operation* union_op;
    };

    // ******************************************************************
    // *                                                                *
    // *                 operation classes for fbuilder                 *
    // *                                                                *
    // ******************************************************************

    struct fbop_union_bool {
        static inline void finalize(const minterm_coll &mc,
            unsigned low, unsigned high, edge_value &cv, node_handle &cp)
        {
            cv.set();
            for (unsigned i=low; i<high; i++) {
                terminal t = mc.at(i).getTerm();
                cp = t.getHandle();
                if (cp) return;
            }
        }

    };
};

// ******************************************************************
// *                                                                *
// *                        minterm  methods                        *
// *                                                                *
// ******************************************************************


MEDDLY::minterm::minterm(const domain* D, set_or_rel sr) : termval(true)
{
    _D = D;
    if (D) {
        num_vars = D->getNumVariables();
    } else {
        num_vars = 0;
    }
    for_relations = sr;
#ifdef USE_VECTOR
    _from.resize(1+num_vars, DONT_CARE);
    if (for_relations) {
        _to.resize(1+num_vars, DONT_CARE);
    }
#else
    _from = new int[1+num_vars];
    if (for_relations) {
        _to = new int[1+num_vars];
    } else {
        _to = nullptr;
    }
#endif
}

MEDDLY::minterm::~minterm()
{
#ifndef USE_VECTOR
    delete[] _from;
    delete[] _to;
#endif
}

void MEDDLY::minterm::show(output &s) const
{
    //
    // Full storage
    //
    s.put("[ top");
    for (unsigned i=num_vars; i; i--) {
        s.put(", ");
        if (_from[i] < 0)   s.put('x');
        else                s.put(_from[i]);
        if (isForSets()) continue;
        s.put("->");
        if (DONT_CARE == _to[i])        s.put('x');
        else if (DONT_CHANGE == _to[i]) s.put('i');
        else                            s.put(_to[i]);
    }
    s.put(", bot]");
}

// ******************************************************************
// *                                                                *
// *                      minterm_coll methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::minterm_coll::minterm_coll(unsigned maxsz, const domain* D,
            set_or_rel sr)
#ifndef USE_VECTOR
    : max_coll_size(maxsz)
#endif
{
    _D = D;
    num_vars = D ? D->getNumVariables() : 0;
    for_relations = sr;

#ifdef USE_VECTOR
    _mtlist.resize(maxsz);
#else
    _mtlist = new minterm*[max_coll_size];
#endif
    for (unsigned i=0; i<maxsize(); i++) {
        _mtlist[i] = new minterm(_D, for_relations);
    }
    clear();
}


MEDDLY::minterm_coll::~minterm_coll()
{
    for (unsigned i=0; i<maxsize(); i++) {
        delete _mtlist[i];
        _mtlist[i] = nullptr;
    }
#ifndef USE_VECTOR
    delete[] _mtlist;
#endif
}

void MEDDLY::minterm_coll::buildFunction(dd_edge &e)
{
    if (!e.getForest()) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if (e.getForest()->getDomain() != _D) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (e.getForest()->isForRelations() != isForRelations())
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    if (e.getForest()->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        fbuilder<fbop_union_bool, EdgeOp_none> fb(e.getForest(), *this, UNION);
        node_handle en;
        fb.compute(num_vars, ~0, 0, first_unused, e.setEdgeValue(), en);
        e.set(en);
        return;
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::minterm_coll::show(output &s,
        const char* pre, const char* post) const
{
    for (unsigned i=0; i<size(); i++) {
        if (pre) {
            s.put(pre);
        } else {
            s.put("  ");
            s.put(i, 4);
            s.put("  ");
        }
        _mtlist[i]->show(s);
        s.put(post);
    }
}

int MEDDLY::minterm_coll::_collect_first(int L, unsigned low,
        unsigned hi, unsigned& lend)
{
    int val = at(low).var(L);
    bool val_dc;
    if ((DONT_CARE == val) && isForRelations() && (L>0)) {
        val_dc = (DONT_CHANGE == at(low).var(-L));
    } else {
        val_dc = false;
    }
    lend = low+1;
    unsigned hptr = hi;

    for (;;)
    {
        //
        // Find smallest element different from val
        //
        for (;;) {
            const int vl = at(lend).var(L);
            if (vl != val) break;
            if ((DONT_CARE == vl) && isForRelations() && (L>0)) {
                const bool vl_dc = (DONT_CHANGE == at(lend).var(-L));
                if (val_dc != vl_dc) break;
            }
            ++lend;
            if (lend >= hptr) {
#ifdef DEBUG_COLLECT_FIRST
                FILE_output out(stdout);
                out << "Collected first on range [" << low << ", "
                    << hi << "):\n";
                for (unsigned i=low; i<hi; i++) {
                    if (i==lend) out << "    ======================\n";
                    out << "    ";
                    out.put(long(i), 2);
                    out << "  ";
                    at(i).show(out);
                    out << "\n";
                }
                out << "val " << val << "\n";
#endif
                return val;
            }
        }

        //
        // Find high element to swap with
        //
        for (;;) {
            const int vh = at(hptr-1).var(L);
            if (vh == val) {
                if ((DONT_CARE == vh) && isForRelations() && (L>0)) {
                    const bool vh_dc = (DONT_CHANGE == at(hptr-1).var(-L));
                    if (val_dc == vh_dc) break;
                } else {
                    break;
                }
            }
            --hptr;
            if (lend >= hptr) {
#ifdef DEBUG_COLLECT_FIRST
                FILE_output out(stdout);
                out << "Collected first on range [" << low << ", "
                    << hi << "):\n";
                for (unsigned i=low; i<hi; i++) {
                    if (i==lend) out << "    ======================\n";
                    out << "    ";
                    out.put(long(i), 2);
                    out << "  ";
                    at(i).show(out);
                    out << "\n";
                }
                out << "val " << val << "\n";
#endif
                return val;
            }
        }

        //
        // Swap elements
        //
        minterm* tmp = _mtlist[lend];
        _mtlist[lend] = _mtlist[hptr-1];
        _mtlist[hptr-1] = tmp;
    }

}


// ******************************************************************
// *                                                                *
// *                        fbuilder methods                        *
// *                                                                *
// ******************************************************************

template <class OP, class EdgeOp>
void MEDDLY::fbuilder<OP,EdgeOp>::compute(int L, unsigned in,
        unsigned low, unsigned high,
        edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(high > low);
    if (0==L) {
        OP::finalize(mtc, low, high, cv, cp);
        return;
    }

    unpacked_node* Cu;

    //
    // Special case: only one minterm left
    //
    if (high - low == 1) {
        OP::finalize(mtc, low, high, cv, cp);
        int k = 0;
        while (k != L) {
            const int prevk = k;
            k = F->isForRelations()
                    ? MXD_levels::upLevel(k)
                    : MDD_levels::upLevel(k)
            ;
            const int vark = mtc.at(low).var(k);

            if (DONT_CHANGE == vark) {
                //
                // Next 2 levels are identity pattern;
                // process them together.
                //
                MEDDLY_DCASSERT(k<0);
                MEDDLY_DCASSERT(DONT_CARE == mtc.at(low).var(-k));

                k = MXD_levels::upLevel(k);
                cp = F->makeIdentitiesTo(cp, k-1, k, ~0);
                continue;
            }

            if (DONT_CARE == vark) {
                //
                // Next level is redundant
                //
                cp = F->makeRedundantsTo(cp, prevk, k);
                continue;
            }

            //
            // Create a singleton node.
            // Except don't if we're at a primed level and
            // the forest is identity reduced.
            //
            if (F->isIdentityReduced() && k<0) {
                //
                // See if we can point to the singleton node
                // we're about to create.
                //
                if (k==L) {
                    if (in == vark) {
                        //
                        // skip this node.
                        // and k==L so we're done
                        //
                        return;
                    }
                } else {
                    if (vark == mtc.at(low).from(-k)) {
                        //
                        // skip this node
                        //
                        continue;
                    }
                }
            }

            Cu = unpacked_node::newSparse(F, k, 1);
            Cu->setSparse(0, vark, cv, cp);
            F->createReducedNode(Cu, cv, cp);
        }
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
    // Allocate unpacked result node
    //
    const int Lnext = F->isForRelations()
                        ? MXD_levels::downLevel(L)
                        : MDD_levels::downLevel(L);
    Cu = unpacked_node::newFull(F, L, F->getLevelSize(L));

    unsigned mid;
    do {
        //
        // Break off next interval
        //
        int val = mtc.collect_first(L, low, high, mid);
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
            if (F->isForRelations() && L>0
                    && (DONT_CHANGE == mtc.at(low).var(-L)))
            {
                //
                // Build the "don't change" portion
                // Note this is down TWO levels, at the next unprimed level
                //
                MEDDLY_DCASSERT(!has_identity);
                compute(L-1, ~0, low, mid, ident_val, ident_node);
                ident_node = F->makeIdentitiesTo(ident_node, L-1, L, ~0);
                has_identity = true;
            } else {
                //
                // Build the "don't care" portion
                //
                MEDDLY_DCASSERT(!has_dont_care);
                compute(Lnext, ~0, low, mid, dnc_val, dnc_node);
                dnc_node = F->makeRedundantsTo(dnc_node, Lnext, L);
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
            compute(Lnext, i, low, mid, tv, tp);
            F->unlinkNode(Cu->down(i));
            Cu->setFull(i, tv, tp);
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "    done val " << val << " on [" << low << ", " << mid << ")\n";
#endif
        low = mid;
    } while (mid < high);

    F->createReducedNode(Cu, cv, cp);
    if (has_dont_care) {
        union_op->compute(L, in, dnc_val, dnc_node, cv, cp, cv, cp);
    }
    if (has_identity) {
        union_op->compute(L, in, ident_val, ident_node, cv, cp, cv, cp);
    }
}

