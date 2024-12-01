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
// #define DEBUG_SORT
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
    template <class OP>
    class fbuilder {
        public:
            fbuilder(forest* f, minterm_coll &mtlist, binary_builtin Union)
                : mtc(mtlist)
            {
                F = f;
                union_op = Union(f, f, f);
                MEDDLY_DCASSERT(union_op);
            }

            inline int unprimed(unsigned i, int L) const {
                return mtc.unprimed(i, L);
            }
            inline void swap(unsigned i, unsigned j) {
                mtc.swap(i, j);
            }
            inline const minterm& element(unsigned i) const {
                return mtc.at(i);
            }

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
                MEDDLY_DCASSERT(L>0);
                MEDDLY_DCASSERT(!m.isForRelations());
                for (int k=1; k<=L; k++) {
                    if (DONT_CARE == m.from(k)) {
                        cp = F->makeRedundantsTo(cp, k-1, k);
                    } else {
                        // make a singleton node
                        unpacked_node* nb = unpacked_node::newSparse(F, k, 1);
                        MEDDLY_DCASSERT(m.from(k) >= 0);
                        nb->setSparse(0, m.from(k), cv, cp);
                        F->createReducedNode(nb, cv, cp);
                    }
                } // for i
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
                MEDDLY_DCASSERT(L>0);
                MEDDLY_DCASSERT(m.isForRelations());
                MEDDLY_DCASSERT(F->isIdentityReduced());
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
                            cp = F->makeRedundantsTo(cp, k-1, k);
                            continue;
                        }
                        cp = F->makeRedundantsTo(cp, k-1, -k);
                    } else {
                        if (m.from(k) != m.to(k)) {
                            //
                            // The singleton node at the primed
                            // level can be pointed at, so
                            // build the node.
                            //
                            unpacked_node* nb
                                = unpacked_node::newSparse(F, k, 1);
                            nb->setSparse(0, m.to(k), cv, cp);
                            F->createReducedNode(nb, cv, cp);
                        }
                    }

                    //
                    // Build node at unprimed level, unless skipped?
                    //
                    if (DONT_CARE != m.from(k)) {
                        unpacked_node* nb
                            = unpacked_node::newSparse(F, k, 1);
                        nb->setSparse(0, m.from(k), cv, cp);
                        F->createReducedNode(nb, cv, cp);
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
                MEDDLY_DCASSERT(L>0);
                MEDDLY_DCASSERT(m.isForRelations());
                MEDDLY_DCASSERT(!F->isIdentityReduced());
                for (int k=1; k<=L; k++) {
                    //
                    // Check for identity pattern at levels (k, k')
                    //
                    if (DONT_CHANGE == m.to(k)) {
                        MEDDLY_DCASSERT(DONT_CARE == m.from(k));
                        cp = F->makeIdentitiesTo(cp, k-1, k, ~0);
                        continue;
                    }

                    //
                    // Build node at primed level
                    //
                    if (DONT_CARE == m.to(k)) {
                        if (DONT_CARE == m.from(k)) {
                            cp = F->makeRedundantsTo(cp, k-1, k);
                            continue;
                        }
                        cp = F->makeRedundantsTo(cp, k-1, -k);
                    } else {
                        unpacked_node* nb
                            = unpacked_node::newSparse(F, k, 1);
                        nb->setSparse(0, m.to(k), cv, cp);
                        F->createReducedNode(nb, cv, cp);
                    }

                    //
                    // Build node at unprimed level
                    //
                    if (DONT_CARE == m.from(k)) {
                        cp = F->makeRedundantsTo(cp, -k, k);
                    } else {
                        unpacked_node* nb
                            = unpacked_node::newSparse(F, k, 1);
                        nb->setSparse(0, m.from(k), cv, cp);
                        F->createReducedNode(nb, cv, cp);
                    }
                } // for k
            }



            void createEdgeSet(int L, unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

            void createEdgeRel(int L, unsigned in,
                    unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

        private:
            /*
            inline void makeRedundantNode(int L, int L_size,
                    edge_value &dv, node_handle &dp)
            {
                if (F->isQuasiReduced()) {
                    // Build a redundant node at level L
                    unpacked_node* ru = unpacked_node::newFull(F, L, L_size);
                    ru->setFull(0, dv, dp);
                    for (unsigned v=1; v<L_size; v++) {
                        ru->setFull(v, dv, F->linkNode(dp));
                    }
                    F->createReducedNode(ru, dv, dp);
                }
            }
            */

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
    _from = new int[1+num_vars];
    if (for_relations) {
        _to = new int[1+num_vars];
    } else {
        _to = nullptr;
    }
}

MEDDLY::minterm::~minterm()
{
    delete[] _from;
    delete[] _to;
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
            set_or_rel sr) : max_coll_size(maxsz)
{
    _D = D;
    num_vars = D ? D->getNumVariables() : 0;
    for_relations = sr;

    _mtlist = new minterm*[max_coll_size];
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
    delete[] _mtlist;
}

void MEDDLY::minterm_coll::moveValueToFront(int L, int &minV,
        unsigned low, unsigned high, unsigned& mid)
{
    const int frontV = minV;
    minV = std::numeric_limits<int>::max();

    mid = low;
    if (L<0) {
        const int absL = -L;
        for (unsigned i=low; i<high; i++) {
            const int pri = primed(i, absL);
            if (frontV == pri) {
                if (mid != i) {
                    swap(mid, i);
                }
                ++mid;
            } else {
                minV = MIN(minV, pri);
            }
        }
    } else {
        for (unsigned i=low; i<high; i++) {
            const int unpr = unprimed(i, L);
            if (frontV == unpr) {
                if (mid != i) {
                    swap(mid, i);
                }
                ++mid;
            } else {
                minV = MIN(minV, unpr);
            }
        }
    }
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

    node_handle en;
    edge_value ev;

    if (e.getForest()->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        fbuilder<fbop_union_bool> fb(e.getForest(), *this, UNION);
        if (isForRelations()) {
            fb.createEdgeRel(num_vars, ~0, 0, first_unused, ev, en);
        } else {
            fb.createEdgeSet(num_vars, 0, first_unused, ev, en);
        }
        e.set(ev, en);
        return;
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

/*
void MEDDLY::minterm_coll::sortOnVariable(int L, unsigned low, unsigned high,
        std::vector<unsigned> &index, std::vector <int> &value)
{
#ifdef DEBUG_SORT
    FILE_output out(stdout);
    out << "Sorting variable " << L << " on range ["
        << low << ", " << high << ")\n";
#endif

    //
    // Determine the smallest and largest values
    //
    const int absL = ABS(L);
    int minV = std::numeric_limits<int>::max();
    int maxV = DONT_CHANGE;
    if (L<0) {
        for (unsigned i=low; i<high; i++) {
            const int pri = primed(i, absL);
            minV = MIN(minV, pri);
            maxV = MAX(maxV, pri);
        }
    } else {
        for (unsigned i=low; i<high; i++) {
            const int unpr = unprimed(i, absL);
            minV = MIN(minV, unpr);
            maxV = MAX(maxV, unpr);
        }
    }

    //
    // Initialize index, value
    //
    index.clear();
    index.push_back(low);
    value.clear();

    //
    // Sort. For each minimum value v, move those values to the front.
    // While doing that, keep track of the next smallest value.
    //
    unsigned indx_next = low;
    for (int v = minV; v<maxV; v = minV) {
        value.push_back(minV);
        minV = std::numeric_limits<int>::max();

        //
        // move anything with value v, to the "new" front
        //
        if (L<0) {
            for (unsigned i=index.back(); i<high; i++) {
                const int pri = primed(i, absL);
                if (v == pri) {
                    if (indx_next != i) {
                        swap(indx_next, i);
                    }
                    ++indx_next;
                } else {
                    minV = MIN(minV, pri);
                }
            }
        } else {
            for (unsigned i=index.back(); i<high; i++) {
                const int unpr = unprimed(i, absL);
                if (v == unpr) {
                    if (indx_next != i) {
                        swap(indx_next, i);
                    }
                    ++indx_next;
                } else {
                    minV = MIN(minV, unpr);
                }
            }
        }
        MEDDLY_DCASSERT(indx_next > index.back());
        index.push_back(indx_next);

    } // for v
    MEDDLY_DCASSERT(index.back() < high);
    index.push_back(high);
    value.push_back(minV);

#ifdef DEBUG_SORT
    out << "Done sort\n";
    out << "Indexes: [" << index[0];
    for (unsigned i=1; i<index.size(); i++) {
        out << ", " << index[i];
    }
    out << "]\n";
    out << "Values: [" << value[0];
    for (unsigned i=1; i<value.size(); i++) {
        out << ", " << value[i];
    }
    out << "]\n";
    for (unsigned i=index[0]; i<index.back(); i++) {
        out << "    ";
        out.put(long(i), 2);
        out << "  ";
        at(i).show(out);
        out << "\n";
    }
#endif
}
*/

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

template <class OP>
void MEDDLY::fbuilder<OP>::createEdgeSet(int L, unsigned low, unsigned high,
        edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(L>=0);
    MEDDLY_DCASSERT(high > low);

    //
    // Terminal case
    //
    if (0==L) {
        OP::finalize(mtc, low, high, cv, cp);
        return;
    }

    //
    // Special case: only one minterm left
    //
    if (high - low == 1) {
        OP::finalize(mtc, low, high, cv, cp);
        setPathToBottom(L, element(low), cv, cp);
        return;
    }

    //
    // Determine values
    //
    unsigned mid = low;
    int minV, maxV;
    mtc.getMinMax(L, low, high, minV, maxV);

    //
    // Special case: all values are DONT_CARE
    //
    if (DONT_CARE == maxV) {
        createEdgeSet(L-1, low, high, cv, cp);
        cp = F->makeRedundantsTo(cp, L-1, L);
        return;
    }

    //
    // Get ready for recursions
    //

    // Don't care stuff
    edge_value dnc_val, cz_val;
    node_handle dnc_node, cz_node;
    bool has_dont_care = false;

    // size of variables at level L
    unsigned L_size = F->getLevelSize(L);
    // sparse node we're going to build
    unpacked_node* Cu = unpacked_node::newSparse(F, L, L_size);
    // number of nonzero edges in our sparse node
    unsigned z = 0;

    //
    // Recursion loop
    //
    while (minV < maxV) {
        const int currV = minV;
        mtc.moveValueToFront(L, minV, low, high, mid);

        if (DONT_CARE == currV) {
            //
            // Build a redundant node at level L
            // but recurse to figure out where it points to.
            //
            createEdgeSet(L-1, low, mid, dnc_val, dnc_node);
            dnc_node = F->makeRedundantsTo(dnc_node, L-1, L);
            has_dont_care = true;

            // Make sure dnc_node isn't recycled yet
            Cu->setTempRoot(dnc_node);

        } else {
            //
            // Recurse and add result to unpacked node Cu
            //
            createEdgeSet(L-1, low, mid, cz_val, cz_node);

            //
            // add to sparse node, unless transparent
            //
            if (!F->isTransparentEdge(cz_val, cz_node)) {
                Cu->setSparse(z, currV, cz_val, cz_node);
                z++;
            }
        }

        low = mid;
    }
    MEDDLY_DCASSERT(maxV>=0);
    MEDDLY_DCASSERT(minV == maxV);
    //
    // Last value
    //
    createEdgeSet(L-1, mid, high, cz_val, cz_node);
    if (!F->isTransparentEdge(cz_val, cz_node)) {
        Cu->setSparse(z, minV, cz_val, cz_node);
        z++;
    }

    //
    // Finish building
    //
    Cu->shrink(z);
    F->createReducedNode(Cu, cv, cp);
    if (has_dont_care) {
        union_op->compute(L, ~0, dnc_val, dnc_node, cv, cp, cv, cp);
    }

}

// ******************************************************************

template <class OP>
void MEDDLY::fbuilder<OP>::createEdgeRel(int L, unsigned in,
        unsigned low, unsigned high, edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(L>=0);
    MEDDLY_DCASSERT(high > low);

    //
    // Terminal case
    //
    if (0==L) {
        OP::finalize(mtc, low, high, cv, cp);
        return;
    }

    //
    // Special case: only one minterm left
    //
    if (high - low == 1) {
        OP::finalize(mtc, low, high, cv, cp);
        if (F->isIdentityReduced()) {
            identityPathToBottom(L, element(low), cv, cp);
        } else {
            relPathToBottom(L, element(low), cv, cp);
        }
        return;
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

