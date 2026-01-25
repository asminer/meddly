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
#include "dd_edge.h"

// #define DEBUG_MOVE_VALUES
// #define DEBUG_MOVE_PAIRS
// #include "operators.h"

// #define DEBUG_CREATE_EDGE_SET
// #define DEBUG_CREATE_EDGE_REL

// ******************************************************************
// *                                                                *
// *                    fbuilder class hierarchy                    *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    class fbuilder_forest {
        public:
            ///
            /// Initialize for forest f, and default value <dv, dp>.
            ///
            ///     @param f    Forest we're building in
            ///     @param def  Default value
            ///
            fbuilder_forest(forest* f, const rangeval &def);
            ~fbuilder_forest();

            /// unlink any edges
            void destroy();

            /// Build a chain of nodes for a single "set" minterm.
            ///     @param L    Last level to build a node.
            ///     @param m    Minterm to build from
            ///     @param cv   Input: bottom value
            ///                 Output: top value
            ///     @param cp   Input: bottom node
            ///                 Output: top node
            ///
            void setPathToBottom(int L, const minterm &m,
                    edge_value &cv, node_handle &cp);

            /// Build a chain of nodes for a single "relation" minterm.
            ///
            ///     @param L    Last level to build a node.
            ///     @param m    Minterm to build from
            ///     @param cv   Input: bottom value
            ///                 Output: top value
            ///     @param cp   Input: bottom node
            ///                 Output: top node
            ///
            void relPathToBottom(int L, const minterm &m,
                    edge_value &cv, node_handle &cp);

        protected: // helpers

            /// Set up a fresh unpacked node, in a set,
            /// with all values equal to the default value.
            ///
            ///     @param  k       Level
            ///     @param  spsz    Max number of non-default values
            ///
            inline unpacked_node* newSetNode(unsigned k, unsigned spsz)
                const
                {
                    if (default_is_zero) {
                        return unpacked_node::newWritable(F, k, spsz,
                                    SPARSE_ONLY);
                    }
                    MEDDLY_DCASSERT(dp_unp);
                    unpacked_node* n =
                        unpacked_node::newWritable(F, k, FULL_ONLY);
                    for (unsigned i=0; i<n->getSize(); i++) {
                        n->setFull(i, dv, F->linkNode(dp_unp[k-1]));
                    }
                    return n;
                }

            /// Set up a fresh unpacked node, in a relation,
            /// at an unprimed level,
            /// with all values equal to the default value.
            ///
            ///     @param  k       Level
            ///     @param  spsz    Max number of non-default values
            ///
            inline unpacked_node* newUnprimedNode(int k, unsigned spsz)
                const
                {
                    if (default_is_zero) {
                        return unpacked_node::newWritable(F, k, spsz,
                                    SPARSE_ONLY);
                    }
                    MEDDLY_DCASSERT(dp_pri);
                    unpacked_node* n =
                        unpacked_node::newWritable(F, k, FULL_ONLY);
                    for (unsigned i=0; i<n->getSize(); i++) {
                        n->setFull(i, dv, F->linkNode(dp_pri[k]));
                    }
                    return n;
                }

            /// Set up a fresh unpacked node, in a relation,
            /// at a primed level,
            /// with all values equal to the default value.
            ///
            ///     @param  k       Level
            ///     @param  spsz    Max number of non-default values
            ///
            inline unpacked_node* newPrimedNode(int k, unsigned spsz)
                const
                {
                    if (default_is_zero) {
                        return unpacked_node::newWritable(F, -k, spsz,
                                    SPARSE_ONLY);
                    }
                    MEDDLY_DCASSERT(dp_unp);
                    unpacked_node* n =
                        unpacked_node::newWritable(F, -k, FULL_ONLY);
                    for (unsigned i=0; i<n->getSize(); i++) {
                        n->setFull(i, dv, F->linkNode(dp_unp[k-1]));
                    }
                    return n;
                }

            /// Add an edge to an unpacked node
            ///     @param  z   Nonzero pointer, if needed
            ///     @param  ndx Index
            ///     @param  cv  Value
            ///     @param  cn  Node
            inline void addToNode(unpacked_node* n, unsigned &z, int ndx,
                    const edge_value &cv, node_handle cn) const
            {
                if (default_is_zero) {
                    // n is sparse
                    if (!F->isTransparentEdge(cv, cn)) {
                        n->setSparse(z, ndx, cv, cn);
                        ++z;
                    }
                } else {
                    // n is full; unlink old pointer first
                    F->unlinkNode(n->down(ndx));
                    n->setFull(ndx, cv, cn);
                }
            }

            /// Build an identity pattern.
            /// Except here we might have default values, not zeroes,
            /// off the diagonals.
            ///     @param  k   Build the pattern for levels (k, k')
            ///     @param  dn  Downward pointer for diagonals
            ///
            inline void identityPattern(int k, edge_value &dv, node_handle &dn)
                const
            {
                if (default_is_zero) {
                    dn = F->makeIdentitiesTo(dn, k-1, k, ~0);
                } else {
                    //
                    // Build an identity pattern except use
                    // the default value off the diagonal.
                    //
                    MEDDLY_DCASSERT(dp_unp);
                    unpacked_node* nu =
                        unpacked_node::newWritable(F, k, FULL_ONLY);
                    for (unsigned i=0; i<nu->getSize(); i++) {
                        unpacked_node* np = newPrimedNode(k, 1);
                        unsigned z = 0;
                        addToNode(np, z, i, dv, F->linkNode(dn));

                        edge_value cv;
                        node_handle cp;
                        F->createReducedNode(np, cv, cp);
                        nu->setFull(i, cv, cp);
                    }
                    F->unlinkNode(dn);
                    F->createReducedNode(nu, dv, dn);
                }
            }

        protected:
            forest* F;

            //
            // Default value.
            //
            rangeval deflt;
            bool default_is_zero;

            //
            // Default value, node at each level.
            // If we're quasi-reduced, these will be actual nodes;
            // otherwise, they will all be terminals.
            //
            // For unprimed levels: use <dv, dp_unp[k]>
            // For   primed levels: use <dv, dp_pri[k]>
            //
            edge_value   dv;
            node_handle* dp_unp;
            node_handle* dp_pri;
    };

    class fbuilder_common : public fbuilder_forest {
        public:
            fbuilder_common(forest* f, const rangeval &def,
                    minterm_coll &mtl, binary_builtin0 Union)
                    : fbuilder_forest(f, def), mtc(mtl)
            {
                union_op = build(Union, f, f, f);
            }

#ifdef ALLOW_OLD_BINARY_0_17_6
            fbuilder_common(forest* f, const rangeval &def,
                    minterm_coll &mtl, binary_builtin2 Union)
                    : fbuilder_forest(f, def), mtc(mtl)
            {
                union_op = Union(f, f, f);
            }
#endif

            /// Get the minimum and maximum values for level L,
            /// on the interval [low, high).
            ///
            inline void getMinMax(int L, unsigned low, unsigned hi,
                int &minV, int &maxV) const
            {
                MEDDLY_DCASSERT(!mtc.isForRelations());
                MEDDLY_DCASSERT(L>0);
                minV = std::numeric_limits<int>::max();
                maxV = DONT_CHANGE;
                for (unsigned i=low; i<hi; i++) {
                    const int unpr = mtc.unprimed(i, L);
                    minV = MIN(minV, unpr);
                    maxV = MAX(maxV, unpr);
                }
            }


            /// Arrange minterms on interval [low, hi)
            /// into two intervals: [low, mid) and [mid, hi)
            /// where [low, mid) all have level L value == minV,
            /// and   [mid, hi)  all have level L value != minv.
            ///
            /// On output, minV is the smallest level L value on [mid, hi).
            ///
            ///     @param  L       Level to sort on.
            ///     @param  minV    Input: Value to sort on.
            ///                     Output: smallest value on [mid, hi).
            ///     @param  low     low end of input interval
            ///     @param  hi      high end of input interval
            ///
            ///     @param  mid     On output: determines where to
            ///                     split [low, hi).
            ///
            void moveValuesToFront(int L, int &minV,
                    unsigned low, unsigned hi, unsigned& mid);


            /// Pair comparison: is a,b < c,d
            static inline bool pairLT(int a, int b, int c, int d)
            {
                if (a<c) return true;
                if (a>c) return false;
                // a == c
                return b<d;
            }
            /// Pair comparison: is a,b <= c,d
            static inline bool pairLE(int a, int b, int c, int d)
            {
                if (a<c) return true;
                if (a>c) return false;
                // a == c
                return b<=d;
            }

            /// Get the minimum and maximum pairs for levels L,L'
            /// on the interval [low, high).
            ///
            inline void getMinMax(int L, unsigned low, unsigned hi,
                    int &minU, int &minP, int &maxU, int &maxP) const
            {
                MEDDLY_DCASSERT(mtc.isForRelations());
                MEDDLY_DCASSERT(L>0);
                minU = std::numeric_limits<int>::max();
                minP = std::numeric_limits<int>::max();
                maxU = -10;
                maxP = -10;
                for (unsigned i=low; i<hi; i++) {
                    const int unp = mtc.unprimed(i, L);
                    const int pri = mtc.primed(i, L);
                    if (pairLT(unp, pri, minU, minP)) {
                        minU = unp;
                        minP = pri;
                    }
                    if (pairLT(maxU, maxP, unp, pri)) {
                        maxU = unp;
                        maxP = pri;
                    }
                }
            }


            /// The relation version of moveValueToFront().
            /// Arrange minterms on interval [low, hi)
            /// into two intervals: [low, mid) and [mid, hi)
            /// where [low, mid) all have level L,L' pairs  = (mU, mP),
            /// and   [mid, hi)  all have level L,L' pairs != (mU, mP).
            ///
            /// On output, (mU, mP) is the smallest level L,L'
            /// pair on [mid, hi).
            ///
            ///     @param  L       Level to sort on.
            ///     @param  mU      Input: Value to sort on, unprimed vars.
            ///                     Output: next value to sort on
            ///     @param  mP      Input: Value to sort on, primed vars.
            ///                     Output: next value to sort on
            ///
            ///     @param  low     low end of input interval
            ///     @param  hi      high end of input interval
            ///
            ///     @param  mid     On output: determines where to
            ///                     split [low, hi).
            ///
            void movePairsToFront(int L, int &mU, int &mP,
                    unsigned low, unsigned hi, unsigned& mid);

        protected:
#ifdef DEBUG_MOVE_VALUES
            static inline output& show_val(output &os, int val) {
                if (val >= 0) {
                    return os << val;
                }
                return os << 'x';
            }
#endif
#ifdef DEBUG_MOVE_PAIRS
            static inline output& show_pair(output &os, int from, int to) {
                if (from >= 0) {
                    os << '(' << from << "->";
                } else {
                    os << "(x->";
                }
                if (to >= 0) {
                    return os << to << ')';
                }
                if (DONT_CARE == to) {
                    return os << "x)";
                }
                return os << "i)";
            }
#endif

        protected:
            // computes <cv, cp> = union(<cv, cp>, <av, ap>)
            inline void accumulate(int L, unsigned in,
                    const edge_value &av, node_handle& ap,
                    edge_value &cv, node_handle& cp)
            {
                if (union_op) {
                    const node_handle orig_cp = cp;
                    union_op->compute(L, in, av, ap, cv, cp, cv, cp);
                    F->unlinkNode(orig_cp);
                    F->unlinkNode(ap);
                    ap = 0;
                } else {
                    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
                }
            }

        protected:
            minterm_coll &mtc;

        private:
            binary_operation* union_op;

    };

    /*
        To build a function for a collection of minterms, struct OP must
        provide the following static method:

            void finalize(const minterm_coll &mc, unsigned low, unsigned high,
                const forest* F, edge_value &cv, node_handle &cp);

            Determine the terminal edge for the minterms in mc with
            indexes in [low, high). Store the result in <cv, cp>.

    */
    template <class OP>
    class fbuilder : public fbuilder_common {
        public:
            fbuilder(forest* f, const rangeval &def,
                    minterm_coll &mtl, binary_builtin0 Union)
                    : fbuilder_common(f, def, mtl, Union) { }

#ifdef ALLOW_OLD_BINARY_0_17_6
            fbuilder(forest* f, const rangeval &def,
                    minterm_coll &mtl, binary_builtin2 Union)
                    : fbuilder_common(f, def, mtl, Union) { }
#endif

            inline void createEdge(int L, unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp)
            {
                if (mtc.isForRelations()) {
                    createEdgeRel(L, low, high, cv, cp);
                } else {
                    createEdgeSet(L, low, high, cv, cp);
                }
            }

        protected:
            void createEdgeSet(int L, unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);

            void createEdgeRel(int L, unsigned low, unsigned high,
                    edge_value &cv, node_handle &cp);
    };

    // ******************************************************************
    // *                                                                *
    // *                 operation classes for fbuilder                 *
    // *                                                                *
    // ******************************************************************

    template <typename T>
    struct fbop_min_tmpl {
        static inline void finalize(
                const minterm_coll &mc, unsigned low, unsigned high,
                const forest* F, edge_value &cv, node_handle &cp)
        {
            rangeval val;
            val = mc.at(low).getValue();
            for (unsigned i=low+1; i<high; i++) {
                const rangeval &mci = mc.at(i).getValue();
                if (mci.isPlusInfinity()) {
                    continue;
                }
                if (val.isPlusInfinity()) {
                    val = mci;
                    continue;
                }
                if (T(mci) < T(val)) {
                    val = mci;
                }
            }
            F->getEdgeForValue(val, cv, cp);
        }
    };
    template <typename T>
    struct fbop_max_tmpl {
        static inline void finalize(
                const minterm_coll &mc, unsigned low, unsigned high,
                const forest* F, edge_value &cv, node_handle &cp)
        {
            rangeval val;
            val = mc.at(low).getValue();
            for (unsigned i=low+1; i<high; i++) {
                if (val.isPlusInfinity()) break;
                const rangeval &mci = mc.at(i).getValue();
                if (mci.isPlusInfinity()) {
                    val = mci;
                    break;
                }
                if (T(mci) > T(val)) {
                    val = T(mci);
                }
            }
            F->getEdgeForValue(val, cv, cp);
        }
    };
};

// ******************************************************************
// *                                                                *
// *                        minterm  methods                        *
// *                                                                *
// ******************************************************************


MEDDLY::minterm::minterm(const domain* D, set_or_rel sr) : value(true)
{
    _D = D;
    for_relations = sr;
    initVectors();
}

MEDDLY::minterm::minterm(const forest* F) : value(true)
{
    if (F) {
        _D = F->getDomain();
        for_relations = F->isForRelations();
    } else {
        _D = nullptr;
        for_relations = false;
    }
    initVectors();
}

MEDDLY::minterm::minterm(const minterm &m) : value(true)
{
    _D = m._D;
    for_relations = m.for_relations;
    initVectors();
    if (for_relations) {
        setAll(m._from, m._to, m.value);
    } else {
        setAll(m._from, m.value);
    }
}


MEDDLY::minterm::~minterm()
{
    delete[] _from;
    delete[] _to;
}

void MEDDLY::minterm::setFrom(const minterm &m)
{
    if (_D != m._D) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (for_relations != m.for_relations) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (for_relations) {
        setAll(m._from, m._to, m.value);
    } else {
        setAll(m._from, m.value);
    }
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
    s.put(", bot]  ");
    value.write(s);
}


void MEDDLY::minterm::buildFunction(rangeval def, dd_edge &e) const
{
    forest* F = e.getForest();
    if (!F) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if (F->getDomain() != _D) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (F->isForRelations() != isForRelations())
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    node_handle cp;
    edge_value cv;

    F->getEdgeForValue(value, cv, cp);

    fbuilder_forest fb(F, def);

    if (isForRelations()) {
        fb.relPathToBottom(num_vars, *this, cv, cp);
    } else {
        fb.setPathToBottom(num_vars, *this, cv, cp);
    }

    e.set(cv, cp);
#ifdef DEVELOPMENT_CODE
    fb.destroy();   // destroy early, for counts
    F->validateIncounts(true, __FILE__, __LINE__, "minterm::buildFunction");
#endif
}

void MEDDLY::minterm::initVectors()
{
    if (!_D) {
        num_vars = 0;
        _from = nullptr;
        _to = nullptr;
        return;
    }

    num_vars = _D->getNumVariables();
    _from = new int[1+num_vars];
    if (for_relations) {
        _to = new int[1+num_vars];
    } else {
        _to = nullptr;
    }
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

MEDDLY::minterm_coll::minterm_coll(unsigned maxsz, const forest* F)
        : max_coll_size(maxsz)
{
    if (!F) {
        throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    }
    _D = F->getDomain();
    num_vars = _D ? _D->getNumVariables() : 0;
    for_relations = F->isForRelations();

    _mtlist = new minterm*[max_coll_size];
    for (unsigned i=0; i<maxsize(); i++) {
        _mtlist[i] = new minterm(F);
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


void MEDDLY::minterm_coll::buildFunctionMax(rangeval def, dd_edge &e)
{
    forest* F = e.getForest();
    if (!F) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if (F->getDomain() != _D) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (F->isForRelations() != isForRelations())
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    //
    // <ev, en> is the function we're building
    //
    node_handle en;
    edge_value ev;

    /*
     * Special case: empty collection
     */
    if (0==first_unused) {
        F->getEdgeForValue(def, ev, en);
        en = F->makeRedundantsTo(en, 0, num_vars);
        e.set(ev, en);
        return;
    }

    /*
     * Determine operation, based on the forest type.
     */
    switch (F->getRangeType()) {
        case range_type::BOOLEAN:
        {
            // fbuilder<fbop_union_bool> fb(F, def, *this, UNION);
            fbuilder< fbop_max_tmpl<bool> > fb(F, def, *this, UNION);
            fb.createEdge(num_vars, 0, first_unused, ev, en);
            break;
        }

        case range_type::INTEGER:
        {
            fbuilder< fbop_max_tmpl<long> > fb(F, def, *this, MAXIMUM);
            fb.createEdge(num_vars, 0, first_unused, ev, en);
            break;
        }

        case range_type::REAL:
        {
            fbuilder< fbop_max_tmpl<double> > fb(F, def, *this, MAXIMUM);
            fb.createEdge(num_vars, 0, first_unused, ev, en);
            break;
        }
    }

    e.set(ev, en);
#ifdef DEVELOPMENT_CODE
    F->validateIncounts(true, __FILE__, __LINE__,
            "minterm_coll::buildFunctionMax");
#endif
}

void MEDDLY::minterm_coll::buildFunctionMin(rangeval def, dd_edge &e)
{
    forest* F = e.getForest();
    if (!F) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if (F->getDomain() != _D) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (F->isForRelations() != isForRelations())
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    //
    // <ev, en> is the function we're building
    //
    node_handle en;
    edge_value ev;

    /*
     * Special case: empty collection
     */
    if (0==first_unused) {
        F->getEdgeForValue(def, ev, en);
        en = F->makeRedundantsTo(en, 0, num_vars);
        e.set(ev, en);
        return;
    }

    /*
     * Determine operation, based on the forest type.
     */
    switch (F->getRangeType()) {
        case range_type::BOOLEAN:
        {
            // fbuilder<fbop_union_bool> fb(F, def, *this, INTERSECTION);
            fbuilder< fbop_min_tmpl <bool> > fb(F, def, *this, INTERSECTION);
            fb.createEdge(num_vars, 0, first_unused, ev, en);
            break;
        }

        case range_type::INTEGER:
        {
            fbuilder< fbop_min_tmpl<long> > fb(F, def, *this, MINIMUM);
            fb.createEdge(num_vars, 0, first_unused, ev, en);
            break;
        }

        case range_type::REAL:
        {
            fbuilder< fbop_min_tmpl<double> > fb(F, def, *this, MINIMUM);
            fb.createEdge(num_vars, 0, first_unused, ev, en);
            break;
        }
    }

    e.set(ev, en);
#ifdef DEVELOPMENT_CODE
    F->validateIncounts(true, __FILE__, __LINE__,
            "minterm_coll::buildFunctionMin");
#endif
}


void MEDDLY::minterm_coll::show(output &s,
        const char* pre, const char* post) const
{
    for (unsigned i=0; i<size(); i++) {
        if (pre) {
            s.put(pre);
        }
        s.put(i, 4);
        s.put("  ");
        _mtlist[i]->show(s);
        s.put(post);
    }
}

void MEDDLY::minterm_coll::_expand()
{
    unsigned new_size;
    if (0==max_coll_size) {
        new_size = 16;
    } else {
        new_size = 2*max_coll_size;
    }
    minterm** new_list = new minterm*[new_size];
    unsigned i;
    for (i=0; i<max_coll_size; i++) {
        new_list[i] = _mtlist[i];
    }
    for (   ; i<new_size; i++) {
        new_list[i] = new minterm(_D, for_relations);
    }

    delete[] _mtlist;
    _mtlist = new_list;
    max_coll_size = new_size;
}

// ******************************************************************
// *                                                                *
// *                    fbuilder_forest  methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::fbuilder_forest::fbuilder_forest(forest* f, const rangeval &def)
{
    F = f;
    deflt = def;

    node_handle dp;
    F->getEdgeForValue(deflt, dv, dp);

    default_is_zero = F->isTransparentEdge(dv, dp);

    if (default_is_zero) {
        dp_unp = nullptr;
        dp_pri = nullptr;
    } else {
        const unsigned L = F->getNumVariables();
        if (F->isForRelations()) {
            dp_unp = new node_handle[L+1];
            dp_pri = new node_handle[L+1];
            dp_unp[0] = dp;
            dp_pri[0] = 0;
            for (int k=1; k<=L; k++) {
                dp_pri[k] = F->makeRedundantsTo(dp_unp[k-1], k-1, -k);
                dp_unp[k] = F->makeRedundantsTo(dp_pri[k], -k, k);
            }
        } else {
            dp_unp = new node_handle[L+1];
            dp_pri = nullptr;
            dp_unp[0] = dp;
            for (unsigned k=1; k<=L; k++) {
                dp_unp[k] = F->makeRedundantsTo(dp_unp[k-1], k-1, k);
            }
        }
    }
}

MEDDLY::fbuilder_forest::~fbuilder_forest()
{
    destroy();
    if (dp_unp) {
        delete[] dp_unp;
    }
    if (dp_pri) {
        delete[] dp_pri;
    }
}

void MEDDLY::fbuilder_forest::destroy()
{
    if (dp_unp) {
        const unsigned L = F->getNumVariables();
        F->unlinkNode(dp_unp[L]);
        dp_unp[L] = 0;
    }
}

void MEDDLY::fbuilder_forest::setPathToBottom(int L, const minterm &m,
        edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(L>0);
    MEDDLY_DCASSERT(!m.isForRelations());

    for (int k=1; k<=L; k++) {
        if (DONT_CARE == m.from(k)) {
            cp = F->makeRedundantsTo(cp, k-1, k);
        } else {
            // make a node with one edge <cv, cp>,
            // all the rest <dv, dp>
            MEDDLY_DCASSERT(m.from(k) >= 0);

            unpacked_node* n = newSetNode(k, 1);
            unsigned z=0;
            addToNode(n, z, m.from(k), cv, cp);

            F->createReducedNode(n, cv, cp);
        }
    } // for i
}

void MEDDLY::fbuilder_forest::relPathToBottom(int L, const minterm &m,
        edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(L>0);
    MEDDLY_DCASSERT(m.isForRelations());

    if (F->isIdentityReduced() && default_is_zero) {
        //
        //
        // Special case: we're identity reduced, and <dv, dp> is the zero edge
        //
        //
        for (int k=1; k<=L; k++) {
            //
            // Check for identity pattern at levels (k, k')
            //
            if (DONT_CHANGE == m.to(k)) {
                MEDDLY_DCASSERT(DONT_CARE == m.from(k));
                continue;
            }

            //
            // Build node at primed level
            //
            if (DONT_CARE == m.to(k)) {
                cp = F->makeRedundantsTo(cp, k-1, -k);
            } else {
                if (m.from(k) != m.to(k)) {
                    unpacked_node* n = newPrimedNode(k, 1);
                    unsigned z=0;
                    addToNode(n, z, m.to(k), cv, cp);
                    F->createReducedNode(n, cv, cp);
                }
            }

            //
            // Build node at unprimed level
            //
            if (DONT_CARE == m.from(k)) {
                cp = F->makeRedundantsTo(cp, -k, k);
            } else {
                unpacked_node* n = newUnprimedNode(k, 1);
                unsigned z=0;
                addToNode(n, z, m.from(k), cv, cp);
                F->createReducedNode(n, cv, cp);
            }
        } // for k
    } else {
        //
        //
        // General case
        //
        //
        for (int k=1; k<=L; k++) {
            //
            // Check for identity pattern at levels (k, k')
            //
            if (DONT_CHANGE == m.to(k)) {
                MEDDLY_DCASSERT(DONT_CARE == m.from(k));

                identityPattern(k, cv, cp);

                continue;
            }

            //
            // Process the primed level only.
            //
            if (DONT_CARE == m.to(k)) {
                cp = F->makeRedundantsTo(cp, k-1, -k);
            } else {
                // make a node with one edge <cv, cp>,
                // all the rest <dv, dp>
                MEDDLY_DCASSERT(m.to(k) >= 0);

                unpacked_node* np = newPrimedNode(k, 1);
                unsigned z=0;
                addToNode(np, z, m.to(k), cv, cp);
                F->createReducedNode(np, cv, cp);
            }

            //
            // Process the unprimed level.
            //
            if (DONT_CARE == m.from(k)) {
                cp = F->makeRedundantsTo(cp, -k, k);
            } else {
                // make a node with one edge <cv, cp>,
                // all the rest <dv, dp>
                MEDDLY_DCASSERT(m.from(k) >= 0);

                unpacked_node* nu = newUnprimedNode(k, 1);
                unsigned z=0;
                addToNode(nu, z, m.from(k), cv,
                    F->redirectSingleton(m.from(k), cp)
                );

                F->createReducedNode(nu, cv, cp);
            }

        } // for k
    }
}


// ******************************************************************
// *                                                                *
// *                    fbuilder_common  methods                    *
// *                                                                *
// ******************************************************************

void MEDDLY::fbuilder_common::moveValuesToFront(int L, int &minV,
        unsigned low, unsigned high, unsigned& mid)
{
    MEDDLY_DCASSERT(!mtc.isForRelations());
    MEDDLY_DCASSERT(L>0);

#ifdef DEBUG_MOVE_VALUES
    FILE_output out(stdout);
    out << "starting moveValuesToFront(" << L << ", ";
    show_val(out, minV) << ", "
        << low << ", " << high << ")\n";
#endif

    const int frontV = minV;
    minV = std::numeric_limits<int>::max();

    mid = low;
    for (unsigned i=low; i<high; i++) {
        const int unpr = mtc.unprimed(i, L);
        if (frontV == unpr) {
            if (mid != i) {
                mtc.swap(mid, i);
            }
            ++mid;
        } else {
            minV = MIN(minV, unpr);
        }
    }

#ifdef DEBUG_MOVE_VALUES
    out << "Done move; mid=" << mid << ", minV=";
    show_val(out, minV) << "\n";
    for (unsigned i=low; i<high; i++) {
        out << "    ";
        out.put(long(i), 2);
        out << "  ";
        mtc.at(i).show(out);
        out << "\n";
    }
#endif

}

void MEDDLY::fbuilder_common::movePairsToFront(int L, int &mU, int &mP,
                    unsigned low, unsigned high, unsigned& mid)
{
    MEDDLY_DCASSERT(mtc.isForRelations());
    MEDDLY_DCASSERT(L>0);

#ifdef DEBUG_MOVE_PAIRS
    FILE_output out(stdout);
    out << "starting movePairsToFront(" << L << ", ";
    show_pair(out, mU, mP) << ", " << low << ", " << high << ")\n";
#endif

    const int frontU = mU;
    const int frontP = mP;
    mU = std::numeric_limits<int>::max();
    mP = std::numeric_limits<int>::max();

    mid = low;
    for (unsigned i=low; i<high; i++) {
        const int unp = mtc.unprimed(i, L);
        const int pri = mtc.primed(i, L);
        if ((frontU == unp) && (frontP == pri)) {
            if (mid != i) {
                mtc.swap(mid, i);
            }
            ++mid;
        } else {
            if (pairLT(unp, pri, mU, mP)) {
                mU = unp;
                mP = pri;
            }
        }
    }

#ifdef DEBUG_MOVE_PAIRS
    out << "Done move; mid=" << mid << ", min=";
    show_pair(out, mU, mP) << "\n";
    for (unsigned i=low; i<high; i++) {
        out << "    ";
        out.put(long(i), 2);
        out << "  ";
        mtc.at(i).show(out);
        out << "\n";
    }
#endif
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
        OP::finalize(mtc, low, high, F, cv, cp);
        return;
    }

    //
    // Special case: only one minterm left
    //
    if (high - low == 1) {
        OP::finalize(mtc, low, high, F, cv, cp);
        setPathToBottom(L, mtc.at(low), cv, cp);
        return;
    }

    //
    // Determine values
    //
    unsigned mid = low;
    int minV, maxV;
    getMinMax(L, low, high, minV, maxV);

    //
    // Special case: all values are DONT_CARE
    //
    if (DONT_CARE == maxV) {
        MEDDLY_DCASSERT(DONT_CARE == minV);
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
    const unsigned L_size = F->getLevelSize(L);
    unsigned z = 0;
    unpacked_node* Cu = newSetNode(L, L_size);

    //
    // Recursion loop
    //
    while (minV <= maxV) {
        const int currV = minV;
        if (minV < maxV) {
            // More interval to split
            moveValuesToFront(L, minV, low, high, mid);
        } else {
            // We're on the last interval
            ++minV;
            mid = high;
        }

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
            // add to unpacked node
            //
            addToNode(Cu, z, currV, cz_val, cz_node);
        }

        low = mid;
    }

    //
    // Finish building
    //
    if (Cu->isSparse()) {
        Cu->shrink(z);
    }
    F->createReducedNode(Cu, cv, cp);
    if (has_dont_care) {
        accumulate(L, ~0, dnc_val, dnc_node, cv, cp);
    }

}

// ******************************************************************

template <class OP>
void MEDDLY::fbuilder<OP>::createEdgeRel(int L, unsigned low, unsigned high,
        edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(L>=0);
    MEDDLY_DCASSERT(high > low);

    //
    // Terminal case
    //
    if (0==L) {
        OP::finalize(mtc, low, high, F, cv, cp);
        return;
    }

    //
    // Special case: only one minterm left
    //
    if (high - low == 1) {
        OP::finalize(mtc, low, high, F, cv, cp);
        relPathToBottom(L, mtc.at(low), cv, cp);
        return;
    }

    //
    // Determine values
    //
    unsigned mid = low;
    int minU, minP, maxU, maxP;
    getMinMax(L, low, high, minU, minP, maxU, maxP);

    //
    // Special case: all pairs are (DONT_CARE, DONT_CARE)
    //
    if ( (DONT_CARE == maxU) && (DONT_CARE == maxP) && (DONT_CARE == minP) )
    {
        MEDDLY_DCASSERT(DONT_CARE == minU);
        createEdgeRel(L-1, low, high, cv, cp);
        cp = F->makeRedundantsTo(cp, L-1, L);
        return;
    }

    //
    // Special case: all pairs are (DONT_CARE, DONT_CHANGE)
    //
    if ( (DONT_CARE == maxU) && (DONT_CHANGE == maxP) )
    {
        MEDDLY_DCASSERT(DONT_CARE == minU);
        createEdgeRel(L-1, low, high, cv, cp);
        identityPattern(L, cv, cp);
        return;
    }

    //
    // Get ready for recursions.
    //
    //  Unprimed level variables:
    //      <ue_v, ue_p>    : unprimed extra function for union
    //      Cu              : unpacked node
    //      zu              : sparse count for Cu
    //
    //  Primed level variables:
    //      <pe_v, pe_p>    : primed extra function for union
    //      Cp              : unpacked node
    //      zp              : sparse count for Cp
    //      cpin            : unprimed value pointing to Cp
    //
    //  Other variables:
    //      <tv, tp>        :   generic temp function
    //

    const unsigned L_size = F->getLevelSize(L);

    edge_value ue_v, pe_v, tv;
    node_handle ue_p, pe_p, tp;
    bool has_unprimed_extra = false;
    bool has_primed_extra = false;

    unpacked_node* Cu = newUnprimedNode(L, L_size);
    unpacked_node* Cp = newPrimedNode(L, L_size);
    unsigned zu = 0;
    unsigned zp = 0;
    int cpin = -10;

    //
    // Recursion loop
    //
    while (pairLE(minU, minP, maxU, maxP)) {
        const int currU = minU;
        const int currP = minP;
        if (pairLT(minU, minP, maxU, maxP)) {
            // More interval to split
            movePairsToFront(L, minU, minP, low, high, mid);
        } else {
            // We're on the last interval
            ++minU;
            mid = high;
        }

        //
        // Start another primed node, or not?
        //
        if (cpin < DONT_CHANGE) {
            // First time through
            cpin = currU;
#ifdef DEBUG_CREATE_EDGE_REL
            printf("Starting Cp level %d in %d\n", -L, cpin);
#endif
        }
        if (cpin != currU) {
            //
            // Done with Cp; check it in and connect it
            //
            if (Cp->isSparse()) {
                Cp->shrink(zp);
            }
            F->createReducedNode(Cp, tv, tp, cpin);

            if (has_primed_extra) {
                MEDDLY_DCASSERT(cpin >= 0);
                accumulate(-L, cpin, pe_v, pe_p, tv, tp);
            }
            if (DONT_CARE == cpin) {
                //
                // Add to the extra unprimed function
                //
                tp = F->makeRedundantsTo(tp, -L, L);
                if (has_unprimed_extra) {
                    accumulate(L, ~0, tv, tp, ue_v, ue_p);
                } else {
                    ue_v = tv;
                    ue_p = tp;
                    has_unprimed_extra = true;
                }
                Cu->setTempRoot(ue_p);
            } else {
                //
                // Connect to cpin in Cu
                //
                addToNode(Cu, zu, cpin, tv, tp);
            }
            //
            // Reset Cp
            //
            Cp = newPrimedNode(L, L_size);
            zp = 0;
            has_primed_extra = false;
            cpin = currU;
#ifdef DEBUG_CREATE_EDGE_REL
            printf("Starting Cp level %d in %d\n", -L, cpin);
#endif
        }

        if (DONT_CARE == currU)
        {
            if (DONT_CHANGE == currP || DONT_CARE == currP)
            {
                //
                // Build an identity (matrix) node at level L
                // or a don't care (matrix) node at level L,
                // but recurse to figure out where it points to.
                //
                createEdgeRel(L-1, low, mid, tv, tp);
                if (DONT_CHANGE == currP) {
                    identityPattern(L, tv, tp);
                } else {
                    tp = F->makeRedundantsTo(tp, L-1, L);
                }
                //
                // Add this to the UNPRIMED extra function
                //
                if (has_unprimed_extra) {
                    accumulate(L, ~0, tv, tp, ue_v, ue_p);
                } else {
                    ue_v = tv;
                    ue_p = tp;
                    has_unprimed_extra = true;
                }

                // Make sure ue_p isn't recycled
                Cu->setTempRoot(ue_p);

            } else {
                //
                // Add to Cp
                //
                createEdgeRel(L-1, low, mid, tv, tp);
                addToNode(Cp, zp, currP, tv, tp);
#ifdef DEBUG_CREATE_EDGE_REL
                printf("Adding to Cp level %d in %d\n", -L, cpin);
                printf("  down[%d] = %d level %d\n",
                        currP, tp, F->getNodeLevel(tp));
#endif
            }

        } else {
            //
            // Unprimed has an actual value.
            // But we still could have DONT_CARE for the primed value.
            //
            MEDDLY_DCASSERT(currP != DONT_CHANGE);
            if (DONT_CARE == currP)
            {
                //
                // Recurse and add to the PRIMED extra function
                //
                createEdgeRel(L-1, low, mid, tv, tp);
                tp = F->makeRedundantsTo(tp, L-1, -L);
                if (has_primed_extra) {
                    accumulate(-L, cpin, tv, tp, pe_v, pe_p);
                } else {
                    pe_v = tv;
                    pe_p = tp;
                    has_primed_extra = true;
                }

                // Make sure ue_p isn't recycled
                Cp->setTempRoot(pe_p);

            } else {
                //
                // Add to Cp
                //
                createEdgeRel(L-1, low, mid, tv, tp);
                addToNode(Cp, zp, currP, tv, tp);
#ifdef DEBUG_CREATE_EDGE_REL
                printf("Adding to Cp level %d in %d\n", -L, cpin);
                printf("  down[%d] = %d level %d\n",
                        currP, tp, F->getNodeLevel(tp));
#endif
            }
        }

        low = mid;
    }

    //
    // Close off the final Cp
    //
    if (Cp->isSparse()) {
        Cp->shrink(zp);
    }
    F->createReducedNode(Cp, tv, tp, cpin);
    if (has_primed_extra) {
        MEDDLY_DCASSERT(cpin >= 0);
        accumulate(-L, cpin, pe_v, pe_p, tv, tp);
    }
    if (DONT_CARE == cpin) {
        //
        // Add this to the UNPRIMED extra function
        //
        tp = F->makeRedundantsTo(tp, -L, L);
        if (has_unprimed_extra) {
            accumulate(L, ~0, tv, tp, ue_v, ue_p);
        } else {
            ue_v = tv;
            ue_p = tp;
            has_unprimed_extra = true;
        }
        Cu->setTempRoot(ue_p);
    } else {
        //
        // Connect to cpin in Cu
        //
        addToNode(Cu, zu, cpin, tv, tp);
    }

    //
    // Finish Cu
    //
    if (Cu->isSparse()) {
        Cu->shrink(zu);
    }
    F->createReducedNode(Cu, cv, cp);
#ifdef DEBUG_CREATE_EDGE_REL
    printf("Reduced Cu level %d -> node %d at level %d\n",
            L, cp, F->getNodeLevel(cp));
#endif
    if (has_unprimed_extra) {
        accumulate(L, ~0, ue_v, ue_p, cv, cp);
    }
#ifdef DEBUG_CREATE_EDGE_REL
    printf("After union level %d -> node %d at level %d\n",
            L, cp, F->getNodeLevel(cp));
#endif
}

