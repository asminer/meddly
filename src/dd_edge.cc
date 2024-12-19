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



#include "defines.h"
#include "dd_edge.h"
#include "forest.h"
#include "minterms.h"
#include "io.h"

#include "node_marker.h"

#include "oper_binary.h"
#include "ops_builtin.h"
#include "operators.h"

#ifdef ALLOW_DEPRECATED_0_17_3
#include "io_dot.h"
#endif

// #define DEBUG_CLEANUP
// #define DEBUG_FIRST
// #define DEBUG_BINSRCH

namespace MEDDLY {
    class iterator_helper;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                     iterator_helper  class                     *
// *                                                                *
// *                                                                *
// ******************************************************************

class MEDDLY::iterator_helper {
    public:
        iterator_helper(dd_edge::iterator &it) : I(it) { }

        inline unpacked_node* U_from(unsigned k) {
            MEDDLY_DCASSERT(I.U_from);
            return I.U_from[k];
        }
        inline unpacked_node* U_to(unsigned k) {
            MEDDLY_DCASSERT(I.U_to);
            return I.U_to[k];
        }
        inline unsigned& Z_from(unsigned k) {
            MEDDLY_DCASSERT(I.Z_from);
            return I.Z_from[k];
        }
        inline unsigned& Z_to(unsigned k) {
            MEDDLY_DCASSERT(I.Z_to);
            return I.Z_to[k];
        }
        inline edge_value& ev_from(unsigned k) {
            MEDDLY_DCASSERT(I.ev_from);
            return I.ev_from[k];
        }
        inline edge_value& ev_to(unsigned k) {
            MEDDLY_DCASSERT(I.ev_to);
            return I.ev_to[k];
        }

        inline unsigned getNumVars() const {
            MEDDLY_DCASSERT(I.M);
            return I.M->getNumVars();
        }
        inline int& M_from(unsigned k) {
            MEDDLY_DCASSERT(I.M);
            return I.M->from(k);
        }
        inline int& M_to(unsigned k) {
            MEDDLY_DCASSERT(I.M);
            return I.M->to(k);
        }
        inline void M_setTerm(const terminal &t) {
            MEDDLY_DCASSERT(I.M);
            I.M->setTerm(t);
        }
        inline bool isForSets() const {
            MEDDLY_DCASSERT(I.M);
            return I.M->isForSets();
        }

        inline int mask_from(unsigned k) const {
            MEDDLY_DCASSERT(I.mask);
            return I.mask->from(k);
        }
        inline int mask_to(unsigned k) const {
            MEDDLY_DCASSERT(I.mask);
            return I.mask->to(k);
        }

        inline void setAtEnd(bool e) {
            I.atEnd = e;
        }

        inline bool isMultiTerminal() const {
            MEDDLY_DCASSERT(I.F);
            return I.F->isMultiTerminal();
        }
        inline terminal_type getTerminalType() const {
            MEDDLY_DCASSERT(I.F);
            return I.F->getTerminalType();
        }
        inline int getNodeLevel(node_handle p) const {
            MEDDLY_DCASSERT(I.F);
            return I.F->getNodeLevel(p);
        }

        inline const forest* F() {
            MEDDLY_DCASSERT(I.F);
            return I.F;
        }

    private:
        dd_edge::iterator &I;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      iterator_templ class                      *
// *                                                                *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <class EOP>
    class iterator_templ : public iterator_helper {
        public:
            iterator_templ(dd_edge::iterator &it)
                : iterator_helper(it) { }

            void next();

            bool first_unpr(unsigned k, node_handle p);
            bool first_pri(unsigned k, node_handle p);
    };
};

// ******************************************************************
// *                     iterator_templ methods                     *
// ******************************************************************

template <class EOP>
void MEDDLY::iterator_templ<EOP>::next()
{
    if (isForSets()) {
        //
        // Set version
        //
        for (unsigned k=1; k<=getNumVars(); k++) {
            // See if we can advance the nonzero ptr at level k
            const unpacked_node* U_p = U_from(k);
            if (U_p)
            {
                // this is a free variable.
                unsigned& z = Z_from(k);
                for (z++; z<U_p->getSize(); z++) {
                    // Update the minterm
                    M_from(k) = U_p->index(z);
                    if (EOP::hasEdgeValues()) {
                        // accumulate edge value to here
                        ev_from(k) =
                            EOP::applyOp(ev_from(k+1), U_p->edgeval(z));
                    }
                    if (first_unpr(k-1, U_p->down(z))) {
                        return;
                    }
                }
            }
        } // for k
    } else {
        //
        // Relation version
        //
        for (unsigned k=1; k<=getNumVars(); k++) {
            //
            // Advance primed variable x'_k
            //
            const unpacked_node* U_p = U_to(k);
            if (U_p)
            {
                // this is a free variable.
                unsigned& z = Z_to(k);
                for (z++; z<U_p->getSize(); z++) {
                    // Update the minterm
                    M_to(k) = U_p->index(z);
                    if (EOP::hasEdgeValues()) {
                        // accumulate edge value to here
                        ev_to(k) = EOP::applyOp(ev_from(k), U_p->edgeval(z));
                    }
                    if (first_unpr(k-1, U_p->down(z))) {
                        return;
                    }
                }
            }
            //
            // Advance unprimed variable x_k
            //
            const unpacked_node* U_u = U_from(k);
            if (U_u)
            {
                // this is a free variable.
                unsigned& z = Z_from(k);
                for (z++; z<U_u->getSize(); z++) {
                    // Update the minterm
                    M_from(k) = U_u->index(z);
                    if (EOP::hasEdgeValues()) {
                        // accumulate edge value to here
                        ev_from(k) = EOP::applyOp(ev_to(k+1), U_u->edgeval(z));
                    }
                    if (first_pri(k, U_u->down(z))) {
                        return;
                    }
                }
            }
        } // for k
    }
    setAtEnd(true);
}

template <class EOP>
bool MEDDLY::iterator_templ<EOP>::first_unpr(unsigned k, node_handle p)
{
    if (0==p) return false;
    if (0==k) {
        // set the terminal value
        if (isMultiTerminal()) {
            M_setTerm(terminal(getTerminalType(), p));
            return true;
        }
        if (isForSets()) {
            M_setTerm( EOP::buildTerm(ev_from(1), p) );
        } else {
            M_setTerm( EOP::buildTerm(ev_to(1), p) );
        }
        return true;
    }

    //
    // Is this a bound variable?
    //
    const int plvl = getNodeLevel(p);
    unpacked_node *U = U_from(k);
    if (!U) {
        //
        // This value is fixed, determine the down pointer and recurse.
        //
        node_handle pdn = p;

        if (EOP::hasEdgeValues()) {
            //
            // There are edge values to update.
            //
            const edge_value& up = isForSets() ? ev_from(k+1) : ev_to(k+1);

            if (int(k) == plvl) {
                F()->getDownPtr(p, mask_from(k), ev_from(k), pdn);
                EOP::accumulateOp(ev_from(k), up);
            } else {
                ev_from(k) = up;
            }
        } else {
            if (int(k) == plvl) {
                pdn = F()->getDownPtr(p, mask_from(k));
            }
        }

        if (isForSets()) {
            return first_unpr(k-1, pdn);
        } else {
            return first_pri(k, pdn);
        }
    }

    //
    // Free variable.
    // Set up the unpacked node.
    //
    if (k != plvl) {
        if (EOP::hasEdgeValues()) {
            edge_value zero;
            EOP::clear(zero);
            U->initRedundant(F(), k, zero, p, SPARSE_ONLY);
        } else {
            U->initRedundant(F(), k, p, SPARSE_ONLY);
        }
    } else {
        F()->unpackNode(U, p, SPARSE_ONLY);
    }

    //
    // Loop until we find a match
    //
    unsigned& z = Z_from(k);
    if (isForSets()) {
        //
        // set version
        //
        for (z=0; z < U->getSize(); ++z)
        {
            // Update minterm
            M_from(k) = U->index(z);
            if (EOP::hasEdgeValues()) {
                // accumulate edge value to here
                ev_from(k) = EOP::applyOp(U->edgeval(z), ev_from(k+1));
            }
            if (first_unpr(k-1, U->down(z))) return true;
        }
    } else {
        //
        // relation version
        //
        for (z=0; z < U->getSize(); ++z)
        {
            // Update minterm
            M_from(k) = U->index(z);
            if (EOP::hasEdgeValues()) {
                // accumulate edge value to here
                ev_from(k) = EOP::applyOp(U->edgeval(z), ev_to(k+1));
            }
            if (first_pri(k, U->down(z))) return true;
        }
    }

    // Still here? couldn't find one
    return false;
}

template <class EOP>
bool MEDDLY::iterator_templ<EOP>::first_pri(unsigned k, node_handle p)
{
    MEDDLY_DCASSERT(!isForSets());
    MEDDLY_DCASSERT(k);

    if (0==p) return false;

    //
    // Is this a bound variable?
    //
    const int plvl = getNodeLevel(p);
    unpacked_node *U = U_to(k);
    if (!U) {
        //
        // This value is fixed, determine the down pointer and recurse.
        //

        //
        // Deal with identity mask
        //
        int i = mask_to(k);
        if (DONT_CHANGE == i) {
            i = M_from(k);
            M_to(k) = i;
        }

        if (EOP::hasEdgeValues()) {
            //
            // There are edge values to update.
            //
            const edge_value& up = ev_from(k);

            if (k == -plvl) {
                //
                // Node at this level.
                //
                node_handle pdn;
                F()->getDownPtr(p, i, ev_to(k), pdn);
                EOP::accumulateOp(ev_to(k), up);
                return first_unpr(k-1, pdn);
            }
            //
            // Recurse for redundant nodes, or identity
            // if from == to.
            //
            ev_to(k) = up;
            if (F()->isFullyReduced() || (M_from(k) == M_to(k))) {
                return first_unpr(k-1, p);
            }
            return false;

        } else {

            //
            // No edge values
            //
            if (k == -plvl) {
                return first_unpr(k-1, F()->getDownPtr(p, i));
            }
            if (F()->isFullyReduced() || (M_from(k) == M_to(k))) {
                return first_unpr(k-1, p);
            }
            return false;
        }
    }

    //
    // Free variable.
    // Set up the unpacked node.
    //
    if (k != -plvl) {
        if (EOP::hasEdgeValues()) {
            edge_value zero;
            EOP::clear(zero);
            if (F()->isFullyReduced()) {
                U->initRedundant(F(), -int(k), zero, p, SPARSE_ONLY);
            } else {
                U->initIdentity(F(), -int(k), M_from(k), zero, p, SPARSE_ONLY);
            }
        } else {
            if (F()->isFullyReduced()) {
                U->initRedundant(F(), -int(k), p, SPARSE_ONLY);
            } else {
                U->initIdentity(F(), -int(k), M_from(k), p, SPARSE_ONLY);
            }
        }
    } else {
        F()->unpackNode(U, p, SPARSE_ONLY);
    }

    //
    // Loop until we find a match
    //
    unsigned& z = Z_to(k);
    for (z=0; z < U->getSize(); ++z)
    {
        // Update minterm
        M_to(k) = U->index(z);
        if (EOP::hasEdgeValues()) {
            // accumulate edge value to here
            ev_to(k) = EOP::applyOp(ev_from(k), U->edgeval(z));
        }
        if (first_unpr(k-1, U->down(z))) return true;
    }

    // Still here? couldn't find one
    return false;
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                   dd_edge::iterator  methods                   *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::dd_edge::iterator::iterator()
{
    // Build an empty, 'end' iterator
    atEnd = true;
    F = nullptr;

    U_from = nullptr;
    U_to = nullptr;
    Z_from = nullptr;
    Z_to = nullptr;
    ev_from = nullptr;
    ev_to = nullptr;
    M = nullptr;
    mask = nullptr;
}

MEDDLY::dd_edge::iterator::iterator(const dd_edge &E, const minterm* _mask)
{
    F = E.getForest();
    M = new minterm(F);

    //
    // Set up array of unpacked nodes
    //

    const unsigned K = M->getNumVars();
    U_from = new unpacked_node* [1+K];
    for (unsigned i=0; i<=K; i++) {
        U_from[i] = nullptr;
    }
    if (M->isForSets()) {
        U_to = nullptr;
    } else {
        U_to = new unpacked_node* [1+K];
        for (unsigned i=0; i<=K; i++) {
            U_to[i] = nullptr;
        }
    }

    //
    // Allocate array of pointers into the unpacked nodes
    //

    Z_from = new unsigned [1+K];
    if (M->isForSets()) {
        Z_to = nullptr;
    } else {
        Z_to = new unsigned [1+K];
    }

    //
    // Allocate array of edge values, if we're an EV forest
    //
    if (F->isMultiTerminal()) {
        ev_from = nullptr;
        ev_to = nullptr;
    } else {
        if (M->isForSets()) {
            ev_from = new edge_value [2+K];
            ev_to = nullptr;
        } else {
            ev_from = new edge_value [1+K];
            ev_to = new edge_value [2+K];
        }
    }

    //
    // Everything is set up except for mask specific stuff.
    //

    restart(E, _mask);
}

MEDDLY::dd_edge::iterator::iterator(iterator &&I)
{
    printf("In iterator's move constructor :)\n");
    // Take everything from I
    U_from = I.U_from;
    U_to   = I.U_to;

    Z_from = I.Z_from;
    Z_to   = I.Z_to;

    ev_from = I.ev_from;
    ev_to   = I.ev_to;

    root_ev   = I.root_ev;
    root_node = I.root_node;

    F = I.F;
    M = I.M;
    mask = I.mask;
    atEnd = I.atEnd;

    // And clear out I, for anything that would be deleted
    I.U_from = nullptr;
    I.U_to = nullptr;
    I.Z_from = nullptr;
    I.Z_to = nullptr;
    I.ev_from = nullptr;
    I.ev_to = nullptr;
    I.M = nullptr;
}

MEDDLY::dd_edge::iterator::~iterator()
{
    if (U_from) {
        MEDDLY_DCASSERT(M);
        for (unsigned i=M->getNumVars(); i; --i) {
            unpacked_node::Recycle(U_from[i]);
        }
        delete[] U_from;
        U_from = nullptr;
    }
    if (U_to) {
        MEDDLY_DCASSERT(M);
        for (unsigned i=M->getNumVars(); i; --i) {
            unpacked_node::Recycle(U_to[i]);
        }
        delete[] U_to;
        U_to = nullptr;
    }
    delete[] Z_from;
    delete[] Z_to;
    delete[] ev_from;
    delete[] ev_to;
    delete M;
}

void MEDDLY::dd_edge::iterator::restart(const dd_edge &E, const minterm* _mask)
{
    if (F != E.getForest()) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    mask = _mask;
    root_ev = E.getEdgeValue();
    root_node = E.getNode();

    //
    // Seed the edge value array if needed
    //
    const unsigned K = M->getNumVars();
    if (!F->isMultiTerminal()) {
        if (M->isForSets()) {
            ev_from[1+K] = root_ev;
        } else {
            ev_to[1+K] = root_ev;
        }
    }

    //
    // Allocate unpacked nodes, but only for free variables.
    // Re-use old unpacked nodes when we can.
    //

    for (unsigned i=K; i; --i) {
        if (_mask && (DONT_CARE != mask->from(i))) {
            if (U_from[i]) {
                unpacked_node::Recycle(U_from[i]);
                U_from[i] = nullptr;
            }
            M->from(i) = mask->from(i);
        } else {
            if (!U_from[i]) {
                U_from[i] = unpacked_node::New(F);
            }
        }
    }
    if (M->isForRelations()) {
        for (unsigned i=K; i; --i) {
            if (_mask && (DONT_CARE != mask->to(i))) {
                if (U_to[i]) {
                    unpacked_node::Recycle(U_to[i]);
                    U_to[i] = nullptr;
                }
                M->to(i) = mask->to(i);
            } else {
                if (!U_to[i]) {
                    U_to[i] = unpacked_node::New(F);
                }
            }
        }
    }

    //
    // Invoke appropriate 'first'
    //
    switch (F->getEdgeLabeling()) {
        case edge_labeling::MULTI_TERMINAL:
            {
                iterator_templ< EdgeOp_none > IT(*this);
                atEnd = ! IT.first_unpr(K, root_node);
                return;
            }

        case edge_labeling::EVTIMES:
            switch (F->getEdgeType()) {
                case edge_type::FLOAT:
                    {
                        iterator_templ< EdgeOp_times<float> > IT(*this);
                        atEnd = ! IT.first_unpr(K, root_node);
                        return;
                    }
                case edge_type::DOUBLE:
                    {
                        iterator_templ< EdgeOp_times<double> > IT(*this);
                        atEnd = ! IT.first_unpr(K, root_node);
                        return;
                    }
                default:
                    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            }

        case edge_labeling::EVPLUS:
        case edge_labeling::INDEX_SET:
            switch (F->getEdgeType()) {
                case edge_type::INT:
                    {
                        iterator_templ< EdgeOp_plus<int> > IT(*this);
                        atEnd = ! IT.first_unpr(K, root_node);
                        return;
                    }
                case edge_type::LONG:
                    {
                        iterator_templ< EdgeOp_plus<long> > IT(*this);
                        atEnd = ! IT.first_unpr(K, root_node);
                        return;
                    }
                default:
                    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            }

        default:
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

    }
}

void MEDDLY::dd_edge::iterator::next()
{
    MEDDLY_DCASSERT(F);

    switch (F->getEdgeLabeling()) {
        case edge_labeling::MULTI_TERMINAL:
            {
                iterator_templ< EdgeOp_none > IT(*this);
                IT.next();
                return;
            }

        case edge_labeling::EVTIMES:
            switch (F->getEdgeType()) {
                case edge_type::FLOAT:
                    {
                        iterator_templ< EdgeOp_times<float> > IT(*this);
                        IT.next();
                        return;
                    }
                case edge_type::DOUBLE:
                    {
                        iterator_templ< EdgeOp_times<double> > IT(*this);
                        IT.next();
                        return;
                    }
                default:
                    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            }

        case edge_labeling::EVPLUS:
        case edge_labeling::INDEX_SET:
            switch (F->getEdgeType()) {
                case edge_type::INT:
                    {
                        iterator_templ< EdgeOp_plus<int> > IT(*this);
                        IT.next();
                        return;
                    }
                case edge_type::LONG:
                    {
                        iterator_templ< EdgeOp_plus<long> > IT(*this);
                        IT.next();
                        return;
                    }
                default:
                    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            }

        default:
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

    }

    /*
    if (M->isForSets()) {
        for (unsigned k=1; k<=M->getNumVars(); k++) {
            // See if we can advance the nonzero ptr at level k
            const unpacked_node* U_p = U_from[k];
            if (U_p)
            {
                // this is a free variable.
                unsigned& z = Z_from[k];
                for (z++; z<U_p->getSize(); z++) {
                    // Update the minterm
                    M->from(k) = U_p->index(z);
                    if (ev_from) {
                        // accumulate edge value to here
                        ev_from[k] = U_p->edgeval(z);
                        if (F->isEVTimes()) {
                            ev_from[k].multiply(ev_from[k+1]);
                        } else {
                            ev_from[k].add(ev_from[k+1]);
                        }
                    }
                    if (first_unprimed(k-1, U_p->down(z))) {
                        return;
                    }
                }
            }
        } // for k
        atEnd = true;
        return;
    }

    //
    // For relations
    //
    for (unsigned k=1; k<=M->getNumVars(); k++) {
        //
        // Advance primed variable x'_k
        //
        const unpacked_node* U_p = U_to[k];
        if (U_p)
        {
            // this is a free variable.
            unsigned& z = Z_to[k];
            for (z++; z<U_p->getSize(); z++) {
                // Update the minterm
                M->to(k) = U_p->index(z);
                if (ev_to) {
                    // accumulate edge value to here
                    ev_to[k] = U_p->edgeval(z);
                    if (F->isEVTimes()) {
                        ev_to[k].multiply(ev_from[k]);
                    } else {
                        ev_to[k].add(ev_from[k]);
                    }
                }
                if (first_unprimed(k-1, U_p->down(z))) {
                    return;
                }
            }
        }
        //
        // Advance unprimed variable x_k
        //
        const unpacked_node* U_u = U_from[k];
        if (U_u)
        {
            // this is a free variable.
            unsigned& z = Z_from[k];
            for (z++; z<U_u->getSize(); z++) {
                // Update the minterm
                M->from(k) = U_u->index(z);
                if (ev_from) {
                    // accumulate edge value to here
                    ev_from[k] = U_u->edgeval(z);
                    if (F->isEVTimes()) {
                        ev_from[k].multiply(ev_to[k+1]);
                    } else {
                        ev_from[k].add(ev_to[k+1]);
                    }
                }
                if (first_primed(k, U_u->down(z))) {
                    return;
                }
            }
        }
    } // for k
    atEnd = true;
    return;
    */
}

/*

bool MEDDLY::dd_edge::iterator::first_unprimed(unsigned k, node_handle p)
{
    MEDDLY_DCASSERT(M);
    MEDDLY_DCASSERT(F);

#ifdef DEBUG_FIRST
    printf("entering first_unprimed(%u, %ld)\n", k, long(p));
#endif

    if (0==p) return false;
    if (0==k) {
        // set the terminal value
        if (F->isMultiTerminal()) {
            M->setTerm(terminal(F->getTerminalType(), p));
#ifdef DEBUG_FIRST
            printf("    completed minterm:\n    ");
            FILE_output out(stdout);
            M->show(out);
            out.put('\n');
#endif
            return true;
        }
        // tbd - need to accumulate edge values
        return true;
    }

    //
    // Is this a bound variable?
    //
    const int plvl = F->getNodeLevel(p);
    unpacked_node *U = U_from[k];
    if (!U) {
#ifdef DEBUG_FIRST
        printf("    value is fixed at %d\n", mask->from(k));
#endif
        //
        // This value is fixed, determine the down pointer and recurse.
        //
        MEDDLY_DCASSERT(mask);

        node_handle pdn = p;
        if (ev_from) {
            //
            // There are edge values to update.
            //
            const edge_value& above =
                M->isForSets() ? ev_from[k+1] : ev_to[k+1];

            if (int(k) == plvl) {
                F->getDownPtr(p, mask->from(k), ev_from[k], pdn);

                if (F->isEVTimes()) {
                    ev_from[k].multiply(above);
                } else {
                    ev_from[k].add(above);
                }
            } else {
                ev_from[k] = above;
            }
        } else {
            //
            // No edge values to accumulate
            //
            if (int(k) == plvl) {
                pdn = F->getDownPtr(p, mask->from(k));
            }
        }

        if (M->isForSets()) {
            return first_unprimed(k-1, pdn);
        } else {
            return first_primed(k, pdn);
        }
    }

    //
    // Free variable.
    // Set up the unpacked node.
    //
    if (k != plvl) {
        U->initRedundant(F, k, p, SPARSE_ONLY);
    } else {
        F->unpackNode(U, p, SPARSE_ONLY);
    }

    //
    // Loop until we find a match
    //
    unsigned& z = Z_from[k];
    if (M->isForSets()) {
        //
        // set version
        //
        for (z=0; z < U->getSize(); ++z)
        {
            // Update minterm
            M->from(k) = U->index(z);
            if (ev_from) {
                // accumulate edge value to here
                ev_from[k] = U->edgeval(z);
                if (F->isEVTimes()) {
                    ev_from[k].multiply(ev_from[k+1]);
                } else {
                    ev_from[k].add(ev_from[k+1]);
                }
            }
            if (first_unprimed(k-1, U->down(z))) return true;
        }
    } else {
        //
        // relation version
        //
        for (z=0; z < U->getSize(); ++z)
        {
            // Update minterm
            M->from(k) = U->index(z);
            if (ev_from) {
                // accumulate edge value to here
                ev_from[k] = U->edgeval(z);
                if (F->isEVTimes()) {
                    ev_from[k].multiply(ev_to[k+1]);
                } else {
                    ev_from[k].add(ev_to[k+1]);
                }
            }
            if (first_primed(k, U->down(z))) return true;
        }
    }

    // Still here? couldn't find one
    return false;
}

bool MEDDLY::dd_edge::iterator::first_primed(unsigned k, node_handle p)
{
    MEDDLY_DCASSERT(M);
    MEDDLY_DCASSERT(F);
    MEDDLY_DCASSERT(M->isForRelations());
    MEDDLY_DCASSERT(k);

#ifdef DEBUG_FIRST
    printf("entering first_unprimed(%u, %ld)\n", k, long(p));
#endif

    if (0==p) return false;

    //
    // Is this a bound variable?
    //
    const int plvl = F->getNodeLevel(p);
    unpacked_node *U = U_to[k];
    if (!U) {
#ifdef DEBUG_FIRST
        printf("    value is fixed at %d\n", mask->to(k));
#endif
        //
        // This value is fixed, determine the down pointer and recurse.
        //
        MEDDLY_DCASSERT(mask);

        //
        // Deal with identity mask
        //
        int i = mask->to(k);
        if (DONT_CHANGE == i) {
            i = M->from(k);
            M->to(k) = i;
        }

        if (k == -plvl)  {
            //
            // Recurse on p[i]
            //
            return first_unprimed(k-1, F->getDownPtr(p, i));
        }

        if (F->isFullyReduced()) {
            //
            // redundant node at this level
            // all down pointers are equal; recurse on p
            //
            return first_unprimed(k-1, p);
        }

        //
        // Identity node at this level
        //
        if (M->from(k) != M->to(k)) return false;
        return first_unprimed(k-1, p);
    }

    //
    // Free variable.
    // Set up the unpacked node.
    //
    if (k != -plvl) {
        if (F->isFullyReduced()) {
            U->initRedundant(F, -int(k), p, SPARSE_ONLY);
        } else {
            U->initIdentity(F, -int(k), M->from(k), p, SPARSE_ONLY);
        }
    } else {
        F->unpackNode(U, p, SPARSE_ONLY);
    }

    //
    // Loop until we find a match
    //
    unsigned& z = Z_to[k];
    for (z=0; z < U->getSize(); ++z)
    {
        // Update minterm
        M->to(k) = U->index(z);
        if (first_unprimed(k-1, U->down(z))) return true;
    }

    // Still here? couldn't find one
    return false;
}

*/

bool MEDDLY::dd_edge::iterator::equals(const iterator &I) const
{
    if (F != I.F) return false;
    if (root_node != I.root_node) return false;
    if (mask != I.mask) return false;

    MEDDLY_DCASSERT(M);
    if (M->isForSets()) {
        for (unsigned i = M->getNumVars(); i; --i) {
            if (Z_from[i] != I.Z_from[i]) return false;
        }
    } else {
        for (unsigned i = M->getNumVars(); i; --i) {
            if (Z_from[i] != I.Z_from[i]) return false;
            if (Z_to[i] != I.Z_to[i]) return false;
        }
    }

    return true;
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        dd_edge  methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::dd_edge::dd_edge(forest* p)
{
#ifdef DEBUG_CLEANUP
    std::cout << "Creating dd_edge" << std::endl;
#endif
    node = 0;
    prev = nullptr;
    next = nullptr;
    if (p)  p->registerEdge(*this);
    else    parentFID = 0;
}

// Copy Constructor.
MEDDLY::dd_edge::dd_edge(const dd_edge& e)
{
#ifdef DEBUG_CLEANUP
    std::cout << "Creating dd_edge via copy" << std::endl;
#endif
    prev = nullptr;
    next = nullptr;
    init(e);
}


// Assignment operator.
MEDDLY::dd_edge& MEDDLY::dd_edge::operator=(const dd_edge& e)
{
    if (equals(e)) return *this;
    detach();
    init(e);
    return *this;
}

// Destructor.  Will notify parent as appropriate.
MEDDLY::dd_edge::~dd_edge()
{
#ifdef DEBUG_CLEANUP
    std::cout << "Deleting dd_edge" << std::endl;
#endif
    detach();
}

void MEDDLY::dd_edge::attach(forest* p)
{
    forest* efp = forest::getForestWithID(parentFID);
    if (p == efp) return;
    if (efp) {
#ifdef DEBUG_CLEANUP
        std::cout << "unlinking " << node << " in dd_edge::attach" << std::endl;
#endif
        efp->unlinkNode(node);
        node = 0;
        efp->unregisterEdge(*this);
    }
    if (p) {
        p->registerEdge(*this);
    } else {
        parentFID = 0;
    }
}

MEDDLY::forest* MEDDLY::dd_edge::getForest() const
{
    return forest::getForestWithID(parentFID);
}

unsigned long MEDDLY::dd_edge::getNodeCount() const
{
    forest* fp = forest::getForestWithID(parentFID);

    if (fp) {
        node_marker M(fp);
        M.mark(node);
        return M.countMarked();
    }
    return 0;
}

unsigned long MEDDLY::dd_edge::getEdgeCount(bool countZeroes) const
{
    forest* fp = forest::getForestWithID(parentFID);

    if (fp) {
        node_marker M(fp);
        M.mark(node);
        return countZeroes ? M.countEdges() : M.countNonzeroEdges();
    }
    return 0;
}

int MEDDLY::dd_edge::getLevel() const
{
    if (0==node) return 0;

    forest* efp = forest::getForestWithID(parentFID);
    MEDDLY_DCASSERT(efp);
    return efp->getNodeLevel(node);
}


void MEDDLY::dd_edge::show(output &s) const
{
    if (!parentFID) {
        s.put("<null edge>");
        return;
    }

    forest* efp = forest::getForestWithID(parentFID);
    if (!efp) {
        s.put("<orphaned edge>");
        return;
    }

    efp->showEdge(s, edgeval, node);
    s.put(" in ");
    if (efp->isMultiTerminal()) {
        s.put("MT");
    }
    if (efp->isEVPlus()) {
        s.put("EV+");
    }
    if (efp->isEVTimes()) {
        s.put("EV*");
    }
    if (efp->isForRelations()) {
        s.put("MxD");
    } else {
        s.put("MDD");
    }
    s.put(" forest ");
    s.put(parentFID);
}

void MEDDLY::dd_edge::showGraph(output &s) const
{
    forest* efp = forest::getForestWithID(parentFID);
    if (!efp) {
        s.put("null graph\n");
        return;
    }
    if (efp->isMultiTerminal()) {
        s.put("MT");
    }
    if (efp->isEVPlus()) {
        s.put("EV+");
    }
    if (efp->isEVTimes()) {
        s.put("EV*");
    }
    if (efp->isForRelations()) {
        s.put("MxD");
    } else {
        s.put("MDD");
    }
    if (efp->isIdentityReduced()) {
        s.put(" (identity reduced)");
    }
    if (efp->isFullyReduced()) {
        s.put(" (fully reduced)");
    }
    if (efp->isQuasiReduced()) {
        s.put(" (quasi reduced)");
    }
    s.put(" rooted at edge ");
    efp->showEdge(s, edgeval, node);
    s.put('\n');

    node_marker M(efp);
    M.mark(node);
    M.showByLevels(s);
}

#ifdef ALLOW_DEPRECATED_0_17_3

void MEDDLY::dd_edge::writePicture(const char* filename,
        const char* extension) const
{
    dot_maker DM(forest::getForestWithID(parentFID), filename);

    DM.addRootEdge(*this);
    DM.doneGraph();
    DM.runDot(extension);
}

#endif


void MEDDLY::dd_edge::write(output &s, const std::vector <unsigned> &map) const
{
    //
    // Edge info
    //
    edgeval.write(s);
    s.put(' ');

    //
    // Node info
    //
    if (node <= 0) {
        //
        // terminal
        //
        forest* efp = forest::getForestWithID(parentFID);
        terminal t;
        t.setFromHandle(efp->getTerminalType(), node);
        t.write(s);
    } else {
        //
        // non-terminal
        //
        s.put("n ");
        s.put(map[unsigned(node)]);
    }
}


void MEDDLY::dd_edge::read(input &s, const std::vector <node_handle> &map)
{
#ifdef DEBUG_READ_DD
    std::cerr << "    in dd_edge::read\n";
#endif
    //
    // Edge info
    //
    edgeval.read(s);
    s.stripWS();

    //
    // Node info
    //
    int c = s.get_char();
    if ('n' == c) {
        //
        // non-terminal
        //
        long d = s.get_integer();
        if (d<0) throw error(error::INVALID_FILE, __FILE__, __LINE__);
        set_and_link(map[unsigned(d)]);
    } else {
        //
        // terminal
        //
        s.unget(c);
        terminal t;
        t.read(s);
        set(t.getHandle());
    }
#ifdef DEBUG_READ_DD
    std::cerr << "    done dd_edge::read\n";
#endif
}


// **********************************************************************
//
// Function evaluation
// Eventually replacing evaluate() in the forest class
//
// **********************************************************************

namespace MEDDLY {

    class evaluator_helper_mt {
            const forest* F;
            const minterm& m;
        public:
            evaluator_helper_mt(const forest* _F, const minterm &_m)
                : F(_F), m(_m)
            {
                MEDDLY_DCASSERT(F);
            }

            /*
            inline void evaluate(node_handle &p)
            {
                if (m.isForRelations()) {
                    if (F->isIdentityReduced()) {
                        ident_rel_eval(p);
                    } else {
                        fully_rel_eval(p);
                    }
                } else {
                    set_eval(p);
                }
            }
            */

            inline void set_eval(node_handle &p)
            {
                MEDDLY_DCASSERT( !m.isForRelations() );
                while (!F->isTerminalNode(p)) {
    	            int level = F->getNodeLevel(p);
                    p = F->getDownPtr(p, m.from(level));
                }
            }

            inline void fully_rel_eval(node_handle &p)
            {
                MEDDLY_DCASSERT( m.isForRelations() );
                while (!F->isTerminalNode(p)) {
    	            int level = F->getNodeLevel(p);
                    if (level > 0) {
                        p = F->getDownPtr(p, m.from(level));
                    } else {
                        p = F->getDownPtr(p, m.to(-level));
                    }
                }
            }

            inline void ident_rel_eval(node_handle &p)
            {
                if (0==p) return;
                MEDDLY_DCASSERT( m.isForRelations() );
                int plvl = F->getNodeLevel(p);
                int L = F->getNumVariables();
                while (L) {
                    //
                    // Unprimed
                    //
                    MEDDLY_DCASSERT(L>0);
                    if (plvl == L) {
                        p = F->getDownPtr(p, m.from(L));
                        if (0==p) return;
                        plvl = F->getNodeLevel(p);
                    }
                    L = MXD_levels::downLevel(L);
                    //
                    // Primed
                    //
                    MEDDLY_DCASSERT(L<0);
                    if (plvl == L) {
                        p = F->getDownPtr(p, m.to(-L));
                        if (0==p) return;
                        plvl = F->getNodeLevel(p);
                    } else {
                        if (m.to(-L) != m.from(-L)) {
                            p = 0;
                            return;
                        }
                    }
                    L = MXD_levels::downLevel(L);
                }
            }
    };

    template <class EOP>
    class evaluator_helper {
            const forest* F;
            const minterm& m;
        public:
            evaluator_helper(const forest* _F, const minterm &_m)
                : F(_F), m(_m)
            {
                MEDDLY_DCASSERT(F);
            }

            inline void evaluate(edge_value &ev, node_handle &p)
            {
                if (m.isForRelations()) {
                    if (F->isIdentityReduced()) {
                        ident_rel_eval(ev, p);
                    } else {
                        fully_rel_eval(ev, p);
                    }
                } else {
                    set_eval(ev, p);
                }
            }

            inline void set_eval(edge_value &ev, node_handle &p)
            {
                MEDDLY_DCASSERT( !m.isForRelations() );
                while (!F->isTerminalNode(p)) {
                    edge_value pv;
    	            int level = F->getNodeLevel(p);
                    MEDDLY_DCASSERT(level>0);
                    F->getDownPtr(p, m.from(level), pv, p);
                    EOP::accumulateOp(ev, pv);
                }
            }

            inline void fully_rel_eval(edge_value &ev, node_handle &p)
            {
                MEDDLY_DCASSERT( m.isForRelations() );
                while (!F->isTerminalNode(p)) {
                    edge_value pv;
                    int level = F->getNodeLevel(p);
                    if (level > 0) {
                        F->getDownPtr(p, m.from(level), pv, p);
                    } else {
                        F->getDownPtr(p, m.to(-level), pv, p);
                    }
                    EOP::accumulateOp(ev, pv);
                }
            }

            inline void ident_rel_eval(edge_value &ev, node_handle &p)
            {
                if (0==p) return;
                MEDDLY_DCASSERT( m.isForRelations() );
                int plvl = F->getNodeLevel(p);
                int L = F->getNumVariables();
                while (L) {
                    edge_value pv;
                    //
                    // Unprimed
                    //
                    MEDDLY_DCASSERT(L>0);
                    if (plvl == L) {
                        F->getDownPtr(p, m.from(L), pv, p);
                        EOP::accumulateOp(ev, pv);
                        if (0==p) return;
                        plvl = F->getNodeLevel(p);
                    }
                    L = MXD_levels::downLevel(L);
                    //
                    // Primed
                    //
                    MEDDLY_DCASSERT(L<0);
                    if (plvl == L) {
                        F->getDownPtr(p, m.to(-L), pv, p);
                        EOP::accumulateOp(ev, pv);
                        if (0==p) return;
                        plvl = F->getNodeLevel(p);
                    } else {
                        if (m.to(-L) != m.from(-L)) {
                            EOP::clear(ev);
                            p = 0;
                            return;
                        }
                    }
                    L = MXD_levels::downLevel(L);
                }
            }

    };
};


// **********************************************************************

void MEDDLY::dd_edge::evaluate(const minterm& m,
            edge_value &ev, node_handle &en) const
{
    forest* fp = forest::getForestWithID(parentFID);
    if (!fp) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if ( fp->isForRelations() != m.isForRelations() ) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if ( fp->getDomain() != m.getDomain() ) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    if ( fp->isMultiTerminal() )
    {
        evaluator_helper_mt EH(fp, m);
        if (m.isForRelations()) {
            if (fp->isIdentityReduced()) {
                EH.ident_rel_eval(en);
            } else {
                EH.fully_rel_eval(en);
            }
        } else {
            EH.set_eval(en);
        }
        return;
    }

    if ( fp->isEVPlus() || fp->isIndexSet() )
    {
        if (edge_type::INT == fp->getEdgeType()) {
            evaluator_helper< EdgeOp_plus <int> > EH(fp, m);
            EH.evaluate(ev, en);
            return;
        }
        if (edge_type::LONG == fp->getEdgeType()) {
            evaluator_helper< EdgeOp_plus <long> > EH(fp, m);
            EH.evaluate(ev, en);
            return;
        }
    }

    if ( fp->isEVTimes() )
    {
        if (edge_type::FLOAT == fp->getEdgeType()) {
            evaluator_helper< EdgeOp_times <float> > EH(fp, m);
            EH.evaluate(ev, en);
            return;
        }
        if (edge_type::DOUBLE == fp->getEdgeType()) {
            evaluator_helper< EdgeOp_times <double> > EH(fp, m);
            EH.evaluate(ev, en);
            return;
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

// **********************************************************************

bool MEDDLY::dd_edge::getElemInt(long index, minterm &m) const
{
    //
    // Sanity checks
    //
    forest* fp = forest::getForestWithID(parentFID);
    if (!fp) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if (!fp->isIndexSet()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    if ( fp->getDomain() != m.getDomain() ) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (fp->isForRelations()) {
        // I don't think we have index set relations?
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    if (m.isForRelations() ) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    MEDDLY_DCASSERT(fp->getEdgeType() == edge_type::LONG);

    if (index < 0) return false;

    node_handle p = node;
    unpacked_node* U = unpacked_node::New(fp);
    for (unsigned k = fp->getNumVariables(); k; --k) {
        //
        // I don't think index sets can skip levels at all
        //
        MEDDLY_DCASSERT(k == fp->getNodeLevel(p));
        fp->unpackNode(U, p, SPARSE_ONLY);

        //
        // Backward Linear search: largest i such that
        // edge value i <= index
        //
        unsigned zmax = U->getSize()-1;
        while (zmax) {
            if (U->edgeval(zmax).getInt() <= index) break;
            --zmax;
        }

        //
        // go down
        //
        m.from(k) = U->index(zmax);
        index -= U->edgeval(zmax).getInt();
        p = U->down(zmax);
    } // for k
    unpacked_node::Recycle(U);
    if (index > 0) return false;
    return true;
}

bool MEDDLY::dd_edge::getElemLong(long index, minterm &m) const
{
    //
    // Sanity checks
    //
    forest* fp = forest::getForestWithID(parentFID);
    if (!fp) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    if (!fp->isIndexSet()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    if ( fp->getDomain() != m.getDomain() ) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    if (fp->isForRelations()) {
        // I don't think we have index set relations?
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    if (m.isForRelations() ) {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }
    MEDDLY_DCASSERT(fp->getEdgeType() == edge_type::LONG);

    if (index < 0) return false;

    node_handle p = node;
    unpacked_node* U = unpacked_node::New(fp);
    for (unsigned k = fp->getNumVariables(); k; --k) {
        //
        // I don't think index sets can skip levels at all
        //
        MEDDLY_DCASSERT(k == fp->getNodeLevel(p));
        fp->unpackNode(U, p, SPARSE_ONLY);

        //
        // Backward Linear search: largest i such that
        // edge value i <= index
        //
        unsigned zmax = U->getSize()-1;
        while (zmax) {
            if (U->edgeval(zmax).getLong() <= index) break;
            --zmax;
        }

        //
        // go down
        //
        m.from(k) = U->index(zmax);
        index -= U->edgeval(zmax).getLong();
        p = U->down(zmax);
    } // for k
    unpacked_node::Recycle(U);
    if (index > 0) return false;
    return true;
}

//
// Private
//

void MEDDLY::dd_edge::init(const dd_edge &e)
{
    forest* efp = forest::getForestWithID(e.parentFID);
    if (efp) {
        efp->registerEdge(*this);
        node = efp->linkNode(e.node);
        edgeval = e.edgeval;
        label = e.label;
    } else {
        parentFID = 0;
        node = 0;
        edgeval.set();
    }
}

void MEDDLY::dd_edge::set(node_handle n)
{
    if (node != n) {
        forest* efp = forest::getForestWithID(parentFID);
        if (efp) {
#ifdef DEBUG_CLEANUP
            std::cout << "unlinking " << node << " in dd_edge::set" << std::endl;
#endif
            efp->unlinkNode(node);
            node = n;
        } else {
            node = 0;
        }
    }
}

void MEDDLY::dd_edge::set_and_link(node_handle n)
{
    if (node != n) {
        forest* efp = forest::getForestWithID(parentFID);
        if (efp) {
#ifdef DEBUG_CLEANUP
            std::cout << "unlinking " << node << " in dd_edge::set" << std::endl;
            std::cout << "linking " << n << " in dd_edge::set" << std::endl;
#endif
            efp->unlinkNode(node);
            node = efp->linkNode(n);
        } else {
            node = 0;
        }
    }
}


