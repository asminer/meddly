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

// ******************************************************************
// *                                                                *
// *                                                                *
// *                   dd_edge::iterator  methods                   *
// *                                                                *
// *                                                                *
// ******************************************************************

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
    root_ev = E.getEdgeValue();
    root_node = E.getNode();

    F = E.getForest();
    M = new minterm(F);
    mask = _mask;

    //
    // Allocate unpacked nodes, but only for free variables.
    //

    U_from = new unpacked_node* [1+M->getNumVars()];
    U_from[0] = nullptr;
    for (unsigned i=M->getNumVars(); i; --i) {
        if (_mask && (DONT_CARE != mask->from(i))) {
            U_from[i] = nullptr;
            M->from(i) = mask->from(i);
        } else {
            U_from[i] = unpacked_node::New(F);
        }
    }
    if (M->isForSets()) {
        U_to = nullptr;
    } else {
        U_to = new unpacked_node* [1+M->getNumVars()];
        U_to[0] = nullptr;
        for (unsigned i=M->getNumVars(); i; --i) {
            if (_mask && (DONT_CARE != mask->to(i))) {
                U_to[i] = nullptr;
                M->to(i) = mask->to(i);
            } else {
                U_to[i] = unpacked_node::New(F);
            }
        }
    }

    //
    // Allocate array of pointers into the unpacked nodes
    //

    Z_from = new unsigned [1+M->getNumVars()];
    if (M->isForSets()) {
        Z_to = nullptr;
    } else {
        Z_to = new unsigned [1+M->getNumVars()];
    }

    //
    // Allocate array of edge values, if we're an EV forest
    //
    if (F->isMultiTerminal()) {
        ev_from = nullptr;
        ev_to = nullptr;
    } else {
        ev_from = new edge_value [1+M->getNumVars()];
        if (M->isForSets()) {
            ev_to = nullptr;
        } else {
            ev_to = new edge_value [1+M->getNumVars()];
        }
    }

    atEnd = ! first_unprimed(M->getNumVars(), root_node);
    // atEnd = ! first(int(M->getNumVars()), root_node);
}

void MEDDLY::dd_edge::iterator::next()
{
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
                if (first_primed(k, U_u->down(z))) {
                    return;
                }
            }
        }
    } // for k
    atEnd = true;
    return;
}

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

        const node_handle pdn = (int(k) == plvl)
            ? F->getDownPtr(p, mask->from(k))
            : p
        ;

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
        for (z=0; z < U->getSize(); ++z)
        {
            // Update minterm
            M->from(k) = U->index(z);
            if (first_unprimed(k-1, U->down(z))) return true;
        }
    } else {
        for (z=0; z < U->getSize(); ++z)
        {
            // Update minterm
            M->from(k) = U->index(z);
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


