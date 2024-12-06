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

            inline node_handle set_eval(node_handle p)
            {
                MEDDLY_DCASSERT( !m.isForRelations() );
                while (!F->isTerminalNode(p)) {
    	            int level = F->getNodeLevel(p);
                    p = F->getDownPtr(p, m.from(level));
                }
                return p;
            }

            inline node_handle fully_rel_eval(node_handle p)
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
                return p;
            }

            inline node_handle ident_rel_eval(node_handle p)
            {
                if (0==p) return 0;
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
                        if (0==p) return 0;
                        plvl = F->getNodeLevel(p);
                    }
                    L = MXD_levels::downLevel(L);
                    //
                    // Primed
                    //
                    MEDDLY_DCASSERT(L<0);
                    if (plvl == L) {
                        p = F->getDownPtr(p, m.to(-L));
                        if (0==p) return 0;
                        plvl = F->getNodeLevel(p);
                    } else {
                        if (m.to(-L) != m.from(-L)) {
                            return 0;
                        }
                    }
                    L = MXD_levels::downLevel(L);
                }
                return p;
            }
    };

    /*

    template <class EOP>
    class evaluator_helper_ev {
            const forest* F;
            const minterm& m;
        public:
            evaluator_helper_ev(const forest* _F, const minterm &_m)
                : F(_F), m(_m)
            {
                MEDDLY_DCASSERT(F);
            }

            inline void set_eval(edge_value &ev, node_handle &p)
            {
                MEDDLY_DCASSERT( !m.isForRelations() );
                while (!F->isTerminalNode(p)) {
                    edge_value pv;
    	            int level = F->getNodeLevel(p);
                    MEDDLY_DCASSERT(level>0);
                    F->getDownPtr(p, m.from(level), pv, p);
                    ev = EOP::accumulate(ev, pv);
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
                    ev = EOP::accumulate(ev, pv);
                }
            }

            inline void ident_rel_eval(edge_value &ev, node_handle &p)
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
                    ev = EOP::accumulate(ev, pv);
                }
            }

    };
    */
};


// **********************************************************************

void MEDDLY::dd_edge::evaluate(const minterm& m,
            edge_value &ev, node_handle &en) const
{
    forest* fp = forest::getForestWithID(parentFID);
    if (!fp) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
    /*
    if (! fp->isRangeType(range_type::BOOLEAN) ) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
    */
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
                en = EH.ident_rel_eval(node);
            } else {
                en = EH.fully_rel_eval(node);
            }
        } else {
            en = EH.set_eval(node);
        }
        ev.set();
        return;
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


