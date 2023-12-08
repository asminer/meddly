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
#include "io.h"

#include "node_marker.h"

#include "opname.h"
#include "oper_binary.h"
#include "ops_builtin.h"
#include "operators.h"

#ifdef ALLOW_DEPRECATED_0_17_3
#include "io_dot.h"
#endif

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
    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
    );
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

    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
    );
    MEDDLY_DCASSERT(efp);
    return efp->getNodeLevel(node);
}


void MEDDLY::dd_edge::show(output &s) const
{
    if (!parentFID) {
        s.put("<null edge>");
        return;
    }

    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
    );
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
    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
    );
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
        expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
        );
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


//
// Private
//

void MEDDLY::dd_edge::init(const dd_edge &e)
{
    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(e.parentFID)
    );
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
        expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
        );
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
        expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
        );
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


