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
#include "encoders.h"
#include "io.h"

#include "opname.h"
#include "oper_binary.h"
#include "ops_builtin.h"
#include "operators.h"

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
    fprintf(stderr, "Creating dd_edge %p\n", this);
#endif
    label = nullptr;
    node = 0;
    edge_int = 0;
    if (p)  p->registerEdge(*this);
    else    parentFID = 0;
}

// Copy Constructor.
MEDDLY::dd_edge::dd_edge(const dd_edge& e)
{
#ifdef DEBUG_CLEANUP
    fprintf(stderr, "Creating dd_edge %p\n", this);
#endif
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
    fprintf(stderr, "Deleting dd_edge %p\n", this);
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
        efp->unlinkNode(node);
        node = 0;
        edge_int = 0;
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

void MEDDLY::dd_edge::setLabel(const char* L)
{
    if (label) free(label);
    label = L ? strdup(L) : nullptr;
}

unsigned long MEDDLY::dd_edge::getNodeCount() const
{
    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
    );
    return efp ? efp->getNodeCount(node) : 0;
}

unsigned long MEDDLY::dd_edge::getEdgeCount(bool countZeroes) const
{
    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
    );
    return efp ? efp->getEdgeCount(node, countZeroes) : 0;
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

    s.put('<');
    if (!efp->isMultiTerminal()) {
        efp->showEdgeValue(s, *this);
        s.put(", ");
    }
    if (efp->isTerminalNode(node)) {
        efp->showTerminal(s, node);
    } else {
        s.put('#');
        s.put(long(node));
    }
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
    s.put('>');
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
    s.put(" rooted at edge <");
    if (!efp->isMultiTerminal()) {
        efp->showEdgeValue(s, *this);
        s.put(", ");
    }
    if (efp->isTerminalNode(node)) {
        efp->showTerminal(s, node);
    } else {
        s.put('#');
        s.put(long(node));
    }
    s.put(">\n");
    efp->showNodeGraph(s, &node, 1);
}

void MEDDLY::dd_edge::writePicture(const char* filename,
        const char* extension) const
{
    expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
    );
    if (efp) efp->writeNodeGraphPicture(filename, extension, &node, &label, 1);
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
        edge_int = e.edge_int;
    } else {
        parentFID = 0;
        node = 0;
        edge_int = 0;
    }
}

void MEDDLY::dd_edge::set(node_handle n)
{
    if (node != n) {
        expert_forest* efp = static_cast <expert_forest*> (
                forest::getForestWithID(parentFID)
        );
        if (efp) {
            efp->unlinkNode(node);
            node = efp->linkNode(n);
        } else {
            node = 0;
        }
    }
}

