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

// #define DEBUG_ITER_BEGIN

#ifdef NEW_DD_EDGES

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
    parentFID = p ? p->FID() : 0;
    node = 0;
    edge_int = 0;
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
    if (efp) {
        node = efp->removeRoot(node);
        edge_int = 0;
    }
    parentFID = p ? p->FID() : 0;
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
        parentFID = e.parentFID;
        node = efp->addRoot(e.node);
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
            efp->removeRoot(node);
            node = efp->addRoot(n);
        } else {
            node = 0;
        }
    }
}

#else

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        dd_edge  methods (old)                  *
// *                                                                *
// *                                                                *
// ******************************************************************

// Helper functions
inline void linkNode(MEDDLY::forest* p, int node)
{
  MEDDLY_DCASSERT(p);
  MEDDLY_DCASSERT(smart_cast<MEDDLY::expert_forest*>(p));
  smart_cast<MEDDLY::expert_forest*>(p)->linkNode(node);
}

inline void unlinkNode(MEDDLY::forest* p, int node)
{
  if (p) {
    MEDDLY_DCASSERT(smart_cast<MEDDLY::expert_forest*>(p));
    smart_cast<MEDDLY::expert_forest*>(p)->unlinkNode(node);
  }
}


MEDDLY::dd_edge::dd_edge()
: parent(0), index(0),
  node(0), raw_value(0),
  opPlus(0), opStar(0), opMinus(0), opDivide(0)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %p\n", this);
#endif
  label = 0;
}

// Constructor.
MEDDLY::dd_edge::dd_edge(forest* p)
: parent(p), index(0),
  node(0), raw_value(0),
  opPlus(0), opStar(0), opMinus(0), opDivide(0)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %p\n", this);
#endif
  MEDDLY_DCASSERT(p != NULL);
  parent->registerEdge(*this);
  MEDDLY_DCASSERT(index);
  label = 0;
}


// Copy Constructor.
MEDDLY::dd_edge::dd_edge(const dd_edge& e)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %p\n", this);
#endif
  init(e);
  MEDDLY_DCASSERT(index);
}


// Assignment operator.
MEDDLY::dd_edge& MEDDLY::dd_edge::operator=(const dd_edge& e)
{
  if (&e != this) {
    destroy();
    init(e);
  }
  MEDDLY_DCASSERT(index);
  return *this;
}

// Destructor.  Will notify parent as appropriate.
MEDDLY::dd_edge::~dd_edge()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Deleting dd_edge %p\n", this);
#endif
  destroy();
}

//
// Initialization helper
//
void MEDDLY::dd_edge::init(const dd_edge &e)
{
  parent = e.parent;
  node = e.node;
  raw_value = e.raw_value;

  linkNode(parent, node);

  opPlus = e.opPlus;
  opStar = e.opStar;
  opMinus = e.opMinus;
  opDivide = e.opDivide;

  if (parent) parent->registerEdge(*this);
  MEDDLY_DCASSERT(index);

  label = e.label ? strdup(e.label) : 0;
}

//
// destruction helper
//
void MEDDLY::dd_edge::destroy()
{
  if (index) {
    // still registered; unregister before discarding
    int old = node;
    node = 0;
    unlinkNode(parent, old);
    if (parent) parent->unregisterEdge(*this);
  }
  parent = 0;
  index = 0;
  free(label);
  label = 0;
}

void MEDDLY::dd_edge::setForest(forest* f)
{
  destroy();
  parent = f;
  if (parent) {
    parent->registerEdge(*this);
    MEDDLY_DCASSERT(index);
  }
}

//void MEDDLY::dd_edge::getEdgeValue(int& ev) const
//{
//  MEDDLY_DCASSERT(parent);
//  MEDDLY_DCASSERT(edge_labeling::MULTI_TERMINAL != parent->getEdgeLabeling());
//  MEDDLY_DCASSERT(range_type::INTEGER == parent->getRangeType());
//  expert_forest::EVencoder<int>::readValue(&raw_value, ev);
//}

void MEDDLY::dd_edge::getEdgeValue(long& ev) const
{
  MEDDLY_DCASSERT(parent);
  MEDDLY_DCASSERT(edge_labeling::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(range_type::INTEGER == parent->getRangeType());
  EVencoder<long>::readValue(&raw_value, ev);
}

void MEDDLY::dd_edge::getEdgeValue(float& ev) const
{
  MEDDLY_DCASSERT(parent);
  MEDDLY_DCASSERT(edge_labeling::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(range_type::REAL == parent->getRangeType());
  EVencoder<float>::readValue(&raw_value, ev);
}

void MEDDLY::dd_edge::set(node_handle n)
{
  MEDDLY_DCASSERT(parent);
  node_handle old = node;
  node = n;
  unlinkNode(parent, old);
}

void MEDDLY::dd_edge::set(node_handle n, int v)
{
  set(n, (long)v);
}

void MEDDLY::dd_edge::set(node_handle n, long v)
{
  set(n);
  setEdgeValue(v);
}

void MEDDLY::dd_edge::set(node_handle n, float v)
{
  set(n);
  setEdgeValue(v);
}

void MEDDLY::dd_edge::setEdgeValue(int value)
{
  setEdgeValue((long)value);
}

void MEDDLY::dd_edge::setEdgeValue(long value)
{
  MEDDLY_DCASSERT(parent);
  MEDDLY_DCASSERT(edge_labeling::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(range_type::INTEGER == parent->getRangeType());
  EVencoder<long>::writeValue(&raw_value, value);
}

void MEDDLY::dd_edge::setEdgeValue(float value)
{
  MEDDLY_DCASSERT(parent);
  MEDDLY_DCASSERT(edge_labeling::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(range_type::REAL == parent->getRangeType());
  EVencoder<float>::writeValue(&raw_value, value);
}

void MEDDLY::dd_edge::setLabel(const char* L)
{
  if (label) free(label);
  label = strdup(L);
}

int MEDDLY::dd_edge::getLevel() const
{
    if (0==node) return 0;
    MEDDLY_DCASSERT(parent);
    const expert_forest* ef = dynamic_cast <expert_forest*>(parent);
    MEDDLY_DCASSERT(ef);
    return ef->getNodeLevel(node);
}


/*
double MEDDLY::dd_edge::getCardinality() const
{
  double c;
  apply(CARDINALITY, *this, c);
  return c;
}
*/

unsigned long MEDDLY::dd_edge::getNodeCount() const
{
  return smart_cast<expert_forest*>(parent)->getNodeCount(node);
}

unsigned long MEDDLY::dd_edge::getEdgeCount(bool countZeroes) const
{
  return smart_cast<expert_forest*>(parent)->getEdgeCount(node, countZeroes);
}

//
// Operator +=
/*
MEDDLY::dd_edge& MEDDLY::dd_edge::operator+=(const dd_edge& e)
{
    if (!opPlus) {
        binary_opname* theop;
        theop = (parent->getRangeType() == range_type::BOOLEAN)
                    ? UNION()
                    : PLUS();
        MEDDLY_DCASSERT(theop);
        opPlus = theop->getOperation(*this, e, *this);
        MEDDLY_DCASSERT(opPlus);
    }
    opPlus->computeTemp(*this, e, *this);
    // apply will call set() which in turn will set updateNeeded to true
    return *this;
}
*/

// Operator *=
/*
MEDDLY::dd_edge& MEDDLY::dd_edge::operator*=(const dd_edge& e)
{
    if (!opStar) {
        binary_opname* theop;
        theop = (parent->getRangeType() == range_type::BOOLEAN)
                    ? INTERSECTION()
                    : MULTIPLY();
        MEDDLY_DCASSERT(theop);
        opStar = theop->getOperation(*this, e, *this);
        MEDDLY_DCASSERT(opStar);
    }
    opStar->computeTemp(*this, e, *this);
    // apply will call set() which in turn will set updateNeeded to true
    return *this;
}
*/


// Operator -=
/*
MEDDLY::dd_edge& MEDDLY::dd_edge::operator-=(const dd_edge& e)
{
    if (!opMinus) {
        binary_opname* theop;
        theop = (parent->getRangeType() == range_type::BOOLEAN)
                    ? DIFFERENCE()
                    : MINUS();
        MEDDLY_DCASSERT(theop);
        opMinus = theop->getOperation(*this, e, *this);
        MEDDLY_DCASSERT(opMinus);
    }
    opMinus->computeTemp(*this, e, *this);
    // apply will call set() which in turn will set updateNeeded to true
    return *this;
}
*/


// Operator /=
/*
MEDDLY::dd_edge& MEDDLY::dd_edge::operator/=(const dd_edge& e)
{
    if (!opDivide) {
        binary_opname* theop = DIVIDE();
        MEDDLY_DCASSERT(theop);
        opDivide = theop->getOperation(*this, e, *this);
        MEDDLY_DCASSERT(opDivide);
    }
    opDivide->computeTemp(*this, e, *this);
    // apply will call set() which in turn will set updateNeeded to true
    return *this;
}
*/

// Display the edge information.
void MEDDLY::dd_edge::show(output &strm, int verbosity) const
{
  expert_forest* eParent = smart_cast<expert_forest*>(parent);

  strm.put("(Forest Addr: ");
  strm.put_hex((unsigned long) parent);
  strm.put(", ");

  strm.put("transparent: ");
  eParent->showTerminal(strm, eParent->getTransparentNode());
  strm.put(", ");

  if (eParent->isTerminalNode(node)) {
    strm.put("node: ");
    eParent->showTerminal(strm, node);
    strm.put("*, ");
  }
  else {
    strm.put("node: ");
    strm.put(long(node));
    strm.put(", ");
  }
  if (!eParent->isMultiTerminal()) {
    strm.put("value: ");
    if (eParent->getRangeType() == range_type::REAL) {
      float ev;
      getEdgeValue(ev);
      strm.put(ev);
    } else {
      long iv = Inf<long>();
      getEdgeValue(iv);
      strm.put(iv);
    }
    strm.put(", ");
  }
  strm.put("level: ");
  strm.put(getLevel());
  if (0 != getLevel()) {
    strm << ", extensible: " << eParent->isExtensible(node);
  }
  strm << ")\n";

  if (verbosity == 2 || verbosity == 3) {
    if (eParent->isMultiTerminal()) {
      strm.put("MT");
    }
    if (eParent->isEVPlus()) {
      strm.put("EV+");
    }
    if (eParent->isEVTimes()) {
      strm.put("EV*");
    }
    if (eParent->isForRelations()) {
      strm.put("MxD");
    } else {
      strm.put("MDD");
    }
    strm.put(" rooted at this node:\n");
    if (eParent->isEVPlus()) {
      long ev = Inf<long>();
      getEdgeValue(ev);
      strm << "Dangling: " << ev << "\n";
    }
    eParent->showNodeGraph(strm, &node, 1);
  }
  if (verbosity == 1 || verbosity == 3) {
    double card;
    apply(CARDINALITY, *this, card);
    strm << "Cardinality of node " << long(node) << ": ";
    strm.put(card, 0, 8, 'e');
    strm.put('\n');
  }
}

void MEDDLY::dd_edge::write(output &s, const node_handle* map) const
{
  expert_forest* eParent = smart_cast<expert_forest*>(parent);

  if (!eParent->isMultiTerminal()) {
    eParent->writeEdgeValue(s, &raw_value);
    s.put(' ');
  }
  if (node > 0) {
    s.put(long(map[node]));
  } else {
    s.put(long(node));
  }
  s.put('\n');
}

void MEDDLY::dd_edge::read(forest* p, input &s, const node_handle* map)
{
  destroy();

  parent = p;
  expert_forest* eParent = smart_cast<expert_forest*>(parent);

  if (!eParent->isMultiTerminal()) {
    s.stripWS();
    eParent->readEdgeValue(s, &raw_value);
  }

  s.stripWS();
  long lnode = s.get_integer();
  if (lnode <= 0) {
    node = lnode;
  } else {
    node = map[lnode];
  }

  linkNode(parent, node);

  opPlus = 0;
  opStar = 0;
  opMinus = 0;
  opDivide = 0;

  if (parent) parent->registerEdge(*this);
  MEDDLY_DCASSERT(index);
}


void MEDDLY::dd_edge::writePicture(const char* filename, const char* extension) const
{
  if (parent) {
    expert_forest* eParent = smart_cast<expert_forest*>(parent);
    eParent->writeNodeGraphPicture(filename, extension, &node, &label, 1);
  }
}

#endif
