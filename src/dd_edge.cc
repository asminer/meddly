
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
#include "operations/mpz_object.h"
#include <cfloat>

// #define DEBUG_CLEANUP

// #define DEBUG_ITER_BEGIN

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

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        dd_edge  methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************

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
//  MEDDLY_DCASSERT(forest::MULTI_TERMINAL != parent->getEdgeLabeling());
//  MEDDLY_DCASSERT(forest::INTEGER == parent->getRangeType());
//  expert_forest::EVencoder<int>::readValue(&raw_value, ev);
//}

void MEDDLY::dd_edge::getEdgeValue(long& ev) const
{
  MEDDLY_DCASSERT(parent);
  MEDDLY_DCASSERT(forest::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(forest::INTEGER == parent->getRangeType());
  expert_forest::EVencoder<long>::readValue(&raw_value, ev);
}

void MEDDLY::dd_edge::getEdgeValue(float& ev) const
{
  MEDDLY_DCASSERT(parent);
  MEDDLY_DCASSERT(forest::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(forest::REAL == parent->getRangeType());
  expert_forest::EVencoder<float>::readValue(&raw_value, ev);
}

int MEDDLY::dd_edge::getLevel() const
{
  if (0==node) return 0;
  MEDDLY_DCASSERT(parent);
  const expert_forest* ef = dynamic_cast <expert_forest*>(parent);
  MEDDLY_DCASSERT(ef);
  return ef->getNodeLevel(node);
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
  MEDDLY_DCASSERT(forest::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(forest::INTEGER == parent->getRangeType());
  expert_forest::EVencoder<long>::writeValue(&raw_value, value);
}

void MEDDLY::dd_edge::setEdgeValue(float value)
{
  MEDDLY_DCASSERT(parent);
  MEDDLY_DCASSERT(forest::MULTI_TERMINAL != parent->getEdgeLabeling());
  MEDDLY_DCASSERT(forest::REAL == parent->getRangeType());
  expert_forest::EVencoder<float>::writeValue(&raw_value, value);
}

void MEDDLY::dd_edge::setLabel(const char* L)
{
  if (label) free(label);
  label = strdup(L);
}

unsigned MEDDLY::dd_edge::getNodeCount() const
{
  return smart_cast<expert_forest*>(parent)->getNodeCount(node);
}

unsigned MEDDLY::dd_edge::getLastHandle() const
{
  return smart_cast<expert_forest*>(parent)->getLastHandle(node);
}

unsigned MEDDLY::dd_edge::getEdgeCount(bool countZeroes) const
{
  return smart_cast<expert_forest*>(parent)->getEdgeCount(node, countZeroes);
}

//
// Operator +=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator+=(const dd_edge& e)
{
  if (opPlus == 0) {
    if (parent->getRangeType() == forest::BOOLEAN)
      opPlus = getOperation(UNION, *this, e, *this);
    else
      opPlus = getOperation(PLUS, *this, e, *this);
    MEDDLY_DCASSERT(opPlus != 0);
  }
  opPlus->computeTemp(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  return *this;
}


// Operator *=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator*=(const dd_edge& e)
{
  if (opStar == 0) {
    if (parent->getRangeType() == forest::BOOLEAN)
      opStar = getOperation(INTERSECTION, *this, e, *this);
    else
      opStar = getOperation(MULTIPLY, *this, e, *this);
    MEDDLY_DCASSERT(opStar != 0);
  }
  opStar->computeTemp(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  return *this;
}


// Operator -=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator-=(const dd_edge& e)
{
  if (opMinus == 0) {
    if (parent->getRangeType() == forest::BOOLEAN)
      opMinus = getOperation(DIFFERENCE, *this, e, *this);
    else
      opMinus = getOperation(MINUS, *this, e, *this);
    MEDDLY_DCASSERT(opMinus != 0);
  }
  opMinus->computeTemp(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  return *this;
}


// Operator /=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator/=(const dd_edge& e)
{
  if (opDivide == 0) {
    opDivide = getOperation(DIVIDE, *this, e, *this);
  }
  opDivide->computeTemp(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  return *this;
}


// Display the edge information.
void MEDDLY::dd_edge::show(output &strm, int verbosity) const
{
  expert_forest* eParent = smart_cast<expert_forest*>(parent);

  strm << "(Forest Addr: ";
  strm.put_hex((unsigned long) parent);
  strm << ", ";
  
  strm << "transparent: ";
  eParent->showTerminal(strm, eParent->getTransparentNode());
  strm << ", ";
  
  if (eParent->isTerminalNode(node)) {
    strm << "node: ";
    eParent->showTerminal(strm, node);
    strm << "*, ";
  }
  else {
    strm << "node: " << long(node) << ", ";
  }
  if (!eParent->isMultiTerminal()) {
    if (eParent->getRangeType() == forest::REAL) {
      float ev;
      getEdgeValue(ev);
      strm << "value: " << ev << ", ";
    } else {
      long iv = Inf<long>();
      getEdgeValue(iv);
      strm << "value: " << iv << ", ";
    }
  }
  strm << "level: " << getLevel();
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
    strm << "Cardinality of node " << long(node) << ": ";
    strm.put(getCardinality(), 0, 8, 'e');
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
