
// $Id$

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
: parent(0),
  node(0), value(0), level(0), index(-1),
  opPlus(0), opStar(0), opMinus(0), opDivide(0)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %x\n", this);
#endif
}

// Constructor.
MEDDLY::dd_edge::dd_edge(forest* p)
: parent(p),
  node(0), value(0), level(0), index(-1),
  opPlus(0), opStar(0), opMinus(0), opDivide(0)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %x\n", this);
#endif
  MEDDLY_DCASSERT(p != NULL);
  parent->registerEdge(*this);
  MEDDLY_DCASSERT(index != -1);
}


// Copy Constructor.
MEDDLY::dd_edge::dd_edge(const dd_edge& e)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating dd_edge %x\n", this);
#endif
  init(e);
  MEDDLY_DCASSERT(index != -1);
}


// Assignment operator.
MEDDLY::dd_edge& MEDDLY::dd_edge::operator=(const dd_edge& e)
{
  if (&e != this) {
    destroy();
    init(e);
  }
  MEDDLY_DCASSERT(index != -1);
  return *this;
}

// Destructor.  Will notify parent as appropriate.
MEDDLY::dd_edge::~dd_edge()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Deleting dd_edge %x\n", this);
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
  value = e.value;
  level = e.level;

  linkNode(parent, node);

  opPlus = e.opPlus;
  opStar = e.opStar;
  opMinus = e.opMinus;
  opDivide = e.opDivide;

  if (parent) parent->registerEdge(*this);
  MEDDLY_DCASSERT(index != -1);
}

//
// destruction helper
//
void MEDDLY::dd_edge::destroy()
{
  if (index != -1) {
    // still registered; unregister before discarding
    int old = node;
    node = 0;
    unlinkNode(parent, old);
    if (parent) parent->unregisterEdge(*this);
  }
}

void MEDDLY::dd_edge::getEdgeValue(float& ev) const
{
  ev = toFloat(value);
}


void MEDDLY::dd_edge::set(int n, int v, int l)
{
  int old = node;
  node = n;
  unlinkNode(parent, old);
  value = v;
  level = l;
}


void MEDDLY::dd_edge::set(int n, float v, int l)
{
  int old = node;
  node = n;
  unlinkNode(parent, old);
  value = toInt(v);
  level = l;
}


unsigned MEDDLY::dd_edge::getNodeCount() const
{
  return smart_cast<expert_forest*>(parent)->getNodeCount(node);
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
  opPlus->compute(*this, e, *this);
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
  opStar->compute(*this, e, *this);
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
  opMinus->compute(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  return *this;
}


// Operator /=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator/=(const dd_edge& e)
{
  if (opDivide == 0) {
    opDivide = getOperation(DIVIDE, *this, e, *this);
  }
  opDivide->compute(*this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  return *this;
}


// Display the edge information.
void MEDDLY::dd_edge::show(FILE* strm, int verbosity) const
{
  expert_forest* eParent = smart_cast<expert_forest*>(parent);
  fprintf(strm, "(Forest Addr: %p, ", parent);
  if (eParent->isTerminalNode(node)) {
    if (eParent->getRangeType() == forest::REAL) {
      fprintf(strm, "node: %f*, ", eParent->getReal(node));
    }
    else if (eParent->getRangeType() == forest::INTEGER) {
      fprintf(strm, "node: %d*, ", eParent->getInteger(node));
    }
    else {
      MEDDLY_DCASSERT(eParent->getRangeType() == forest::BOOLEAN);
      fprintf(strm, "node: %s*, ",(eParent->getBoolean(node)? "true": "false"));
    }
  }
  else {
    fprintf(strm, "node: %d, ", node);
  }
  if (eParent->getEdgeLabeling() == forest::EVTIMES) {
    // fprintf(strm, "value: %f, level: %d)\n", toFloat(value), level);
    float ev = 0;
    getEdgeValue(ev);
    fprintf(strm, "value: %f, level: %d)\n", ev, level);
  } else {
    fprintf(strm, "value: %d, level: %d)\n", value, level);
  }
  if (verbosity == 2 || verbosity == 3) {
    fprintf(strm, "MDD rooted at this node:\n");
    smart_cast<expert_forest*>(parent)->showNodeGraph(strm, node);
  }
  if (verbosity == 1 || verbosity == 3) {
    fprintf(strm, "Cardinality of node %d: %0.8e\n", node, getCardinality());
  }
}

