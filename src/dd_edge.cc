
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



#include "../include/mddexpert.h"


// Constructor.
dd_edge::dd_edge(forest* p)
: parent(p),
  node(0), value(0), level(0), index(-1),
  opPlus(0), opStar(0), opMinus(0), opDivide(0)
{
  assert(p != NULL);
  smart_cast<expert_forest*>(parent)->registerEdge(*this);
  assert(index != -1);
}


// Destructor.  Will notify parent as appropriate.
dd_edge::~dd_edge()
{
  if (index != -1) {
    // still registered; unregister before discarding
    smart_cast<expert_forest*>(parent)->unlinkNode(node);
    smart_cast<expert_forest*>(parent)->unregisterEdge(*this);
  }
}


// Copy Constructor.
dd_edge::dd_edge(const dd_edge& e)
: parent(e.parent), node(e.node), value(e.value), level(e.level), index(-1),
  opPlus(e.opPlus), opStar(e.opStar),
  opMinus(e.opMinus), opDivide(e.opDivide)
{
  smart_cast<expert_forest*>(parent)->registerEdge(*this);
  smart_cast<expert_forest*>(parent)->linkNode(node);
}


// Assignment operator.
dd_edge& dd_edge::operator=(const dd_edge& e)
{
  if (this != &e) {
    smart_cast<expert_forest*>(parent)->unlinkNode(node);
    if (parent != e.parent) {
      smart_cast<expert_forest*>(parent)->unregisterEdge(*this);
      assert(index == -1);
      parent = e.parent;
      smart_cast<expert_forest*>(parent)->registerEdge(*this);
      opPlus = e.opPlus;
      opStar = e.opStar;
      opMinus = e.opMinus;
      opDivide = e.opDivide;
    }
    smart_cast<expert_forest*>(parent)->linkNode(e.node);
    node = e.node;
    value = e.value;
    level = e.level;
    // operation pointers remain the same since its the same parent
  }
  return *this;
}


void dd_edge::set(int n, int v, int l)
{
  smart_cast<expert_forest*>(parent)->unlinkNode(node);
  node = n;
  value = v;
  level = l;
}


// Operator +=
dd_edge& dd_edge::operator+=(const dd_edge& e)
{
  if (opPlus == 0) {
    const int nOperands = 3;
    forest* forests[nOperands] = {parent, parent, parent};
    compute_manager::op_code opCode =
      (parent->getRangeType() == forest::BOOLEAN)
      ? compute_manager::UNION
      : compute_manager::PLUS;
    opPlus =
      smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
      getOpInfo(opCode, forests, nOperands);
    assert(opPlus != 0);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
    apply(opPlus, *this, e, *this);
  return *this;
}


// Operator *=
dd_edge& dd_edge::operator*=(const dd_edge& e)
{
  if (opStar == 0) {
    const int nOperands = 3;
    forest* forests[nOperands] = {parent, parent, parent};
    compute_manager::op_code opCode =
      (parent->getRangeType() == forest::BOOLEAN)
      ? compute_manager::INTERSECTION
      : compute_manager::MULTIPLY;
    opStar =
      smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
      getOpInfo(opCode, forests, nOperands);
    assert(opStar != 0);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
    apply(opStar, *this, e, *this);
  return *this;
}


// Operator -=
dd_edge& dd_edge::operator-=(const dd_edge& e)
{
  if (opMinus == 0) {
    const int nOperands = 3;
    forest* forests[nOperands] = {parent, parent, parent};
    compute_manager::op_code opCode =
      (parent->getRangeType() == forest::BOOLEAN)
      ? compute_manager::DIFFERENCE
      : compute_manager::MINUS;
    opMinus =
      smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
      getOpInfo(opCode, forests, nOperands);
    assert(opMinus != 0);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
    apply(opMinus, *this, e, *this);
  return *this;
}


// Operator /=
dd_edge& dd_edge::operator/=(const dd_edge& e)
{
  if (opDivide == 0) {
    const int nOperands = 3;
    forest* forests[nOperands] = {parent, parent, parent};
    assert(parent->getRangeType() != forest::BOOLEAN);
    opDivide =
      smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
      getOpInfo(compute_manager::DIVIDE, forests, nOperands);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MDDLIB_getComputeManager())->
    apply(opDivide, *this, e, *this);
  return *this;
}


// Display the edge information.
void dd_edge::show(FILE* strm, int verbosity) const
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
      assert(eParent->getRangeType() == forest::BOOLEAN);
      fprintf(strm, "node: %s*, ",(eParent->getBoolean(node)? "true": "false"));
    }
  }
  else {
    fprintf(strm, "node: %d, ", node);
  }
  fprintf(strm, "value: %d, level: %d)\n", value, level);
  if (verbosity == 2 || verbosity == 3) {
    fprintf(strm, "MDD rooted at this node:\n");
    smart_cast<expert_forest*>(parent)->showNodeGraph(strm, node);
  }
  if (verbosity == 1 || verbosity == 3) {
    fprintf(strm, "Cardinality of node %d: %0.8e\n", node,
        smart_cast<expert_forest*>(parent)->getCardinality(node));
  }
}


