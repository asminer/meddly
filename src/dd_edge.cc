
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

// #define DEBUG_ITER_BEGIN

// Constructor.
MEDDLY::dd_edge::dd_edge(forest* p)
: parent(p),
  node(0), value(0), level(0), index(-1),
  opPlus(0), opStar(0), opMinus(0), opDivide(0),
  updateNeeded(true), beginIterator(0)
{
  DCASSERT(p != NULL);
  smart_cast<expert_forest*>(parent)->registerEdge(*this);
  DCASSERT(index != -1);
}


// Destructor.  Will notify parent as appropriate.
MEDDLY::dd_edge::~dd_edge()
{
  if (index != -1) {
    // still registered; unregister before discarding
    smart_cast<expert_forest*>(parent)->unlinkNode(node);
    smart_cast<expert_forest*>(parent)->unregisterEdge(*this);
  }
  if (beginIterator != 0) delete beginIterator;
}


// Copy Constructor.
MEDDLY::dd_edge::dd_edge(const dd_edge& e)
: parent(e.parent), node(e.node), value(e.value), level(e.level), index(-1),
  opPlus(e.opPlus), opStar(e.opStar),
  opMinus(e.opMinus), opDivide(e.opDivide),
  updateNeeded(e.updateNeeded), beginIterator(0)
{
  smart_cast<expert_forest*>(parent)->registerEdge(*this);
  smart_cast<expert_forest*>(parent)->linkNode(node);
  if (!updateNeeded) {
    beginIterator = new const_iterator(*(e.beginIterator));
  }
}


// Assignment operator.
MEDDLY::dd_edge& MEDDLY::dd_edge::operator=(const dd_edge& e)
{
  if (this != &e) {
    smart_cast<expert_forest*>(parent)->unlinkNode(node);
    if (parent != e.parent) {
      smart_cast<expert_forest*>(parent)->unregisterEdge(*this);
      DCASSERT(index == -1);
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

    updateNeeded = e.updateNeeded;
    if (updateNeeded) {
      beginIterator = 0;
    } else {
      beginIterator = new const_iterator(*(e.beginIterator));
    }
  }
  return *this;
}


void MEDDLY::dd_edge::getEdgeValue(float& ev) const
{
  ev = toFloat(value);
}


void MEDDLY::dd_edge::set(int n, int v, int l)
{
  if (node != n) { updateNeeded = true; }
  smart_cast<expert_forest*>(parent)->unlinkNode(node);
  node = n;
  value = v;
  level = l;
}


void MEDDLY::dd_edge::set(int n, float v, int l)
{
  if (node != n) { updateNeeded = true; }
  smart_cast<expert_forest*>(parent)->unlinkNode(node);
  node = n;
  value = toInt(v);
  level = l;
}


double MEDDLY::dd_edge::getCardinality() const
{
#if 0
  return smart_cast<expert_forest*>(parent)->getCardinality(node);
#else
  double cardinality = 0;
  getComputeManager()->apply(compute_manager::CARDINALITY,
      *this, cardinality);
  return cardinality;
#endif
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
    const int nOperands = 3;
    op_param plist[nOperands];
    plist[0].set(parent);
    plist[1].set(parent);
    plist[2].set(parent);
    compute_manager::op_code opCode =
      (parent->getRangeType() == forest::BOOLEAN &&
        parent->getEdgeLabeling() == forest::MULTI_TERMINAL)
      ? compute_manager::UNION
      : compute_manager::PLUS;
    opPlus =
      smart_cast<expert_compute_manager*>(getComputeManager())->
      getOpInfo(opCode, plist, nOperands);
    DCASSERT(opPlus != 0);
  }
  DCASSERT(e.parent == parent);
  smart_cast<expert_compute_manager*>(getComputeManager())->
    apply(opPlus, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
  return *this;
}


// Operator *=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator*=(const dd_edge& e)
{
  if (opStar == 0) {
    const int nOperands = 3;
    op_param plist[nOperands];
    plist[0].set(parent);
    plist[1].set(parent);
    plist[2].set(parent);
    compute_manager::op_code opCode =
      (parent->getRangeType() == forest::BOOLEAN &&
        parent->getEdgeLabeling() == forest::MULTI_TERMINAL)
      ? compute_manager::INTERSECTION
      : compute_manager::MULTIPLY;
    opStar =
      smart_cast<expert_compute_manager*>(getComputeManager())->
      getOpInfo(opCode, plist, nOperands);
    DCASSERT(opStar != 0);
  }
  DCASSERT(e.parent == parent);
  smart_cast<expert_compute_manager*>(getComputeManager())->
    apply(opStar, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
  return *this;
}


// Operator -=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator-=(const dd_edge& e)
{
  if (opMinus == 0) {
    const int nOperands = 3;
    op_param plist[nOperands];
    plist[0].set(parent);
    plist[1].set(parent);
    plist[2].set(parent);
    compute_manager::op_code opCode =
      (parent->getRangeType() == forest::BOOLEAN &&
        parent->getEdgeLabeling() == forest::MULTI_TERMINAL)
      ? compute_manager::DIFFERENCE
      : compute_manager::MINUS;
    opMinus =
      smart_cast<expert_compute_manager*>(getComputeManager())->
      getOpInfo(opCode, plist, nOperands);
    DCASSERT(opMinus != 0);
  }
  DCASSERT(e.parent == parent);
  smart_cast<expert_compute_manager*>(getComputeManager())->
    apply(opMinus, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
  return *this;
}


// Operator /=
MEDDLY::dd_edge& MEDDLY::dd_edge::operator/=(const dd_edge& e)
{
  if (opDivide == 0) {
    const int nOperands = 3;
    op_param plist[nOperands];
    plist[0].set(parent);
    plist[1].set(parent);
    plist[2].set(parent);
    DCASSERT(!(parent->getRangeType() == forest::BOOLEAN &&
        parent->getEdgeLabeling() == forest::MULTI_TERMINAL));
    opDivide =
      smart_cast<expert_compute_manager*>(getComputeManager())->
      getOpInfo(compute_manager::DIVIDE, plist, nOperands);
  }
  DCASSERT(e.parent == parent);
  smart_cast<expert_compute_manager*>(getComputeManager())->
    apply(opDivide, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
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
      DCASSERT(eParent->getRangeType() == forest::BOOLEAN);
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


MEDDLY::dd_edge::iterator MEDDLY::dd_edge::begin()
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  DCASSERT(beginIterator != 0);
  return *beginIterator;
}


MEDDLY::dd_edge::iterator MEDDLY::dd_edge::beginRow(const int* minterm)
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  if (this->parent->isForRelations())
    // return iterator(this, true, minterm);
    return iterator(this, iterator::ROW, minterm);
  else
    return iterator();
}


MEDDLY::dd_edge::iterator MEDDLY::dd_edge::beginColumn(const int* minterm)
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  // return iterator(this, false, minterm);
  return iterator(this, iterator::COLUMN, minterm);
}


void MEDDLY::dd_edge::updateIterators()
{
  // update beginIterator
  if (beginIterator != 0) delete beginIterator;

  beginIterator = new iterator(this, iterator::DEFAULT, 0);
}


MEDDLY::dd_edge::iterator::iterator()
: e(0), size(0), element(0), nodes(0), pelement(0), pnodes(0), type(DEFAULT)
{
  foundPathAtLevel = 0;
}


MEDDLY::dd_edge::iterator::iterator(dd_edge* e, iter_type t, const int* minterm)
: e(e), size(0), element(0), nodes(0), pelement(0), pnodes(0), type(t)
{
  if (e == 0) return;
  if (type == ROW || type == COLUMN) {
    if (!e->parent->isForRelations()) return;
  }

  foundPathAtLevel = e->parent->getDomain()->getNumVariables();
  size = e->parent->getDomain()->getNumVariables() + 1;
  element = (int*) malloc(size * sizeof(int));
  nodes = (int*) malloc(size * sizeof(int));
  memset(element, 0, size * sizeof(int));
  memset(nodes, 0, size * sizeof(int));

  if (e->parent->isForRelations()) {
    pelement = (int*) malloc(size * sizeof(int));
    pnodes = (int*) malloc(size * sizeof(int));
    memset(pelement, 0, size * sizeof(int));
    memset(pnodes, 0, size * sizeof(int));
  }

#ifdef DEBUG_ITER_BEGIN
  printf("Begin (nodes, element): [");
  for (int i = size-1; i > -1; --i)
  {
    printf("%d:%d ", nodes[i], element[i]);
  }
  if (e->parent->isForRelations()) {
    printf("] --> [");
    for (int i = size-1; i > -1; --i)
    {
      printf("%d:%d ", pnodes[i], pelement[i]);
    }
  }
  printf("]\n");
#endif

  switch (type) {
    case DEFAULT:
      // do findFirstElement() and set element[] and nodes[]
      ++(*this);
      break;

    case ROW:
      memcpy(element, minterm, size * sizeof(int));
      findFirstColumn(e->parent->getDomain()->getNumVariables(), e->node);
      break;

    case COLUMN:
      memcpy(pelement, minterm, size * sizeof(int));
      findFirstRow(e->parent->getDomain()->getNumVariables(), e->node);
      break;

    default:
      break;
  }

#ifdef DEBUG_ITER_BEGIN
  printf("Begin (nodes, element): [");
  for (int i = size-1; i > -1; --i)
  {
    printf("%d:%d ", nodes[i], element[i]);
  }
  if (e->parent->isForRelations()) {
    printf("] --> [");
    for (int i = size-1; i > -1; --i)
    {
      printf("%d:%d ", pnodes[i], pelement[i]);
    }
  }
  printf("]\n");
#endif

}


MEDDLY::dd_edge::iterator::~iterator()
{
  if (e != 0) {
    free(element);
    free(nodes);
    if (pelement != 0) free(pelement);
    if (pnodes != 0) free (pnodes);
  }
}


MEDDLY::dd_edge::iterator::iterator(const iterator& iter)
: e(iter.e), size(iter.size),
  element(0), nodes(0), pelement(0), pnodes(0), type(iter.type)
{
  if (e != 0) {
    element = (int*) malloc(size * sizeof(int));
    nodes = (int*) malloc(size * sizeof(int));
    memcpy(element, iter.element, size * sizeof(int));
    memcpy(nodes, iter.nodes, size * sizeof(int));

    if (e->parent->isForRelations()) {
      pelement = (int*) malloc(size * sizeof(int));
      pnodes = (int*) malloc(size * sizeof(int));
      memcpy(pelement, iter.pelement, size * sizeof(int));
      memcpy(pnodes, iter.pnodes, size * sizeof(int));
    }

    foundPathAtLevel = e->parent->getDomain()->getNumVariables();
  }
}


MEDDLY::dd_edge::iterator& MEDDLY::dd_edge::iterator::operator=(const iterator& iter)
{
  if (this != &iter) {
    if (e != 0) {
      free(element);
      free(nodes);
      if (pelement != 0) free(pelement);
      if (pnodes != 0) free(pnodes);
      e = 0;
      size = 0;
      element = nodes = pelement = pnodes = 0;
      type = DEFAULT;
    }
    if (iter.e != 0) {
      e = iter.e;
      size = iter.size;
      type = iter.type;

      element = (int*) malloc(size * sizeof(int));
      nodes = (int*) malloc(size * sizeof(int));
      memcpy(element, iter.element, size * sizeof(int));
      memcpy(nodes, iter.nodes, size * sizeof(int));

      if (e->parent->isForRelations()) {
        pelement = (int*) malloc(size * sizeof(int));
        pnodes = (int*) malloc(size * sizeof(int));
        memcpy(pelement, iter.pelement, size * sizeof(int));
        memcpy(pnodes, iter.pnodes, size * sizeof(int));
      }

      foundPathAtLevel = e->parent->getDomain()->getNumVariables();
    }
  }
  return *this;
}


// PRE: minterm[] is stored in pelement[]
bool MEDDLY::dd_edge::iterator::findNextRow(int height)
{
  // See if you advance at a lower height
  // If yes, return true.
  // Otherwise, pick next index at this height, and find next path.
  // If no such path exists, try next index until you run out of indexes.
  // If you run our of indexes return false.
  // Else, return true.

  if (height == 0) {
    nodes[0] = 0;
    return false;
  }

  DCASSERT(e != 0);
  DCASSERT(type == COLUMN);
  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  DCASSERT(f->getReductionRule() == forest::IDENTITY_REDUCED);

  int nextHeight = height - 1;

  if (findNextRow(nextHeight)) return true;

  // Find next available index at the unprimed level
  // and call findNextRow() until it returns
  // true or till you run out of indexes.

  // pelement[height] stays constant
  // modify only element[height]

  if (nodes[height] == 0) {
    // Previously skipped level.
    // Since pelement[height] cannot be changed,
    // and since this is a skipped level, element[i] must equal pelement[i].
    // Therefore, no other indexes are available at this level.
    DCASSERT(pnodes[height] == 0);
    element[height] = 0;
    return false;
  }

  DCASSERT(!f->isTerminalNode(nodes[height]));
  DCASSERT(!f->isTerminalNode(pnodes[height]));

  bool found = false;

  const int* dptrs = 0;
  assert(f->getDownPtrs(nodes[height], dptrs));

  // Select next index at the primed level
  if (f->isFullNode(nodes[height])) {
    const int sz = f->getFullNodeSize(nodes[height]);
    for (int i = element[height] + 1; i < sz; i++)
    {
      if (dptrs[i] != 0) {
        // explore path
        // dptrs[i] must be a primed node
        DCASSERT(f->getNodeLevel(dptrs[i]) == -height);
        int unprimedNode = f->getDownPtr(dptrs[i], pelement[height]);
        if (findFirstRow(nextHeight, unprimedNode)) {
          DCASSERT(nodes[f->getNodeLevel(unprimedNode)] == unprimedNode);
          pnodes[height] = dptrs[i];
          element[height] = i;
          found = true;
          break;
        }
      }
    }
  }
  else {
    DCASSERT(f->isSparseNode(nodes[height]));
    const int sz = f->getSparseNodeSize(nodes[height]);
    const int* iptrs = 0;
    assert(f->getSparseNodeIndexes(nodes[height], iptrs));
    int i = 0;
    for ( ; i < sz && iptrs[i] <= element[height]; i++);
    for ( ; i < sz; i++)
    {
      DCASSERT(dptrs[i] != 0);
      // Explore path.
      // dptrs[i] must be primed.
      DCASSERT(f->getNodeLevel(dptrs[i]) == -height);
      int unprimedNode = f->getDownPtr(dptrs[i], pelement[height]);
      if (findFirstRow(nextHeight, unprimedNode)) {
        DCASSERT(nodes[f->getNodeLevel(unprimedNode)] == unprimedNode);
        pnodes[height] = dptrs[i];
        element[height] = iptrs[i];
        found = true;
        break;
      }
    }
  }

  if (!found) {
    element[height] = 0;
    pnodes[height] = 0;
    nodes[height] = 0;
  } else {
    foundPathAtLevel = height;
  }

  return found;
}


// PRE: minterm[] is stored in element[]
bool MEDDLY::dd_edge::iterator::findNextColumn(int height)
{
  // See if you advance at a lower height
  // If yes, return true.
  // Otherwise, pick next index at this height, and find next path.
  // If no such path exists, try next index until you run out of indexes.
  // If you run our of indexes return false.
  // Else, return true.

  if (height == 0) {
    nodes[0] = 0;
    return false;
  }

  DCASSERT(e != 0);
  DCASSERT(type == ROW);
  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  DCASSERT(f->getReductionRule() == forest::IDENTITY_REDUCED);

  int nextHeight = height - 1;

  if (findNextColumn(nextHeight)) return true;

  // Find next available index at the primed level
  // and call findNextColumn() until it returns
  // true or till you run out of indexes.

  // element[height] stays constant
  // modify only pelement[height]

  if (nodes[height] == 0) {
    // Previously skipped level.
    // Since element[height] cannot be changed,
    // and since this is a skipped level, pelement[i] must equal element[i].
    // Therefore, no other indexes are available at this level.
    DCASSERT(pnodes[height] == 0);
    pelement[height] = 0;
    return false;
  }

  DCASSERT(!f->isTerminalNode(nodes[height]));
  DCASSERT(!f->isTerminalNode(pnodes[height]));

  bool found = false;

  const int* dptrs = 0;
  assert(f->getDownPtrs(pnodes[height], dptrs));
  
  // Select next index at the primed level
  if (f->isFullNode(pnodes[height])) {
    const int sz = f->getFullNodeSize(pnodes[height]);
    for (int i = pelement[height] + 1; i < sz; i++)
    {
      if (dptrs[i] != 0) {
        // explore path
        // dptrs[i] must be an unprimed node
        DCASSERT(f->getNodeLevel(dptrs[i]) >= 0);
        if (findFirstColumn(nextHeight, dptrs[i])) {
          DCASSERT(nodes[f->getNodeLevel(dptrs[i])] == dptrs[i]);
          pelement[height] = i;
          found = true;
          break;
        }
      }
    }
  }
  else {
    DCASSERT(f->isSparseNode(pnodes[height]));
    const int sz = f->getSparseNodeSize(pnodes[height]);
    const int* iptrs = 0;
    assert(f->getSparseNodeIndexes(pnodes[height], iptrs));
    int i = 0;
    for ( ; i < sz && iptrs[i] <= pelement[height]; i++);
    for ( ; i < sz; i++)
    {
      DCASSERT(dptrs[i] != 0);
      // Explore path.
      // dptrs[i] must be unprimed.
      DCASSERT(f->getNodeLevel(dptrs[i]) >= 0);
      if (findFirstColumn(nextHeight, dptrs[i])) {
        DCASSERT(nodes[f->getNodeLevel(dptrs[i])] == dptrs[i]);
        pelement[height] = iptrs[i];
        found = true;
        break;
      }
    }
  }

  if (!found) {
    pelement[height] = 0;
    pnodes[height] = 0;
    nodes[height] = 0;
  } else {
    foundPathAtLevel = -height;
  }

  return found;
}


// Return true, if a row has been found starting at given level, and
// the element and node vectors have been filled (at and below given level).
// PRE: minterm[] have been copied into pelement[]
bool MEDDLY::dd_edge::iterator::findFirstRow(int height, int node)
{
  DCASSERT(e != 0);
  DCASSERT(type == COLUMN);
  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  DCASSERT(f->getReductionRule() == forest::IDENTITY_REDUCED);

  int nodeHeight = f->isTerminalNode(node)? 0: f->getNodeHeight(node);

  // Fill skipped levels
  for (int i = height; i != nodeHeight; i--) {
    element[i] = pelement[i];
    nodes[i] = pnodes[i] = 0;
  }

  if (f->isTerminalNode(node)) {
    // If terminal zero, then return false.
    // Fill up the node and element vectors at the terminal level.
    int i = nodeHeight;
    element[i] = pelement[i] = pnodes[i] = 0;
    nodes[i] = node;
    return node != 0;
  }

  DCASSERT(height > 0);

  // node must be an unprimed node (based on identity-reduction rule).
  DCASSERT(f->getNodeLevel(node) > 0);

  // Find the first i such that node[i][pelement[nodeHeight]] != 0

  int nextHeight = height - 1;
  bool found = false;

  const int* dptrs = 0;
  if (!f->getDownPtrs(node, dptrs)) return false;

  if (f->isFullNode(node)) {
    const int sz = f->getFullNodeSize(node);
    for (int i = 0; i < sz; i++) {
      if (dptrs[i] != 0) {
        // explore path
        int unprimedNode = f->getDownPtr(dptrs[i], pelement[nodeHeight]);
        if (findFirstRow(nextHeight, unprimedNode)) {
          // found a path, break
          nodes[nodeHeight] = node;
          pnodes[nodeHeight] = dptrs[i];
          element[nodeHeight] = i;
          found = true;
          break;
        }
      }
    }
  }
  else {
    DCASSERT(f->isSparseNode(node));
    const int sz = f->getSparseNodeSize(node);
    for (int i = 0; i < sz; i++) {
      DCASSERT(dptrs[i] != 0);
      // explore path
      int unprimedNode = f->getDownPtr(dptrs[i], pelement[nodeHeight]);
      if (findFirstRow(nextHeight, unprimedNode)) {
        // found a path, break
        nodes[nodeHeight] = node;
        pnodes[nodeHeight] = dptrs[i];
        element[nodeHeight] = f->getSparseNodeIndex(node, i);
        found = true;
        break;
      }
    }
  }

  return found;
}


// Return true, if a column has been found starting at given level, and
// the element and node vectors have been filled (at and below given level).
// PRE: minterm[] have been copied into element[]
bool MEDDLY::dd_edge::iterator::findFirstColumn(int height, int node)
{
  DCASSERT(e != 0);
  DCASSERT(type == ROW);
  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  DCASSERT(f->getReductionRule() == forest::IDENTITY_REDUCED);

  int nodeHeight = f->isTerminalNode(node)? 0: f->getNodeHeight(node);

  // Fill skipped levels
  for (int i = height; i != nodeHeight; i--) {
    pelement[i] = element[i];
    nodes[i] = pnodes[i] = 0;
  }

  if (f->isTerminalNode(node)) {
    // If terminal zero, then return false.
    // Fill up the node and element vectors at the terminal level.
    int i = nodeHeight;
    element[i] = pelement[i] = pnodes[i] = 0;
    nodes[i] = node;
    return node != 0;
  }

  DCASSERT(height > 0);

  // node must be an unprimed node (based on identity-reduction rule).
  DCASSERT(f->getNodeLevel(node) > 0);

  // Find the first i such that node[element[nodeHeight]][i] != 0

  int nextHeight = height - 1;
  bool found = false;
  int primedNode = f->getDownPtr(node, element[nodeHeight]);

  const int* dptrs = 0;
  if (!f->getDownPtrs(primedNode, dptrs)) return false;

  if (f->isFullNode(primedNode)) {
    // Find first i such that primedNode[i] does not lead to a 0.
    const int sz = f->getFullNodeSize(primedNode);
    for (int i = 0; i < sz; i++)
    {
      if (dptrs[i] != 0) {
        // explore path
        if (findFirstColumn(nextHeight, dptrs[i])) {
          // found a path, break
          pnodes[nodeHeight] = primedNode;
          nodes[nodeHeight] = node;
          pelement[nodeHeight] = i;
          found = true;
          break;
        }
      }
    }
  }
  else {
    DCASSERT(f->isSparseNode(primedNode));
    // Find first i such that primedNode[i] does not lead to a 0.
    const int sz = f->getSparseNodeSize(primedNode);
    for (int i = 0; i < sz; i++)
    {
      DCASSERT(dptrs[i] != 0);
      // explore path
      if (findFirstColumn(nextHeight, dptrs[i])) {
        // found a path, break
        pnodes[nodeHeight] = primedNode;
        nodes[nodeHeight] = node;
        pelement[nodeHeight] = f->getSparseNodeIndex(primedNode, i);
        found = true;
        break;
      }
    }
  }

  return found;
}


void MEDDLY::dd_edge::iterator::incrNonRelation()
{
  DCASSERT(e != 0);
  DCASSERT(e->node != 0);
  DCASSERT(type == DEFAULT);

  expert_domain* d = smart_cast<expert_domain*>(e->parent->useDomain());
  expert_forest* f = smart_cast<expert_forest*>(e->parent);

  int currLevel = d->getNumVariables();
  foundPathAtLevel = currLevel;

  if (nodes[0] == 0) {
    memset(nodes, 0, size * sizeof(int));
    nodes[f->getNodeLevel(e->node)] = e->node;
  }
  else {
    // Start from the bottom level and see if you can find the next
    // valid edge. Move up a level if no more edges exist from this level.

    int lastSkipped = e->node;
    currLevel = 1;
    bool found = false;

    while (currLevel <= d->getNumVariables())
    {
      int node = nodes[currLevel];
      if (node == 0) {
        // level was skipped previously
        // see if next index is available, otherwise move up one level
        if (element[currLevel] + 1 < f->getLevelSize(currLevel)) {
          // index is available, use it and break out of loop
          element[currLevel]++;
          nodes[f->getNodeLevel(lastSkipped)] = lastSkipped;
          found = true;
          break;
        }
      }
      else if (f->isFullNode(node)) {
        // find next valid index
        for (int i = element[currLevel] + 1, j = f->getFullNodeSize(node);
            i < j; ++i)
        {
          if (f->getFullNodeDownPtr(node, i) != 0) {
            // found new path
            element[currLevel] = i;
            node = f->getFullNodeDownPtr(node, i);
            nodes[f->getNodeLevel(node)] = node;
            found = true;
            break;
          }
        }
        if (found) break;
      }
      else {
        DCASSERT(f->isSparseNode(node));
        // find sparse index corresponding to element[currLevel]
        int begin = 0;
        int end = f->getSparseNodeSize(node);
        if (element[currLevel] < f->getSparseNodeIndex(node, end-1)) {
          while (begin + 1 < end) {
            int mid = (begin + end) / 2;
            DCASSERT(mid > begin && mid < end);
            if (f->getSparseNodeIndex(node, mid) > element[currLevel]) {
              end = mid;
            } else {
              begin = mid;
            }
            //int mid = 
            //if (f->getSparseNodeDownPtr(node, begin)
          }
          // range is [begin, end), therefore
          DCASSERT(f->getSparseNodeIndex(node, begin) == element[currLevel]);
          DCASSERT(begin + 1 < f->getSparseNodeSize(node));
          element[currLevel] = f->getSparseNodeIndex(node, begin + 1);
          node = f->getSparseNodeDownPtr(node, begin + 1);
          nodes[f->getNodeLevel(node)] = node;
          found = true;
          break;
        }
      }
      // no new path from current level so reset nodes[currLevel]
      if (nodes[currLevel] != 0) lastSkipped = nodes[currLevel];
      nodes[currLevel] = 0;
      element[currLevel] = 0;
      currLevel++;
    }

    if (!found) {
      // no more paths
      nodes[0] = 0;
      return;
    }

    // found new path
    // select the first path starting from level below currLevel
    foundPathAtLevel = currLevel;
    currLevel--;
  }

  // select the first path starting at currLevel
  for ( ; currLevel; currLevel--)
  {
    int node = nodes[currLevel];
    if (node == 0) {
      element[currLevel] = 0;
    } else if (f->isFullNode(node)) {
      for (int i = 0, j = f->getFullNodeSize(node); i < j; ++i)
      {
        int n = f->getFullNodeDownPtr(node, i);
        if (n != 0) {
          element[currLevel] = i;
          nodes[f->getNodeLevel(n)] = n;
          break;
        }
      }
    } else {
      DCASSERT(f->isSparseNode(node));
      element[currLevel] = f->getSparseNodeIndex(node, 0);
      int n = f->getSparseNodeDownPtr(node, 0);
      nodes[f->getNodeLevel(n)] = n;
    }
  }
  DCASSERT(nodes[0] != 0);
#ifdef DEBUG_ITER_BEGIN
  printf("nodes[]: [");
  for (int i = size - 1; i > 0; --i)
  {
    printf("%d ", nodes[i]);
  }
  printf("]\n");
#endif
}


void MEDDLY::dd_edge::iterator::incrNonIdentRelation()
{
  DCASSERT(e != 0);
  DCASSERT(e->node != 0);
  DCASSERT(type == DEFAULT);

  expert_domain* d = smart_cast<expert_domain*>(e->parent->useDomain());
  expert_forest* f = smart_cast<expert_forest*>(e->parent);

  DCASSERT(f->getReductionRule() == forest::FULLY_REDUCED
      || f->getReductionRule() == forest::QUASI_REDUCED);

  int currLevel = d->getNumVariables();
  bool isCurrLevelPrime = false;
  foundPathAtLevel = currLevel;

  if (nodes[0] == 0) {
    memset(nodes, 0, size * sizeof(int));
    memset(pnodes, 0, size * sizeof(int));
    int nodeLevel = f->getNodeLevel(e->node);
    if (nodeLevel < 0) {
      pnodes[-nodeLevel] = e->node;
    } else {
      nodes[nodeLevel] = e->node;
    }
  }
  else {
    // Start from the bottom level and see if you can find the next
    // valid edge. Move up a level if no more edges exist from this level.

    int lastSkipped = e->node;
    currLevel = 1;
    isCurrLevelPrime = true;
    bool found = false;

    while (currLevel <= d->getNumVariables())
    {
      int node = isCurrLevelPrime? pnodes[currLevel]: nodes[currLevel];
      int& currElement = isCurrLevelPrime?
          pelement[currLevel]: element[currLevel];
      if (node == 0) {
        // level was skipped previously
        // see if next index is available, otherwise move up one level.
        if (currElement + 1 <
            d->getVariableBound(currLevel, isCurrLevelPrime)) {
          // index is available, use it and break out of loop
          currElement++;
          node = lastSkipped;
          found = true;
        }
      }
      else if (f->isFullNode(node)) {
        // find next valid index
        int i = 1 + (isCurrLevelPrime? pelement[currLevel]: element[currLevel]);
        int j = f->getFullNodeSize(node);
        for ( ; i < j; ++i)
        {
          if (f->getFullNodeDownPtr(node, i) != 0) {
            // found new path
            currElement = i;
            node = f->getFullNodeDownPtr(node, i);
            found = true;
            break;
          }
        }
      }
      else {
        DCASSERT(f->isSparseNode(node));
        // find sparse index corresponding to element[currLevel]
        int begin = 0;
        int end = f->getSparseNodeSize(node);
        int searchFor =
            isCurrLevelPrime? pelement[currLevel]: element[currLevel];
        if (searchFor < f->getSparseNodeIndex(node, end-1)) {
          while (begin + 1 < end) {
            int mid = (begin + end) / 2;
            DCASSERT(mid > begin && mid < end);
            if (f->getSparseNodeIndex(node, mid) > searchFor) {
              end = mid;
            } else {
              begin = mid;
            }
          }
          // range is [begin, end), therefore
          DCASSERT(f->getSparseNodeIndex(node, begin) == searchFor);
          DCASSERT(begin + 1 < f->getSparseNodeSize(node));
          currElement = f->getSparseNodeIndex(node, begin + 1);
          node = f->getSparseNodeDownPtr(node, begin + 1);
          found = true;
        }
      }

      if (found) {
        int nodeLevel = f->getNodeLevel(node);
        if (nodeLevel < 0) {
          pnodes[-nodeLevel] = node;
        } else {
          nodes[nodeLevel] = node;
        }
        break;
      }

      // no new path from current level so reset nodes[currLevel]
      currElement = 0;
      if (isCurrLevelPrime) {
        if (pnodes[currLevel] != 0) {
          lastSkipped = pnodes[currLevel];
          pnodes[currLevel] = 0;
        }
        isCurrLevelPrime = false;
      } else {
        if (nodes[currLevel] != 0) {
          lastSkipped = nodes[currLevel];
          nodes[currLevel] = 0;
        }
        currLevel++;
        isCurrLevelPrime = true;
      }
    }

    if (!found) {
      // no more paths
      nodes[0] = 0;
      return;
    }

    // found new path
    // select the first path starting from level below currLevel
    if (isCurrLevelPrime) {
      foundPathAtLevel = -currLevel;
      isCurrLevelPrime = false;
      currLevel--;
    } else {
      foundPathAtLevel = currLevel;
      isCurrLevelPrime = true;
    }
  }

  // select the first path starting at currLevel
  for ( ; currLevel; )
  {
    int node = isCurrLevelPrime? pnodes[currLevel]: nodes[currLevel];
    int& currElement = isCurrLevelPrime
        ? pelement[currLevel]: element[currLevel];
    if (node == 0) {
      currElement = 0;
    } else {
      int n = -1;
      int index = -1;
      if (f->isFullNode(node)) {
        for (int i = 0, j = f->getFullNodeSize(node); i < j; ++i)
        {
          n = f->getFullNodeDownPtr(node, i);
          if (n != 0) { index = i; break; }
        }
      } else {
        DCASSERT(f->isSparseNode(node));
        index = f->getSparseNodeIndex(node, 0);
        n = f->getSparseNodeDownPtr(node, 0);
      }
      DCASSERT(index != -1);
      currElement = index;
      int nodeLevel = f->getNodeLevel(n);
      if (nodeLevel < 0) {
        pnodes[-nodeLevel] = n;
      } else {
        nodes[nodeLevel] = n;
      }
    }

    if (isCurrLevelPrime) {
      isCurrLevelPrime = false;
      currLevel--;
    } else {
      isCurrLevelPrime = true;
    }
  }
  DCASSERT(nodes[0] != 0);
#ifdef DEBUG_ITER_BEGIN
  printf("nodes[]: [");
  for (int i = size - 1; i > 0; --i)
  {
    printf("%d ", nodes[i]);
  }
  printf("]\n");
#endif
}


void MEDDLY::dd_edge::iterator::incrRelation()
{
  DCASSERT(e != 0);
  DCASSERT(e->node != 0);
  DCASSERT(type == DEFAULT);

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  if (f->getReductionRule() != forest::IDENTITY_REDUCED) {
    incrNonIdentRelation();
    return;
  }

  expert_domain* d = smart_cast<expert_domain*>(e->parent->useDomain());
  int currLevel = d->getNumVariables();
  bool isCurrLevelPrime = false;
  foundPathAtLevel = currLevel;

  if (nodes[0] == 0) {
    memset(nodes, 0, size * sizeof(int));
    memset(pnodes, 0, size * sizeof(int));
    int nodeLevel = f->getNodeLevel(e->node);
    assert(nodeLevel >= 0);
    nodes[nodeLevel] = e->node;
    //TODO: test mtmxd iterator
  }
  else {
    // Start from the bottom level and see if you can find the next
    // valid edge. Move up a level if no more edges exist from this level.

    int lastSkipped = e->node;
    currLevel = 1;
    isCurrLevelPrime = true;
    bool found = false;

    while (currLevel <= d->getNumVariables())
    {
      int node = isCurrLevelPrime? pnodes[currLevel]: nodes[currLevel];
      if (node == 0) {
        // level was skipped previously
        // see if next index is available, otherwise move up one level
        // since this forest uses identity-reduction, look for an
        // available index only at the prime level (the index at the
        // unprime level is the same).
        if (isCurrLevelPrime) {
          if (element[currLevel] + 1 < f->getLevelSize(currLevel)) {
            // index is available, use it and break out of loop
            element[currLevel]++;
            pelement[currLevel] = element[currLevel];
            int nodeLevel = f->getNodeLevel(lastSkipped);
            if (nodeLevel < 0) {
              pnodes[-nodeLevel] = lastSkipped;
            } else {
              nodes[nodeLevel] = lastSkipped;
            }
            found = true;
          }
        }
      }
      else if (f->isFullNode(node)) {
        // find next valid index
        int i = 1 + (isCurrLevelPrime? pelement[currLevel]: element[currLevel]);
        int j = f->getFullNodeSize(node);
        for ( ; i < j; ++i)
        {
          if (f->getFullNodeDownPtr(node, i) != 0) {
            // found new path
            if (isCurrLevelPrime) {
              pelement[currLevel] = i;
            } else {
              element[currLevel] = i;
            }
            node = f->getFullNodeDownPtr(node, i);
            int nodeLevel = f->getNodeLevel(node);
            if (nodeLevel < 0) {
              pnodes[-nodeLevel] = node;
            } else {
              nodes[nodeLevel] = node;
            }
            found = true;
            break;
          }
        }
      }
      else {
        DCASSERT(f->isSparseNode(node));
        // find sparse index corresponding to element[currLevel]
        int begin = 0;
        int end = f->getSparseNodeSize(node);
        int searchFor =
            isCurrLevelPrime? pelement[currLevel]: element[currLevel];
        if (searchFor < f->getSparseNodeIndex(node, end-1)) {
          while (begin + 1 < end) {
            int mid = (begin + end) / 2;
            DCASSERT(mid > begin && mid < end);
            if (f->getSparseNodeIndex(node, mid) > searchFor) {
              end = mid;
            } else {
              begin = mid;
            }
          }
          // range is [begin, end), therefore
          DCASSERT(f->getSparseNodeIndex(node, begin) == searchFor);
          DCASSERT(begin + 1 < f->getSparseNodeSize(node));
          if (isCurrLevelPrime) {
            pelement[currLevel] = f->getSparseNodeIndex(node, begin + 1);
          } else {
            element[currLevel] = f->getSparseNodeIndex(node, begin + 1);
          }
          node = f->getSparseNodeDownPtr(node, begin + 1);
          int nodeLevel = f->getNodeLevel(node);
          if (nodeLevel < 0) {
            pnodes[-nodeLevel] = node;
          } else {
            nodes[nodeLevel] = node;
          }
          found = true;
        }
      }

      if (found) break;

      // no new path from current level so reset nodes[currLevel]
      if (isCurrLevelPrime) {
        if (pnodes[currLevel] != 0) lastSkipped = pnodes[currLevel];
        pnodes[currLevel] = 0;
        pelement[currLevel] = 0;
        isCurrLevelPrime = false;
      } else {
        if (nodes[currLevel] != 0) lastSkipped = nodes[currLevel];
        nodes[currLevel] = 0;
        element[currLevel] = 0;
        currLevel++;
        isCurrLevelPrime = true;
      }
    }

    if (!found) {
      // no more paths
      nodes[0] = 0;
      return;
    }

    // found new path
    // select the first path starting from level below currLevel
    if (isCurrLevelPrime) {
      foundPathAtLevel = -currLevel;
      isCurrLevelPrime = false;
      currLevel--;
    } else {
      foundPathAtLevel = currLevel;
      isCurrLevelPrime = true;
    }
  }

  // select the first path starting at currLevel
  for ( ; currLevel; )
  {
    int node = isCurrLevelPrime? pnodes[currLevel]: nodes[currLevel];
    if (node == 0) {
      DCASSERT(!isCurrLevelPrime);
      element[currLevel] = 0;
      pelement[currLevel] = 0;
      // jump over the prime level
      isCurrLevelPrime = true;
    } else {
      int n = -1;
      int index = -1;
      if (f->isFullNode(node)) {
        for (int i = 0, j = f->getFullNodeSize(node); i < j; ++i)
        {
          n = f->getFullNodeDownPtr(node, i);
          if (n != 0) { index = i; break; }
        }
      } else {
        DCASSERT(f->isSparseNode(node));
        index = f->getSparseNodeIndex(node, 0);
        n = f->getSparseNodeDownPtr(node, 0);
      }
      DCASSERT(index != -1);
      if (isCurrLevelPrime) {
        pelement[currLevel] = index;
      } else {
        element[currLevel] = index;
      }
      int nodeLevel = f->getNodeLevel(n);
      if (nodeLevel < 0) {
        pnodes[-nodeLevel] = n;
      } else {
        nodes[nodeLevel] = n;
      }
    }

    if (isCurrLevelPrime) {
      isCurrLevelPrime = false;
      currLevel--;
    } else {
      isCurrLevelPrime = true;
    }
  }
  DCASSERT(nodes[0] != 0);
#ifdef DEBUG_ITER_BEGIN
  printf("nodes[]: [");
  for (int i = size - 1; i > 0; --i)
  {
    printf("%d ", nodes[i]);
  }
  printf("]\n");
#endif
}


void MEDDLY::dd_edge::iterator::operator++()
{
  // find next
  // set element to next

  if (e != 0 && e->node != 0) {
    if (e->parent->isForRelations()) {
      switch (type) {
        case DEFAULT:
          incrRelation();
          break;
        case ROW:
          if (!findNextColumn(e->parent->getDomain()->getNumVariables()))
            foundPathAtLevel = -(e->parent->getDomain()->getNumVariables());
          break;
        case COLUMN:
          if (!findNextRow(e->parent->getDomain()->getNumVariables()))
            foundPathAtLevel = e->parent->getDomain()->getNumVariables();
          break;
        default:
          break;
      }
    } else {
      DCASSERT(type == DEFAULT);
      incrNonRelation();
    }
  }
}


void MEDDLY::dd_edge::iterator::operator--()
{
  // find prev
  // set element to prev

  if (e == 0) return;
  if (e->node == 0) return;

  if (nodes[0] == 0) {
    // reached end, find last element in dd_edge
    assert(false);
  } else {
    // find prev element
    assert(false);
  }
}


bool MEDDLY::dd_edge::iterator::operator!=(const iterator& iter) const
{
#if 0
  DCASSERT((e != iter.e) || (size == iter.size));

  // if terminals are different, return true.
  // else if terminals are 0, return false.
  // else return (element != iter.element)
  return e == 0
          ? true
          : e != iter.e
            ? true
            : type != iter.type
              ? true
              : (nodes[0] != iter.nodes[0])
                ? true
                : (nodes[0] == 0)
                  ? false
                  : e->parent->isForRelations()
                    ? (0 != memcmp(element, iter.element, size) ||
                      0 != memcmp(pelement, iter.pelement, size))
                    : 0 != memcmp(element, iter.element, size);
#else
  return !(this->operator==(iter));
#endif
}


bool MEDDLY::dd_edge::iterator::operator==(const iterator& iter) const
{
  DCASSERT((e != iter.e) || (size == iter.size));

#if 0
  // if both terminals are 0 return true.
  // if only one terminal is 0 return false.
  // if neither terminal is 0, return (element == iter.element)
  return !(*this != iter);
#else
  if (e != iter.e) return false;
  if (type != iter.type) return false;
  if (e == 0) return true;
  if (e->parent->isForRelations()) {
    if (0 != memcmp(pelement, iter.pelement, size)) return false;
  }
  return 0 == memcmp(element, iter.element, size);
#endif
}


void MEDDLY::dd_edge::iterator::getValue(int& val) const
{
  if (e == 0 || nodes[0] == 0) return;
  if (e->parent->getEdgeLabeling() == forest::MULTI_TERMINAL)
    val = nodes[0];
  else {
    // return the value using evaluate
    e->parent->evaluate(*e, element, val);
  }
}


void MEDDLY::dd_edge::iterator::getValue(float& val) const
{
  if (e == 0 || nodes[0] == 0) return;
  if (e->parent->getEdgeLabeling() == forest::MULTI_TERMINAL)
    val = toFloat(nodes[0]);
  else {
    // return the value using evaluate
    e->parent->evaluate(*e, element, val);
  }
}


const int* MEDDLY::dd_edge::iterator::getAssignments() const
{
  return e == 0
          ? 0
          : nodes[0] == 0? 0: element;
}


const int* MEDDLY::dd_edge::iterator::getPrimedAssignments() const
{
  return e == 0
          ? 0
          : nodes[0] == 0? 0: e->parent->isForRelations()? pelement: 0;
}

