
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
dd_edge::dd_edge(forest* p)
: parent(p),
  node(0), value(0), level(0), index(-1),
  opPlus(0), opStar(0), opMinus(0), opDivide(0),
  updateNeeded(true), beginIterator(0), endIterator(0)
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
  if (beginIterator != 0) delete beginIterator;
  if (endIterator != 0) delete endIterator;
}


// Copy Constructor.
dd_edge::dd_edge(const dd_edge& e)
: parent(e.parent), node(e.node), value(e.value), level(e.level), index(-1),
  opPlus(e.opPlus), opStar(e.opStar),
  opMinus(e.opMinus), opDivide(e.opDivide),
  updateNeeded(e.updateNeeded), beginIterator(0), endIterator(0)
{
  smart_cast<expert_forest*>(parent)->registerEdge(*this);
  smart_cast<expert_forest*>(parent)->linkNode(node);
  if (!updateNeeded) {
    beginIterator = new const_iterator(*(e.beginIterator));
    endIterator = new const_iterator(*(e.endIterator));
  }
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

    updateNeeded = e.updateNeeded;
    if (updateNeeded) {
      beginIterator = 0; endIterator = 0;
    } else {
      beginIterator = new const_iterator(*(e.beginIterator));
      endIterator = new const_iterator(*(e.endIterator));
    }
  }
  return *this;
}


void dd_edge::getEdgeValue(float& ev) const
{
  ev = toFloat(value);
}


void dd_edge::set(int n, int v, int l)
{
  if (node != n) { updateNeeded = true; }
  smart_cast<expert_forest*>(parent)->unlinkNode(node);
  node = n;
  value = v;
  level = l;
}


void dd_edge::set(int n, float v, int l)
{
  if (node != n) { updateNeeded = true; }
  smart_cast<expert_forest*>(parent)->unlinkNode(node);
  node = n;
  value = toInt(v);
  level = l;
}


double dd_edge::getCardinality() const
{
  return smart_cast<expert_forest*>(parent)->getCardinality(node);
}

unsigned dd_edge::getNodeCount() const
{
  return smart_cast<expert_forest*>(parent)->getNodeCount(node);
}

unsigned dd_edge::getEdgeCount(bool countZeroes) const
{
  return smart_cast<expert_forest*>(parent)->getEdgeCount(node, countZeroes);
}

//
// Operator +=
dd_edge& dd_edge::operator+=(const dd_edge& e)
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
      smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
      getOpInfo(opCode, plist, nOperands);
    assert(opPlus != 0);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
    apply(opPlus, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
  return *this;
}


// Operator *=
dd_edge& dd_edge::operator*=(const dd_edge& e)
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
      smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
      getOpInfo(opCode, plist, nOperands);
    assert(opStar != 0);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
    apply(opStar, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
  return *this;
}


// Operator -=
dd_edge& dd_edge::operator-=(const dd_edge& e)
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
      smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
      getOpInfo(opCode, plist, nOperands);
    assert(opMinus != 0);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
    apply(opMinus, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
  return *this;
}


// Operator /=
dd_edge& dd_edge::operator/=(const dd_edge& e)
{
  if (opDivide == 0) {
    const int nOperands = 3;
    op_param plist[nOperands];
    plist[0].set(parent);
    plist[1].set(parent);
    plist[2].set(parent);
    assert(!(parent->getRangeType() == forest::BOOLEAN &&
        parent->getEdgeLabeling() == forest::MULTI_TERMINAL));
    opDivide =
      smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
      getOpInfo(compute_manager::DIVIDE, plist, nOperands);
  }
  assert(e.parent == parent);
  smart_cast<expert_compute_manager*>(MEDDLY_getComputeManager())->
    apply(opDivide, *this, e, *this);
  // apply will call set() which in turn will set updateNeeded to true
  DCASSERT(updateNeeded == true);
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
    fprintf(strm, "Cardinality of node %d: %0.8e\n", node,
        smart_cast<expert_forest*>(parent)->getCardinality(node));
  }
}


dd_edge::iterator dd_edge::begin()
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  DCASSERT(beginIterator != 0);
  return *beginIterator;
}


dd_edge::iterator dd_edge::end()
{
  if (updateNeeded) {
    updateIterators();
    updateNeeded = false;
  }
  DCASSERT(endIterator != 0);
  return *endIterator;
}


void dd_edge::updateIterators()
{
  // update beginIterator and endIterator
  if (beginIterator != 0) delete beginIterator;
  if (endIterator != 0) delete endIterator;

  beginIterator = new iterator(this, true);
  endIterator = new iterator(this, false);
}


dd_edge::iterator::iterator()
: e(0), size(0), element(0), nodes(0), pelement(0), pnodes(0)
{ }


dd_edge::iterator::iterator(dd_edge* e, bool begin)
: e(e), size(0), element(0), nodes(0), pelement(0), pnodes(0)
{
  if (e == 0) return;
  
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

  if (begin) {
    // do findFirstElement() and set element[] and nodes[]
#ifdef DEBUG_ITER_BEGIN
    printf("Begin: [");
    for (int i = size-1; i > -1; --i)
    {
      printf("%d:%d ", nodes[i], element[i]);
    }
    printf("]\n");
#endif
    ++(*this);
#ifdef DEBUG_ITER_BEGIN
    printf("Begin: [");
    for (int i = size-1; i > -1; --i)
    {
      printf("%d:%d ", nodes[i], element[i]);
    }
    printf("]\n");
#endif
  }
}


dd_edge::iterator::~iterator()
{
  if (e != 0) {
    free(element);
    free(nodes);
    if (pelement != 0) free(pelement);
    if (pnodes != 0) free (pnodes);
  }
}


dd_edge::iterator::iterator(const iterator& iter)
: e(iter.e), size(iter.size),
  element(0), nodes(0), pelement(0), pnodes(0)
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
  }
}


dd_edge::iterator& dd_edge::iterator::operator=(const iterator& iter)
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
    }
    if (iter.e != 0) {
      e = iter.e;
      size = iter.size;

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
    }
  }
  return *this;
}


void dd_edge::iterator::operator--()
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


void dd_edge::iterator::incrNonRelation()
{
  DCASSERT(e != 0);
  DCASSERT(e->node != 0);

  expert_domain* d = smart_cast<expert_domain*>(e->parent->useDomain());
  expert_forest* f = smart_cast<expert_forest*>(e->parent);

  int currLevel = d->getTopVariable();

  if (nodes[0] == 0) {
    memset(nodes, 0, size * sizeof(int));
    nodes[f->getNodeLevel(e->node)] = e->node;
  }
  else {
    // Start from the bottom level and see if you can find the next
    // valid edge. Move up a level if no more edges exist from this level.

    int lastSkipped = e->node;
    currLevel = d->getVariableAbove(domain::TERMINALS);
    bool found = false;

    while (currLevel != -1)
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
      currLevel = d->getVariableAbove(currLevel);
    }

    if (!found) {
      // no more paths
      nodes[0] = 0;
      return;
    }

    // found new path
    // select the first path starting from level below currLevel
    currLevel = d->getVariableBelow(currLevel);
  }

  // select the first path starting at currLevel
  for ( ; currLevel != domain::TERMINALS;
      currLevel = d->getVariableBelow(currLevel))
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
  assert(nodes[0] != 0);
#ifdef DEBUG_ITER_BEGIN
  printf("nodes[]: [");
  for (int i = size - 1; i > 0; --i)
  {
    printf("%d ", nodes[i]);
  }
  printf("]\n");
#endif
}


void dd_edge::iterator::incrNonIdentRelation()
{
  DCASSERT(e != 0);
  DCASSERT(e->node != 0);

  expert_domain* d = smart_cast<expert_domain*>(e->parent->useDomain());
  expert_forest* f = smart_cast<expert_forest*>(e->parent);

  DCASSERT(f->getReductionRule() == forest::FULLY_REDUCED
      || f->getReductionRule() == forest::QUASI_REDUCED);

  int currLevel = d->getTopVariable();
  bool isCurrLevelPrime = false;

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
    currLevel = d->getVariableAbove(domain::TERMINALS);
    isCurrLevelPrime = true;
    bool found = false;

    while (currLevel != -1)
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
        currLevel = d->getVariableAbove(currLevel);
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
      isCurrLevelPrime = false;
      currLevel = d->getVariableBelow(currLevel);
    } else {
      isCurrLevelPrime = true;
    }
  }

  // select the first path starting at currLevel
  for ( ; currLevel != domain::TERMINALS; )
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
      currLevel = d->getVariableBelow(currLevel);
    } else {
      isCurrLevelPrime = true;
    }
  }
  assert(nodes[0] != 0);
#ifdef DEBUG_ITER_BEGIN
  printf("nodes[]: [");
  for (int i = size - 1; i > 0; --i)
  {
    printf("%d ", nodes[i]);
  }
  printf("]\n");
#endif
}


void dd_edge::iterator::incrRelation()
{
  DCASSERT(e != 0);
  DCASSERT(e->node != 0);

  expert_forest* f = smart_cast<expert_forest*>(e->parent);
  if (f->getReductionRule() != forest::IDENTITY_REDUCED) {
    incrNonIdentRelation();
    return;
  }

  expert_domain* d = smart_cast<expert_domain*>(e->parent->useDomain());
  int currLevel = d->getTopVariable();
  bool isCurrLevelPrime = false;

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
    currLevel = d->getVariableAbove(domain::TERMINALS);
    isCurrLevelPrime = true;
    bool found = false;

    while (currLevel != -1)
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
        currLevel = d->getVariableAbove(currLevel);
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
      isCurrLevelPrime = false;
      currLevel = d->getVariableBelow(currLevel);
    } else {
      isCurrLevelPrime = true;
    }
  }

  // select the first path starting at currLevel
  for ( ; currLevel != domain::TERMINALS; )
  {
    int node = isCurrLevelPrime? pnodes[currLevel]: nodes[currLevel];
    if (node == 0) {
      assert(!isCurrLevelPrime);
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
      assert(index != -1);
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
      currLevel = d->getVariableBelow(currLevel);
    } else {
      isCurrLevelPrime = true;
    }
  }
  assert(nodes[0] != 0);
#ifdef DEBUG_ITER_BEGIN
  printf("nodes[]: [");
  for (int i = size - 1; i > 0; --i)
  {
    printf("%d ", nodes[i]);
  }
  printf("]\n");
#endif
}


void dd_edge::iterator::operator++()
{
  // find next
  // set element to next

  if (e != 0 && e->node != 0) {
    if (e->parent->isForRelations()) {
      incrRelation();
    } else {
      incrNonRelation();
    }
  }
}


bool dd_edge::iterator::operator!=(const iterator& iter) const
{
  DCASSERT((e != iter.e) || (size == iter.size));

  // if terminals are different, return true.
  // else if terminals are 0, return false.
  // else return (element != iter.element)
  return e == 0
          ? true
          : e != iter.e
            ? true
            : (nodes[0] != iter.nodes[0])
              ? true
              : (nodes[0] == 0)
                ? false
                : e->parent->isForRelations()
                  ? (0 != memcmp(element, iter.element, size) ||
                    0 != memcmp(pelement, iter.pelement, size))
                  : 0 != memcmp(element, iter.element, size);
}


bool dd_edge::iterator::operator==(const iterator& iter) const
{
  DCASSERT((e != iter.e) || (size == iter.size));

  // if both terminals are 0 return true.
  // if only one terminal is 0 return false.
  // if neither terminal is 0, return (element == iter.element)

  return !(*this != iter);
}


const int* dd_edge::iterator::getAssignments() const
{
  return e == 0
          ? 0
          : nodes[0] == 0? 0: element;
}


const int* dd_edge::iterator::getPrimedAssignments() const
{
  return e == 0
          ? 0
          : nodes[0] == 0? 0: e->parent->isForRelations()? pelement: 0;
}
