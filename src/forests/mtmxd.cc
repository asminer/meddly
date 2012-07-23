
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


#include "mtmxd.h"
#include "../unique_table.h"

MEDDLY::mtmxd_forest::mtmxd_forest(int dsl, domain *d,
    bool relation, forest::range_type t,
    forest::edge_labeling e, const policies &p)
: mt_forest(dsl, d, relation, t, e, p)
{
  unpList = 0;
  pList = 0;
  tList = 0;
  listSize = 0;
  count = 0;
  slot = 0;
  countSize = 0;
}


MEDDLY::mtmxd_forest::~mtmxd_forest()
{
  if (unpList) free(unpList);
  if (pList) free(pList);
  if (tList) free(tList);
  if (count) free(count);
  if (slot) free(slot);
}


void MEDDLY::mtmxd_forest::expandCountAndSlotArrays(int size)
{
  if (size <= countSize) return;

  int newCountSize = countSize == 0? 8: countSize;
  while (newCountSize < size) { newCountSize *= 2; }

  count = (int*) realloc(count, newCountSize * sizeof(int));
  slot = (int*) realloc(slot, newCountSize * sizeof(int));
  memset(count + countSize, 0, (newCountSize - countSize) * sizeof(int));
  memset(slot + countSize, 0, (newCountSize - countSize) * sizeof(int));
  countSize = newCountSize;
}

int MEDDLY::mtmxd_forest::createNode(int k, int index, int dptr)
{
  MEDDLY_DCASSERT(index >= 0 && index < getLevelSize(k) && isValidNodeIndex(dptr));

  if (dptr == 0) return 0;

  node_builder& nb = useSparseBuilder(k, 1);
  nb.d(0) = dptr;
  nb.i(0) = index;
  return createReducedNode(-1, nb);
}

int MEDDLY::mtmxd_forest::createNode(int k, int index1, int index2, int dptr)
{
  MEDDLY_DCASSERT((index1 >= 0 && index2 >= 0) ||
      (index1 >= -1 && index2 >= -1) ||
      (index1 >= -2 && index2 >= -2 && index1 == index2));

  int result = 0;

  if (isIdentityReduced()) {

    if (index1 == -2) {
      // "don't change"
      result = dptr;
    }
    else {
      int p = 0;
      if (index2 == -1) {
        // represents "don't care"
        insertRedundantNode(-k, dptr);
        p = dptr;
      } else {
        p = createNode(-k, index2, dptr);
      }
      if (index1 == -1) {
        // represents "don't care"
        insertRedundantNode(k, p);
        result = p;
      } else {
        result = createNode(k, index1, p);
      }
    }

  }
  else if (index1 == -2) {

    // "don't change"
    MEDDLY_DCASSERT(isQuasiReduced() || isFullyReduced());
    int sz = getLevelSize(k);
    node_builder &nb = useNodeBuilder(k, sz);
    for (int i=0; i<sz; i++) {
      nb.d(i) = createNode(-k, i, dptr);
    }
    result = createReducedNode(-1, nb);
  }
  else if (isQuasiReduced()) {

    int p = 0;
    if (index2 == -1) {
      // represents "don't care"
      insertRedundantNode(-k, dptr);
      p = dptr;
    } else {
      p = createNode(-k, index2, dptr);
    }
    if (index1 == -1) {
      // represents "don't care"
      insertRedundantNode(k, p);
      result = p;
    } else {
      result = createNode(k, index1, p);
    }
  }
  else {

    // deal with "don't care" for primed level
    int p = (index2 == -1) ? dptr : createNode(-k, index2, dptr);
    // deal with "don't care" for unprimed level
    result = (index1 == -1) ? p : createNode(k, index1, p);

  }

  return result;
}


int MEDDLY::mtmxd_forest::createNode(const int* v, const int* vp, int term,
    int startAtHeight, bool primedLevel)
{
  // construct the edge bottom-up
  for (int k=1; k<startAtHeight; k++) {
    term = createNode(k, v[k], vp[k], term);
  } // for k

  // deal with height == startAtHeight
  // handle primed level first
  if (primedLevel) {
    // only primed level to be handled at this height
    term = createNode(-startAtHeight, vp[startAtHeight], term);
  } else {
    // both primed and unprimed levels to be handled at this height
    term = createNode(startAtHeight, v[startAtHeight], vp[startAtHeight], term);
  }

  return term;
}


void MEDDLY::mtmxd_forest::createEdge(const int* v, const int* vp, int term,
    int startAtHeight, bool primedLevel, dd_edge& e)
{
  term = createNode(v, vp, term, startAtHeight, primedLevel);
  e.set(term, 0, getNodeLevel(term));
}


void MEDDLY::mtmxd_forest::createEdge(const int* v, const int* vp, int term,
    dd_edge& e)
{
  createEdge(v, vp, term, getExpertDomain()->getNumVariables(), false, e);
}


int MEDDLY::mtmxd_forest::createEdgeTo(int dptr)
{
  MEDDLY_DCASSERT(isTerminalNode(dptr));
  if (dptr == 0) return 0;

  if (isFullyReduced()) return linkNode(dptr);

  // construct the edge bottom-up
  int curr = dptr;
  for (int i=1; i<=getExpertDomain()->getNumVariables(); i++) {
    curr = createNode(i, -1, -1, curr);
  }
  return curr;
}


int MEDDLY::mtmxd_forest::getTerminalNodeForEdge(int n, const int* vlist,
    const int* vplist) const
{
  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contain a single element.

  if (isIdentityReduced()) {

    while (!isTerminalNode(n)) {
      int nLevel = getNodeLevel(n);
      if (nLevel < 0) {
        // Primed Node
        int next = getDownPtr(n, vplist[-nLevel]);
        MEDDLY_DCASSERT(isTerminalNode(next) || isUnprimedNode(next));
        int currHeight = getNodeHeight(n) - 1;
        int nextHeight = getNodeHeight(next);
        if (nextHeight < currHeight) {
          // skipped levels
          for ( ; nextHeight != currHeight; --currHeight)
          {
            if (vlist[currHeight] != vplist[currHeight]) {
              next = 0;
              break;
            }
          }
        }
        n = next;
      }
      else {
        // Unprimed Node
        MEDDLY_DCASSERT(getDownPtr(n, vlist[nLevel]) == 0
            || -nLevel == getNodeLevel(getDownPtr(n, vlist[nLevel])));
        n = getDownPtr(n, vlist[nLevel]);
      }
    }

  }
  else {

    while (!isTerminalNode(n)) {
      n = isPrimedNode(n)
        ? getDownPtr(n, vplist[-getNodeLevel(n)])
        : getDownPtr(n, vlist[getNodeLevel(n)]);
    }

  }
  return n;
}

