
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


// TODO: implement changes made to mt_forest interface
//
// TODO: HERE: go through every function in mdds.h and mdds.cc


#include "mt.h"
#include "../unique_table.h"

#define ENABLE_CACHE_COUNTING 0
#define ENABLE_IN_COUNTING 0

#ifdef NEW_MT

// ******************************************************************
// *                                                                *
// *                     mt_base_forest methods                     *
// *                                                                *
// ******************************************************************


MEDDLY::mt_base_forest::mt_base_forest(int dsl, domain *d, bool rel,
  range_type t, const policies &p)
: expert_forest(dsl, d, rel, t, MULTI_TERMINAL, p)
{
  unionOp = 0;
}

bool MEDDLY::mt_base_forest::isRedundant(const node_builder &nb) const
{
  if (isQuasiReduced()) return false;
  if (nb.getLevel() < 0 && isIdentityReduced()) return false;
  int common = nb.d(0);
  for (int i=1; i<nb.rawSize(); i++) 
    if (nb.d(i) != common) return false;
  return true;
}

bool MEDDLY::mt_base_forest::isIdentityEdge(const node_builder &nb, int i) const
{
  if (nb.getLevel() > 0) return false;
  if (!isIdentityReduced()) return false;
  if (i<0) return false;
  return nb.d(i) != 0;
}


#else 

// ******************************************************************
// *                                                                *
// *                         public methods                         *
// *                                                                *
// ******************************************************************


MEDDLY::mt_forest::mt_forest(int dsl, domain *d, bool rel, range_type t,
  edge_labeling ev, const policies &p)
: expert_forest(dsl, d, rel, t, ev, p)
{
  counting = false;

  dptrsSize = 0;
  dptrs = 0;
}


MEDDLY::mt_forest::~mt_forest()
{
  if (dptrsSize > 0) {
    free(dptrs);
    dptrsSize = 0;
    dptrs = 0;
  }
}


// ******************************************************************
// *                                                                *
// *            virtual  and overriding default behavior            *
// *                                                                *
// ******************************************************************


void MEDDLY::mt_forest::createEdgeForVar(int vh, bool primedLevel,
    const bool* terms, dd_edge& result)
{
  if (getRangeType() != forest::BOOLEAN) 
      throw error(error::INVALID_OPERATION);
  if (getLevelSize(vh) != 2) 
      throw error(error::INVALID_OPERATION);

  edgeForVarInternal(vh, primedLevel, terms, result);
}


void MEDDLY::mt_forest::createEdgeForVar(int vh, bool primedLevel,
    const int* terms, dd_edge& result)
{
  if (getRangeType() != forest::INTEGER) 
      throw error(error::INVALID_OPERATION);

  edgeForVarInternal(vh, primedLevel, terms, result);
}


void MEDDLY::mt_forest::createEdgeForVar(int vh, bool primedLevel,
    const float* terms, dd_edge& result)
{
  if (getRangeType() != forest::REAL) 
      throw error(error::INVALID_OPERATION);

  edgeForVarInternal(vh, primedLevel, terms, result);
}


bool MEDDLY::mt_forest::isRedundant(const node_builder &nb) const
{
  if (isQuasiReduced()) return false;
  if (nb.getLevel() < 0 && isIdentityReduced()) return false;
  int common = nb.d(0);
  for (int i=1; i<nb.rawSize(); i++) 
    if (nb.d(i) != common) return false;
  return true;
}

bool MEDDLY::mt_forest::isIdentityEdge(const node_builder &nb, int i) const
{
  if (nb.getLevel() > 0) return false;
  if (!isIdentityReduced()) return false;
  if (i<0) return false;
  return nb.d(i) != 0;
}

// ******************************************************************
// *                                                                *
// *                      disorganized methods                      *
// *                                                                *
// ******************************************************************


MEDDLY::node_handle 
MEDDLY::mt_forest::buildLevelNodeHelper(int lh, node_handle* dptrs, int sz)
{
  MEDDLY_DCASSERT(dptrs != 0);
  MEDDLY_DCASSERT(sz > 0);

  const int absLh = (lh < 0)? -lh: lh;

  if (isForRelations()) {
    // build from bottom up
    // for i = 1 to lh-1
    //    for j = 0 to sz
    //      dptrs[j] = node at level i with all downpointers to prev dptrs[j]
    //      do for primed first and then unprimed

    if (!isFullyReduced()) {
      for (int k = 1; k < absLh; ++k)
      {
        for (int j = 0; j < sz; ++j)
        {
          // primed
          insertRedundantNode(-k, dptrs[j]);
          // unprimed
          insertRedundantNode(k, dptrs[j]);
        } // for j
      } // for k

      // Finally, deal with lh level
      // if lh is unprimed, need to create nodes at primed level

      if (lh > 0) {
        // create nodes at level -lh
        for (int j = 0; j < sz; ++j)
        {
          // primed
          insertRedundantNode(-lh, dptrs[j]);
        }
      }
    }
  }
  else if (isQuasiReduced()) {
    MEDDLY_DCASSERT(!isForRelations());
    // build from bottom up
    // for i = 1 to lh-1
    //    for j = 0 to sz
    //      dptrs[j] = node at level i with all downpointers to prev dptrs[j]

    for (int i = 1; i < absLh; ++i)
    {
      for (int j = 0; j < sz; ++j)
      {
        insertRedundantNode(i, dptrs[j]);
      }
    }
  }

  // Now, deal with lh level
  node_builder& nb = useNodeBuilder(lh, sz);
  for (int i=0; i<sz; i++) {
    nb.d(i) = dptrs[i];
  }
  node_handle node = createReducedNode(-1, nb);

  // now build the levels above this node
  if (isForRelations()) {
    if (!isFullyReduced()) {
      // build additional node at lh level if necessary
      if (lh < 0) {
        // build unprimed node at level ABS(lh)
        insertRedundantNode(absLh, node);
      }

      // build primed and unprimed nodes for levels lh+1 to topLevel
      int topHeight = getDomain()->getNumVariables();
      for (int i = absLh + 1; i <= topHeight; ++i)
      {
        // primed
        insertRedundantNode(-i, node);
        // unprimed
        insertRedundantNode(i, node);
      }
    }
    // done building node for Relations
  }
  else if (isQuasiReduced()) {
    MEDDLY_DCASSERT(!isForRelations());
    // build nodes for levels above lh
    int topHeight = getDomain()->getNumVariables();
    for (int i = absLh + 1; i <= topHeight; ++i)
    {
      insertRedundantNode(i, node);
    }
    // done building node for (MT)MDDs
  }

  // MEDDLY_DCASSERT(isReducedNode(node));
  return node;
}

MEDDLY::node_handle* 
MEDDLY::mt_forest::getTerminalNodes(int n, const bool* terms)
{
  MEDDLY_DCASSERT(n == 2);
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    stats.incMemAlloc((n - dptrsSize) * sizeof(node_handle));
    dptrsSize = n;
    dptrs = (node_handle *) realloc(dptrs, dptrsSize * sizeof(node_handle));
    MEDDLY_DCASSERT(NULL != dptrs);
  }

  // fill array with terminal nodes
  if (terms) {
    for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(terms[i]);
  } else {
    dptrs[0] = getTerminalNode(false);
    dptrs[1] = getTerminalNode(true);
  }
  return dptrs;
}


MEDDLY::node_handle* 
MEDDLY::mt_forest::getTerminalNodes(int n, const int* terms)
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    stats.incMemAlloc((n - dptrsSize) * sizeof(node_handle));
    dptrsSize = n;
    dptrs = (node_handle *) realloc(dptrs, dptrsSize * sizeof(node_handle));
    MEDDLY_DCASSERT(NULL != dptrs);
  }

  // fill array with terminal nodes
  if (terms) {
    for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(terms[i]);
  } else {
    for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(i);
  }
  return dptrs;
}


MEDDLY::node_handle* 
MEDDLY::mt_forest::getTerminalNodes(int n, const float* terms)
{
  MEDDLY_DCASSERT(getRangeType() == forest::REAL);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    stats.incMemAlloc((n - dptrsSize) * sizeof(node_handle));
    dptrsSize = n;
    dptrs = (node_handle *) realloc(dptrs, dptrsSize * sizeof(node_handle));
    MEDDLY_DCASSERT(NULL != dptrs);
  }
  // fill array with terminal nodes
  if (terms) {
    for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(terms[i]);
  } else {
    for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(float(i));
  }
  return dptrs;
}

#endif

