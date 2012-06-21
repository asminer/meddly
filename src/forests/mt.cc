
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

// #include <map>    // for getCardinality()
// #include <queue>  // for showNodeGraph
// #include <vector> // for showNodeGraph
// #include <set>
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#include <limits.h>

#define DONT_USE_FULL 0

#define ENABLE_GC 1
// #define DEBUG_GC
#define ENABLE_CACHE_COUNTING 0
#define ENABLE_IN_COUNTING 0

#define DEBUG_DELETE_NM 0

// #define DEBUG_CARD

//#define TRACK_DELETIONS

const int add_size = 1024;
// const int l_add_size = 24;

// ******************************************************************
// *                                                                *
// *                         public methods                         *
// *                                                                *
// ******************************************************************


MEDDLY::mt_forest::mt_forest(int dsl, domain *d, bool rel, range_type t,
  edge_labeling ev, const policies &p, int dataHeaderSize)
: expert_forest(dsl, d, rel, t, ev, p)
{
  this->dataHeaderSize = dataHeaderSize;

  unique = new mdd_hash_table<mt_forest> (this);

  delete_terminal_nodes = false;

  counting = false;

  dptrsSize = 0;
  dptrs = 0;

  nodeA = 0;
  nodeB = 0;
}

// ******************************************************************
// *                                                                *
// *                       protected  methods                       *
// *                                                                *
// ******************************************************************

void MEDDLY::mt_forest::dumpUniqueTable(FILE* s) const
{
  fprintf(s, "Uniqueness table:\n");
  unique->show(s);
}

// ******************************************************************
// *                                                                *
// *                      disorganized methods                      *
// *                                                                *
// ******************************************************************

#if 1


int MEDDLY::mt_forest::buildLevelNodeHelper(int lh, int* dptrs, int sz)
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
      for (int i = 1; i < absLh; ++i)
      {
        for (int j = 0; j < sz; ++j)
        {
          // primed
          int temp = createTempNodeMaxSize(-i, false);
          setAllDownPtrsWoUnlink(temp, dptrs[j]);
          unlinkNode(dptrs[j]);
          dptrs[j] = reduceNode(temp);
          // unprimed
          temp = createTempNodeMaxSize(i, false);
          setAllDownPtrsWoUnlink(temp, dptrs[j]);
          unlinkNode(dptrs[j]);
          dptrs[j] = reduceNode(temp);
        }
      }

      // Finally, deal with lh level
      // if lh is unprimed, need to create nodes at primed level

      if (lh > 0) {
        // create nodes at level -lh
        for (int j = 0; j < sz; ++j)
        {
          // primed
          int temp = createTempNodeMaxSize(-lh, false);
          setAllDownPtrsWoUnlink(temp, dptrs[j]);
          unlinkNode(dptrs[j]);
          dptrs[j] = reduceNode(temp);
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
        int temp = createTempNodeMaxSize(i, false);
        setAllDownPtrsWoUnlink(temp, dptrs[j]);
        unlinkNode(dptrs[j]);
        dptrs[j] = reduceNode(temp);
      }
    }
  }

  // Now, deal with lh level
  int node = createTempNode(lh, sz, false);
  int* curr = getFullNodeDownPtrs(node);
  int* stop = curr + getFullNodeSize(node);
  // no need to link/unlink nodes since we pass the link
  // from dptrs[] to curr[]
  for ( ; curr != stop; ++curr, ++dptrs) { *curr = *dptrs; }
  node = reduceNode(node);

  // now build the levels above this node
  if (isForRelations()) {
    if (!isFullyReduced()) {
      // build additional node at lh level if necessary
      if (lh < 0) {
        // build unprimed node at level ABS(lh)
        int temp = createTempNodeMaxSize(absLh, false);
        setAllDownPtrsWoUnlink(temp, node);
        unlinkNode(node);
        node = reduceNode(temp);
      }

      // build primed and unprimed nodes for levels lh+1 to topLevel
      int topHeight = getDomain()->getNumVariables();
      for (int i = absLh + 1; i <= topHeight; ++i)
      {
        // primed
        int temp = createTempNodeMaxSize(-i, false);
        setAllDownPtrsWoUnlink(temp, node);
        unlinkNode(node);
        node = reduceNode(temp);
        // unprimed
        temp = createTempNodeMaxSize(i, false);
        setAllDownPtrsWoUnlink(temp, node);
        unlinkNode(node);
        node = reduceNode(temp);
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
      int temp = createTempNodeMaxSize(i, false);
      setAllDownPtrsWoUnlink(temp, node);
      unlinkNode(node);
      node = reduceNode(temp);
    }
    // done building node for (MT)MDDs
  }

  MEDDLY_DCASSERT(isReducedNode(node));
  return node;
}


void MEDDLY::mt_forest::buildLevelNode(int k, int* dptrs, int sz)
{
  MEDDLY_DCASSERT(getLevelNode(k) == 0);
  MEDDLY_DCASSERT(dptrs != 0);
  MEDDLY_DCASSERT(sz > 0);

  levels[k].levelNode = buildLevelNodeHelper(k, dptrs, sz);
  MEDDLY_DCASSERT(getLevelNode(k) != 0 && isReducedNode(getLevelNode(k)) ||
      getLevelNode(k) == 0 && sz == 1 && dptrs[0] == 0);
}


void MEDDLY::mt_forest::clearLevelNode(int k)
{
  unlinkNode(levels[k].levelNode);
  levels[k].levelNode = 0;
}


void MEDDLY::mt_forest::clearLevelNodes()
{
  // for each level, unlink the level node
  for (int i=getMinLevelIndex(); i<=getNumVariables(); i++) {
    clearLevelNode(i);
  }
}


int* MEDDLY::mt_forest::getTerminalNodes(int n, bool* terms)
{
  MEDDLY_DCASSERT(n == 2);
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    stats.incMemAlloc((n - dptrsSize) * sizeof(int));
    dptrsSize = n;
    dptrs = (int *) realloc(dptrs, dptrsSize * sizeof(int));
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


int* MEDDLY::mt_forest::getTerminalNodes(int n, int* terms)
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    stats.incMemAlloc((n - dptrsSize) * sizeof(int));
    dptrsSize = n;
    dptrs = (int *) realloc(dptrs, dptrsSize * sizeof(int));
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


int* MEDDLY::mt_forest::getTerminalNodes(int n, float* terms)
{
  MEDDLY_DCASSERT(getRangeType() == forest::REAL);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    stats.incMemAlloc((n - dptrsSize) * sizeof(int));
    dptrsSize = n;
    dptrs = (int *) realloc(dptrs, dptrsSize * sizeof(int));
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


void MEDDLY::mt_forest::createEdgeForVar(int vh, bool primedLevel,
    bool* terms, dd_edge& result)
{
  if (!isValidVariable(vh)) 
    throw error(error::INVALID_VARIABLE);
  if (result.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (getRangeType() != forest::BOOLEAN) 
    throw error(error::INVALID_OPERATION);
  if (getLevelSize(vh) != 2) 
    throw error(error::INVALID_OPERATION);

  int k = primedLevel? -vh: vh;
  MEDDLY_DCASSERT(isValidLevel(k));

  if (!isForRelations() && primedLevel) 
    throw error(error::INVALID_ASSIGNMENT);
  if (getEdgeLabeling() != MULTI_TERMINAL)
    throw error(error::INVALID_OPERATION);
  int *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
  int node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

  result.set(node, 0, getNodeLevel(node));
}


void MEDDLY::mt_forest::createSubMatrix(const dd_edge& rows,
    const dd_edge& cols, const dd_edge& a, dd_edge& result)
{
  throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::mt_forest::createSubMatrix(const bool* const* vlist,
    const bool* const* vplist, const dd_edge a, dd_edge& b)
{
  if (a.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (b.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (!isMxd()) 
    throw error(error::INVALID_OPERATION);

  // Build Mask: go bottom up
  // When mask for level i is done, create node at level i+1 (higher)
  // such that all downpointers with vlist[i+1]==1 point to mask.
  int mask = getTerminalNode(true);
  int nVars = getExpertDomain()->getNumVariables();
  for (int level = 1; level <= nVars; level++)
  {
    // create node at prime level
    int nodeSize = getExpertDomain()->getVariableBound(level, true);
    int node = createTempNode(-level, nodeSize, false);
    for (int i = 0; i < nodeSize; i++)
    {
      setDownPtrWoUnlink(node, i, (vplist[level][i]? mask: 0));
    }
    unlinkNode(mask);
    mask = reduceNode(node);

    // create node at unprime level
    nodeSize = getExpertDomain()->getVariableBound(level, false);
    node = createTempNode(level, nodeSize, false);
    for (int i = 0; i < nodeSize; i++)
    {
      setDownPtrWoUnlink(node, i, (vlist[level][i]? mask: 0));
    }
    unlinkNode(mask);
    mask = reduceNode(node);
  }

  b.set(mask, 0, getNodeLevel(mask));
#if 0
  b.show(stdout, 3);
#endif
  b *= a;
}


void MEDDLY::mt_forest::getElement(const dd_edge& a,
    int index, int* e)
{
  throw error(error::INVALID_OPERATION);
}


void MEDDLY::mt_forest::createEdgeForVar(int vh, bool primedLevel,
    int* terms, dd_edge& result)
{
  if (!isValidVariable(vh)) 
    throw error(error::INVALID_VARIABLE);
  if (result.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (getRangeType() != forest::INTEGER) 
    throw error(error::INVALID_OPERATION);

  int k = primedLevel? -vh: vh;
  MEDDLY_DCASSERT(isValidLevel(k));

  if (!isForRelations() && primedLevel) 
    throw error(error::INVALID_ASSIGNMENT);
  if (getEdgeLabeling() != MULTI_TERMINAL)
    throw error(error::INVALID_OPERATION);
  int *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
  int node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

  result.set(node, 0, getNodeLevel(node));
}


void MEDDLY::mt_forest::createEdgeForVar(int vh, bool primedLevel,
    float* terms, dd_edge& result)
{
  if (!isValidVariable(vh)) 
    throw error(error::INVALID_VARIABLE);
  if (result.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (getRangeType() != forest::REAL) 
    throw error(error::INVALID_OPERATION);

  int k = primedLevel? -vh: vh;
  MEDDLY_DCASSERT(isValidLevel(k));

  if (!isForRelations() && primedLevel) 
    throw error(error::INVALID_ASSIGNMENT);
  if (getEdgeLabeling() != MULTI_TERMINAL)
    throw error(error::INVALID_OPERATION);
  int *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
  int node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

  result.set(node, 0, getNodeLevel(node));
}


#endif


void MEDDLY::mt_forest::clearAllNodes()
{
  int level = getExpertDomain()->getNumVariables();
  while (level > 0 && stats.active_nodes > 0)
  {
    // find all nodes at curr level and make them into orphans
    for (int i = 1; i < a_last; i++)
    {
      if (isActiveNode(i) && getNodeLevel(i) == level && getInCount(i) > 0) {
        getInCount(i) = 1;
        unlinkNode(i);
      }
    }

    if (stats.active_nodes > 0 && isForRelations()) {
      level = -level;
      // find all nodes at curr level and make them into orphans
      for (int i = 1; i < a_last; i++)
      {
        if (isActiveNode(i) && getNodeLevel(i) == level && getInCount(i) > 0) {
          getInCount(i) = 1;
          unlinkNode(i);
        }
      }
    }

    level--;
  }
}


MEDDLY::mt_forest::~mt_forest()
{
  if (nodeA) delete nodeA;
  if (nodeB) delete nodeB;

  // setting delete_terminal_nodes will ensure that all nodes in compute
  // tables are removed during garbage collection.
  delete_terminal_nodes = true;
  // unlink all nodes that are being stored at the respective levels
  clearLevelNodes();

  /*

  // All these things should now be unnecessary...

  // remove all disconnected nodes
  gc(true);

  if (stats.active_nodes != 0) {
    printf("MEDDLY ALERT: In %s, active_nodes > 0.\n", __func__);
    printf("This usually means that your application is still referring\n");
    printf("to one or more MDD nodes. ");
    printf("Fixing this may benefit your application.\n");
#ifdef DEBUG_GC
    printf("%p: active %ld, zombie %ld, orphan %ld\n",
        this, active_nodes, zombie_nodes, orphan_nodes);
#endif
    clearAllNodes();
    gc(true);
  }
  */

  if (dptrsSize > 0) {
    free(dptrs);
    dptrsSize = 0;
    dptrs = 0;
  }

  delete unique;
#if DEBUG_DELETE_NM
  printf("Deleted unique table\n");
  fflush(stdout);
#endif
}

// *********************************************************************
// TODO: test this out
int MEDDLY::mt_forest::buildQuasiReducedNodeAtLevel(int k, int p)
{
  MEDDLY_DCASSERT(isQuasiReduced());
  int curr = p;
  int p_level = getNodeLevel(p);
  for (int i = p_level + 1; i <= k; i++)
  {
    int n = createTempNodeMaxSize(i);
    setAllDownPtrs(n, curr);
    curr = reduceNode(n);
  }
  return curr;
}
// *********************************************************************

int MEDDLY::mt_forest::getMddLevelMaxBound(int k) const
{
  // Go through each node in this level and check its size - in terms of 
  // # of downpointers. If it has more downpointers that the current max_bound
  // update max_bound. After going throug all nodes, return max_bound.
  if (k<getMinLevelIndex()) return 0;
  if (k>getNumVariables()) return 0;
  if (0 == levels[k].data) return 0;
  int* data = levels[k].data;

  int max_bound = 0;
  for (int a=1; a<levels[k].last; ) {
    if (data[a] < 0) {
      // hole; skip ahead
      // hole size is stored as a negative number so subtract
      a -= data[a]; 
    } else {
      // proper node - full or sparse node
      if (data[a+2] > 0) {
        // full
        if (max_bound < data[a+2]) max_bound = data[a+2];
        a += 4 + data[a+2];
      } else {
        // sparse
        int max_index = data[a + 3 - data[a+2] - 1];
        if (max_bound < (max_index + 1)) max_bound = max_index + 1;
        a += 4 - 2 * data[a+2];
      }
    }
  }
  return max_bound;
}

int MEDDLY::mt_forest::getMxdLevelMaxBound(int k) const
{
  return MAX(getMddLevelMaxBound(k), getMddLevelMaxBound(-k));
}

int MEDDLY::mt_forest::getLevelMaxBound(int k) const
{
  return isForRelations()?
            getMxdLevelMaxBound(k):
            getMddLevelMaxBound(k);
}

void MEDDLY::mt_forest::showAll() const
{
  dumpInternal(stderr);
}

void MEDDLY::mt_forest::showNodeGraph(FILE *s, int p) const
{
  int* list = markNodesInSubgraph(p, true);
  if (0==list) return;

  // Print by levels
  for (int k = getNumVariables(); k; )
  {
    bool printed = false;
    for (int i=0; list[i]; i++) {
      if (getNodeLevel(list[i]) != k) continue;

      if (!printed) {
        const variable* v = getDomain()->getVar(ABS(k));
        if (v->getName()) {
          fprintf(s, "Level: %s \n", v->getName());
        } else {
          fprintf(s, "Level: %d \n", ABS(k));
        }
        printed = true;
      }

      fprintf(s, "  ");
      showNode(s, list[i]);
      fprintf(s, "\n");
    }
    
    // next level
    k *= -1;
    if (k>0) k--;
  } // for k

  free(list);
}


int ifTermGetInt(const MEDDLY::mt_forest *nm, int node)
{
  return nm->isTerminalNode(node) ? nm->getInteger(node) : node;
}

void MEDDLY::mt_forest::showNode(FILE *s, int p, int verbose) const
{
  if (isTerminalNode(p)) {
    fprintf(s, "(terminal)");
    return;
  }
  if (isDeletedNode(p)) {
    fprintf(s, "DELETED");
    return;
  }
  if (isZombieNode(p)) {
    fprintf(s, "Zombie cc: %d", -address[p].cache_count);
    return;
  }
  int a = getNodeOffset(p);
  int l = getNodeLevel(p);
#if 0
  int p_width = digits(a_last);
  int l_width = digits(l_size);
#endif
  int* data = levels[l].data;
  if (verbose) {
    const variable* v = getDomain()->getVar(ABS(l));
    if (v->getName()) {
      fprintf(s, " level: %s", v->getName());
    } else {
      fprintf(s, " level: %d", ABS(l));
    }
    if (getNodeLevel(p) < 0)
      fprintf(s, "'");
    else
      fprintf(s, " ");
    fprintf(s, " in: %d", data[a]);
    fprintf(s, " cc: %d", address[p].cache_count);
  } else {
    fprintf(s, "%snode: %d", (isReducedNode(p)? " ": "+"), p);
  }
  if (isSparseNode(p)) {
    // sparse
    if (verbose)
      fprintf(s, " nnz : %d", getSparseNodeSize(p));
    fprintf(s, " down: (");
    for (int z=0; z<getSparseNodeSize(p); z++) {
      if (z) fprintf(s, ", ");
      if (isEVPlus()) {
        int e = 0;
        getSparseNodeEdgeValue(p, z, e);
        if (e == INF) {
          fprintf(s, "%d:<INF,%d>",
              getSparseNodeIndex(p, z),
              getSparseNodeDownPtr(p, z));
        } else {
          fprintf(s, "%d:<%d,%d>",
              getSparseNodeIndex(p, z),
              e,
              getSparseNodeDownPtr(p, z));
        }
      } else if (isEVTimes()) {
        float e = 0;
        getSparseNodeEdgeValue(p, z, e);
        fprintf(s, "%d:<%f,%d>",
            getSparseNodeIndex(p, z),
            e,
            getSparseNodeDownPtr(p, z));
      } else {
        if (isTerminalNode(getSparseNodeDownPtr(p, z))) {
          fprintf(s, "%d:", getSparseNodeIndex(p, z));
          if (getRangeType() == forest::REAL) {
            fprintf(s, "%f", getReal(getSparseNodeDownPtr(p, z)));
          } else if (getRangeType() == forest::INTEGER) {
            fprintf(s, "%d", getInteger(getSparseNodeDownPtr(p, z)));
          } else {
            MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);
            fprintf(s, "%s",
                (getBoolean(getSparseNodeDownPtr(p, z))? "T": "F"));
          }
          fprintf(s, "*");
        } else {
          fprintf(s, "%d:%d",
              getSparseNodeIndex(p, z),
              getSparseNodeDownPtr(p, z));
        }
      }
    }
    fprintf(s, ")");
  } else {
    int size = *(getNodeAddress(p) + 2);
    // fprintf(s, " size: %d down: [", p_width, getFullNodeSize(p));
    if (verbose) fprintf(s, " size: %d", size);
    fprintf(s, " down: [");
    for (int i=0; i<getFullNodeSize(p); i++) {
      if (i) fprintf(s, "|");
      if (isEVPlus()) {
        int e = 0;
        getFullNodeEdgeValue(p, i, e);
        if  (e == INF) {
          fprintf(s, "<INF,%d>",
              getFullNodeDownPtr(p, i));
        } else {
          fprintf(s, "<%d,%d>", e,
              getFullNodeDownPtr(p, i));
        }
      } else if (isEVTimes()) {
        float e = 0;
        getFullNodeEdgeValue(p, i, e);
        fprintf(s, "<%f,%d>", e,
            getFullNodeDownPtr(p, i));
      } else {
        if (isTerminalNode(getFullNodeDownPtr(p, i))) {
          if (getRangeType() == forest::REAL) {
            fprintf(s, "%f", getReal(getFullNodeDownPtr(p, i)));
          } else if (getRangeType() == forest::INTEGER) {
            fprintf(s, "%d", getInteger(getFullNodeDownPtr(p, i)));
          } else {
            fprintf(s, "%s",
                (getBoolean(getFullNodeDownPtr(p, i))? "T": "F"));
          }
          fprintf(s, "*");
        } else {
          fprintf(s, "%d", getFullNodeDownPtr(p, i));
        }
      }
    }
    fprintf(s, "]");
  }
}

void MEDDLY::mt_forest::showNode(int p) const
{
  MEDDLY_DCASSERT(isEVPlus() || isEVTimes());
  if (isTerminalNode(p)) {
    fprintf(stderr, "(terminal)");
    return;
  }
  if (isDeletedNode(p)) {
    fprintf(stderr, "DELETED");
    return;
  }
  int a = getNodeOffset(p);
  int l = getNodeLevel(p);
  int p_width = digits(a_last);
  int l_width = digits(getNumVariables());
  int* data = levels[l].data;
  fprintf(stderr, "node: %*d level: %*d", p_width, p, l_width, ABS(l));
  if (getNodeLevel(p) < 0)
    fprintf(stderr, "'");
  else
    fprintf(stderr, " ");
  fprintf(stderr, " in: %*d", p_width, data[a]);
  if (isSparseNode(p)) {
    // sparse
    fprintf(stderr, " nnz : %*d down: (", p_width, getSparseNodeSize(p));
    for (int z=0; z<getSparseNodeSize(p); z++) {
      if (z) fprintf(stderr, ", ");
      fprintf(stderr, "%*d:%*d",
          p_width, getSparseNodeIndex(p, z),
          p_width, getSparseNodeDownPtr(p, z));
    }
    fprintf(stderr, ")");
  } else {
    int size = *(getNodeAddress(p) + 2);
    // fprintf(stderr, " size: %*d down: [", p_width, getFullNodeSize(p));
    fprintf(stderr, " size: %*d down: [", p_width, size);
    for (int i=0; i<getFullNodeSize(p); i++) {
      if (i) fprintf(stderr, "|");
      fprintf(stderr, "%*d", p_width, getFullNodeDownPtr(p, i));
    }
    fprintf(stderr, "]");
  }
}

void MEDDLY::mt_forest::compactAllLevels()
{
  for (int i=getMinLevelIndex(); i<=getNumVariables(); i++) {
    levels[i].compactLevel = true;
    levels[i].compact(address);
  }
}


// ------------------------------------------------------------------
//  For uniqueness table
// ------------------------------------------------------------------

/*
 * Bob Jenkin's Hash
 * Free to use for educational or commerical purposes
 * http://burtleburtle.net/bob/hash/doobs.html
 */
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
#define mix(a,b,c) \
  { \
      a -= c;  a ^= rot(c, 4);  c += b; \
      b -= a;  b ^= rot(a, 6);  a += c; \
      c -= b;  c ^= rot(b, 8);  b += a; \
      a -= c;  a ^= rot(c,16);  c += b; \
      b -= a;  b ^= rot(a,19);  a += c; \
      c -= b;  c ^= rot(b, 4);  b += a; \
  }
#define final(a,b,c) \
  { \
      c ^= b; c -= rot(b,14); \
      a ^= c; a -= rot(c,11); \
      b ^= a; b -= rot(a,25); \
      c ^= b; c -= rot(b,16); \
      a ^= c; a -= rot(c,4);  \
      b ^= a; b -= rot(a,14); \
      c ^= b; c -= rot(b,24); \
  }

#if 1

unsigned MEDDLY::mt_forest::hash(int h) const 
{
  int* k = getNodeAddress(h);
  int length = k[2];
  MEDDLY_DCASSERT(length != 0);
  k += 3;

  uint32_t a[] = { uint32_t(getNodeLevel(h)), 0, 0xdeadbeef };

  if (length > 0) {
    // Full node
    int* ptr = k;
    int* stop = k + length;
    unsigned nnz = 1;
    for (int i = 0; ptr != stop; ) {
      if (*ptr == 0) { ++ptr; ++i; continue; }
      a[nnz++] += i++;
      if (nnz == 3) {
        mix(a[0], a[1], a[2]);
        nnz = 0;
      }
      a[nnz++] += *ptr++;
      if (nnz == 3) {
        mix(a[0], a[1], a[2]);
        nnz = 0;
      }
    }
  }
  else {
    // Sparse node
    int* indexes = k;
    int* ptr = k - length;
    int* stop = ptr - length;
    unsigned nnz = 1;
    for ( ; ptr != stop; ) {
      a[nnz++] += *indexes++;
      if (nnz == 3) {
        mix(a[0], a[1], a[2]);
        nnz = 0;
      }
      a[nnz++] += *ptr++;
      if (nnz == 3) {
        mix(a[0], a[1], a[2]);
        nnz = 0;
      }
    }
  }

  //if (nnz > 0) {
    final(a[0], a[1], a[2]);
  //}

  // report the result
  return a[2];
}

#else

unsigned MEDDLY::mt_forest::hash(int h) const 
{
  int* k = getNodeAddress(h);
  int length = k[2];
  MEDDLY_DCASSERT(length != 0);

  uint32_t a = 0;
  uint32_t b = 0;
  uint32_t c = 0;

  if (length > 0) {
    // full node
    // move forward 3 slots since thats the start of the pointers
    k += 3;

    // find nnz
    unsigned nnz = 0;
    {
      int *ptr = k;
      int *stop = k + length;
      while (ptr != stop)
      {
        if (*ptr) nnz++;
        ++ptr;
      }
    }

    // Set up the internal state
#if 1
    a = 0xdeadbeef;
    b = uint32_t(nnz)<<16;
    c = uint32_t(getNodeLevel(h));
#endif

    // handle most of the key
    while (nnz > 3)
    {
      while (*k == 0) { ++k; }
      a += *k;
      ++k;
      while (*k == 0) { ++k; }
      b += *k;
      ++k;
      while (*k == 0) { ++k; }
      c += *k;
      ++k;
      mix(a,b,c);
      nnz -= 3;
    }

    // handle the last 3 uint32_t's
    switch(nnz)
    { 
      case 0:
        // nothing left to add
        final(a,b,c);
        break;
      case 1:
        while (*k == 0) { ++k; }
        a += *k;
        final(a,b,c);
        break;
      case 2:
        while (*k == 0) { ++k; }
        a += *k;
        ++k;
        while (*k == 0) { ++k; }
        b += *k;
        final(a,b,c);
        break;
      case 3:
        while (*k == 0) { ++k; }
        a += *k;
        ++k;
        while (*k == 0) { ++k; }
        b += *k;
        ++k;
        while (*k == 0) { ++k; }
        c += *k;
        final(a,b,c);
        break;
    }
  }
  else {
    // sparse node
    // -ve length indicates sparse nodes; make it +ve
    length = -length;
    // move forward (3 + length) slots since thats the start of the pointers
    k += 3 + length;
    // Set up the internal state;
    // a = b = c = 0xdeadbeef;
    // a = b = c = 0xdeadbeef + uint32_t(length)<<2;
#if 1
    a = 0xdeadbeef;
    b = uint32_t(length)<<16;
    c = uint32_t(getNodeLevel(h));
#endif
    // handle most of the key
    while (length > 3)
    {
      a += k[0];
      b += k[1];
      c += k[2];
      mix(a,b,c);
      length -= 3;
      k += 3;
    }

    // handle the last 3 uint32_t's
    switch(length)
    { 
      // all the case statements fall through
      case 3: c += k[2];
      case 2: b += k[1];
      case 1: a += k[0];
              //final(a,b,c);
      case 0: // nothing left to add
              final(a,b,c);
              break;
    }
  }

  // report the result
  return c;
}

#endif


bool MEDDLY::mt_forest::equalsFF(int h1, int h2) const
{
  MEDDLY_DCASSERT(isFullNode(h1));
  MEDDLY_DCASSERT(isFullNode(h2));

  int *ptr1 = getNodeAddress(h1) + 2;
  int *ptr2 = getNodeAddress(h2) + 2;
  int sz1 = *ptr1++;
  int sz2 = *ptr2++;

  int* h1Stop = ptr1 + sz1;
  int* h2Stop = ptr2 + sz2;

  if (sz1 > sz2) {
    while (ptr2 != h2Stop) { if (*ptr1++ != *ptr2++) return false; }
    while (ptr1 != h1Stop) { if (*ptr1++ != 0) return false; }
  }
  else {
    MEDDLY_DCASSERT(sz1 <= sz2);
    while (ptr1 != h1Stop) { if (*ptr1++ != *ptr2++) return false; }
    while (ptr2 != h2Stop) { if (*ptr2++ != 0) return false; }
  }

  if (isMultiTerminal()) return true;

  // Check edge-values
  MEDDLY_DCASSERT(ptr1 == h1Stop);
  MEDDLY_DCASSERT(ptr2 == h2Stop);

  h1Stop += MIN(sz1, sz2);
  if (isEVPlus()) {
    while (ptr1 != h1Stop) { if (*ptr1++ != *ptr2++) return false; }
  } else {
    MEDDLY_DCASSERT(isEVTimes());
    while (ptr1 != h1Stop) {
      if (!isAlmostEqual(*ptr1++, *ptr2++)) return false;
    }
  }
  return true;
}


bool MEDDLY::mt_forest::equalsSS(int h1, int h2) const
{
  MEDDLY_DCASSERT(isSparseNode(h1));
  MEDDLY_DCASSERT(isSparseNode(h2));

  int *ptr1 = getNodeAddress(h1) + 2;
  int *ptr2 = getNodeAddress(h2) + 2;
  int sz1 = -(*ptr1++);
  int sz2 = -(*ptr2++);

  if (sz1 != sz2) return false;

  int* h1Stop = ptr1 + sz1 + sz1;
  while (ptr1 != h1Stop) { if (*ptr1++ != *ptr2++) return false; }

  if (isMultiTerminal()) return true;

  // Check edge-values
  MEDDLY_DCASSERT(ptr1 == h1Stop);
  MEDDLY_DCASSERT(ptr2 == (getNodeAddress(h2) + 3 + sz1 + sz1));

  h1Stop += sz1;
  if (isEVPlus()) {
    while (ptr1 != h1Stop) { if (*ptr1++ != *ptr2++) return false; }
  } else {
    MEDDLY_DCASSERT(isEVTimes());
    while (ptr1 != h1Stop) {
      if (!isAlmostEqual(*ptr1++, *ptr2++)) return false;
    }
  }
  return true;
}


bool MEDDLY::mt_forest::equalsFS(int h1, int h2) const
{
  MEDDLY_DCASSERT(isFullNode(h1));
  MEDDLY_DCASSERT(isSparseNode(h2));

  int *ptr1 = getNodeAddress(h1) + 2;
  int *ptr2 = getNodeAddress(h2) + 2;
  int sz1 = *ptr1++;
  int sz2 = -(*ptr2++);

  int* h1Start = ptr1;
  int* h1Stop = ptr1 + sz1;
  int* h2Stop = ptr2 + sz2;
  int* down2 = h2Stop;

  // If the last index in h2 does not exist in h1, return false.
  // Otherwise, h1 is either the same "size" as h2 or larger than h2.

  if (h2Stop[-1] >= sz1) {
    // Last index of h2 does not exist in h1.
    return false;
  }

  while (ptr2 != h2Stop) {
    int index = *ptr2++;
    MEDDLY_DCASSERT(index < sz1);
    int* stop = h1Start + index;
    while (ptr1 != stop) { if (*ptr1++ != 0) return false; }
    if (*ptr1++ != *down2++) return false;
  }

  while (ptr1 != h1Stop) {
    if (*ptr1++ != 0) return false;
  }

  if (isMultiTerminal()) return true;

  // Check edge-values
  MEDDLY_DCASSERT(ptr1 == h1Stop);
  MEDDLY_DCASSERT(ptr2 == h2Stop);
  MEDDLY_DCASSERT(down2 == h2Stop + sz2);

  // ptr1 and down2 are pointing at the start of edge-values
  // Reset the index pointer for h2 (sparse node).
  ptr2 -= sz2;
  if (isEVPlus()) {
    while (ptr2 != h2Stop) {
      if (ptr1[*ptr2++] != *down2++) return false;
    }
  } else {
    MEDDLY_DCASSERT(isEVTimes());
    while (ptr2 != h2Stop) {
      if (!isAlmostEqual(ptr1[*ptr2++], *down2++)) return false;
    }
  }
  return true;
}


bool MEDDLY::mt_forest::equals(int h1, int h2) const 
{
  MEDDLY_DCASSERT(h1);	
  MEDDLY_DCASSERT(h2);
  MEDDLY_DCASSERT(isActiveNode(h1));
  MEDDLY_DCASSERT(isActiveNode(h2));
  MEDDLY_DCASSERT(!isTerminalNode(h1));
  MEDDLY_DCASSERT(!isTerminalNode(h2));

  if (getNodeLevel(h1) != getNodeLevel(h2)) { return false; }

  return
    isFullNode(h1)
    ? (isFullNode(h2) ? equalsFF(h1, h2) : equalsFS(h1, h2))
    : (isFullNode(h2) ? equalsFS(h2, h1) : equalsSS(h1, h2));
}


// ------------------------------------------------------------------
//  Protected methods
// ------------------------------------------------------------------

void MEDDLY::mt_forest::deleteNode(int p)
{
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(getInCount(p) == 0);
  MEDDLY_DCASSERT(isActiveNode(p));

#if 0
  validateIncounts();
#endif

  int* foo = getNodeAddress(p);
  int k = getNodeLevel(p);

  // remove from unique table (only applicable to reduced nodes)
  if (isReducedNode(p)) {
#ifdef DEVELOPMENT_CODE
    MEDDLY_DCASSERT(unique->find(p) == p);
#endif

#ifdef TRACK_DELETIONS
    showNode(stdout, p);
#endif

#ifdef DEVELOPMENT_CODE
    int x = unique->remove(p);
#else
    unique->remove(p);
#endif

#ifdef TRACK_DELETIONS
    printf("%s: p = %d, unique->remove(p) = %d\n", __func__, p, x);
    fflush(stdout);
#endif

    MEDDLY_DCASSERT(x != -1);
    MEDDLY_DCASSERT(p == x);
    MEDDLY_DCASSERT(address[p].cache_count == 0);
  }
  else {
    // Temporary node
    // TODO:
    // clear cache of corresponding temporary node?
    decrTempNodeCount(k);
  }

  // unlink children
  const int nDptrs = ABS(foo[2]);
  int* downptr = foo + 3 + (foo[2] < 0? nDptrs: 0);
  int* stop = downptr + nDptrs;
#if ENABLE_IN_COUNTING
  while (downptr < stop) {
    int temp = *downptr;
    *downptr++ = 0;
    unlinkNode(temp);
  }
#else
  while (downptr < stop) {
    unlinkNode(*downptr++);
  }
#endif

  // Recycle node memory
  levels[k].recycleNode(getNodeOffset(p));
  /*
  levels[k].makeHole(getNodeOffset(p), getDataHeaderSize() +
      nDptrs * ((foo[2] < 0)
        ? (isMultiTerminal() ? 2: 3)
        : (isMultiTerminal() ? 1: 2)
        ));
  */

  // recycle the index
  freeNode(p);

  if (levels[k].compactLevel) levels[k].compact(address);

#if 0
  validateIncounts();
#endif

}

void MEDDLY::mt_forest::zombifyNode(int p)
{
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(isReducedNode(p));
  MEDDLY_DCASSERT(getCacheCount(p) > 0);  // otherwise this node should be deleted
  MEDDLY_DCASSERT(getInCount(p) == 0);
  MEDDLY_DCASSERT(address[p].cache_count > 0);

  stats.zombie_nodes++;
  levels[getNodeLevel(p)].zombie_nodes++;
  stats.decActive(1);

  // mark node as zombie
  address[p].cache_count = -address[p].cache_count;

  MEDDLY_DCASSERT(unique->find(p) == p);
#ifdef DEVELOPMENT_CODE 
  int x = unique->remove(p);
  MEDDLY_DCASSERT(x==p);
#else
  unique->remove(p);
#endif

  int node_level = getNodeLevel(p);
  int node_offset = getNodeOffset(p);
  int* foo = getNodeAddress(p);

  address[p].offset = 0;

  // unlinkNode children
  if (foo[2] < 0) {
    // Sparse encoding
    int* downptr = foo + 3 - foo[2];
    int* stop = downptr - foo[2];
    for (; downptr < stop; ++downptr) {
#if ENABLE_IN_COUNTING
      int temp = *downptr;
      *downptr = 0;
      unlinkNode(temp);
#else
      unlinkNode(*downptr);
#endif
    }
    // Recycle node memory
    levels[node_level].recycleNode(node_offset);
    /*
    levels[node_level].makeHole(node_offset, getDataHeaderSize()
        - (isMultiTerminal() ? 2: 3) * foo[2]);  
    */
  } else {
    // Full encoding
    int* downptr = foo + 3;
    int* stop = downptr + foo[2];
    for (; downptr < stop; ++downptr) {
#if ENABLE_IN_COUNTING
      int temp = *downptr;
      *downptr = 0;
      unlinkNode(temp);
#else
      unlinkNode(*downptr);
#endif
    }
    // Recycle node memory
    levels[node_level].recycleNode(node_offset);
    /*
    levels[node_level].makeHole(node_offset, getDataHeaderSize()
        + (isMultiTerminal() ? 1: 2) * foo[2]);  
    */
  }
}


void MEDDLY::mt_forest::garbageCollect() {
  gc();
}


bool MEDDLY::mt_forest::gc(bool zombifyOrphanNodes) {
#if ENABLE_GC
  // change isStale such that all nodes with (incount == 0)
  // are considered to be stale
  if (performing_gc) return false;

  performing_gc = true;

  bool freed_some = false;
  nodes_activated_since_gc = 0;

#ifdef DEBUG_GC
  printf("Garbage collection in progress... \n");
  fflush(stdout);
#endif

  if (isPessimistic()) {
#ifdef DEBUG_GC
    printf("Zombie nodes: %ld\n", zombie_nodes);
#endif
    // remove the stale nodes entries from caches
    removeStaleComputeTableEntries();
#ifdef DEBUG_GC
    printf("Zombie nodes: %ld\n", zombie_nodes);
#endif
#ifdef DEVELOPMENT_CODE
    if (stats.zombie_nodes != 0) {
      showInfo(stderr, 1);
      showComputeTable(stderr, 6);
    }
    MEDDLY_DCASSERT(stats.zombie_nodes == 0);
#endif
    freed_some = true;
  } else {

    if (zombifyOrphanNodes) {
#ifdef DEBUG_GC
      fprintf(stderr, "Active: %ld, Zombie: %ld, Orphan: %ld\n",
          stats.active_nodes, stats.zombie_nodes, stats.orphan_nodes);
#endif
      // convert orphans to zombies
      // nodeDeletionPolicy = forest::PESSIMISTIC_DELETION;
      stats.orphan_nodes = 0;
      MEDDLY_DCASSERT(stats.zombie_nodes == 0);
      for (int i = 1; i <= a_last; i++) {
        MEDDLY_DCASSERT(!isTerminalNode(i));
        if (isActiveNode(i) && getInCount(i) == 0) {
          MEDDLY_DCASSERT(getCacheCount(i) > 0);
          zombifyNode(i);
        }
      }
#ifdef DEBUG_GC
      fprintf(stderr, "Active: %ld, Zombie: %ld, Orphan: %ld\n",
          stats.active_nodes, zombie_nodes, orphan_nodes);
#endif
      // remove the stale nodes entries from caches
      removeStaleComputeTableEntries();
#ifdef DEVELOPMENT_CODE
      if (stats.zombie_nodes != 0) {
        showInfo(stderr, 1);
        showComputeTable(stderr, 5);
      }
      // TBD: better error message here
#endif
      MEDDLY_DCASSERT(stats.zombie_nodes == 0);
      // nodeDeletionPolicy = forest::OPTIMISTIC_DELETION;
    } else {

      // remove the stale nodes entries from caches
      removeStaleComputeTableEntries();

    }

    freed_some = true;
  }

#ifdef DEBUG_GC
  printf("Compacting levels...\n");
  fflush(stdout);
#endif

  compactAllLevels();

#ifdef DEBUG_GC
  printf("  done compacting.\n");
  fflush(stdout);
#endif

  performing_gc = false;

  return freed_some;

#else // ENABLE_GC

  return false;

#endif // ENABLE_GC
}


void MEDDLY::mt_forest::removeZombies(int max_zombies) {
#if 1
  return;
#else
  // too many zombies? kill em!
  if (zombie_nodes > max_zombies && stats.active_nodes/zombie_nodes < 3) {
#if 0
    for (int i = 1; i <= a_last; i++) {
      if (isActiveNode(i) && isZombieNode(i)) {
        showNode(stdout, i); printf("\n");
      }
    }
#endif
    // remove the stale nodes entries from caches
    removeStaleComputeTableEntries();
#if 0
    if (zombie_nodes > 0) {
      for (int i = 1; i <= a_last; i++) {
        if (isActiveNode(i) && isZombieNode(i)) {
          showNode(stdout, i); printf("\n");
        }
      }
    }
#endif
    MEDDLY_DCASSERT(zombie_nodes == 0);
  }
#endif
}

int MEDDLY::mt_forest::getFreeNode(int k)
{
  if (a_unused) {
    // grab a recycled index
    int p = a_unused;
    MEDDLY_DCASSERT(address[p].offset < 1);
    a_unused = -address[p].offset;
    stats.incActive(1);
    return p;
  }
  // new index
  // TODO:
  if (a_last + 1 >= a_size) {
    // compactLevel(k);
    // if (a_last + 1 < a_size) return getFreeNode(k);
    // int new_a_size = (a_size > 16384)? a_size + 16384: a_size * 2;

#if 0
    int new_a_size = (a_size > 1024)? a_size + 1024: a_size * 2;
#else
    // increase size by 50%
    int min_size = (a_last + 1) * 0.375;
    int new_a_size = min_size * 4;
#endif

    mdd_node_data *temp = (mdd_node_data*) realloc(address,
        new_a_size * sizeof(mdd_node_data));
    // MEDDLY_DCASSERT(NULL != temp);
    if (NULL == temp) {
      fprintf(stderr, "Memory allocation error while allocating MDD nodes.\n");
      exit(1);
    }
    stats.incMemAlloc((new_a_size - a_size) * sizeof(mdd_node_data));
    address = temp;
    memset(address + a_size, 0, (new_a_size - a_size) * sizeof(mdd_node_data));
    a_size = new_a_size;
  }
  a_last++;
  stats.incActive(1);
  if (getCurrentNumNodes() > peak_nodes) peak_nodes = getCurrentNumNodes();
  return a_last;
}

void MEDDLY::mt_forest::freeZombieNode(int p)
{
  MEDDLY_DCASSERT(address[p].level != 0);
  MEDDLY_DCASSERT(address[p].cache_count == 0);
  stats.zombie_nodes--;
  levels[address[p].level].zombie_nodes--;
  address[p].level = 0;
  address[p].cache_count = 0;
  if (p == a_last) { 
    // special case
    address[p].offset = 0;
    a_last--;
  } else {
    address[p].offset = -a_unused;
    a_unused = p;
  }
#ifdef TRACK_DELETIONS
  printf("reclaimed zombie %d\n", p);
#endif
}

void MEDDLY::mt_forest::freeNode(int p)
{
#ifdef TRACE_REDUCE
  printf("%s: p = %d, a_last = %d\n", __func__, p, a_last);
#endif

  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
  MEDDLY_DCASSERT(address[p].cache_count == 0);

  stats.decActive(1);

  address[p].level = 0;
  address[p].cache_count = 0;
  if (p == a_last) { 
    // special case
    address[p].offset = 0;
    a_last--;
    if (a_size > add_size && a_last < a_size/2) {
      address = (mdd_node_data *)
          realloc(address, a_size/2 * sizeof(mdd_node_data));
      if (NULL == address) throw MEDDLY::error(MEDDLY::error::INSUFFICIENT_MEMORY);
      a_size /= 2;
      stats.decMemAlloc(a_size * sizeof(mdd_node_data));
#ifdef MEMORY_TRACE
      printf("Reduced node[] by a factor of 2. New size: %d.\n", a_size);
#endif
    }
  } else {
    address[p].offset = -a_unused;
    a_unused = p;
  }
}

void MEDDLY::mt_forest::reportMemoryUsage(FILE * s, const char filler) {
  fprintf(s, "%cPeak Nodes:             %ld\n", filler, getPeakNumNodes());
  fprintf(s, "%cActive Nodes:           %ld\n", filler, getCurrentNumNodes());
#if 0
  unsigned count = 0;
  for (int i = 1; i <= a_last; ++i) if (isActiveNode(i)) ++count;
  fprintf(s, "%cActive Nodes (manual):\t\t%d\n", filler, count);
  fprintf(s, "%c%cZombie Nodes:\t\t%d\n", filler, filler,
      getZombieNodeCount());
  fprintf(s, "%c%cTemp Nodes:\t\t%d\n", filler, filler, getTempNodeCount());
  fprintf(s, "%c%cOrphan Nodes:\t\t%d\n", filler, filler,
      getOrphanNodeCount());
#endif
  fprintf(s, "%cReclaimed Nodes:        %ld\n", filler, stats.reclaimed_nodes);
  fprintf(s, "%cMem Used:               %ld\n", filler,
      getCurrentMemoryUsed());
  fprintf(s, "%cPeak Mem Used:          %ld\n", filler, getPeakMemoryUsed());
  fprintf(s, "%cMem Allocated:          %ld\n", filler,
      getCurrentMemoryAllocated());
  fprintf(s, "%cPeak Mem Allocated:     %ld\n",
      filler, getPeakMemoryAllocated());
  fprintf(s, "%cUnique Tbl Mem Used:    %ld\n", filler,
      getUniqueTableMemoryUsed());
  fprintf(s, "%cCompactions:            %ld\n", filler, stats.num_compactions);
#if 0
  fprintf(s, "%cHole Memory Usage:\t%d\n", filler, getHoleMemoryUsage());
  fprintf(s, "%cMax Hole Chain:\t%d\n", filler, getMaxHoleChain());
  fprintf(s, "%cCompactions:\t\t%d\n", filler, getCompactionsCount());
  // compareCacheCounts();
#endif

#if 1
  // Print hole-recyling info
  // Compute chain lengths
  std::map<int, int> chainLengths;

  for (int k=getMinLevelIndex(); k<=getNumVariables(); k++) 
  {
    levels[k].addToChainCounts(chainLengths);
  }

  fprintf(s, "Hole Chains (size, count):\n");
  for (std::map<int, int>::iterator iter = chainLengths.begin();
    iter != chainLengths.end(); ++iter)
  {
    fprintf(s, "\t%d: %d\n", iter->first, iter->second);
  }
#endif
}

void MEDDLY::mt_forest::compareCacheCounts(int p)
{
#if ENABLE_CACHE_COUNTING
  counting = true;
  if (p == -1) {
    // get cache counts
    unsigned sz = a_last + 1;
    unsigned count[sz];
    memset(count, 0, sizeof(unsigned) * sz);
    for (unsigned i = 0; i < nm_users.size(); i++) {
      if (nm_users[i] != NULL) nm_users[i]->getCacheCounts(count, sz);
    }
    // verify counts
    if (isPessimistic()) {
      assert(false);
    } else {
      // active nodes count should match and inactive nodes' count must be 0.
      for (int i = 0; i < (a_last + 1); ++i) {
        if (isActiveNode(i)) {
          if (!isTerminalNode(i)) {
#if USE_MDD_HASH_TABLE
            assert(count[i] == getCacheCount(i));
#else
            assert(count[i] == getCacheCount(i) - 1);
#endif
          }
        } else {
          assert(count[i] == 0);
        }
      }
    }
  } else {
    // get cache counts
    unsigned count = 0;
    for (unsigned i = 0; i < nm_users.size(); i++) {
      if (nm_users[i] != NULL) count += nm_users[i]->getCacheCount(p);
    }
    // verify counts
    if (isPessimistic()) {
      assert(false);
    } else {
      // active nodes count should match and inactive nodes' count must be 0.
      assert(count == getCacheCount(p));
   }
  }
  counting = false;
#endif
}


void MEDDLY::mt_forest::validateIncounts()
{
#if ENABLE_IN_COUNTING

  // Inspect every active node's down pointers to determine
  // the incoming count for every active node.
  
  unsigned sz = a_last + 1;
  unsigned in_count[sz];
  memset(in_count, 0, sizeof(unsigned) * sz);
  const int *dptrs = NULL;
  for (int i = 0; i < sz; ++i) {
    if (isActiveNode(i) && !(isTerminalNode(i))) {
      // count down-pointers
      if (isFullNode(i)) {
        dptrs = getFullNodeDownPtrsReadOnly(i);
        for (int j = getFullNodeSize(i) - 1; j >=0 ; --j) {
          if (isTerminalNode(dptrs[j])) continue;
          MEDDLY_CHECK_RANGE(0, dptrs[j], sz);
          in_count[dptrs[j]]++;
        }
      } else {
        MEDDLY_DCASSERT(isSparseNode(i));
        dptrs = getSparseNodeDownPtrs(i);
        for (int j = getSparseNodeSize(i) - 1; j >= 0; --j) {
          if (isTerminalNode(dptrs[j])) continue;
          MEDDLY_CHECK_RANGE(0, dptrs[j], sz);
          in_count[dptrs[j]]++;
        }
      }
    }
  }

  // Validate the incoming count stored with each active node using the
  // in_count array computed above
  for (int i = 0; i < sz; ++i) {
    if (isActiveNode(i) && !(isTerminalNode(i))) {
      assert(in_count[i] <= getInCount(i));
    }
  }
  
#endif
}


void MEDDLY::mt_forest::showLevel(FILE *s, int k) const {
  dumpInternalLevel(s, k);
}
void MEDDLY::mt_forest::showAll(FILE *s, int verb) const { 
  if (0==verb)  return;
  if (1==verb)  dump(s);
  else          dumpInternal(s); 
}

void MEDDLY::mt_forest::show(FILE *s, int h) const { fprintf(s, "%d", h); }

long MEDDLY::mt_forest::getHoleMemoryUsage() const {
  long sum = 0;
  for (int i=getMinLevelIndex(); i<=getNumVariables(); i++)
    sum += levels[i].getHoleSlots();
  return sum * sizeof(int); 
}

void MEDDLY::mt_forest::validateDownPointers(int p, bool recursive)
{
  if (isTerminalNode(p)) return;

  int nodeHeight = getNodeHeight(p);
  int nodeLevel = getNodeLevel(p);
  int nodeSize = isFullNode(p)? getFullNodeSize(p): getSparseNodeSize(p);
  const int* ptr = isFullNode(p)? getFullNodeDownPtrsReadOnly(p):
    getSparseNodeDownPtrs(p);

  switch (getReductionRule()) {
    case policies::FULLY_REDUCED:
      if (isUnprimedNode(p)) {
        // unprimed node
        for (int i = 0; i < nodeSize; ++i) {
          //assert(isReducedNode(ptr[i]));
          assert(!isForRelations() ||
              isTerminalNode(ptr[i]) ||
              getNodeHeight(ptr[i]) < nodeHeight ||
              getNodeLevel(ptr[i]) == -nodeLevel);
        }
      } else {
        // primed node
        for (int i = 0; i < nodeSize; ++i) {
          //assert(isReducedNode(ptr[i]));
          assert(isTerminalNode(ptr[i]) ||
              getNodeHeight(ptr[i]) < nodeHeight);
        }
      }
      break;

    case policies::QUASI_REDUCED:
      if (isUnprimedNode(p)) {
        // unprimed node
        for (int i = 0; i < nodeSize; ++i) {
          //assert(isReducedNode(ptr[i]));
          assert(!isForRelations() ||
              isTerminalNode(ptr[i]) ||
              getNodeLevel(ptr[i]) == -nodeLevel);
        }
      } else {
        // primed node
        for (int i = 0; i < nodeSize; ++i) {
          //assert(isReducedNode(ptr[i]));
          assert(isTerminalNode(ptr[i]) ||
              (getNodeHeight(ptr[i]) == (nodeHeight - 1) &&
               isUnprimedNode(ptr[i])));
        }
      }
      break;

    case policies::IDENTITY_REDUCED:
      assert(isForRelations());
      if (isUnprimedNode(p)) {
        // unprimed node
        for (int i = 0; i < nodeSize; ++i) {
          //assert(isReducedNode(ptr[i]));
          assert(ptr[i] == 0 || (getNodeLevel(ptr[i]) == -nodeLevel));
        }
      } else {
        // primed node
        for (int i = 0; i < nodeSize; ++i) {
          //assert(isReducedNode(ptr[i]));
          assert(getNodeHeight(ptr[i]) < nodeHeight);
          assert(isTerminalNode(ptr[i]) || isUnprimedNode(ptr[i]));
        }
      }
      break;

    default:
      break;
  }

  if (recursive) {
    for (int i = 0; i < nodeSize; ++i) {
      validateDownPointers(ptr[i], true);
    }
  }
}


int MEDDLY::mt_forest::addReducedNodes(int a, int b)
{
  MEDDLY_DCASSERT(isReducedNode(a));
  MEDDLY_DCASSERT(isReducedNode(b));

  // Neither a nor b is a terminal node.
  // Compute result using dd_edge::operator+=.
  if (nodeA == 0) {
    nodeA = new dd_edge(this);
    MEDDLY_DCASSERT(nodeB == 0);
    nodeB = new dd_edge(this); 
  }

  dd_edge& A = *nodeA;
  dd_edge& B = *nodeB;

  linkNode(a);
  linkNode(b);
  A.set(a, 0, getNodeLevel(a));
  B.set(b, 0, getNodeLevel(b));

  A += B;
  int result = sharedCopy(A.getNode());
  A.clear();
  B.clear();

  return result;
}


// Creates a temporary node as a copy of node a.
// The size variable can be used to create a node a size larger than a.
int MEDDLY::mt_forest::makeACopy(int a, int size)
{
  MEDDLY_DCASSERT(isMultiTerminal());
  int result = 0;
  int* rDptrs = 0;

  if (isFullNode(a)) {
    int aSize = getFullNodeSize(a);
    int newSize = MAX ( size ,  aSize ) ;
    result = createTempNode(getNodeLevel(a), newSize, false);
    assert(getDownPtrs(result, rDptrs));
    const int* aDptrs = 0;
    assert(getDownPtrs(a, aDptrs));
    for (const int* end = aDptrs + aSize; aDptrs != end; ) {
      *rDptrs++ = sharedCopy(*aDptrs++);
    }
    for (const int* end = rDptrs + (newSize - aSize); rDptrs != end; ) {
      *rDptrs++ = 0;
    }
  }
  else {
    MEDDLY_DCASSERT(isSparseNode(a));
    int nDptrs = getSparseNodeSize(a);
    int aSize = 1 + getSparseNodeIndex(a, nDptrs - 1);
    int newSize = MAX ( size, aSize ) ;
    result = createTempNode(getNodeLevel(a), newSize, true);
    assert(getDownPtrs(result, rDptrs));
    const int* aDptrs = 0;
    assert(getDownPtrs(a, aDptrs));
    const int* aIndexes = 0;
    assert(getSparseNodeIndexes(a, aIndexes));
    for (const int* end = aDptrs + nDptrs; aDptrs != end; ) {
      rDptrs[*aIndexes++] = sharedCopy(*aDptrs++);
    }
  }
  return result;
}


// For all i, a[i] += b
int MEDDLY::mt_forest::accumulateExpandA(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(!isIdentityReduced());

  bool needsToMakeACopy = cBM;
  int savedTempNode = a;

  int aSize = getFullNodeSize(a);
  int aLevel = getNodeLevel(a);
  int levelSize = getLevelSize(aLevel);

  if (aSize < levelSize) {
    if (needsToMakeACopy) {
      a = makeACopy(a, levelSize);
      needsToMakeACopy = false;
    } else {
      resizeNode(a, levelSize);
    }
    aSize = getFullNodeSize(a);
  }

  MEDDLY_DCASSERT(aSize == levelSize);

  for (int i = 0; i < aSize; i++) {
    int dptr = getFullNodeDownPtr(a, i);
    int result = accumulateMdd(dptr, b, cBM);
    if (result != dptr) {
      if (needsToMakeACopy) {
        a = makeACopy(a);
        needsToMakeACopy = false;
      }
      setDownPtr(a, i, result);
    }
    unlinkNode(result);
  }

  return savedTempNode == a? sharedCopy(a): a;
}


int MEDDLY::mt_forest::accumulateMdd(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(!isIdentityReduced());
  MEDDLY_DCASSERT(isReducedNode(b));

  // Terminal nodes
  if (a == 0 || b == 0) { return sharedCopy(a + b); }
  if (a == -1 || b == -1) { return sharedCopy(-1); }

  MEDDLY_DCASSERT(!isTerminalNode(a) && !isTerminalNode(b));

  // a is a reduced node
  if (isReducedNode(a)) {
    return addReducedNodes(a, b);
  }

  // a is a temporary node
  int aHeight = getNodeHeight(a);
  int bHeight = getNodeHeight(b);

  if (getInCount(a) > 1) cBM = true;

  if (aHeight > bHeight) {
    // b's levels were skipped.
    // only quasi- and fully- reduced Mdds.
    // a[i] += b
    return accumulateExpandA(a, b, cBM);
  }

  bool needsToMakeACopy = cBM;
  int savedTempNode = a;

  if (aHeight < bHeight) {
    // Build node c at the same level as b.
    // set all c[i] = a;
    int temp = a;
    a = createTempNodeMaxSize(getNodeLevel(b), false);
    setAllDownPtrsWoUnlink(a, temp);
    needsToMakeACopy = false;
  }

  // Expand both nodes. a is a full node, b can be either sparse or full.

  // Node b is Truncated-Full
  if (isFullNode(b)) {
    int size = getFullNodeSize(b);
    // Resize a.
    if (getFullNodeSize(a) < size) {
      if (needsToMakeACopy) {
        a = makeACopy(a, size);
        needsToMakeACopy = false;
      } else {
        resizeNode(a, size);
      }
      MEDDLY_DCASSERT(getFullNodeSize(a) == size);
    }
    // Accumulate into a.
    for (int i = 0; i < size; ++i) {
      int dptr = getFullNodeDownPtr(a, i);
      int result = accumulateMdd(dptr, getFullNodeDownPtr(b, i), cBM);
      if (result != dptr) {
        if (needsToMakeACopy) {
          a = makeACopy(a);
          needsToMakeACopy = false;
        }
        setDownPtr(a, i, result);
      }
      unlinkNode(result);
    }
  }
  // Node b is Sparse
  else {
    MEDDLY_DCASSERT(isSparseNode(b));
    int nDptrs = getSparseNodeSize(b);
    int size = 1 + getSparseNodeIndex(b, nDptrs - 1);
    // Resize a.
    if (getFullNodeSize(a) < size) {
      if (needsToMakeACopy) {
        a = makeACopy(a, size);
        needsToMakeACopy = false;
      } else {
        resizeNode(a, size);
      }
      MEDDLY_DCASSERT(getFullNodeSize(a) == size);
    }
    // Accumulate into a.
    for (int i = 0; i < nDptrs; ++i) {
      int index = getSparseNodeIndex(b, i);
      int dptr = getFullNodeDownPtr(a, index);
      int result = accumulateMdd(dptr, getSparseNodeDownPtr(b, i), cBM);
      if (result != dptr) {
        if (needsToMakeACopy) {
          a = makeACopy(a);
          needsToMakeACopy = false;
        }
        setDownPtr(a, index, result);
      }
      unlinkNode(result);
    }
  }

  return savedTempNode == a? sharedCopy(a): a;
}


void MEDDLY::mt_forest::accumulate(int& a, int b)
{
  if (isActiveNode(a) && isActiveNode(b)) {
    int result = accumulateMdd(a, b, false);
    unlinkNode(a);
    a = result;
    return;
  }
  throw error(error::INVALID_OPERATION);
}


// Add an element to a temporary edge
// Start this recursion at the top level in the domain.
// Use expert_domain::getNumVariables() to obtain the topmost level in
// the domain.
// cBM: copy before modifying.
int MEDDLY::mt_forest::accumulate(int tempNode, bool cBM,
    int* element, int level)
{
  MEDDLY_DCASSERT(isMdd());

  if (tempNode == -1) return -1;
  if (level == 0) {
    accumulateMintermAddedElement = true;
    return -1;
  }

  int index = element[level];
  int nodeLevel = getNodeLevel(tempNode);
  int nextLevel = level-1;

  int dptr = 0;
  int newDptr = 0;
  int inCount = 0;

  if (level == nodeLevel) {
    inCount = getInCount(tempNode);
    dptr = getDownPtr(tempNode, index);
  }
  else {
    // Levels have been skipped.
    // We are only dealing with MDDs here.
    inCount = getLevelSize(level);
    dptr = tempNode;
  }

  // An incount > 1 indicates a need to duplicate the node before
  // modifying.
  if (inCount > 1) cBM = true;

  newDptr = accumulate(dptr, cBM, element, nextLevel);

  if (newDptr == dptr) {
    // Element got absorbed into dptr
    return tempNode;
  }

  // If tempNode is 0, create a temporary node.
  // If tempNode is a reduced node or if its incount > 1,
  //    create a copy (which is a temporary node).
  // Otherwise, use tempNode (should be a temporary node with incount == 1).
  int newNode = 0;
  if (tempNode == 0) {
    newNode = createTempNode(level, index + 1, true);
  } else if (level != nodeLevel) {
    newNode = createTempNodeMaxSize(level, false);
    setAllDownPtrsWoUnlink(newNode, dptr);
  } else if (isReducedNode(tempNode)) {
    newNode = makeACopy(tempNode, index + 1);
  } else if (cBM) {
    newNode = makeACopy(tempNode, index + 1);
  } else {
    newNode = tempNode;
  }

  MEDDLY_DCASSERT(!isReducedNode(newNode));
  if (getFullNodeSize(newNode) < (index + 1)) {
    resizeNode(newNode, index + 1);
  }
  setDownPtr(newNode, index, newDptr);
  unlinkNode(newDptr);

  return newNode;
}


// Add an element to a temporary edge
bool MEDDLY::mt_forest::accumulate(int& tempNode, int* element)
{
  assert(isActiveNode(tempNode));
  assert(element != 0);

  // Enlarge variable bounds if necessary
  for (int level=1; level<=getExpertDomain()->getNumVariables(); level++) {
    int sz = element[level] + 1;
    if (sz > getExpertDomain()->getVariableBound(level)) {
      useExpertDomain()->enlargeVariableBound(level, false, sz);
    }
  }

  accumulateMintermAddedElement = false;
  int result = accumulate(tempNode, false,
      element, getExpertDomain()->getNumVariables());
  if (tempNode != result) {
    // tempNode had to be copied into another node by accumulate().
    // This could be either because tempNode was a reduced node,
    // or because tempNode had incount > 1.
    unlinkNode(tempNode);
    tempNode = result;
  }
  // Note: tempNode == result indicates that the element was added
  // to the existing temporary node. Therefore, there is no need to
  // change incounts.
  return accumulateMintermAddedElement;
}


int MEDDLY::mt_forest::recursiveReduceNode(std::map<int, int>& cache,
    int root)
{
  MEDDLY_DCASSERT(!isReducedNode(root));
  MEDDLY_DCASSERT(isFullNode(root));

  // Check cache for result
  std::map<int, int>::iterator iter = cache.find(root);
  if (iter != cache.end()) {
    // Cache hit
    linkNode(iter->second);
    unlinkNode(root);
    return iter->second;
  }

  int size = getFullNodeSize(root);

  for (int i = 0; i < size; ++i)
  {
    int dptr = getFullNodeDownPtr(root, i);

    // Ignore terminal nodes and reduced nodes
    if (isReducedNode(dptr)) continue;

    int temp = recursiveReduceNode(cache, dptr);
    // At this point, temp's incount has been increased and
    // dptr's incount has been decreased by 1.

    MEDDLY_DCASSERT(isReducedNode(temp));

    // Using WoUnlink, since dptr's has already been decreased
    setDownPtrWoUnlink(root, i, temp);
    // Unlinking temp because setDownPtr increases temp's incount by 1
    unlinkNode(temp);
  }

  int result = reduceNode(root);

  // Save result in cache
  if (isActiveNode(root)) cache[root] = result;

  return result;
}


int MEDDLY::mt_forest::recursiveReduceNode(int tempNode, bool clearCache)
{
  // Recursive procedure:
  // Start at root (i.e. tempNode).
  // Build a temporary node for root.
  // -- Build reduced nodes for each child(root).
  // Keep track of duplicates via a compute cache which maps
  //   each temporary node to a reduced node.

  MEDDLY_DCASSERT(!isReducedNode(tempNode));

  std::map<int, int>& cache = recursiveReduceCache;
  if (clearCache) cache.clear();
  return recursiveReduceNode(cache, tempNode);
}


int MEDDLY::mt_forest::createTempNode(int k, int sz, bool clear)
{
  MEDDLY_DCASSERT(k != 0);

  if (isTimeToGc()) {
    fprintf(stderr, "Started forest garbage collector.\n");
    gc();
    fprintf(stderr, "Stopped forest garbage collector.\n");
  }

  MEDDLY_DCASSERT(k);
  MEDDLY_DCASSERT(isValidLevel(k));
  MEDDLY_CHECK_RANGE(1, sz, getLevelSize(k) + 1);

  // get a location in address[] to store the node
  int p = getFreeNode(k);

#ifdef DEBUG_MDD_H
  printf("%s: k: %d, sz: %d, new p: %d\n", __func__, k, sz, p);
  fflush(stdout);
#endif

  // fill in the location with p's address info
  MEDDLY_DCASSERT(isMultiTerminal());
  address[p].level = k;
  // address[p].offset = levels[k].getHole(4 + sz, true);
  address[p].offset = levels[k].allocNode(sz, p, clear);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_H
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

  /*
  int* foo = levels[k].data + address[p].offset;
  foo[0] = 1;                   // #incoming
  foo[1] = temp_node;
  foo[2] = sz;                  // size
  foo[3 + sz] = p;              // pointer to this node in the address array

  // initialize
  if (clear) initDownPtrs(p);
  */

#ifdef TRACK_DELETIONS
  cout << "Creating node " << p << "\n";
  cout.flush();
#endif

  incrTempNodeCount(k);
  nodes_activated_since_gc++;

  return p;
}


void MEDDLY::mt_forest::handleNewOrphanNode(int p) {
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(getInCount(p) == 0);

  // insted of orphan_nodes++ here; do it only when the orphan is not going
  // to get deleted or converted into a zombie

  // Two possible scenarios:
  // (1) a reduced node, or
  // (2) a temporary node ready to be deleted.
  MEDDLY_DCASSERT(isReducedNode(p) || getCacheCount(p) == 0);

  if (getCacheCount(p) == 0) {
    // delete node
    // this should take care of the temporary nodes also
#ifdef TRACK_DELETIONS
    cout << "Deleting node " << p << " from unlinkNode\t";
    showNode(stdout, p);
    cout << "\n";
    cout.flush();
#endif
    deleteNode(p);
  }
  else if (isPessimistic()) {
    // zombify node
    zombifyNode(p);
  }
  else {
    stats.orphan_nodes++;
  }

#if 0
  if (getOrphanNodeCount() > 100000)
    smart_cast<expert_compute_manager*>(MEDDLY::getComputeManager())
      ->removeStales(this);
#endif
}


// ********************* utils ************************

bool MEDDLY::mt_forest::singleNonZeroAt(int p, int val, int index) const
{
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(!isZombieNode(p));
  MEDDLY_DCASSERT(val != 0);
  if (isFullNode(p)) {
    const int* dptr = getFullNodeDownPtrsReadOnly(p);
    const int sz = getFullNodeSize(p);
    if (index >= sz || dptr[index] != val) return false;
    int i = 0;
    for ( ; i < index; ++i) { if (dptr[i] != 0) return false; }
    for (i = index + 1 ; i < sz; ++i) { if (dptr[i] != 0) return false; }
  } else {
    if (getSparseNodeSize(p) != 1) return false;
    if (getSparseNodeIndex(p, 0) != index) return false;
    if (getSparseNodeDownPtr(p, 0) != val) return false;
  }
  return true;
}


bool MEDDLY::mt_forest::checkForReductions(int p, int nnz, int& result)
{
  if (isQuasiReduced()) return false;
  if (nnz != getLevelSize(getNodeLevel(p))) return false;

  const int* ptr = getFullNodeDownPtrs(p);
  int size = getFullNodeSize(p);

  switch (getReductionRule()) {

    case policies::FULLY_REDUCED:
      result = ptr[0];
      for (int i = 1; i < size; ++i) {
        if (ptr[i] != result) return false;
      }
      break;

    case policies::IDENTITY_REDUCED:
      if (isForRelations()) {
        if (isPrimedNode(p)) return false;
        if (isFullNode(ptr[0])) {
          result = getFullNodeDownPtr(ptr[0], 0);
          if (result == 0) return false;
        } else {
          int index = getSparseNodeIndex(ptr[0], 0);
          if (index != 0) return false;
          result = getSparseNodeDownPtr(ptr[0], 0);
          MEDDLY_DCASSERT(result != 0);
        }
        for (int i = 0; i < size; i++) {
          if (!singleNonZeroAt(ptr[i], result, i)) return false;
        }
      }
      else {
        printf("Identity-Reduction is valid only for forests that ");
        printf("store relations.\n");
        printf("Either change reduction rule for forest %p or enable\n", this);
        printf("relations for it.\n");
        printf("Terminating.\n");
        exit(1);
      }
      break;

    default:
      return false;
  }

  return true;
}


