
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

// #define MEMORY_TRACE

const int add_size = 1024;
// const int l_add_size = 24;

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

#ifdef ACCUMULATE_ON
  nodeA = 0;
  nodeB = 0;
#endif
}


MEDDLY::mt_forest::~mt_forest()
{
#ifdef ACCUMULATE_ON
  if (nodeA) delete nodeA;
  if (nodeB) delete nodeB;
#endif

  clearLevelNodes();

  if (dptrsSize > 0) {
    free(dptrs);
    dptrsSize = 0;
    dptrs = 0;
  }

  // delete unique;
#if DEBUG_DELETE_NM
  printf("Deleted unique table\n");
  fflush(stdout);
#endif
}


// ******************************************************************
// *                                                                *
// *            virtual  and overriding default behavior            *
// *                                                                *
// ******************************************************************


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

bool MEDDLY::mt_forest::areDuplicates(int p, const node_builder &nb) const
{
  return areDupsInternal(p, nb);
}

bool MEDDLY::mt_forest::areDuplicates(int p, const node_reader &nr) const
{
  return areDupsInternal(p, nr);
}

/*
int MEDDLY::mt_forest
::createReducedHelper(int in, const node_builder &nb, bool &u)
{
#ifdef DEVELOPMENT_CODE
  validateDownPointers(nb);
#endif

  u = true;

  // get sparse, truncated full sizes and check
  // for redundant / identity reductions.
  int nnz;
  int truncsize;
  if (nb.isSparse()) {
    // Reductions for sparse nodes
    truncsize = -1;
    nnz = nb.getNNZs();
    for (int z=0; z<nnz; z++) {
      MEDDLY_DCASSERT(nb.d(z));
      truncsize = MAX(truncsize, nb.i(z));
    } // for z

    // Is this an identity node, and should we eliminate it?
    if (1==nnz && nb.getLevel()<0 && in==nb.i(0)) {
      if (isIdentityReduced()) {
        return linkNode(nb.d(0));
      }
    }

  } else {
    // Reductions for full nodes
    bool redundant = true;
    int common = nb.d(0);
    if (common) {
      nnz = 1;
      truncsize = 0;
    } else {
      nnz = 0;
      truncsize = -1;
    }
    for (int i=1; i<nb.getSize(); i++) {
      if (redundant) {
        redundant = (nb.d(i) == common);
      }
      if (nb.d(i)) {
        nnz++;
        truncsize = i;
      }
    } // for i

    // Is this a redundant node, and should we eliminate it?
    if (redundant) {
      if (isFullyReduced() || (isIdentityReduced() && nb.getLevel()>0))
        return linkNode(common);
    }

    // Is this an identity node, and should we eliminate it?
    if (isIdentityReduced()) {
      if (in>=0 && 1==nnz && nb.getLevel()<0 && nb.d(in)) 
        return linkNode(nb.d(in));
    }
  }
  truncsize++;

  // Is this a zero node?
  if (0==nnz) {
    MEDDLY_DCASSERT(0==truncsize);  // sanity check
    return 0;
  }

  // check for duplicates in unique table
  int q = unique->find(nb);
  if (q) {
    return linkNode(q);
  }

  // 
  // Not eliminated by reduction rule.
  // Not a duplicate.
  //
  // We need to create a new node for this.
  u = false;
  int p = getFreeNodeHandle();
  address[p].level = nb.getLevel();
  MEDDLY_DCASSERT(0 == address[p].cache_count);
  level_data &ld = levels[nb.getLevel()];

  // First, determine if it should be full or sparse
  if (ld.slotsForNode(-nnz) < ld.slotsForNode(truncsize)) { 
    // sparse node wins
    address[p].offset = ld.allocNode(-nnz, p, false);
    MEDDLY_DCASSERT(1==ld.countOf(address[p].offset));
    int* index = ld.sparseIndexesOf(address[p].offset);
    int* down  = ld.sparseDownOf(address[p].offset);
    nb.copyIntoSparse(down, index, nnz);
  } else {
    // full node wins
    address[p].offset = ld.allocNode(truncsize, p, false);
    MEDDLY_DCASSERT(1==ld.countOf(address[p].offset));
    int* down = ld.fullDownOf(address[p].offset);
    nb.copyIntoFull(down, truncsize);
  } // if 

  // add to UT 
  unique->add(nb.hash(), p);
  
#ifdef DEVELOPMENT_CODE
  node_finder key(this, p);
  MEDDLY_DCASSERT(key.hash() == nb.hash());
  MEDDLY_DCASSERT(unique->find(key) == p);
#endif
  return p;
}
*/

// ******************************************************************
// *                                                                *
// *                         Helper methods                         *
// *                                                                *
// ******************************************************************

/*
void MEDDLY::mt_forest::validateDownPointers(const node_builder &nb) const
{
  int nextLevel;
  switch (getReductionRule()) {
    case policies::IDENTITY_REDUCED:
    case policies::FULLY_REDUCED:
      if (nb.isSparse()) {
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), getNodeLevel(nb.d(z))));
        } 
      } else {
        for (int i=0; i<nb.getSize(); i++) {
          MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), getNodeLevel(nb.d(i))));
        }
      }
      break;

    case policies::QUASI_REDUCED:
      if (isForRelations()) 
        nextLevel = (nb.getLevel()<0) ? -(nb.getLevel()+1) : -nb.getLevel();
      else
        nextLevel = nb.getLevel()-1;
      if (nb.isSparse()) {
        for (int z=0; z<nb.getNNZs(); z++) {
          MEDDLY_DCASSERT(getNodeLevel(nb.d(z)) == nextLevel);
        } 
      } else {
        for (int i=0; i<nb.getSize(); i++) {
          MEDDLY_DCASSERT(getNodeLevel(nb.d(i)) == nextLevel);
        }
      }
      break;

    default:
      throw error(error::NOT_IMPLEMENTED);
  }

}
*/

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
  int node = createReducedNode(-1, nb);

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


/*
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
    nodeBuilder& nb = useNodeBuilder(-level, nodeSize);
    for (int i=0; i<nodeSize; i++) {
      nb.d(i) = vplist[level][i] ? linkNode(mask) : 0;
    }
    unlinkNode(mask);
    mask = createReducedNode(-1, nb);

    // create node at unprime level
    nodeSize = getExpertDomain()->getVariableBound(level, false);
    nb = useNodeBuilder(level, nodeSize);
    for (int i=0; i<nodeSize; i++) {
      nb.d(i) = vlist[level][i] ? linkNode(mask) : 0;
    }
    unlinkNode(mask);
    mask = createReducedNode(-1, nb);
  }

  b.set(mask, 0, getNodeLevel(mask));
#if 0
  b.show(stdout, 3);
#endif
  b *= a;
}
*/

#endif


/*
void MEDDLY::mt_forest::clearAllNodes()
{
  int level = getExpertDomain()->getNumVariables();
  while (level > 0 && stats.active_nodes > 0)
  {
    // find all nodes at curr level and make them into orphans
    for (int i = 1; i < getLastNode(); i++)
    {
      if (isActiveNode(i) && getNodeLevel(i) == level && getInCount(i) > 0) {
        getInCount(i) = 1;
        unlinkNode(i);
      }
    }

    if (stats.active_nodes > 0 && isForRelations()) {
      level = -level;
      // find all nodes at curr level and make them into orphans
      for (int i = 1; i < getLastNode(); i++)
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
*/

// *********************************************************************
// TODO: test this out
int MEDDLY::mt_forest::buildQuasiReducedNodeAtLevel(int k, int p)
{
  MEDDLY_DCASSERT(isQuasiReduced());
  int curr = p;
  int p_level = getNodeLevel(p);
  for (int i = p_level + 1; i <= k; i++)
  {
    insertRedundantNode(i, curr);
  }
  return curr;
}
// *********************************************************************

#if 0
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
#endif

/*
int ifTermGetInt(const MEDDLY::mt_forest *nm, int node)
{
  return nm->isTerminalNode(node) ? nm->getInteger(node) : node;
}
*/

/*
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
  int p_width = digits(getLastNode());
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
*/

// ------------------------------------------------------------------
//  Protected methods
// ------------------------------------------------------------------

/*
void MEDDLY::mt_forest::removeZombies(int max_zombies) {
#if 1
  return;
#else
  // too many zombies? kill em!
  if (zombie_nodes > max_zombies && stats.active_nodes/zombie_nodes < 3) {
#if 0
    for (int i = 1; i <= getLastNode(); i++) {
      if (isActiveNode(i) && isZombieNode(i)) {
        showNode(stdout, i); printf("\n");
      }
    }
#endif
    // remove the stale nodes entries from caches
    removeStaleComputeTableEntries();
#if 0
    if (zombie_nodes > 0) {
      for (int i = 1; i <= getLastNode(); i++) {
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
*/

void MEDDLY::mt_forest::compareCacheCounts(int p)
{
#if ENABLE_CACHE_COUNTING
  counting = true;
  if (p == -1) {
    // get cache counts
    unsigned sz = getLastNode() + 1;
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
      for (int i = 0; i < (getLastNode() + 1); ++i) {
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


long MEDDLY::mt_forest::getHoleMemoryUsage() const {
  long sum = 0;
  for (int i=getMinLevelIndex(); i<=getNumVariables(); i++)
    sum += levels[i].getHoleSlots();
  return sum * sizeof(int); 
}


#ifdef USE_OLD_TEMPNODES

/*
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
*/

#endif

#ifdef ACCUMULATE_ON

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
  int result = linkNode(A.getNode());
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
      *rDptrs++ = linkNode(*aDptrs++);
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
      rDptrs[*aIndexes++] = linkNode(*aDptrs++);
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

  return savedTempNode == a? linkNode(a): a;
}


int MEDDLY::mt_forest::accumulateMdd(int a, int b, bool cBM)
{
  MEDDLY_DCASSERT(!isIdentityReduced());
  MEDDLY_DCASSERT(isReducedNode(b));

  // Terminal nodes
  if (a == 0 || b == 0) { return linkNode(a + b); }
  if (a == -1 || b == -1) { return linkNode(-1); }

  MEDDLY_DCASSERT(!isTerminalNode(a) && !isTerminalNode(b));

  // a is a reduced node
  if (isReducedNode(a)) {
    return addReducedNodes(a, b);
  }

  // a is a temporary node
  int aHeight = getNodeHeight(a);
  int bHeight = getNodeHeight(b);

  if (readInCount(a) > 1) cBM = true;

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

  return savedTempNode == a? linkNode(a): a;
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
    inCount = readInCount(tempNode);
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

int MEDDLY::mt_forest::createTempNode(int k, int sz, bool clear)
{
  MEDDLY_DCASSERT(k != 0);

  if (isTimeToGc()) {
    fprintf(stderr, "Started forest garbage collector.\n");
    garbageCollect();
    fprintf(stderr, "Stopped forest garbage collector.\n");
  }

  MEDDLY_DCASSERT(k);
  MEDDLY_DCASSERT(isValidLevel(k));

  // get a location in address[] to store the node
  int p = getFreeNodeHandle();  // getFreeNode(k);

#ifdef DEBUG_MDD_H
  printf("%s: k: %d, sz: %d, new p: %d\n", __func__, k, sz, p);
  fflush(stdout);
#endif

  // fill in the location with p's address info
  MEDDLY_DCASSERT(isMultiTerminal());
  address[p].level = k;
  address[p].offset = levels[k].allocNode(sz, p, clear);
  address[p].cache_count = 0;

#ifdef DEBUG_MDD_H
  printf("%s: offset: %d\n", __func__, address[p].offset);
  fflush(stdout);
#endif

#ifdef TRACK_DELETIONS
  cout << "Creating node " << p << "\n";
  cout.flush();
#endif

  incrTempNodeCount(k);
  nodes_activated_since_gc++;

  return p;
}

#endif // ACCUMULATE_ON




// ********************* utils ************************

/*
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
*/
