
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


// TODO: implement changes made to node_manager interface
//
// TODO: HERE: go through every function in mdds.h and mdds.cc


#include "mdds.h"

#include <map>    // for getCardinality()
#include <queue>  // for showNodeGraph
#include <vector> // for showNodeGraph
#include <set>
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#define DONT_USE_FULL 0

#define ENABLE_GC 1
#define DEBUG_GC
// #define ALT_ORPHAN_GC
#define ENABLE_CACHE_COUNTING 0
#define ENABLE_IN_COUNTING 0

#define DEBUG_DELETE_NM 0

const int add_size = 1024;
const int l_add_size = 24;

#if 0
double cardinality(const node_manager* const nm, int p, int k,
    std::map<int, double>& visited);
#endif

node_manager::node_manager(domain *d, bool rel, range_type t,
  edge_labeling ev, reduction_rule r, node_storage s, node_deletion_policy nd)
: expert_forest(d, rel, t, ev, r, s, nd)
{
  DCASSERT(d != NULL);
  expertDomain = smart_cast<expert_domain*>(d);
  DCASSERT(expertDomain != NULL);

  curr_mem_alloc = max_mem_alloc = 0;
  a_size = add_size;
  address = (mdd_node_data *) malloc(a_size * sizeof(mdd_node_data));
  assert(NULL != address);
  updateMemoryAllocated(a_size * sizeof(mdd_node_data));
  memset(address, 0, a_size * sizeof(mdd_node_data));
  a_last = peak_nodes = a_unused = 0;
  
  l_size = l_add_size;
  level = (mdd_level_data *) malloc(l_size * sizeof(mdd_level_data));
  assert(NULL != level);
  updateMemoryAllocated(l_size * sizeof(mdd_level_data));
  memset(level, 0, l_size * sizeof(mdd_level_data));

  unique = new mdd_hash_table<node_manager> (this);
  curr_slots = max_slots = max_hole_chain = num_compactions = 0;
  active_nodes = zombie_nodes = orphan_nodes = temp_nodes = reclaimed_nodes = 0;

  // default policies
  performing_gc = false;
  nodes_activated_since_gc = 0;
#if 0
  enable_garbageCollection = false;  // default
#else
  enable_garbageCollection = true;
#endif
#if 1
  holeRecycling = true;
#else
  holeRecycling = false;
#endif

#if 1
  // atmost 100% of the data can be holes, effectively infinite
#if 0
  setCompactionThreshold(50);
#else
  setCompactionThreshold(40);
#endif

  // set level sizes
  assert(0 == setLevelBoundsAndHeights());
#endif

  counting = false;

  dptrsSize = 0;
  dptrs = 0;
}


#if 1


int node_manager::buildLevelNodeHelper(int lh, int* dptrs, int sz)
{
  DCASSERT(dptrs != 0);
  DCASSERT(sz > 0);

  int node = 0;

  if (isForRelations()) {
    // build from bottom up
    // for i = 1 to lh-1
    //    for j = 0 to sz
    //      dptrs[j] = node at level i with all downpointers to prev dptrs[j]
    //      do for primed first and then unprimed

    const int absLh = (lh < 0)? -lh: lh;
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
  else if (getReductionRule() == forest::QUASI_REDUCED) {
    DCASSERT(!isForRelations());
    // build from bottom up
    // for i = 1 to lh-1
    //    for j = 0 to sz
    //      dptrs[j] = node at level i with all downpointers to prev dptrs[j]

    for (int i = 1; i < lh; ++i)
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
  node = createTempNode(lh, sz, false);
  int* curr = getFullNodeDownPtrs(node);
  int* stop = curr + getFullNodeSize(node);
  // no need to link/unlink nodes since we pass the link
  // from dptrs[] to curr[]
  for ( ; curr != stop; ++curr, ++dptrs) { *curr = *dptrs; }
  node = reduceNode(node);

  // now build the levels above this node
  if (isForRelations()) {
    // build additional node at lh level if necessary
    if (lh < 0) {
      // build unprimed node at level ABS(lh) (same as -lh since lh < 0)
      int temp = createTempNodeMaxSize(-lh, false);
      setAllDownPtrsWoUnlink(temp, node);
      unlinkNode(node);
      node = reduceNode(temp);
    }

    // build primed and unprimed nodes for levels lh+1 to topLevel
    int topLevel = d->getTopVariable();
    for (int i = ABS(lh) + 1; i <= topLevel; ++i)
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
    // done building node for Relations
  }
  else if (getReductionRule() == forest::QUASI_REDUCED) {
    DCASSERT(!isForRelations());
    // build nodes for levels lh+1 to topLevel
    int topLevel = d->getTopVariable();
    for (int i = ABS(lh) + 1; i <= topLevel; ++i)
    {
      int temp = createTempNodeMaxSize(i, false);
      setAllDownPtrsWoUnlink(temp, node);
      unlinkNode(node);
      node = reduceNode(temp);
    }
    // done building node for (MT)MDDs
  }

  DCASSERT(isReducedNode(node));
  return node;
}


void node_manager::buildLevelNode(int lh, int* dptrs, int sz)
{
  DCASSERT(getLevelNode(lh) == 0);
  DCASSERT(dptrs != 0);
  DCASSERT(sz > 0);

  level[mapLevel(lh)].levelNode = buildLevelNodeHelper(lh, dptrs, sz);
  DCASSERT(getLevelNode(lh) != 0 && isReducedNode(getLevelNode(lh)));
}


int* node_manager::getTerminalNodes(int n)
{
  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    updateMemoryAllocated((n - dptrsSize) * sizeof(int));
    dptrsSize = n;
    dptrs = (int *) realloc(dptrs, dptrsSize * sizeof(int));
    DCASSERT(NULL != dptrs);
  }

  // store the terminals in the corresponding indexes
  switch (getRangeType()) {
    case BOOLEAN:
      DCASSERT(n == 2);
      dptrs[0] = getTerminalNode(false);
      dptrs[1] = getTerminalNode(true);
      break;
    case INTEGER:
      for (int i = 0; i < n; ++i)
        dptrs[i] = getTerminalNode(int(i));
      break;
    case REAL:
      for (int i = 0; i < n; ++i)
        dptrs[i] = getTerminalNode(float(i));
      break;
  }

  return dptrs;
}


int* node_manager::getTerminalNodes(int n, bool* terms)
{
  DCASSERT(n == 2);
  DCASSERT(getRangeType() == forest::BOOLEAN);
  DCASSERT(terms != 0);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    updateMemoryAllocated((n - dptrsSize) * sizeof(int));
    dptrsSize = n;
    dptrs = (int *) realloc(dptrs, dptrsSize * sizeof(int));
    DCASSERT(NULL != dptrs);
  }
  // fill array with terminal nodes
  for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(terms[i]);
  return dptrs;
}


int* node_manager::getTerminalNodes(int n, int* terms)
{
  DCASSERT(getRangeType() == forest::INTEGER);
  DCASSERT(terms != 0);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    updateMemoryAllocated((n - dptrsSize) * sizeof(int));
    dptrsSize = n;
    dptrs = (int *) realloc(dptrs, dptrsSize * sizeof(int));
    DCASSERT(NULL != dptrs);
  }
  // fill array with terminal nodes
  for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(terms[i]);
  return dptrs;
}


int* node_manager::getTerminalNodes(int n, float* terms)
{
  DCASSERT(getRangeType() == forest::REAL);
  DCASSERT(terms != 0);

  // use the array that comes with object (saves having to alloc/dealloc)
  if (dptrsSize < n) {
    // array not large enough, expand
    updateMemoryAllocated((n - dptrsSize) * sizeof(int));
    dptrsSize = n;
    dptrs = (int *) realloc(dptrs, dptrsSize * sizeof(int));
    DCASSERT(NULL != dptrs);
  }
  // fill array with terminal nodes
  for (int i = 0; i < n; ++i) dptrs[i] = getTerminalNode(terms[i]);
  return dptrs;
}


forest::error node_manager::createEdgeForVar(int vh, bool primedLevel,
    dd_edge& result)
{
  if (!isValidVariable(vh)) return forest::INVALID_VARIABLE;
  if (result.getForest() != this) return forest::INVALID_OPERATION;

  int k = primedLevel? -vh: vh;
  DCASSERT(isValidLevel(k));
  int node = getLevelNode(k);

  if (node == 0) {
    if (!isForRelations() && primedLevel) return forest::INVALID_ASSIGNMENT;
    if (getRangeType() == forest::BOOLEAN && getLevelSize(vh) > 2)
      return forest::INVALID_OPERATION;
    if (getEdgeLabeling() != forest::MULTI_TERMINAL)
      return forest::INVALID_OPERATION;
    int *terminalNodes = getTerminalNodes(getLevelSize(vh));
    buildLevelNode(k, terminalNodes, getLevelSize(vh));
    node = getLevelNode(k);
    DCASSERT(node != 0);
  }

  linkNode(node);
  result.set(node, 0, getNodeLevel(node));
  return forest::SUCCESS;
}


forest::error node_manager::createEdgeForVar(int vh, bool primedLevel,
    bool* terms, dd_edge& result)
{
  if (!isValidVariable(vh)) return forest::INVALID_VARIABLE;
  if (result.getForest() != this) return forest::INVALID_OPERATION;
  if (getRangeType() != forest::BOOLEAN) return forest::INVALID_OPERATION;
  if (getLevelSize(vh) != 2) return forest::INVALID_OPERATION;

  int k = primedLevel? -vh: vh;
  DCASSERT(isValidLevel(k));

  if (!isForRelations() && primedLevel) return forest::INVALID_ASSIGNMENT;
  if (getEdgeLabeling() != forest::MULTI_TERMINAL)
    return forest::INVALID_OPERATION;
  int *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
  int node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

  linkNode(node);
  result.set(node, 0, getNodeLevel(node));
  return forest::SUCCESS;
}


forest::error node_manager::createEdgeForVar(int vh, bool primedLevel,
    int* terms, dd_edge& result)
{
  if (!isValidVariable(vh)) return forest::INVALID_VARIABLE;
  if (result.getForest() != this) return forest::INVALID_OPERATION;
  if (getRangeType() != forest::INTEGER) return forest::INVALID_OPERATION;

  int k = primedLevel? -vh: vh;
  DCASSERT(isValidLevel(k));

  if (!isForRelations() && primedLevel) return forest::INVALID_ASSIGNMENT;
  if (getEdgeLabeling() != forest::MULTI_TERMINAL)
    return forest::INVALID_OPERATION;
  int *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
  int node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

  linkNode(node);
  result.set(node, 0, getNodeLevel(node));
  return forest::SUCCESS;
}


forest::error node_manager::createEdgeForVar(int vh, bool primedLevel,
    float* terms, dd_edge& result)
{
  if (!isValidVariable(vh)) return forest::INVALID_VARIABLE;
  if (result.getForest() != this) return forest::INVALID_OPERATION;
  if (getRangeType() != forest::REAL) return forest::INVALID_OPERATION;

  int k = primedLevel? -vh: vh;
  DCASSERT(isValidLevel(k));

  if (!isForRelations() && primedLevel) return forest::INVALID_ASSIGNMENT;
  if (getEdgeLabeling() != forest::MULTI_TERMINAL)
    return forest::INVALID_OPERATION;
  int *terminalNodes = getTerminalNodes(getLevelSize(vh), terms);
  int node = buildLevelNodeHelper(k, terminalNodes, getLevelSize(vh));

  linkNode(node);
  result.set(node, 0, getNodeLevel(node));
  return forest::SUCCESS;
}


#endif

int node_manager::setLevelBoundsAndHeights()
{
  // find biggest level handle
  const int* h2lMap = expertDomain->getHeightsToLevelsMap();
  const int* bounds = expertDomain->getLevelBounds();

  for (int i = expertDomain->getNumVariables(); i >= 1; --i) {
    if (setLevelBoundAndHeight(h2lMap[i], bounds[h2lMap[i]], i) < 0)
      return -1;
    if (isForRelations()) {
      // primed level
      if (setLevelBoundAndHeight(-h2lMap[i], bounds[h2lMap[i]], i) < 0)
        return -1;
    }
  }
  // success
  return 0;
}

int node_manager::setLevelBoundAndHeight(int k, int sz, int h)
{
  DCASSERT(k != 0);
  int mapped_k = mapLevel(k);
  if (mapped_k >= l_size) {
    // level doesn't exist, add additional levels
    int old_l_size = l_size;
    // l_size = l_size * 2;
    l_size = mapped_k + 2;
    level = (mdd_level_data *) realloc(level, l_size * sizeof(mdd_level_data));
    if (level == NULL) {
      fprintf(stderr, "Memory allocation error while allocating new level.\n");
      exit(1);
    }
    updateMemoryAllocated((l_size - old_l_size) * sizeof(mdd_level_data));
    // wipe new level data
    memset(level + old_l_size, 0,
        (l_size - old_l_size) * sizeof(mdd_level_data));
  }
  // level already defined
  if (level[mapped_k].data != NULL) return -1;
  
  DCASSERT(mapped_k < l_size);
  DCASSERT(level[mapped_k].data == NULL);
  
  level[mapped_k].size = add_size;
  level[mapped_k].data = (int *) malloc(level[mapped_k].size * sizeof(int));
  DCASSERT(NULL != level[mapped_k].data);
  if (level[mapped_k].data == NULL) return -1;
  updateMemoryAllocated(level[mapped_k].size * sizeof(int));
  memset(level[mapped_k].data, 0, level[mapped_k].size * sizeof(int));
  level[mapped_k].holes_top = level[mapped_k].holes_bottom =
    level[mapped_k].hole_slots =
    level[mapped_k].max_hole_chain = level[mapped_k].num_compactions = 0;
  level[mapped_k].last = 0;

  level[mapped_k].height = h;
  level[mapped_k].temp_nodes = 0;
  level[mapped_k].compactLevel = false;

  level[mapped_k].levelNode = 0;

  // success
  return 0;
}


void node_manager::setHoleRecycling(bool policy)
{
  if (policy == holeRecycling) return;
  if (policy) {
    // if we don't compact, some holes will not be tracked.
    compactAllLevels();
  } else {
    // trash hole tracking mechanism
    for (int i=0; i<l_size; i++) {
      level[i].holes_top = level[i].holes_bottom = 0;
    }
  }
  holeRecycling = policy;
}


node_manager::~node_manager()
{
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
  free(address);
#if DEBUG_DELETE_NM
  printf("Deleted address[]\n");
  fflush(stdout);
#endif
  for (int i = 0; i < l_size; i++) { free(level[i].data); }
#if DEBUG_DELETE_NM
  printf("Deleted level[i].data\n");
  fflush(stdout);
#endif
  free(level);
#if DEBUG_DELETE_NM
  printf("Deleted level\n");
  fflush(stdout);
#endif
}

// *********************************************************************
// TODO: test this out
int node_manager::buildQuasiReducedNodeAtLevel(int k, int p)
{
  DCASSERT(reductionRule == forest::QUASI_REDUCED);
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


int node_manager::getMddLevelMaxBound(int k) const
{
  // Go through each node in this level and check its size - in terms of 
  // # of downpointers. If it has more downpointers that the current max_bound
  // update max_bound. After going throug all nodes, return max_bound.
  int mapped_k = mapLevel(k);
  if (mapped_k <= 0 || mapped_k >= l_size) return 0;
  if (level[mapped_k].data == NULL) return 0;
  mdd_level_data* l_info = &level[mapped_k];
  int* data = l_info->data;

  int max_bound = 0;
  for (int a=1; a<l_info->last; ) {
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

int node_manager::getMxdLevelMaxBound(int k) const
{
  return MAX(getMddLevelMaxBound(k), getMddLevelMaxBound(-k));
}

int node_manager::getLevelMaxBound(int k) const
{
  return isForRelations()?
            getMxdLevelMaxBound(k):
            getMddLevelMaxBound(k);
}

int digits(int a) 
{
  int d = 1;
  while (a) { d++; a/=10; }
  return d;
}

void node_manager::dump(FILE *s) const
{
  int nwidth = digits(a_last);
  for (int p=0; p<=a_last; p++) {
    fprintf(s, "%*d\t", nwidth, p);
    showNode(s, p, 1);
    fprintf(s, "\n");
    fflush(s);
  }
}

void node_manager::showAll() const
{
  dumpInternal(stderr);
}

void node_manager::dumpInternal(FILE *s) const
{
  fprintf(s, "Internal forest storage\n");
  fprintf(s, "First unused node index: %d\n", a_unused);
  int awidth = digits(a_last);
  fprintf(s, " Node# :  ");
  for (int p=1; p<=a_last; p++) {
    if (p) fprintf(s, " ");
    fprintf(s, "%*d", awidth, p);
  }
  fprintf(s, "\nLevel  : [");
  for (int p=1; p<=a_last; p++) {
    if (p) fprintf(s, "|");
    fprintf(s, "%*d", awidth, address[p].level);
  }
  fprintf(s, "]\n");
  fprintf(s, "\nOffset : [");
  for (int p=1; p<=a_last; p++) {
    if (p) fprintf(s, "|");
    fprintf(s, "%*d", awidth, address[p].offset);
  }
  fprintf(s, "]\n");
  fprintf(s, "\nCache  : [");
  for (int p=1; p<=a_last; p++) {
    if (p) fprintf(s, "|");
    fprintf(s, "%*d", awidth, address[p].cache_count);
  }
  fprintf(s, "]\n\n");

  for (int i=1; i<l_size; i++) {
    dumpInternalLevel(s, unmapLevel(i));
  }
  
  fprintf(s, "Uniqueness table:\n");
  unique->show(s); 
  fflush(s);
}

void node_manager::dumpInternalLevel(FILE *s, int k) const
{
  int p_level = mapLevel(k);
  mdd_level_data* l_info = &level[p_level];
  int* data = l_info->data;

  if (data == NULL) return; // nothing to display
  
  fprintf(s, "Level %d: ", k);
  fprintf(s, "Height %d: ", l_info->height);
  fprintf(s, "Last slot used: %d\n", l_info->last);
  fprintf(s, "Grid: top = %d bottom = %d\n",
      l_info->holes_top, l_info->holes_bottom);

  fprintf(s, "Data array by record: \n");
  int awidth = digits(a_last);
  int a;
  for (a=1; a<=l_info->last; ) {
    fflush(s);
    fprintf(s, "%*d : [%d", awidth, a, data[a]);
    for (int i=1; i<3; i++) {
      fprintf(s, "|%d", data[a+i]);
    }
    if (data[a]<0) { 
      // hole
      fprintf(s, "| ... ");
      a -= data[a];  
    } else {
      // proper node
      if (data[a+2]>0) {
        // Full storage
        if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
          for (int i=0; i<2*data[a+2]; i++) {
            fprintf(s, "|%d", data[a+3+i]);
          }
          a += 4 + 2 * data[a + 2];
        } else {
          for (int i=0; i<data[a+2]; i++) {
            fprintf(s, "|%d", data[a+3+i]);
          }
          a += 4 + data[a + 2];
        }
      } else {
        // sparse storage
        if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
          for (int i=0; i<-3*data[a+2]; i++) {
            fprintf(s, "|%d", data[a+3+i]);
          }
          a += 4 - 3 * data[a + 2];
        } else {
          for (int i=0; i<-2*data[a+2]; i++) {
            fprintf(s, "|%d", data[a+3+i]);
          }
          a += 4 - 2 * data[a + 2];
        }
      }
    }
    fprintf(s, "|%d]\n", data[a-1]);
  } // for a
  fprintf(s, "%*d : free slots\n", awidth, a);
  fflush(s);
  DCASSERT(a == ((l_info->last)+1));
}

double cardinality(const node_manager* const nm, int p, int k,
    std::map<int, double>& visited)
{
  int pLevel = nm->getNodeLevel(p);
  // fprintf(stderr, "p: %d, k: %d, pLevel: %d\n", p, k, pLevel);
  DCASSERT(pLevel <= k);

  if (pLevel < k)
    return nm->getLevelSize(k) * cardinality(nm, p, k-1, visited);

  if (nm->isTerminalNode(p)) {
    return p == 0? 0: 1;
  }

  // check in cache
  std::map<int, double>::iterator curr = visited.find(p);
  if (curr != visited.end())
    return curr->second;

  // not found in cache; compute by exanding
  double result = 0;
  if (nm->isFullNode(p)) {
    const int sz = nm->getFullNodeSize(p);
    for (int i = 0; i < sz; ++i)
      result += cardinality(nm, nm->getFullNodeDownPtr(p, i), k-1, visited);
  }
  else {
    const int sz = nm->getSparseNodeSize(p);
    for (int i = 0; i < sz; ++i)
      result += cardinality(nm, nm->getSparseNodeDownPtr(p, i), k-1, visited);
  }

  // save result and return
  visited[p] = result;
  // fprintf(stderr, "p: %d, k: %d, pLevel: %d, card: %f\n",
  //    p, k, pLevel, result);
  return result;
}

double node_manager::getCardinality(int node) const
{
  if (!(isMdd() || isMtMdd())) return 0;
  std::map<int, double> visited;
#if 0
  return getCardinality(node, d->getTopVariable(), visited);
#else
  return cardinality(this, node, d->getTopVariable(), visited);
#endif
}

#if 0
double node_manager::getCardinality(int p, int k, std::map<int, double>& visited)
  const
{
  int pLevel = getNodeLevel(p);
  // fprintf(stderr, "p: %d, k: %d, pLevel: %d\n", p, k, pLevel);
  DCASSERT(pLevel <= k);

  if (pLevel < k)
    return getLevelSize(k) * getCardinality(p, k-1, visited);

  if (isTerminalNode(p)) {
    return p == 0? 0: 1;
  }

  // check in cache
  std::map<int, double>::iterator curr = visited.find(p);
  if (curr != visited.end())
    return curr->second;

  // not found in cache; compute by exanding
  double result = 0;
  if (isFullNode(p)) {
    const int sz = getFullNodeSize(p);
    for (int i = 0; i < sz; ++i)
      result += getCardinality(getFullNodeDownPtr(p, i), k-1, visited);
  }
  else {
    const int sz = getSparseNodeSize(p);
    for (int i = 0; i < sz; ++i)
      result += getCardinality(getSparseNodeDownPtr(p, i), k-1, visited);
  }

  // save result and return
  visited[p] = result;
  // fprintf(stderr, "p: %d, k: %d, pLevel: %d, card: %f\n",
  //    p, k, pLevel, result);
  return result;
}
#endif

void node_manager::showNodeGraph(FILE *s, int p) const
{
  std::vector< std::set<int> >
    discovered(mapLevel(expertDomain->getTopVariable()) + 1);
  std::queue<int> toExpand;

  toExpand.push(p);
  discovered[getNodeLevelMapping(p)].insert(p);

  // expand the front of toExpand;
  // add newly discovered ones to discovered and toExpand

  while (!toExpand.empty()) {
    int p = toExpand.front();
    toExpand.pop();
    if (isTerminalNode(p)) continue;
    // expand
    if (isFullNode(p)) {
      const int sz = getFullNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getFullNodeDownPtr(p, i);
        int k = getNodeLevelMapping(temp);
        // insert into discovered and toExpand if new
        if (discovered[k].find(temp) == discovered[k].end()) {
          toExpand.push(temp);
          discovered[k].insert(temp);
        }
      }
    }
    else {
      const int sz = getSparseNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getSparseNodeDownPtr(p, i);
        int k = getNodeLevelMapping(temp);
        // insert into discovered and toExpand if new
        if (discovered[k].find(temp) == discovered[k].end()) {
          toExpand.push(temp);
          discovered[k].insert(temp);
        }
      }
    }
  }

  // iterate through discovered and print
  for (unsigned i = discovered.size() - 1; i > 0u; i--)
  {
    if (discovered[i].empty()) continue;
    int k = unmapLevel(i);
    fprintf(s, "Level: %d%s\n", ABS(k), (k < 0? "'": " "));
    for (std::set<int>::iterator iter = discovered[i].begin();
        iter != discovered[i].end(); iter++)
    {
      fprintf(s, "  ");
      showNode(s, *iter);
      fprintf(s, "\n");
    }
  }
}


int ifTermGetInt(const node_manager *nm, int node)
{
  return nm->isTerminalNode(node) ? nm->getInteger(node) : node;
}

void node_manager::showNode(FILE *s, int p, int verbose) const
{
  if (isTerminalNode(p)) {
    fprintf(s, "(terminal)");
    return;
  }
  if (isDeletedNode(p)) {
    fprintf(s, "DELETED");
    return;
  }
  int a = getNodeOffset(p);
  int l = getNodeLevelMapping(p);
#if 0
  int p_width = digits(a_last);
  int l_width = digits(l_size);
#endif
  int* data = level[l].data;
  if (verbose) {
    fprintf(s, " level: %d", ABS(unmapLevel(l)));
    if (getNodeLevel(p) < 0)
      fprintf(s, "'");
    else
      fprintf(s, " ");
    fprintf(s, " in: %d", data[a]);
    fprintf(s, " cc: %d", address[p].cache_count);
  } else {
    fprintf(s, "node: %d", p);
  }
  if (isSparseNode(p)) {
    // sparse
    if (verbose)
      fprintf(s, " nnz : %d", getSparseNodeSize(p));
    fprintf(s, " down: (");
    for (int z=0; z<getSparseNodeSize(p); z++) {
      if (z) fprintf(s, ", ");
      if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
        int e = getSparseNodeEdgeValue(p, z);
        if (e == INF) {
          fprintf(s, "%d:INF:%d",
              getSparseNodeIndex(p, z),
              getSparseNodeDownPtr(p, z));
        } else {
          fprintf(s, "%d:%d:%d",
              getSparseNodeIndex(p, z),
              e,
              getSparseNodeDownPtr(p, z));
        }
      } else {
        if (isTerminalNode(getSparseNodeDownPtr(p, z))) {
          fprintf(s, "%d:", getSparseNodeIndex(p, z));
          if (getRangeType() == forest::REAL) {
            fprintf(s, "%f", getReal(getSparseNodeDownPtr(p, z)));
          } else if (getRangeType() == forest::INTEGER) {
            fprintf(s, "%d", getInteger(getSparseNodeDownPtr(p, z)));
          } else {
            assert(getRangeType() == forest::BOOLEAN);
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
      if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
        int e = getFullNodeEdgeValue(p, i);
        if  (e == INF) {
          fprintf(s, "INF:%d",
              getFullNodeDownPtr(p, i));
        } else {
          fprintf(s, "%d:%d", e,
              getFullNodeDownPtr(p, i));
        }
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

void node_manager::showNode(int p) const
{
  assert(edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES);
  if (isTerminalNode(p)) {
    fprintf(stderr, "(terminal)");
    return;
  }
  if (isDeletedNode(p)) {
    fprintf(stderr, "DELETED");
    return;
  }
  int a = getNodeOffset(p);
  int l = getNodeLevelMapping(p);
  int p_width = digits(a_last);
  int l_width = digits(l_size);
  int* data = level[l].data;
  fprintf(stderr, "node: %*d level: %*d", p_width, p, l_width, ABS(unmapLevel(l)));
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

void node_manager::compactLevel(int k)
{
  if (k == 0) { level[0].compactLevel = false; return; }
  // mapped p_level
  int p_level = mapLevel(k);
  CHECK_RANGE(0, p_level, l_size);
  if (0 == level[p_level].hole_slots ||  // Already compact
      !doesLevelNeedCompaction(k)) {  // Level is compact enough!
    level[p_level].compactLevel = false;
#if 0
    printf("%s: level %d... compact enough\n", __func__, k);
#endif
    return;
  }

  if (0 < level[p_level].temp_nodes) return;   // Temp nodes; do not compact
#if 0
  printf("%s: level %d\n", __func__, k);
#endif

#if 0
  printf("Before compaction:\n");
  dumpInternalLevel(stdout, k);
  printf("\n");
#endif

  // alternate algorithm -- since we now have the node ids in the node data
  int *node_ptr = level[p_level].data + 1;  // since we leave [0] empty
  int *end_ptr = level[p_level].data + level[p_level].last + 1;
  int *curr_ptr = node_ptr;
  int node_size = 0;
  int curr_node = 0;

  if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
    while (node_ptr != end_ptr) {
      // find new node
      if (*node_ptr < 0) {
        // found a hole, advance
        assert(node_ptr[0] == node_ptr[-(*node_ptr)-1]);
        node_size = -(*node_ptr);
        memset(node_ptr, 0, node_size * sizeof(int));
      } else {
        // found an existing node
        if (isPessimistic()) assert(*node_ptr != 0);

        node_size = *(node_ptr + 2);  // [2] = size
        assert (node_size != 0);      // assuming zombies have been deleted
        if (node_size < 0) {
          // sparse
          node_size = 4 - (3 * node_size);
        } else {
          // full
          node_size = 4 + (2 * node_size);
        }
        curr_node = node_ptr[node_size - 1];
        assert(getNodeOffset(curr_node) == (node_ptr - level[p_level].data));
        if (node_ptr != curr_ptr) {
#if 1
          for (int i = 0; i < node_size; ++i) {
            curr_ptr[i] = node_ptr[i];
            node_ptr[i] = 0;
          }
#else
          // copy node_ptr to curr_ptr
          memmove(curr_ptr, node_ptr, node_size * sizeof(int));
#endif
          // change node offset
          address[curr_node].offset = (curr_ptr - level[p_level].data);
        }
        assert(getNodeOffset(curr_node) == (curr_ptr - level[p_level].data));
        curr_ptr += node_size;
      }
      node_ptr += node_size;
    }
  } else {
    while (node_ptr != end_ptr) {
      // find new node
      if (*node_ptr < 0) {
        // found a hole, advance
        assert(node_ptr[0] == node_ptr[-(*node_ptr)-1]);
        node_size = -(*node_ptr);
        memset(node_ptr, 0, node_size * sizeof(int));
      } else {
        // found an existing node
        if (isPessimistic()) assert(*node_ptr != 0);

        node_size = *(node_ptr + 2);  // [2] = size
        assert (node_size != 0);      // assuming zombies have been deleted
        if (node_size < 0) {
          // sparse
          node_size = 4 - (2 * node_size);
        } else {
          // full
          node_size = 4 + node_size;
        }
        curr_node = node_ptr[node_size - 1];
        assert(getNodeOffset(curr_node) == (node_ptr - level[p_level].data));
        if (node_ptr != curr_ptr) {
#if 1
          for (int i = 0; i < node_size; ++i) {
            curr_ptr[i] = node_ptr[i];
            node_ptr[i] = 0;
          }
#else
          // copy node_ptr to curr_ptr
          memmove(curr_ptr, node_ptr, node_size * sizeof(int));
#endif
          // change node offset
          address[curr_node].offset = (curr_ptr - level[p_level].data);
        }
        assert(getNodeOffset(curr_node) == (curr_ptr - level[p_level].data));
        curr_ptr += node_size;
      }
      node_ptr += node_size;
    }
  }
  curr_slots -= level[p_level].last;
  level[p_level].last = (curr_ptr - 1 - level[p_level].data);
  curr_slots += level[p_level].last;
  if (max_slots < curr_slots) max_slots = curr_slots;

  // set up hole pointers and such
  level[p_level].holes_top = level[p_level].holes_bottom = 0;
  level[p_level].hole_slots = 0;

  num_compactions++;
  level[p_level].num_compactions++;
  level[p_level].compactLevel = false;

  if (level[p_level].size > add_size &&
      level[p_level].last < level[p_level].size/2) {
    int new_size = level[p_level].size/2;
    while (new_size > add_size && new_size > level[p_level].last * 3)
    { new_size /= 2; }
    updateMemoryAllocated((new_size - level[p_level].size) * sizeof(int));
    level[p_level].data = (int *)
      realloc(level[p_level].data, new_size * sizeof(int));
    assert(NULL != level[p_level].data);
    level[p_level].size = new_size;
#ifdef MEMORY_TRACE
    printf("Reduced data[] by a factor of 2. New size: %d, Last: %d.\n",
        level[p_level].size, level[p_level].last);
#endif
  }

#if 0
  printf("After compaction:\n");
  dumpInternalLevel(stdout, k);
  printf("\n");
#endif
}

void node_manager::compactAllLevels()
{
  for (int i=0; i<l_size; i++) {
    level[i].compactLevel = true;
    compactLevel(unmapLevel(i));
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
#if ALT_HASH_CALL
unsigned node_manager::hash(int h) const 
#else
unsigned node_manager::hash(int h, unsigned n) const 
// when enabling this make sure you enable the %n in this function
#endif
{
  uint32_t a, b, c;
  int* k = getNodeAddress(h);
  int length = k[2];
  DCASSERT(length != 0);

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
    // a = b = c = 0xdeadbeef;
    // a = b = c = 0xdeadbeef + uint32_t(nnz)<<2;
    // a = b = c = 0xdeadbeef + uint32_t(nnz)<<2 + uint32_t(getNodeLevel(h));
    // a = b = c = 0xdeadbeef + uint32_t(nnz)<<16 + uint32_t(getNodeLevel(h));
#if 1
    a = 0xdeadbeef;
    b = uint32_t(nnz)<<16;
    c = uint32_t(getNodeLevel(h));
#else
    a = b = c = 0;
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
#else
    a = b = c = 0;
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
              final(a,b,c);
      case 0: // nothing left to add
              break;
    }
  }

  // report the result
#if ALT_HASH_CALL
  return c;
#else
  return c % n;
#endif
}


bool node_manager::equals(int h1, int h2) const 
{
  DCASSERT(h1);	
  DCASSERT(h2);
  DCASSERT(isActiveNode(h1));
  DCASSERT(isActiveNode(h2));
  DCASSERT(!isTerminalNode(h1));
  DCASSERT(!isTerminalNode(h2));

  if (getNodeLevel(h1) != getNodeLevel(h2)) return false;
  int *ptr1 = getNodeAddress(h1);
  int *ptr2 = getNodeAddress(h2);
  int sz1 = ptr1[2];
  int sz2 = ptr2[2];
  ptr1 += 3;
  ptr2 += 3;
  if (sz1<0 && sz2<0) {
    // both sparse
    if (sz1 != sz2) return false;
    if (edgeLabel == forest::MULTI_TERMINAL) {
      return 0==memcmp(ptr1, ptr2, -2*sz1*sizeof(int));
    } else {
      return 0==memcmp(ptr1, ptr2, -3*sz1*sizeof(int));
    }
  }
  if (sz1>0 && sz2>0) {
    // both full
    int ms = MIN(sz1, sz2);
    if (memcmp(ptr1, ptr2, ms*sizeof(int))) return false;
    // tails must be zero.  Only one loop will go.
    for (int i=ms; i<sz1; i++)
      if (ptr1[i]) return false; 
    for (int i=ms; i<sz2; i++)
      if (ptr2[i]) return false;
    if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
      // check edge-values
      if (memcmp(ptr1+sz1, ptr2+sz2, ms*sizeof(int))) return false;
    }
    return true;
  }
  if (sz1>0) {
    // node1 is full, node2 is sparse
    int* down2 = ptr2 - sz2;
    int i = 0;
    if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
      int* edge2 = down2 - sz2;
      for (int z=0; z<-sz2; z++) {
        for (; i<ptr2[z]; i++)  // node1 must be zeroes in between
          if (ptr1[i]) return false;
        if (ptr1[i] != down2[z]) return false;
        if (ptr1[i+sz1] != edge2[z]) return false;
        i++;
      }
    } else {
      for (int z=0; z<-sz2; z++) {
        for (; i<ptr2[z]; i++)  // node1 must be zeroes in between
          if (ptr1[i]) return false;
        if (ptr1[i] != down2[z]) return false;
        i++;
      }
    }
    // tail of node1
    for (; i<sz1; i++)
      if (ptr1[i]) return false;
    return true;
  }
  // node2 is full, node1 is sparse
  int* down1 = ptr1 - sz1;
  int i = 0;
  if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
    int* edge1 = down1 - sz1;
    for (int z=0; z<-sz1; z++) {
      for (; i<ptr1[z]; i++)  // node2 must be zeroes in between
        if (ptr2[i]) return false;
      if (ptr2[i] != down1[z]) return false;
      if (ptr2[i+sz2] != edge1[z]) return false;
      i++;
    }
  } else {
    for (int z=0; z<-sz1; z++) {
      for (; i<ptr1[z]; i++)  // node2 must be zeroes in between
        if (ptr2[i]) return false;
      if (ptr2[i] != down1[z]) return false;
      i++;
    }
  }
  // tail of node2
  for (; i<sz2; i++)
    if (ptr2[i]) return false;
  return true;
}

// ------------------------------------------------------------------
//  Protected methods
// ------------------------------------------------------------------

void node_manager::deleteTempNode(int p)
{

#if 0
  validateIncounts();
#endif

  int *p_data = getNodeAddress(p);
  int k = address[p].level;

  // 0: incount
  // 1: temp_node
  // 2: size (full_node)
  // 3..: downpointers (children)
  // last: p

  DCASSERT(address[p].cache_count == 0);
  DCASSERT(p_data[1] == temp_node);
  DCASSERT(p_data[0] == 1);
  DCASSERT(p_data[2] > 0);  // full node

  // unlinkNode children

  for (int *curr = p_data + 3; curr < p_data + 3 + p_data[2]; ++curr) {
#if ENABLE_IN_COUNTING
    int temp = *curr;
    *curr = 0;
    unlinkNode(temp);
#else
    unlinkNode(*curr);
#endif
  }
  // make hole
  if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
    makeHole(getNodeLevel(p), getNodeOffset(p), 4 + 2 * p_data[2]);
  } else {
    makeHole(getNodeLevel(p), getNodeOffset(p), 4 + p_data[2]);
  }
  // reclaim node id
  freeNode(p);

  if (level[mapLevel(k)].compactLevel) compactLevel(k);

#if 0
  validateIncounts();
#endif

}


void node_manager::deleteNode(int p)
{
  DCASSERT(!isTerminalNode(p));
  DCASSERT(getNext(p) != temp_node);
  DCASSERT(getInCount(p) == 0);
  DCASSERT(isReducedNode(p));   // it's in the unique table

#if 0
  validateIncounts();
#endif

  int* foo = getNodeAddress(p);
  int k = getNodeLevel(p);

  // remove from unique table
  assert(unique->find(p) == p);
#ifdef TRACK_DELETIONS
  showNode(stdout, p);
  int x = unique->remove(p);
  printf("%s: p = %d, unique->remove(p) = %d\n", __func__, p, x);
  fflush(stdout);
  DCASSERT(x==p);
#else
  int x = unique->remove(p);
  assert(x != -1);
  assert(p == x);
#endif
  assert(address[p].cache_count == 0);

  // unlinkNode children
  if (foo[2]<0) {
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
    if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
      makeHole(getNodeLevel(p), getNodeOffset(p), 4 - 3 * foo[2]);  
    } else {
      makeHole(getNodeLevel(p), getNodeOffset(p), 4 - 2 * foo[2]);  
    }
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
    if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
      makeHole(getNodeLevel(p), getNodeOffset(p), 4 + 2 * foo[2]);  
    } else {
      makeHole(getNodeLevel(p), getNodeOffset(p), 4 + foo[2]);  
    }
  }

  // recycle the index
  freeNode(p);

  if (level[mapLevel(k)].compactLevel) compactLevel(k);

#if 0
  validateIncounts();
#endif

}

#if 1

void node_manager::zombifyNode(int p)
{
  DCASSERT(isActiveNode(p));
  DCASSERT(!isTerminalNode(p));
  DCASSERT(isReducedNode(p));
  DCASSERT(getCacheCount(p) > 0);  // otherwise this node should be deleted
  DCASSERT(getInCount(p) == 0);

  assert(address[p].cache_count > 0);
# if 0
  assert(zombie_nodes == 0);
#endif

  zombie_nodes++;
  level[getNodeLevelMapping(p)].zombie_nodes++;
  active_nodes--;

  // mark node as zombie
  assert(unique->find(p) == p);
  // tag p to be a zombie node
  address[p].cache_count = -address[p].cache_count;

#ifdef DEVELOPMENT_CODE 
  int x = unique->remove(p);
  DCASSERT(x==p);
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
    if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
      makeHole(node_level, node_offset, 4 - 3 * foo[2]);  
    } else {
      makeHole(node_level, node_offset, 4 - 2 * foo[2]);  
    }
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
    if (edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES) {
      makeHole(node_level, node_offset, 4 + 2 * foo[2]);  
    } else {
      makeHole(node_level, node_offset, 4 + foo[2]);  
    }
  }
}

#endif  // zombifyNode


forest::error node_manager::garbageCollect() {
  gc();
  return forest::SUCCESS;
}


bool node_manager::gc() {
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
    printf("Zombie nodes: %d\n", zombie_nodes);
#endif
    // remove the stale nodes entries from caches
    smart_cast<expert_compute_manager *>(MEDDLY_getComputeManager())->
      removeStales();
#ifdef DEBUG_GC
    printf("Zombie nodes: %d\n", zombie_nodes);
#endif
    assert(zombie_nodes == 0);
    freed_some = true;
  } else {
#ifdef ALT_ORPHAN_GC
    int org_orphan_nodes = orphan_nodes;
    int reps = 0;

    bool saved_gc_state = enable_garbageCollection;
    enable_garbageCollection = true;

    do {
      reps++;

#ifdef DEBUG_GC
      // count the number of orphaned nodes
      fprintf(stderr, "Active: %d, Orphan: %d\n", active_nodes, orphan_nodes);
#endif
      // remove the stale nodes entries from caches
      smart_cast<expert_compute_manager *>(MEDDLY_getComputeManager())->
        removeStales();
    } while (isTimeToGc());

    // return isStale to usual policy
    enable_garbageCollection = saved_gc_state;
    if (orphan_nodes < org_orphan_nodes) freed_some = true;

#ifdef DEBUG_GC
    fprintf(stderr,"  done GC, loop repeated %d times.\n", reps);
    fflush(stderr);
#endif

#else

#ifdef DEBUG_GC
    fprintf(stderr, "Active: %d, Zombie: %d, Orphan: %d\n",
        active_nodes, zombie_nodes, orphan_nodes);
#endif

    // convert orphans to zombies
    nodeDeletionPolicy = forest::PESSIMISTIC_DELETION;
    orphan_nodes = 0;
    for (int i = 1; i <= a_last; i++) {
      DCASSERT(!isTerminalNode(i));
      if (isActiveNode(i) && getInCount(i) == 0) {
        assert(getCacheCount(i) > 0);
        zombifyNode(i);
      }
    }
#ifdef DEBUG_GC
    fprintf(stderr, "Active: %d, Zombie: %d, Orphan: %d\n",
        active_nodes, zombie_nodes, orphan_nodes);
#endif
    // remove the stale nodes entries from caches
    smart_cast<expert_compute_manager *>(MEDDLY_getComputeManager())->
      removeStales();
    assert(zombie_nodes == 0);
    nodeDeletionPolicy = forest::OPTIMISTIC_DELETION;
    freed_some = true;
#endif
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


void node_manager::removeZombies(int max_zombies) {
#if 1
  return;
#else
  // too many zombies? kill em!
  if (zombie_nodes > max_zombies && active_nodes/zombie_nodes < 3) {
#if 0
    for (int i = 1; i <= a_last; i++) {
      if (isActiveNode(i) && isZombieNode(i)) {
        showNode(stdout, i); printf("\n");
      }
    }
#endif
    // remove the stale nodes entries from caches
    smart_cast<expert_compute_manager *>(MEDDLY_getComputeManager())->
      removeStales();
#if 0
    if (zombie_nodes > 0) {
      for (int i = 1; i <= a_last; i++) {
        if (isActiveNode(i) && isZombieNode(i)) {
          showNode(stdout, i); printf("\n");
        }
      }
    }
#endif
    assert(zombie_nodes == 0);
  }
#endif
}

int node_manager::getFreeNode(int k)
{
  if (a_unused) {
    // grab a recycled index
    int p = a_unused;
    DCASSERT(address[p].offset < 1);
    a_unused = -address[p].offset;
    active_nodes++;
    return p;
  }
  // new index
  // TODO:
  if (a_last + 1 >= a_size) {
    // compactLevel(k);
    // if (a_last + 1 < a_size) return getFreeNode(k);
    // int new_a_size = (a_size > 16384)? a_size + 16384: a_size * 2;
    int new_a_size = (a_size > 1024)? a_size + 1024: a_size * 2;
    mdd_node_data *temp = (mdd_node_data*) realloc(address,
        new_a_size * sizeof(mdd_node_data));
    // assert(NULL != temp);
    if (NULL == temp) {
      fprintf(stderr, "Memory allocation error while allocating MDD nodes.\n");
      exit(1);
    }
    updateMemoryAllocated((new_a_size - a_size) * sizeof(mdd_node_data));
    address = temp;
    memset(address + a_size, 0, (new_a_size - a_size) * sizeof(mdd_node_data));
    a_size = new_a_size;
  }
  a_last++;
  active_nodes++;
  if (getCurrentNumNodes() > peak_nodes) peak_nodes = getCurrentNumNodes();
  return a_last;
}

void node_manager::freeZombieNode(int p)
{
  DCASSERT(address[p].level != 0);
  DCASSERT(address[p].cache_count == 0);
  zombie_nodes--;
  level[mapLevel(address[p].level)].zombie_nodes--;
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
}

void node_manager::freeNode(int p)
{
#ifdef TRACE_REDUCE
  printf("%s: p = %d, a_last = %d\n", __func__, p, a_last);
#endif

  DCASSERT(!isTerminalNode(p));
  DCASSERT(!isPessimistic() || !isZombieNode(p));
  DCASSERT(address[p].cache_count == 0);

  active_nodes--;

  address[p].level = 0;
  address[p].cache_count = 0;
  if (p == a_last) { 
    // special case
    address[p].offset = 0;
    a_last--;
    if (a_size > add_size && a_last < a_size/2) {
      address = (mdd_node_data *)
          realloc(address, a_size/2 * sizeof(mdd_node_data));
      assert(NULL != address);
      a_size /= 2;
      updateMemoryAllocated(-a_size * sizeof(mdd_node_data));
#ifdef MEMORY_TRACE
      printf("Reduced node[] by a factor of 2. New size: %d.\n", a_size);
#endif
    }
  } else {
    address[p].offset = -a_unused;
    a_unused = p;
  }
}

void node_manager::gridInsert(int k, int p_offset)
{
  // sanity check to make sure that the first and last slots in this hole
  // have the same value, i.e. -(# of slots in the hole)
  int p_level = mapLevel(k);
  mdd_level_data* l_info = &level[p_level];
  int* l_data = l_info->data;
  DCASSERT(l_data[p_offset] == l_data[p_offset - l_data[p_offset] - 1]);
  // special case: empty
  if (0 == l_info->holes_bottom) {
    // index hole
    l_data[p_offset + 1] = l_data[p_offset + 2] = l_data[p_offset + 3] = 0;
    l_info->holes_top = l_info->holes_bottom = p_offset;
    return;
  }
  // special case: at top
  if (l_data[p_offset] < l_data[l_info->holes_top]) {
    // index hole
    l_data[p_offset + 1] = l_data[p_offset + 3] = 0;
    l_data[p_offset + 2] = l_info->holes_top;
    l_data[l_info->holes_top + 1] = p_offset;
    l_info->holes_top = p_offset;
    return;
  }
  int chain = 0;
  int above = l_info->holes_bottom;
  int below = 0;
  while (l_data[p_offset] < l_data[above]) {
    below = above;
    above = l_data[below + 1];
    chain++;
    DCASSERT(l_data[above + 2] == below);
    DCASSERT(above);  
  }
  l_info->max_hole_chain = MAX(l_info->max_hole_chain, chain);
  max_hole_chain = MAX(max_hole_chain, chain);
  if (l_data[p_offset] == l_data[above]) {
    // Found, add this to chain
    // making a non-index hole
    int right = l_data[above + 3];
    l_data[p_offset + 1] = non_index_hole;
    l_data[p_offset + 2] = above;
    l_data[p_offset + 3] = right;
    if (right) l_data[right + 2] = p_offset;
    l_data[above + 3] = p_offset;
    return; 
  }
  // we should have above < p_offset < below  (remember, -sizes)
  // create an index hole since there were no holes of this size
  l_data[p_offset + 1] = above;
  l_data[p_offset + 2] = below;
  l_data[p_offset + 3] = 0;
  l_data[above + 2] = p_offset;
  if (below) {
    l_data[below + 1] = p_offset;
  } else {
    DCASSERT(above == l_info->holes_bottom);
    l_info->holes_bottom = p_offset;
  }
}

void node_manager::indexRemove(int k, int p_offset)
{
#ifdef MEMORY_TRACE
  cout << __func__ << "(" << k << ", " << p_offset << ")\n";
#endif

  int p_level = mapLevel(k);
  DCASSERT(p_level >= 0);
  DCASSERT(!isHoleNonIndex(k, p_offset));
  int above = level[p_level].data[p_offset + 1];
  int below = level[p_level].data[p_offset + 2];
  int right = level[p_level].data[p_offset + 3];

  if (right >= 1) {
    // there are nodes to the right!
    DCASSERT(level[p_level].data[right + 1] < 0);
    level[p_level].data[right + 1] = above;
    level[p_level].data[right + 2] = below;

    // update the pointers of the holes (index) above and below it
    if (above) {
      level[p_level].data[above + 2] = right;
    } else {
      level[p_level].holes_top = right;
    }

    if (below) {
      level[p_level].data[below + 1] = right;
    } else {
      level[p_level].holes_bottom = right;
    }
    
  } else {
    // there are no non-index nodes
    DCASSERT(right < 1);

    // this was the last node of its size
    // update the pointers of the holes (index) above and below it
    if (above) {
      level[p_level].data[above + 2] = below;
    } else {
      level[p_level].holes_top = below;
    }

    if (below) {
      level[p_level].data[below + 1] = above;
    } else {
      level[p_level].holes_bottom = above;
    }
  }
}

int node_manager::getHole(int k, int slots, bool search_holes)
{
  const int min_node_size = 5;
  int p_level = mapLevel(k);
  DCASSERT(p_level >= 0);
#if 0
  DCASSERT(p_level < 50);
#endif
  DCASSERT(level[p_level].data != NULL);
  DCASSERT(slots >= min_node_size);
  DCASSERT((slots - min_node_size) + 1 <= (getLevelSize(k) *
    ((edgeLabel == forest::EVPLUS || edgeLabel == forest::EVTIMES)? 2: 1)));
  DCASSERT(0 < (slots - min_node_size) + 1);
  CHECK_RANGE(0, p_level, l_size);

  curr_slots += slots;
  if (max_slots < curr_slots) max_slots = curr_slots;

#ifdef DEBUG_MDD_SET
  printf("%s: p_level: %d, slots: %d\n", __func__, p_level, slots);
  fflush(stdout);
#endif

#if 0
  // Do compaction here
  if (doesLevelNeedCompaction(k)) compactLevel(k);
#endif

#ifdef DEBUG_MDD_SET
  printf("%s: done compacting\n", __func__);
  fflush(stdout);
#endif

#if 1
  if (search_holes && areHolesRecycled()) {
#else
  if (areHolesRecycled) {
#endif
    // First, try for a hole exactly of this size
    // by traversing the index nodes in the hole grid
    int chain = 0;
    int curr = level[p_level].holes_bottom;
    while (curr) {
      if (slots == -(level[p_level].data[curr])) break;
      if (slots < -(level[p_level].data[curr])) {
        // no exact match possible
        curr = 0;
        break;
      }
      // move up the hole grid
      curr = level[p_level].data[curr+1];
      chain++;
    }

    // update max hole chain for the level and the entire mdd
    level[p_level].max_hole_chain = MAX(level[p_level].max_hole_chain, chain);
    max_hole_chain = MAX(max_hole_chain, chain);

    if (curr) {
      // perfect fit
      level[p_level].hole_slots -= slots;
      // try to not remove the "index" node
      int next = level[p_level].data[curr + 3];
      if (next) {
        midRemove(k, next);
#ifdef MEMORY_TRACE
        cout << "Removed Non-Index Hole " << next << "\n";
        dumpInternal(stdout);
#endif
        return next;
      }
      indexRemove(k, curr);
#ifdef MEMORY_TRACE
      cout << "Removed Index Hole " << curr << "\n";
      dumpInternal(stdout);
#endif
      return curr;
    }

    // No hole with exact size, try the largest hole
    curr = level[p_level].holes_top;
    if (slots < -(level[p_level].data[curr]) - min_node_size) {
      // we have a hole large enough
      level[p_level].hole_slots -= slots;
      if (level[p_level].data[curr + 3]) {
        // remove middle node
        curr = level[p_level].data[curr + 3];
        midRemove(k, curr);
      } else {
        // remove index node
        indexRemove(k, curr);
      }
      // create a hole for the leftovers
      int newhole = curr + slots;
      int newsize = -(level[p_level].data[curr]) - slots;
      level[p_level].data[newhole] = -newsize;
      level[p_level].data[newhole + newsize - 1] = -newsize;
      gridInsert(k, newhole); 
#ifdef MEMORY_TRACE
      // level[p_level].data[curr] = -slots;  // only necessary for display
      cout << "Removed part of hole " << curr << "\n";
      dumpInternal(stdout);
#endif
      return curr;
    }
  }
  // can't recycle; grab from the end
  if (level[p_level].last + slots >= level[p_level].size) {
#if 0
    printf("Expand level %d: old size=%d", p_level, level[p_level].size);
#endif
    // not enough space, extend
    int *old_data = level[p_level].data;
    int old_size = level[p_level].size;
    while (level[p_level].last + slots >= level[p_level].size) {
      // level[p_level].size *= 2;
      level[p_level].size = (level[p_level].size > 1024)?
        level[p_level].size + 1024: level[p_level].size * 2;
    }
#if 0
    printf(" new size=%d slots=%d last=%d\n",
        level[p_level].size, slots, level[p_level].last);
    fflush(stdout);
#endif
    level[p_level].data =
      (int*) realloc(level[p_level].data, level[p_level].size * sizeof(int));
    // assert(NULL != level[p_level].data);
    if (NULL == level[p_level].data) {
      // garbage collect and try again
      level[p_level].data = old_data;
      level[p_level].size = old_size;
      reportMemoryUsage(stdout);
      fprintf(stderr, "Memory allocation error while expand MDD level.\n");
      exit(1);
    } else {
      updateMemoryAllocated((level[p_level].size - old_size) * sizeof(int));
      memset(level[p_level].data + old_size, 0,
          (level[p_level].size - old_size) * sizeof(int));
    }
  }
  int h = level[p_level].last + 1;
  level[p_level].last += slots;
  // level[p_level].max_slots = MAX(level[p_level].max_slots, level[p_level].last);
  return h;
}


void node_manager::makeHole(int k, int addr, int slots)
{
  // need to map level
  int mapped_k = mapLevel(k);
#ifdef MEMORY_TRACE
  cout << "Calling makeHole(" << k << ", " << addr << ", " << slots << ")\n";
#endif

  curr_slots -= slots;

  int* data = level[mapped_k].data;
  level[mapped_k].hole_slots += slots;
  data[addr] = data[addr+slots-1] = -slots;

  if (!areHolesRecycled()) return;

  bool control_flag = true;
  bool enable_merging = (slots < getLevelSize(k)) || control_flag;
  // Check for a hole to the left
  if (enable_merging && data[addr-1] < 0) {
    // Merge!
    int lefthole = addr + data[addr-1];
    DCASSERT(data[lefthole] == data[addr-1]);
    if (data[lefthole+1] == non_index_hole) midRemove(k, lefthole);
    else indexRemove(k, lefthole);
    slots += (-data[lefthole]);
    addr = lefthole;
    data[addr] = data[addr+slots-1] = -slots;
  }
  
  // if addr is the last hole, absorb into free part of array
  assert(addr + slots - 1 <= level[mapped_k].last);
  if (addr+slots-1 == level[mapped_k].last) {
    level[mapped_k].last -= slots;
    level[mapped_k].hole_slots -= slots;
    if (level[mapped_k].size > add_size &&
        level[mapped_k].last < level[mapped_k].size/2) {
      int new_size = level[mapped_k].size/2;
      while (new_size > add_size && new_size > level[mapped_k].last * 3)
      { new_size /= 2; }
      updateMemoryAllocated((new_size - level[mapped_k].size) * sizeof(int));
      level[mapped_k].data = (int *)
        realloc(level[mapped_k].data, new_size * sizeof(int));
      assert(NULL != level[mapped_k].data);
      level[mapped_k].size = new_size;
#ifdef MEMORY_TRACE
      printf("Reduced data[] by a factor of 2. New size: %d, Last: %d.\n",
          level[mapped_k].size, level[mapped_k].last);
#endif
    }
#ifdef MEMORY_TRACE
    cout << "Level " << k << ", Made Last Hole " << addr << "\n";
    dumpInternal(stdout);
#endif
    return;
  }

  enable_merging = (slots < getLevelSize(k)) || control_flag;
  // Check for a hole to the right
  if (enable_merging && data[addr+slots]<0) {
    // Merge!
    int righthole = addr+slots;
    if (data[righthole+1] == non_index_hole) midRemove(k, righthole);
    else indexRemove(k, righthole);
    slots += (-data[righthole]);
    data[addr] = data[addr+slots-1] = -slots;
  }

  // Add hole to grid
  gridInsert(k, addr); 

#ifdef MEMORY_TRACE
  cout << "Level " << k << ", Made Last Hole " << addr << "\n";
  dumpInternal(stdout);
#endif
}

void node_manager::reportMemoryUsage(FILE * s, const char filler) {
  fprintf(s, "%cPeak Nodes:             %d\n", filler, getPeakNumNodes());
  fprintf(s, "%cActive Nodes:           %d\n", filler, getCurrentNumNodes());
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
  fprintf(s, "%cReclaimed Nodes:        %d\n", filler, reclaimed_nodes);
  fprintf(s, "%cMem Used:               %d\n", filler,
      getCurrentMemoryUsed());
  fprintf(s, "%cPeak Mem Used:          %d\n", filler, getPeakMemoryUsed());
  fprintf(s, "%cMem Allocated:          %d\n", filler,
      getCurrentMemoryAllocated());
  fprintf(s, "%cPeak Mem Allocated:     %d\n",
      filler, getPeakMemoryAllocated());
  fprintf(s, "%cUnique Tbl Mem Used:    %d\n", filler,
      getUniqueTableMemoryUsed());
  fprintf(s, "%cCompactions:            %d\n", filler, getCompactionsCount());
#if 0
  fprintf(s, "%cHole Memory Usage:\t%d\n", filler, getHoleMemoryUsage());
  fprintf(s, "%cMax Hole Chain:\t%d\n", filler, getMaxHoleChain());
  fprintf(s, "%cCompactions:\t\t%d\n", filler, getCompactionsCount());
  // compareCacheCounts();
#endif
}

void node_manager::compareCacheCounts(int p)
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

void node_manager::validateIncounts()
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
          CHECK_RANGE(0, dptrs[j], sz);
          in_count[dptrs[j]]++;
        }
      } else {
        DCASSERT(isSparseNode(i));
        dptrs = getSparseNodeDownPtrs(i);
        for (int j = getSparseNodeSize(i) - 1; j >= 0; --j) {
          CHECK_RANGE(0, dptrs[j], sz);
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


void node_manager::showLevel(FILE *s, int k) const {
  dumpInternalLevel(s, k);
}
void node_manager::showAll(FILE *s, int verb) const { 
  if (0==verb)  return;
  if (1==verb)  dump(s);
  else          dumpInternal(s); 
}

void node_manager::show(FILE *s, int h) const { fprintf(s, "%d", h); }

int node_manager::getHoleMemoryUsage() const {
  int sum = 0;
  for(int i=0; i<l_size; i++) sum += level[i].hole_slots;
  return sum * sizeof(int); 
}

int node_manager::getMaxHoleChain() const { return max_hole_chain; }
int node_manager::getCompactionsCount() const { return num_compactions; }

int node_manager::getCurrentMemoryUsed() const { 
#if 0
  int sum = 0;
  for(int i=0; i<l_size; i++) sum += level[i].last - level[i].hole_slots;
  return sum * sizeof(int); 
#else
  return curr_slots * sizeof(int);
#endif
}

