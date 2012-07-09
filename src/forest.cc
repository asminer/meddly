
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



// TODO: Testing

#include "defines.h"
#include "unique_table.h"
#include "hash_stream.h"
#include <set>
#include <queue>
#include <vector>

// #define DEBUG_CLEANUP
// #define MERGE_RIGHT
// #define MERGE_LEFT
// #define ENABLE_BREAKING_UP_HOLES

// #define DEBUG_SLOW
// #define DEBUG_ADDRESS_RESIZE
// #define DEBUG_GC

// #define MEMORY_TRACE

const int a_min_size = 1024;

// ******************************************************************
// *                                                                *
// *                    forest::statset  methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::forest::statset::statset()
{
  reclaimed_nodes = 0;
  num_compactions = 0;
  zombie_nodes = 0;
  orphan_nodes = 0;
  active_nodes = 0;
  peak_active = 0;
  memory_used = 0;
  memory_alloc = 0;
  peak_memory_used = 0;
  peak_memory_alloc = 0;
  memory_UT = 0;
  peak_memory_UT = 0;
  max_UT_chain = 0;
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                         forest methods                         *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::forest
::forest(int ds, domain* _d, bool rel, range_type t, edge_labeling ev, 
  const policies &p) : deflt(p)
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Creating forest #%d in domain #%d\n", ds, _d->ID());
#endif
  d_slot = ds;
  is_marked_for_deletion = false;
  d = _d;
  isRelation = rel;
  rangeType = t;
  edgeLabel = ev;
  // check policies
  if (!isRelation) {
    if (policies::IDENTITY_REDUCED == deflt.reduction)
      throw error(error::INVALID_POLICY);
  } 
  //
  // Initialize array of operations
  //
  opCount = 0;
  szOpCount = 0;
  //
  // Initialize list of registered dd_edges
  //
  firstHole = -1; // firstHole < 0 indicates no holes.
  firstFree = 0;
  sz = 256;

  // Create an array to store pointers to dd_edges.
  edge = (edge_data *) malloc(sz * sizeof(edge_data));
  for (unsigned i = 0; i < sz; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
}

MEDDLY::forest::~forest()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Deleting forest #%d in domain #%d\n", d_slot, d->ID());
#endif
  // operations are deleted elsewhere...
  free(opCount);
  // Make SURE our edges are orphaned
  for (unsigned i = 0; i < firstFree; ++i) {
    if (edge[i].edge) edge[i].edge->orphan();
  }
  free(edge);
  // NOTE: since the user is provided with the dd_edges instances (as opposed
  // to a pointer), the user program will automatically call the
  // destructor for each dd_edge when the corresponding variable goes out of
  // scope. Therefore there is no need to destruct dd_edges from here.
  d->unlinkForest(this, d_slot);
}

void MEDDLY::forest::markForDeletion()
{
#ifdef DEBUG_CLEANUP
  fprintf(stderr, "Marking forest #%d for deletion in domain #%d\n", d_slot, d->ID());
#endif
  is_marked_for_deletion = true;
  // deal with operations associated with this forest
  for (int i=0; i<szOpCount; i++) 
    if (opCount[i]) {
      operation* op = operation::getOpWithIndex(i);
      op->markForDeletion();
    }
  unregisterDDEdges();
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, bool* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, int* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdgeForVar(int vh, bool vp, float* terms, dd_edge& a)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(const int* const* vlist, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vl, const int* const* vpl, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* const* vplist,
      const int* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::createEdge(const int* const* vlist, const int* const* vplist,
      const float* terms, int N, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(bool val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(int val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::createEdge(float val, dd_edge &e)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, bool &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, int &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::evaluate(const dd_edge &f, const int* vl, float &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, bool &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, int &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::evaluate(const dd_edge& f, const int* vl, const int* vpl, float &t) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest::getElement(const dd_edge& a, int index, int* e)
{
  throw error(error::INVALID_OPERATION);
}

void MEDDLY::forest::findFirstElement(const dd_edge& f, int* vlist) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::forest
::findFirstElement(const dd_edge& f, int* vlist, int* vplist) const
{
  throw error(error::TYPE_MISMATCH);
}




void MEDDLY::forest::removeStaleComputeTableEntries()
{
  if (operation::useMonolithicComputeTable()) {
    operation::removeStalesFromMonolithic();
  } else {
    for (int i=0; i<szOpCount; i++) 
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->removeStaleComputeTableEntries();
      }
  }
}

void MEDDLY::forest::removeAllComputeTableEntries()
{
  if (is_marked_for_deletion) return;
  if (operation::useMonolithicComputeTable()) {
    is_marked_for_deletion = true;
    operation::removeStalesFromMonolithic();
    is_marked_for_deletion = false;
  } else {
    for (int i=0; i<szOpCount; i++) 
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->removeAllComputeTableEntries();
      }
  }
}

void MEDDLY::forest::showComputeTable(FILE* s, int verbLevel) const
{
  if (operation::useMonolithicComputeTable()) {
    operation::showMonolithicComputeTable(s, verbLevel);
  } else {
    for (int i=0; i<szOpCount; i++) 
      if (opCount[i]) {
        operation* op = operation::getOpWithIndex(i);
        op->showComputeTable(s, verbLevel);
      }
  }
}

void MEDDLY::forest::registerOperation(const operation* op)
{
  MEDDLY_DCASSERT(op->getIndex() >= 0);
  if (op->getIndex() >= szOpCount) {
    // need to expand
    int newSize = ((op->getIndex() / 16) +1 )*16; // expand in chunks of 16
    int* tmp = (int*) realloc(opCount, newSize * sizeof(int));
    if (0==tmp) throw error(error::INSUFFICIENT_MEMORY);
    for ( ; szOpCount < newSize; szOpCount++) {
      tmp[szOpCount] = 0;
    }
    opCount = tmp;
  }
  opCount[op->getIndex()] ++;
}

void MEDDLY::forest::unregisterOperation(const operation* op)
{
  MEDDLY_DCASSERT(op->getIndex() >= 0);
  MEDDLY_DCASSERT(szOpCount > op->getIndex());
  MEDDLY_DCASSERT(opCount[op->getIndex()]>0);
  opCount[op->getIndex()] --;
}

void MEDDLY::forest::registerEdge(dd_edge& e) 
{
  // add to collection of edges for this forest.
  // change e.index to help find this edge at a later time.
  if (firstHole >= 0) {
    // hole available; fill it up
    int index = firstHole;
    firstHole = edge[firstHole].nextHole;
    edge[index].edge = &e;
    edge[index].nextHole = -1;
    e.setIndex(index);
  } else {
    // no holes available, add to end of array
    if (firstFree >= sz) {
      // expand edge[]
      int new_sz = sz * 2;
      edge_data* new_edge =
          (edge_data*) realloc(edge, new_sz * sizeof(edge_data));
      if (0 == new_edge) throw error(error::INSUFFICIENT_MEMORY);
      edge = new_edge;
      for (int i = sz; i < new_sz; ++i)
      {
        edge[i].nextHole = -1;
        edge[i].edge = 0;
      }
      sz = new_sz;
    }
    MEDDLY_DCASSERT(firstFree < sz);
    edge[firstFree].nextHole = -1;
    edge[firstFree].edge = &e;
    e.setIndex(firstFree);
    ++firstFree;
  }
}


void MEDDLY::forest::unregisterEdge(dd_edge& e) 
{
  // remove this edge from the collection of edges for this forest.
  // change e.index to -1.
  MEDDLY_DCASSERT(e.getIndex() >= 0);
  int index = e.getIndex();
  MEDDLY_DCASSERT(edge[index].edge == &e);
  edge[index].edge = 0;
  edge[index].nextHole = firstHole;
  firstHole = index;
  e.setIndex(-1);
}


void MEDDLY::forest::unregisterDDEdges() 
{
  // Go through the list of valid edges (value > 0), and set
  // the e.index to -1 (indicating unregistered edge).

  // ignore the NULLs; release the rest
  for (unsigned i = 0; i < firstFree; ++i) {
    if (edge[i].edge != 0) {
      MEDDLY_DCASSERT(edge[i].nextHole == -1);
      edge[i].edge->set(0, 0, 0);
      edge[i].edge->setIndex(-1);
    }
  }

  // firstHole < 0 indicates no holes.
  for (unsigned i = 0; i < firstFree; ++i) {
    edge[i].nextHole = -1;
    edge[i].edge = 0;
  }
  firstHole = -1;
  firstFree = 0;
}



// ******************************************************************
// *                                                                *
// *               expert_forest::nodeFinder  methods               *
// *                                                                *
// ******************************************************************

MEDDLY::expert_forest::nodeFinder::nodeFinder(const expert_forest* p, int n)
{
  parent = p;
  node = n;
  nodeLevel = parent->getNodeLevel(node);
  h = parent->hashNode(node);
  MEDDLY_DCASSERT(node);
  MEDDLY_DCASSERT(!parent->isTerminalNode(node));
  MEDDLY_DCASSERT(parent->isActiveNode(node));
}

bool MEDDLY::expert_forest::nodeFinder::equalsFF(int h1, int h2) const
{
  MEDDLY_DCASSERT(parent->isFullNode(h1));
  MEDDLY_DCASSERT(parent->isFullNode(h2));

  int *ptr1 = parent->getNodeAddress(h1) + 2;
  int *ptr2 = parent->getNodeAddress(h2) + 2;
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

  if (parent->isMultiTerminal()) return true;

  // Check edge-values
  MEDDLY_DCASSERT(ptr1 == h1Stop);
  MEDDLY_DCASSERT(ptr2 == h2Stop);

  h1Stop += MIN(sz1, sz2);
  if (parent->isEVPlus()) {
    while (ptr1 != h1Stop) { if (*ptr1++ != *ptr2++) return false; }
  } else {
    MEDDLY_DCASSERT(parent->isEVTimes());
    while (ptr1 != h1Stop) {
      if (!isAlmostEqual(*ptr1++, *ptr2++)) return false;
    }
  }
  return true;
}


bool MEDDLY::expert_forest::nodeFinder::equalsSS(int h1, int h2) const
{
  MEDDLY_DCASSERT(parent->isSparseNode(h1));
  MEDDLY_DCASSERT(parent->isSparseNode(h2));

  int *ptr1 = parent->getNodeAddress(h1) + 2;
  int *ptr2 = parent->getNodeAddress(h2) + 2;
  int sz1 = -(*ptr1++);
  int sz2 = -(*ptr2++);

  if (sz1 != sz2) return false;

  int* h1Stop = ptr1 + sz1 + sz1;
  while (ptr1 != h1Stop) { if (*ptr1++ != *ptr2++) return false; }

  if (parent->isMultiTerminal()) return true;

  // Check edge-values
  MEDDLY_DCASSERT(ptr1 == h1Stop);
  MEDDLY_DCASSERT(ptr2 == (parent->getNodeAddress(h2) + 3 + sz1 + sz1));

  h1Stop += sz1;
  if (parent->isEVPlus()) {
    while (ptr1 != h1Stop) { if (*ptr1++ != *ptr2++) return false; }
  } else {
    MEDDLY_DCASSERT(parent->isEVTimes());
    while (ptr1 != h1Stop) {
      if (!isAlmostEqual(*ptr1++, *ptr2++)) return false;
    }
  }
  return true;
}


bool MEDDLY::expert_forest::nodeFinder::equalsFS(int h1, int h2) const
{
  MEDDLY_DCASSERT(parent->isFullNode(h1));
  MEDDLY_DCASSERT(parent->isSparseNode(h2));

  int *ptr1 = parent->getNodeAddress(h1) + 2;
  int *ptr2 = parent->getNodeAddress(h2) + 2;
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

  if (parent->isMultiTerminal()) return true;

  // Check edge-values
  MEDDLY_DCASSERT(ptr1 == h1Stop);
  MEDDLY_DCASSERT(ptr2 == h2Stop);
  MEDDLY_DCASSERT(down2 == h2Stop + sz2);

  // ptr1 and down2 are pointing at the start of edge-values
  // Reset the index pointer for h2 (sparse node).
  ptr2 -= sz2;
  if (parent->isEVPlus()) {
    while (ptr2 != h2Stop) {
      if (ptr1[*ptr2++] != *down2++) return false;
    }
  } else {
    MEDDLY_DCASSERT(parent->isEVTimes());
    while (ptr2 != h2Stop) {
      if (!isAlmostEqual(ptr1[*ptr2++], *down2++)) return false;
    }
  }
  return true;
}




// ******************************************************************
// *                                                                *
// *               expert_forest::level_data  methods               *
// *                                                                *
// ******************************************************************

/*

  Hole management.
  ==============================================================
  There are two kinds of holes depending on their location in the "hole grid":
  Index Holes and Non-index Holes.

  The hole grid structure:
  ------------------------
  (holes_bottom)
  holes_of_size_0 (index) -- holes_of_size_0 (non_index) -- (non_index) -- NULL
  |
  holes_of_size_1 -- ""
  |
  :
  :
  (holes_top)

  Index holes are represented as follows:
  ---------------------------------------
  [0] -size (number of slots in hole)     
  [1] up
  [2] down 
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size

  Non-index holes are represented as follows:
  [0] -size (number of slots in hole)     
  [1] flag (<0, indicates non-index node)
  [2] prev pointer (nodes of same size)
  [3] next pointer (nodes of same size)
  [4..size-2] Unused
  :
  :
  [size-1] -size
*/

MEDDLY::expert_forest::level_data::level_data()
{
  parent = 0;
}

MEDDLY::expert_forest::level_data::~level_data()
{
  // printf("Destroying level; max hole chain: %d\n", max_hole_chain);
  if (parent) parent->stats.decMemAlloc(size*sizeof(int));
  free(data);
}

void MEDDLY::expert_forest::level_data
::init(expert_forest* p, char eS, char uH, char hH)
{
  MEDDLY_DCASSERT(0==parent);
  parent = p;
  data = 0;
  size = 0;
  last = 0;
  holes_top = 0;
  holes_bottom = 0;
  hole_slots = 0;
  max_hole_chain = 0;
  zombie_nodes = 0;
  temp_nodes = 0;
  num_compactions = 0;
  levelNode = 0;
  edgeSize = eS;
  unhashedHeader = uH;
  hashedHeader = hH;
  compactLevel = false;
}


int MEDDLY::expert_forest::level_data::allocNode(int sz, int tail, bool clear)
{
  MEDDLY_DCASSERT(parent);
  int slots = slotsForNode(sz);
  MEDDLY_DCASSERT(slots >= commonExtra + unhashedHeader + hashedHeader);

  int off = getHole(slots);
  if (clear) memset(data+off, 0, slots*sizeof(int));
  countOf(off) = 1;                       // #incoming
  nextOf(off) = temp_node_value;          // mark as a temp node
  sizeOf(off) = sz;                       // size
  data[off+slots-1] = tail;               // tail entry
#ifdef MEMORY_TRACE
  printf("Allocated new node size %d, position %d\n", sz, off);
  dumpInternal(stdout, parent->levels - this);
#endif
  return off;
}

void MEDDLY::expert_forest::level_data::compact(node_data* address)
{
  if (0==data) {
    compactLevel = false;
    return;
  }
  if (0 == hole_slots || !needsCompaction()) {  // Level is compact enough!
    compactLevel = false;
    return;
  }

  if (0 < temp_nodes) return;   // Temp nodes; do not compact

#ifdef DEBUG_SLOW
  fprintf(stderr, "Compacting forest level\n");
#endif

  // alternate algorithm -- since we now have the node ids in the node data
  int *node_ptr = data + 1;  // since we leave [0] empty
  int *end_ptr = data + last + 1;
  int *curr_ptr = node_ptr;
  int node_size = 0;
  int curr_node = 0;

  while (node_ptr != end_ptr) {
    // find new node
    if (*node_ptr < 0) {
      // found a hole, advance
      MEDDLY_DCASSERT(node_ptr[0] == node_ptr[-(*node_ptr)-1]);
      node_size = -(*node_ptr);
      memset(node_ptr, 0, node_size * sizeof(int));
    } else {
      // found an existing node
      MEDDLY_DCASSERT(!parent->isPessimistic() || *node_ptr != 0);

      node_size = activeNodeActualSlots(node_ptr - data); 
      curr_node = node_ptr[node_size - 1];
      MEDDLY_DCASSERT(parent->getNodeOffset(curr_node) == (node_ptr - data));
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
        address[curr_node].offset = (curr_ptr - data);
      }
      MEDDLY_DCASSERT(parent->getNodeOffset(curr_node) == (curr_ptr - data));
      curr_ptr += node_size;
    }
    node_ptr += node_size;
  }

  last = (curr_ptr - 1 - data);

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  hole_slots = 0;

  parent->stats.num_compactions++;
  num_compactions++;
  compactLevel = false;

  if (size > add_size && last < size/2) {
    int new_size = size/2;
    while (new_size > add_size && new_size > last * 3) { new_size /= 2; }
    parent->stats.decMemAlloc((size - new_size) * sizeof(int));
    data = (int *) realloc(data, new_size * sizeof(int));
    if (0 == data) throw error(error::INSUFFICIENT_MEMORY);
    size = new_size;
#ifdef MEMORY_TRACE
    printf("Reduced data[] by a factor of 2. New size: %d, Last: %d.\n",
        size, last);
#endif
  }

#if 0
  printf("After compaction:\n");
  dumpInternalLevel(stdout, k);
  printf("\n");
#endif
}

void MEDDLY::expert_forest::level_data::dumpInternal(FILE *s) const
{
  if (0==data) return; // nothing to display
  
  // super clever hack:
  fprintf(s, "Level %ld: ", long(this - parent->levels));
  fprintf(s, "Last slot used: %d\n", last);
  fprintf(s, "Grid: top = %d bottom = %d\n", holes_top, holes_bottom);

  fprintf(s, "Data array by record: \n");
  int awidth = digits(parent->a_last);
  int a;
  for (a=1; a<=last; ) {
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
      int nElements = activeNodeActualSlots(a);
      for (int i=3; i<nElements-1; i++) {
        fprintf(s, "|%d", data[a+i]);
      }
      a += activeNodeActualSlots(a);
    }
    fprintf(s, "|%d]\n", data[a-1]);
  } // for a
  fprintf(s, "%*d : free slots\n", awidth, a);
  fflush(s);
  MEDDLY_DCASSERT(a == (last)+1);
}

void MEDDLY::expert_forest::level_data
::addToChainCounts(std::map<int, int> &chainLengths) const
{
  for (int curr = holes_bottom; curr; curr = data[curr + 1]) {
    int currHoleOffset = curr;
    int count = 0;
    // traverse this chain to get its length
    for (count = 0; currHoleOffset; count++) {
      currHoleOffset = data[currHoleOffset + 3];
    }
    int currHoleSize = -data[curr];
    chainLengths[currHoleSize] += count;
  }
}

//
//
// Private
//
//

int MEDDLY::expert_forest::level_data::getHole(int slots)
{
  MEDDLY_DCASSERT(parent);
  parent->stats.incMemUsed(slots * sizeof(int));

  // First, try for a hole exactly of this size
  // by traversing the index nodes in the hole grid
  int chain = 0;
  int curr = holes_bottom;
  while (curr) {
    if (slots == -(data[curr])) break;
    if (slots < -(data[curr])) {
      // no exact match possible
      curr = 0;
      break;
    }
    // move up the hole grid
    curr = data[curr+1];
    chain++;
  }

  // update max hole chain for the level and the entire mdd
  max_hole_chain = MAX(max_hole_chain, chain);

  if (curr) {
    // perfect fit
    hole_slots -= slots;
    // try to not remove the "index" node
    int next = data[curr + 3];
    if (next) {
      midRemove(next);
#ifdef MEMORY_TRACE
      cout << "Removed Non-Index Hole " << next << "\n";
      dumpInternal(stdout);
#endif
      return next;
    }
    indexRemove(curr);
#ifdef MEMORY_TRACE
    cout << "Removed Index Hole " << curr << "\n";
    dumpInternal(stdout);
#endif
    return curr;
  }

#ifdef ENABLE_BREAKING_UP_HOLES
  // No hole with exact size, try the largest hole
  const int min_node_size = slotsForNode(0);

  curr = holes_top;
  if (slots < -(data[curr]) - min_node_size) {
    // we have a hole large enough
    hole_slots -= slots;
    if (data[curr + 3]) {
      // remove middle node
      curr = data[curr + 3];
      midRemove(k, curr);
    } else {
      // remove index node
      indexRemove(k, curr);
    }
    // create a hole for the leftovers
    int newhole = curr + slots;
    int newsize = -(data[curr]) - slots;
    data[newhole] = -newsize;
    data[newhole + newsize - 1] = -newsize;
    gridInsert(k, newhole); 
#ifdef MEMORY_TRACE
    // data[curr] = -slots;  // only necessary for display
    cout << "Removed part of hole " << curr << "\n";
    dumpInternal(stdout);
#endif
    return curr;
  }
#endif

  // can't recycle; grab from the end
  if (last + slots >= size) {
    // not enough space, extend
    int old_size = size;

    // new size is 50% more than previous (37.5% * 4 = 1.5 => 50% growth)
    size = MAX( old_size, last + slots ) * 1.5;

    data = (int*) realloc(data, size * sizeof(int));
    if (0 == data) {
      // TBD: garbage collect and try again
      throw error(error::INSUFFICIENT_MEMORY);
    } else {
      parent->stats.incMemAlloc((size - old_size) * sizeof(int));
      memset(data + old_size, 0, (size - old_size) * sizeof(int));
    }
  }
  int h = last + 1;
  last += slots;
  return h;
}

void MEDDLY::expert_forest::level_data::makeHole(int addr, int slots)
{
#ifdef MEMORY_TRACE
  printf("Calling makeHole(%d, %d)\n", addr, slots);
#endif
  MEDDLY_DCASSERT(parent);

  parent->stats.decMemUsed(slots * sizeof(int));

  hole_slots += slots;
  data[addr] = data[addr+slots-1] = -slots;

  if (!parent->areHolesRecycled()) return;

  // Check for a hole to the left
#ifdef MERGE_LEFT
  if (data[addr-1] < 0) {
    // Merge!
    int lefthole = addr + data[addr-1];
    MEDDLY_DCASSERT(data[lefthole] == data[addr-1]);
    if (data[lefthole+1] == non_index_hole) midRemove(k, lefthole);
    else indexRemove(k, lefthole);
    slots += (-data[lefthole]);
    addr = lefthole;
    data[addr] = data[addr+slots-1] = -slots;
  }
#endif

  // if addr is the last hole, absorb into free part of array
  MEDDLY_DCASSERT(addr + slots - 1 <= last);
  if (addr+slots-1 == last) {
    last -= slots;
    hole_slots -= slots;
    if (size > add_size && (last + 1) < size/2) {
      int new_size = size/2;
      while (new_size > (last + 1) * 2) new_size /= 2;
      if (new_size < add_size) new_size = add_size;
      parent->stats.incMemAlloc((new_size - size) * sizeof(int));
      data = (int *) realloc(data, new_size * sizeof(int));
      if (0 == data) throw error(error::INSUFFICIENT_MEMORY);
      size = new_size;
#ifdef MEMORY_TRACE
      printf("Reduced data[]. New size: %d, Last: %d.\n", size, last);
#endif
    }
#ifdef MEMORY_TRACE
    printf("Made Last Hole %d\n", addr);
    dumpInternal(stdout, parent->levels - this);
#endif
    return;
  }

#ifdef MERGE_RIGHT
  // Check for a hole to the right
  if (data[addr+slots]<0) {
    // Merge!
    int righthole = addr+slots;
    if (data[righthole+1] == non_index_hole) midRemove(k, righthole);
    else indexRemove(k, righthole);
    slots += (-data[righthole]);
    data[addr] = data[addr+slots-1] = -slots;
  }
#endif

  // Add hole to grid
  gridInsert(addr); 

#ifdef MEMORY_TRACE
  printf("Made Last Hole %d\n", addr);
  dumpInternal(stdout, parent->levels - this);
#endif
}

void MEDDLY::expert_forest::level_data::gridInsert(int p_offset)
{
  // sanity check to make sure that the first and last slots in this hole
  // have the same value, i.e. -(# of slots in the hole)
  MEDDLY_DCASSERT(data[p_offset] == data[p_offset - data[p_offset] - 1]);
  // special case: empty
  if (0 == holes_bottom) {
    // index hole
    data[p_offset + 1] = data[p_offset + 2] = data[p_offset + 3] = 0;
    holes_top = holes_bottom = p_offset;
    return;
  }
  // special case: at top
  if (data[p_offset] < data[holes_top]) {
    // index hole
    data[p_offset + 1] = data[p_offset + 3] = 0;
    data[p_offset + 2] = holes_top;
    data[holes_top + 1] = p_offset;
    holes_top = p_offset;
    return;
  }
  int above = holes_bottom;
  int below = 0;
  while (data[p_offset] < data[above]) {
    below = above;
    above = data[below + 1];
    MEDDLY_DCASSERT(data[above + 2] == below);
    MEDDLY_DCASSERT(above);  
  }
  if (data[p_offset] == data[above]) {
    // Found, add this to chain
    // making a non-index hole
    int right = data[above + 3];
    data[p_offset + 1] = non_index_hole;
    data[p_offset + 2] = above;
    data[p_offset + 3] = right;
    if (right) data[right + 2] = p_offset;
    data[above + 3] = p_offset;
    return; 
  }
  // we should have above < p_offset < below  (remember, -sizes)
  // create an index hole since there were no holes of this size
  data[p_offset + 1] = above;
  data[p_offset + 2] = below;
  data[p_offset + 3] = 0;
  data[above + 2] = p_offset;
  if (below) {
    data[below + 1] = p_offset;
  } else {
    MEDDLY_DCASSERT(above == holes_bottom);
    holes_bottom = p_offset;
  }
}

void MEDDLY::expert_forest::level_data::midRemove(int p_offset)
{
  MEDDLY_DCASSERT(isHoleNonIndex(p_offset));
  int left = data[p_offset+2];
  MEDDLY_DCASSERT(left);
  int right = data[p_offset+3];

  data[left + 3] = right;
  if (right) data[right + 2] = left;
}

void MEDDLY::expert_forest::level_data::indexRemove(int p_offset)
{
#ifdef MEMORY_TRACE
  printf("indexRemove(%d)\n", p_offset);
#endif

  MEDDLY_DCASSERT(!isHoleNonIndex(p_offset));
  int above = data[p_offset + 1];
  int below = data[p_offset + 2];
  int right = data[p_offset + 3];

  if (right >= 1) {
    // there are nodes to the right!
    MEDDLY_DCASSERT(data[right + 1] < 0);
    data[right + 1] = above;
    data[right + 2] = below;

    // update the pointers of the holes (index) above and below it
    if (above) {
      data[above + 2] = right;
    } else {
      holes_top = right;
    }

    if (below) {
      data[below + 1] = right;
    } else {
      holes_bottom = right;
    }
    
  } else {
    // there are no non-index nodes
    MEDDLY_DCASSERT(right < 1);

    // this was the last node of its size
    // update the pointers of the holes (index) above and below it
    if (above) {
      data[above + 2] = below;
    } else {
      holes_top = below;
    }

    if (below) {
      data[below + 1] = above;
    } else {
      holes_bottom = above;
    }
  }
}

// ******************************************************************
// *                                                                *
// *               expert_forest::nodeReader  methods               *
// *                                                                *
// ******************************************************************

MEDDLY::expert_forest::nodeReader::nodeReader(int k)
{
  buffer = 0;
  alloc = 0;
  size = 0;
  level = k;
}

MEDDLY::expert_forest::nodeReader::~nodeReader()
{
  free(buffer);
}

void MEDDLY::expert_forest::nodeReader::dump(FILE* s) const
{
  fprintf(s, "[%d", buffer[0]);
  for (int i=1; i<size; i++)
    fprintf(s, ", %d", buffer[i]);
  fprintf(s, "]");
}

void MEDDLY::expert_forest::nodeReader::resize(int ns)
{
  size = ns;
  if (size <= alloc) return;
  int nalloc = ((ns/8)+1)*8;
  MEDDLY_DCASSERT(nalloc > ns);
  MEDDLY_DCASSERT(nalloc>0);
  MEDDLY_DCASSERT(nalloc>alloc);
  buffer = (int*) realloc(buffer, nalloc*sizeof(int));
  if (0==buffer) throw error(error::INSUFFICIENT_MEMORY);
  alloc = nalloc;
}

// ******************************************************************
// *                                                                *
// *               expert_forest::nodeBuilder methods               *
// *                                                                *
// ******************************************************************

MEDDLY::expert_forest::nodeBuilder::nodeBuilder()
{
  parent = 0;
  down = 0;
  indexes = 0;
  edges = 0;
}

MEDDLY::expert_forest::nodeBuilder::~nodeBuilder()
{
  free(down);
  free(indexes);
  free(edges);
}

void MEDDLY::expert_forest::nodeBuilder::init(int k, const level_data* ld)
{
  MEDDLY_DCASSERT(ld);
  parent = ld;
  level = k;
  size = 0;
  alloc = 0;
  lock = false;
}

bool MEDDLY::expert_forest::nodeBuilder::equals(int p) const
{
  const node_data &node = parent->parent->getNode(p);
  if (node.level != level) return false;

  MEDDLY_DCASSERT(0==edges);  // edge value stuff not implemented yet

  // p is full:
  if (parent->isFull(node.offset)) {
    int fs = parent->fullSizeOf(node.offset);
    const int* pd = parent->fullDownOf(node.offset);
    if (is_sparse) {
      int i = 0;
      for (int z=0; z<size; z++) {
        if (indexes[z] >= fs) return false;
        for (; i<indexes[z]; i++) if (pd[i]) return false;
        if (down[z] != pd[i]) return false;
        i++;
      } // for z
      for (; i<fs; i++) if (pd[i]) return false;
      return true;
    }
    // we're also full
    if (fs > size) return false;
    if (memcmp(down, pd, fs * sizeof(int))) return false;
    for (int i=fs; i<size; i++) if (down[i]) return false;
    return true;
  }

  // p is sparse:
  int nnz = parent->sparseSizeOf(node.offset);
  const int* pd = parent->sparseDownOf(node.offset);
  const int* pi = parent->sparseIndexesOf(node.offset);

  if (is_sparse) {
    if (nnz != size) return false;
    if (memcmp(down, pd, nnz * sizeof(int))) return false;
    if (memcmp(indexes, pi, nnz * sizeof(int))) return false;
    return true;
  }

  // We must be full
  int i = 0;
  for (int z=0; z<nnz; z++) {
    for (; i<pi[z]; i++) if (down[i]) return false;
    if (down[i] != pd[z]) return false;
    i++;
  }
  for (; i<size; i++) if (down[i]) return false;
  return true;
}

void MEDDLY::expert_forest::nodeBuilder::computeHash()
{
  MEDDLY_DCASSERT(!has_hash);
  
  hash_stream s;
  s.start(level);
  
  if (is_sparse) {
    for (int z=0; z<size; z++) {
      MEDDLY_DCASSERT(down[z]);
      s.push(indexes[z], down[z]);
    }
  } else {
    for (int i=0; i<size; i++) {
      if (0==down[i]) continue;
      s.push(i, down[i]);
    }
  }
  h = s.finish();
  has_hash = true;
}

void MEDDLY::expert_forest::nodeBuilder::copyIntoFull(int* d, int N) const
{
  if (is_sparse) {
    memset(d, 0, N * sizeof(int));
    for (int z=0; z<size; z++) {
      if (indexes[z] < N) d[indexes[z]] = down[z];
    }
  } else {
    memcpy(d, down, N * sizeof(int));
  }
}

void MEDDLY::expert_forest::nodeBuilder
::copyIntoSparse(int* d, int* ix, int Z) const
{
  if (is_sparse) {
    memcpy(d, down, Z * sizeof(int));
    memcpy(ix, indexes, Z * sizeof(int));
  } else {
    int z=0;
    MEDDLY_DCASSERT(Z>0);
    for (int i=0; i<size; i++)
      if (down[i]) {
        d[z] = down[i];
        ix[z] = i;
        z++;
        if (z>=Z) return;
      }
  }
}

void MEDDLY::expert_forest::nodeBuilder::enlarge()
{
  if (size <= alloc) return;
  alloc = ((size / 8)+1) * 8;
  MEDDLY_DCASSERT(alloc > size);
  down = (int*) realloc(down, alloc * sizeof(int));
  if (0==down) throw error(error::INSUFFICIENT_MEMORY);
  indexes = (int*) realloc(indexes, alloc * sizeof(int));
  if (0==indexes) throw error(error::INSUFFICIENT_MEMORY);
  if (parent->edgeSize>0) {
    edges = (int*) realloc(down, alloc * parent->edgeSize * sizeof(int));
    if (0==edges) throw error(error::INSUFFICIENT_MEMORY);
  }
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     expert_forest  methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************



MEDDLY::expert_forest::expert_forest(int ds, domain *d, bool rel, range_type t,
  edge_labeling ev, const policies &p)
: forest(ds, d, rel, t, ev, p)
{
  //
  // Inltialize address array
  //
  a_size = a_min_size;
  address = (node_data *) malloc(a_size * sizeof(node_data));
  if (0 == address) throw error(error::INSUFFICIENT_MEMORY);
  stats.incMemAlloc(a_size * sizeof(node_data));
  memset(address, 0, a_size * sizeof(node_data));
  a_last = a_unused = a_next_shrink = 0;
  // peak_nodes 0;
  

  //
  // Initialize level array
  //
  int N = getNumVariables();
  if (rel) {
    raw_levels = new level_data[2*N+1];
    levels = raw_levels + N;
    stats.incMemAlloc(2*N+1 * sizeof(level_data));
  } else {
    raw_levels = new level_data[N+1];
    levels = raw_levels;
    stats.incMemAlloc(N+1 * sizeof(level_data));
  }

  //
  // Initialize builders array
  //
  if (rel) {
    raw_builders = new nodeBuilder[2*N+1];
    builders = raw_builders + N;
  } else {
    raw_builders = new nodeBuilder[N+1];
    builders = raw_builders;
  }
  for (int k=getMinLevelIndex(); k<=N; k++) {
    builders[k].init(k, levels+k);
  }

  //
  // Initialize node readers array
  //
  if (rel) {
    raw_readers = new nodeReader*[2*N+1];
    free_reader = raw_readers + N;
  } else {
    raw_readers = new nodeReader*[N+1];
    free_reader = raw_readers;
  }
  for (int k=getMinLevelIndex(); k<=N; k++) {
    free_reader[k] = 0;
  }

  //
  // Initialize misc. private data
  //
  unique = new unique_table(this);
  performing_gc = false;
}


MEDDLY::expert_forest::~expert_forest() 
{
  // Address array
  free(address);

  // Level array
  delete[] raw_levels;

  // builders array
  delete[] raw_builders;

  // node readers
  int N = getNumVariables();
  for (int k=getMinLevelIndex(); k<=N; k++) {
    while (free_reader[k]) {
      nodeReader* nr = free_reader[k];
      free_reader[k] = nr->next;
      delete nr;
    }
  }
  delete[] raw_readers;

  // unique table
  delete unique;
}


int* MEDDLY::expert_forest::markNodesInSubgraph(int root, bool sort) const
{
  if (isTerminalNode(root)) return 0;

  // initialize lists
  bool* inList = new bool[a_last];
  for (int i=0; i<a_last; i++) inList[i] = false;
  inList--;
  int mlen = 0;
  int msize = 1024;
  int* marked = (int*) malloc(msize * sizeof(int));

  // Initialize search
  marked[mlen] = root;
  mlen++;
  inList[root] = true;

  // Breadth-first search
  for (int mexpl=0; mexpl<mlen; mexpl++) {
    // explore node marked[mexpl]
    if (isFullNode(marked[mexpl])) {
      const int sz = getFullNodeSize(marked[mexpl]);
      for (int i = 0; i < sz; ++i)
      {
        int dn = getFullNodeDownPtr(marked[mexpl], i);
        if (isTerminalNode(dn)) continue;
        MEDDLY_CHECK_RANGE(0, dn-1, a_last);
        if (inList[dn]) continue;
        // add dn to list
        if (mlen+1 >= msize) { 
          // expand.  Note we're leaving an extra slot
          // at the end, for the terminal 0.
          msize += 1024;
          marked = (int*) realloc(marked, msize*sizeof(int));
          if (0==marked) throw error(error::INSUFFICIENT_MEMORY);
        }
        inList[dn] = true;
        marked[mlen] = dn;
        mlen++;
      } // for i
    }
    else { 
      const int sz = getSparseNodeSize(marked[mexpl]);
      for (int i = 0; i < sz; ++i)
      {
        int dn = getSparseNodeDownPtr(marked[mexpl], i);
        if (isTerminalNode(dn)) continue;
        MEDDLY_CHECK_RANGE(0, dn-1, a_last);
        if (inList[dn]) continue;
        // add dn to list
        if (mlen+1 >= msize) { 
          // expand.  Note we're leaving an extra slot
          // at the end, for the terminal 0.
          msize += 1024;
          marked = (int*) realloc(marked, msize*sizeof(int));
          if (0==marked) throw error(error::INSUFFICIENT_MEMORY);
        }
        inList[dn] = true;
        marked[mlen] = dn;
        mlen++;
      } // for i
    }
  } // for mexpl

  // sort
  if (sort && mlen>0) {
    mlen = 0;
    for (int i=1; i<=a_last; i++) if (inList[i]) {
      marked[mlen] = i;
      mlen++;
    }
  }

  // cleanup
  inList++;
  delete[] inList;
  if (0==mlen) {
    free(marked);
    return 0;
  }
  // add 0 to the list
  marked[mlen] = 0;
  return marked;
}

void MEDDLY::expert_forest::dump(FILE *s) const
{
  int nwidth = digits(a_last);
  for (int p=0; p<=a_last; p++) {
    fprintf(s, "%*d\t", nwidth, p);
    showNode(s, p, 1);
    fprintf(s, "\n");
    fflush(s);
  }
}

void MEDDLY::expert_forest::dumpInternal(FILE *s) const
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

  for (int i=1; i<=getNumVariables(); i++) {
    dumpInternalLevel(s, i);
    if (isForRelations()) dumpInternalLevel(s, -i);
  }
 
  unique->show(s);
  fflush(s);
}

void MEDDLY::expert_forest::dumpInternalLevel(FILE *s, int k) const
{
  levels[k].dumpInternal(s);
}

void MEDDLY::expert_forest::dumpUniqueTable(FILE *s) const
{
  unique->show(s);
}

void MEDDLY::expert_forest::showNodeGraph(FILE *s, int p) const
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
        char primed = (k>0) ? ' ' : '\'';
        if (v->getName()) {
          fprintf(s, "Level: %s%c\n", v->getName(), primed);
        } else {
          fprintf(s, "Level: %d%c\n", ABS(k), primed);
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

void MEDDLY::expert_forest
::reportMemoryUsage(FILE * s, const char* pad, int verb)
{
  if (verb>0) {
    fprintf(s, "%sPeak Nodes:             %ld\n", pad, getPeakNumNodes());
    fprintf(s, "%sActive Nodes:           %ld\n", pad, getCurrentNumNodes());
  }
#if 0
  unsigned count = 0;
  for (int i = 1; i <= getLastNode(); ++i) if (isActiveNode(i)) ++count;
  fprintf(s, "%cActive Nodes (manual):\t\t%d\n", pad, count);
  fprintf(s, "%c%cZombie Nodes:\t\t%d\n", pad, filler,
      getZombieNodeCount());
  fprintf(s, "%c%cTemp Nodes:\t\t%d\n", pad, filler, getTempNodeCount());
  fprintf(s, "%c%cOrphan Nodes:\t\t%d\n", pad, filler,
      getOrphanNodeCount());
#endif
  if (verb>2) {
    fprintf(s, "%sReclaimed Nodes:        %ld\n", pad, stats.reclaimed_nodes);
  }
  fprintf(s, "%sMem Used:               %ld\n", pad, getCurrentMemoryUsed());
  fprintf(s, "%sPeak Mem Used:          %ld\n", pad, getPeakMemoryUsed());
  if (verb>1) {
    fprintf(s, "%sMem Allocated:          %ld\n", pad,
      getCurrentMemoryAllocated());
    fprintf(s, "%sPeak Mem Allocated:     %ld\n",
      pad, getPeakMemoryAllocated());
  }
  if (verb>3) {
    fprintf(s, "%sUnique Tbl Mem Used:    %u\n", pad, 
      unique->getMemUsed());
  }
  if (verb>5) {
    fprintf(s, "%sCompactions:            %ld\n", pad, stats.num_compactions);
  }
  if (verb>7) {
    long holemem = 0;
    for (int i=getMinLevelIndex(); i<=getNumVariables(); i++) {
      holemem += levels[i].getHoleSlots();
    }
    holemem *= sizeof(int);
    fprintf(s, "%sHole Memory Usage:\t%ld\n", pad, holemem);
  }

  if (verb>6) {
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
  }
}

unsigned MEDDLY::expert_forest::hashNode(int h) const 
{
  int* d = getNodeAddress(h);
  int length = d[2];
  MEDDLY_DCASSERT(length != 0);
  d += 3;

  hash_stream s;
  s.start(getNodeLevel(h));

  if (length > 0) {
    // Full node
    for (int i=0; i<length; i++) {
      if (0==d[i]) continue;
      s.push(i, d[i]);
    } 
  } else {
    int* indexes = d - length;
    for (int i=0; i<-length; i++) {
      s.push(d[i], indexes[i]); // currently, nodes are stored "indexes first"
      // s.push(indexes[i], d[i]);
    }
  }
  return s.finish();
}


MEDDLY::expert_forest::nodeReader*
MEDDLY::expert_forest::initNodeReader(int node)
{
  const node_data& n = getNode(node);
  int nsize = getLevelSize(n.level);
  nodeReader* nr;
  if (free_reader[n.level]) {
    nr = free_reader[n.level];
    free_reader[n.level] = nr->next;
    nr->next = 0;
  } else {
    nr = new nodeReader(n.level);
  }
  nr->resize(nsize);
  level_data& ld = levels[n.level];
  if (ld.isFull(n.offset)) {
    int i;
    int stop = ld.fullSizeOf(n.offset);
    const int* dn = ld.fullDownOf(n.offset);
    for (i=0; i<stop; i++) {
      nr->buffer[i] = dn[i];
    }
    for (; i<nsize; i++) {
      nr->buffer[i] = 0;
    }
  } else {
    int i = 0;
    int nnz = ld.sparseSizeOf(n.offset);
    const int* dn = ld.sparseDownOf(n.offset);
    const int* ix = ld.sparseIndexesOf(n.offset);
    for (int z=0; z<nnz; z++) {
      for (; i<ix[z]; i++) nr->buffer[i] = 0;
      nr->buffer[i] = dn[z];
      i++;
    }
    for (; i<nsize; i++) nr->buffer[i] = 0;
  }
  return nr;
}

MEDDLY::expert_forest::nodeReader*
MEDDLY::expert_forest::initRedundantReader(int k, int node)
{
  int nsize = getLevelSize(k);
  nodeReader* nr;
  if (free_reader[k]) {
    nr = free_reader[k];
    free_reader[k] = nr->next;
    nr->next = 0;
  } else {
    nr = new nodeReader(k);
  }
  nr->resize(nsize);
  for (int i=0; i<nsize; i++) 
    nr->buffer[i] = node;
  return nr;
}

MEDDLY::expert_forest::nodeReader*
MEDDLY::expert_forest::initIdentityReader(int k, int i, int node)
{
  int nsize = getLevelSize(k);
  nodeReader* nr;
  if (free_reader[k]) {
    nr = free_reader[k];
    free_reader[k] = nr->next;
    nr->next = 0;
  } else {
    nr = new nodeReader(k);
  }
  nr->resize(nsize);
  for (int j=0; j<i; j++) nr->buffer[j] = 0;
  nr->buffer[i] = node;
  for (int j=i+1; j<nsize; j++) nr->buffer[j] = 0;
  return nr;
}

void MEDDLY::expert_forest::recycle(nodeReader *r)
{
  MEDDLY_DCASSERT(r);
  r->next = free_reader[r->level];
  free_reader[r->level] = r;
}

int MEDDLY::expert_forest::getSingletonIndex(int n, int &down) const
{
  const node_data& node = getNode(n);
  MEDDLY_DCASSERT(node.level);
  const level_data &ld = levels[node.level];
  if (ld.isFull(node.offset)) {
    // full node
    int fs = ld.fullSizeOf(node.offset);
    const int* dn = ld.fullDownOf(node.offset);
    for (int i=fs-2; i>=0; i--) if (dn[i]) return -1;
    down = dn[fs-1];
    return fs-1;
  } else {
    // sparse node --- easy
    if (ld.sparseSizeOf(node.offset) != 1) return -1;
    down = ld.sparseDownOf(node.offset)[0];
    return ld.sparseIndexesOf(node.offset)[0];
  }
}

void MEDDLY::expert_forest::garbageCollect()
{
  if (performing_gc) return;
  performing_gc = true;

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
  } else {
    // remove the stale nodes entries from caches
    removeStaleComputeTableEntries();
  }

#ifdef DEBUG_GC
  printf("Compacting levels...\n");
  fflush(stdout);
#endif

  compactMemory(); // might want to remove this

#ifdef DEBUG_GC
  printf("  done compacting.\n");
  fflush(stdout);
#endif

  performing_gc = false;
}

void MEDDLY::expert_forest::compactMemory()
{
  for (int i=getMinLevelIndex(); i<=getNumVariables(); i++) {
    levels[i].compactLevel = true;
    levels[i].compact(address);
  }
}

void MEDDLY::expert_forest::showInfo(FILE* s, int verb)
{
  // Show forest with appropriate level of detail
  if (1==verb)  dump(s);
  else          dumpInternal(s); 
  fprintf(s, "DD stats:\n");
  reportMemoryUsage(s, "    ", verb);
  fprintf(s, "Unique table stats:\n");
  fprintf(s, "\t%-24s%u\n", "Current size:", unique->getSize());
  fprintf(s, "\t%-24s%u\n", "Current entries:", unique->getNumEntries());
}

int MEDDLY::expert_forest::reduceNode(int node)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::normalizeAndReduceNode(int& node, int& ev)
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::normalizeAndReduceNode(int& node, float& ev)
{
  throw error(error::TYPE_MISMATCH);
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                       Protected  methods                       '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


//
// ------------------------------------------------------------------
//  Still need to be organized:
// ------------------------------------------------------------------
//

// TODO: make use of pointers to speed this up.
int MEDDLY::expert_forest::getDownPtr(int p, int i) const {
  MEDDLY_DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return p;
  MEDDLY_DCASSERT(i >= 0);
  if (isFullNode(p)) {
    // full or trunc-full node
    if (getFullNodeSize(p) > i) return getFullNodeDownPtr(p, i);

    // full or trunc-full node, but i lies after the last non-zero downpointer
    return 0;
  } else {
    // sparse node
    // binary search to find the index ptr corresponding to i
    int start = 0;
    int stop = getSparseNodeSize(p) - 1;
    
    // the index ptr corresponding to i has been compressed i.e. its a zero.
    if (getSparseNodeIndex(p, start) > i) return 0;
    if (getSparseNodeIndex(p, stop) < i) return 0;

    int mid = (start + stop)/2;
    while(start < stop) {
      if (getSparseNodeIndex(p, mid) == i)
        return getSparseNodeDownPtr(p, mid);
      if (getSparseNodeIndex(p, mid) > i)
        stop = mid - 1;
      else
        start = mid + 1;
      mid = (start + stop)/2;
    }
    if (getSparseNodeIndex(p, mid) == i)
      return getSparseNodeDownPtr(p, mid);
   
    // index ptr not found - it was compressed because it was a zero.
    return 0;
  }
}


bool MEDDLY::expert_forest::getDownPtrs(int p, std::vector<int>& dptrs) const {
  if (!isActiveNode(p) || isTerminalNode(p)) return false;

  const int* ptrs = 0;
  assert(getDownPtrs(p, ptrs));

  if (isFullNode(p)) {
    int size = getFullNodeSize(p);
    if (dptrs.size() < unsigned(size)) dptrs.resize(size, 0);
    const int* end = ptrs + size;
    std::vector<int>::iterator iter = dptrs.begin();
    while (ptrs != end) { *iter++ = *ptrs++; }
  }
  else {
    int nnz = getSparseNodeSize(p);
    const int* index = 0;
    assert(getSparseNodeIndexes(p, index));
    const int* end = ptrs + nnz;

    int size = index[nnz-1] + 1;
    if (dptrs.size() < unsigned(size)) dptrs.resize(size, 0);

    while (ptrs != end) { dptrs[*index++] = *ptrs++; }
  }
  return true;
}

bool MEDDLY::expert_forest::isStale(int h) const {
  return
    isMarkedForDeletion() || (
      isTerminalNode(h)
      ? discardTemporaryNodesFromComputeCache()
      : isPessimistic()
        ? isZombieNode(h)
        : (getInCount(h) == 0)
    );
}

#ifndef INLINED_REALS

float MEDDLY::expert_forest::getReal(int term) const
{
  return (term == 0)? 0.0: *((float*)&(term <<= 1));
}

int MEDDLY::expert_forest::getTerminalNode(float a) const
{
  return (a == 0.0)? 0: (*((int*)((void*)&a)) >> 1) | 0x80000000;
}

#endif


unsigned MEDDLY::expert_forest::getNodeCount(int p) const
{
  std::set<int> discovered;
  std::queue<int> toExpand;

  if (p != 0 && p != -1) {
    toExpand.push(p);
    discovered.insert(p);
  }

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
        if (temp == 0 || temp == -1) continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
    else {
      const int sz = getSparseNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getSparseNodeDownPtr(p, i);
        if (temp == 0 || temp == -1) continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
  }

  // Add 2 to discovered.size() for terminals 0 and -1.
  return discovered.size() + 2;
}


unsigned MEDDLY::expert_forest::getEdgeCount(int p, bool countZeroes) const
{
  std::set<int> discovered;
  std::queue<int> toExpand;

  toExpand.push(p);
  discovered.insert(p);

  unsigned count = 0;

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
        if (countZeroes) count++;
        else if (temp != 0) count++;
        if (temp == 0 || temp == -1)  continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
    else {
      const int sz = getSparseNodeSize(p);
      for (int i = 0; i < sz; ++i)
      {
        int temp = getSparseNodeDownPtr(p, i);
        if (countZeroes) count++;
        else if (temp != 0) count++;
        if (temp == 0 || temp == -1)  continue;
        // insert into discovered and toExpand if new
        if (discovered.find(temp) == discovered.end()) {
          toExpand.push(temp);
          discovered.insert(temp);
        }
      }
    }
  }

  return count;
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '          Methods for managing  available node handles          '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest::expandHandleList()
{
  // increase size by 50%
  int delta = a_size / 2;
  MEDDLY_DCASSERT(delta>=0);
  address = (node_data*) realloc(address, (a_size+delta) * sizeof(node_data));
  if (0==address) {
    throw error(error::INSUFFICIENT_MEMORY);
  }
  stats.incMemAlloc(delta * sizeof(node_data));
  memset(address + a_size, 0, delta * sizeof(node_data));
  a_size += delta;
  a_next_shrink = a_size / 2;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %d\n", a_size);
#endif
}

void MEDDLY::expert_forest::shrinkHandleList()
{
  // Determine new size
  int new_size = a_min_size;
  while (a_last > new_size) new_size += new_size/2;
  int delta = a_size - new_size;
  if (0==delta) {
    a_next_shrink = 0;
    return;
  }
  // shrink the array
  MEDDLY_DCASSERT(delta>=0);
  MEDDLY_DCASSERT(a_size-delta>=a_min_size);
  address = (node_data*) realloc(address, new_size * sizeof(node_data));
  if (0==address) {
    throw error(error::INSUFFICIENT_MEMORY);
  }
  stats.decMemAlloc(delta * sizeof(node_data));
  a_size -= delta;
  a_next_shrink = a_size / 2;
  MEDDLY_DCASSERT(a_last < a_size);
  // rebuild the free list
  a_unused = 0;
  for (int i=a_last; i; i--) {
    if (address[i].isDeleted()) {
      address[i].setNextDeleted(a_unused);
      a_unused = i;
    }
  }
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank address array, new size %d\n", a_size);
#endif
}

int MEDDLY::expert_forest
::createReducedHelper(int in, const nodeBuilder &nb, bool &u)
{
  throw error(error::TYPE_MISMATCH);
}

