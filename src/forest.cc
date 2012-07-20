
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
// #define TRACK_DELETIONS

// #define DEBUG_CREATE_REDUCED

// Thoroughly check reference counts.
// Very slow.  Use only for debugging.
// #define VALIDATE_INCOUNTS
// #define VALIDATE_INCOUNTS_ON_DELETE

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

MEDDLY::forest::edge_visitor::edge_visitor()
{
}

MEDDLY::forest::edge_visitor::~edge_visitor()
{
}

#if 0
// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_reader  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_reader::node_reader()
{
  down = 0;
  index = 0;
  edge = 0;
  alloc = 0;
  ealloc = 0;
  size = 0;
  nnzs = 0;
  level = 0;
}

MEDDLY::node_reader::~node_reader()
{
  clear();
}

void MEDDLY::node_reader::clear()
{
  free(down);
  free(index);
  free(edge);
  down = 0;
  index = 0;
  edge = 0;
  alloc = 0;
  ealloc = 0;
  size = 0;
  nnzs = 0;
  level = 0;
}

/*
void MEDDLY::node_reader::dump(FILE* s) const
{
  if (is_full) {
    fprintf(s, "[%d", down[0]);
    for (int i=1; i<size; i++)
      fprintf(s, ", %d", down[i]);
    fprintf(s, "]");
  } else {
    fprintf(s, "(%d:%d", index[0], down[0]);
    for (int z=1; z<nnzs; z++) 
      fprintf(s, ", %d:%d", index[z], down[z]);
    fprintf(s, ")");
  }
}
*/

void MEDDLY::node_reader::resize(int k, int ns, char es, bool full)
{
  level = k;
  is_full = full;
  size = ns;
  edge_bytes = es;
  if (size > alloc) {
    int nalloc = ((ns/8)+1)*8;
    MEDDLY_DCASSERT(nalloc > ns);
    MEDDLY_DCASSERT(nalloc>0);
    MEDDLY_DCASSERT(nalloc>alloc);
    down = (int*) realloc(down, nalloc*sizeof(int));
    if (0==down) throw error(error::INSUFFICIENT_MEMORY);
    index = (int*) realloc(index, nalloc*sizeof(int));
    if (0==index) throw error(error::INSUFFICIENT_MEMORY);
    alloc = nalloc;
  }
  if (edge_bytes * size > ealloc) {
    int nalloc = ((edge_bytes * size)/8+1)*8;
    MEDDLY_DCASSERT(nalloc>0);
    MEDDLY_DCASSERT(nalloc>ealloc);
    edge = realloc(edge, nalloc);
    if (0==edge) throw error(error::INSUFFICIENT_MEMORY);
    ealloc = nalloc;
  }
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
  parent->initNodeReader(thisnode, node, false);
  h = parent->hashNode(node);
}

MEDDLY::expert_forest::nodeFinder::~nodeFinder()
{
}

bool MEDDLY::expert_forest::nodeFinder::equals(int p)
{
  parent->initNodeReader(compare, p, false);

  if (thisnode.getLevel() != compare.getLevel()) return false;
  if (thisnode.getNNZs() != compare.getNNZs()) return false;

  for (int z=0; z<thisnode.getNNZs(); z++) {
    if (thisnode.d(z) != compare.d(z)) return false;
    if (thisnode.i(z) != compare.i(z)) return false;
  }
  // TBD : edge values

  return true;
}
#endif

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
::init(expert_forest* p, int k)
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
  edgeSize = parent->edgeSize(k);
  unhashedHeader = parent->unhashedHeaderSize(k);
  hashedHeader = parent->hashedHeaderSize(k);
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

void MEDDLY::expert_forest::level_data::compact(nodeData* address)
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
// *               expert_forest::nodecounter methods               *
// *                                                                *
// ******************************************************************

MEDDLY::expert_forest::nodecounter::nodecounter(expert_forest *p, int* c)
 : edge_visitor()
{
  parent = p;
  counts = c;
}

MEDDLY::expert_forest::nodecounter::~nodecounter()
{
  // DO NOT delete counts.
}

void MEDDLY::expert_forest::nodecounter::visit(dd_edge &e)
{
  int n = e.getNode();
  if (parent->isTerminalNode(n)) return;
  MEDDLY_DCASSERT(n>0);
  MEDDLY_DCASSERT(n<=parent->getLastNode());
  counts[n]++;
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
  address = (nodeData *) malloc(a_size * sizeof(nodeData));
  if (0 == address) throw error(error::INSUFFICIENT_MEMORY);
  stats.incMemAlloc(a_size * sizeof(nodeData));
  memset(address, 0, a_size * sizeof(nodeData));
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
    raw_builders = new node_builder[2*N+1];
    builders = raw_builders + N;
  } else {
    raw_builders = new node_builder[N+1];
    builders = raw_builders;
  }

  //
  // Initialize misc. protected data
  //
  terminalNodesAreStale = false;

  //
  // Initialize misc. private data
  //
  unique = new unique_table(this);
  performing_gc = false;
  in_validate = 0;
  in_val_size = 0;
}


MEDDLY::expert_forest::~expert_forest() 
{
  // Address array
  free(address);

  // Level array
  delete[] raw_levels;

  // builders array
  delete[] raw_builders;

  // unique table
  delete unique;

  // Misc. private data
  free(in_validate);
}

void MEDDLY::expert_forest::initializeForest()
{
  //
  // Initialize level array
  //
  for (int k=getMinLevelIndex(); k<=getNumVariables(); k++) {
    levels[k].init(this, k);
  }

  //
  // Initialize builders array
  //
  for (int k=getMinLevelIndex(); k<=getNumVariables(); k++) {
    builders[k].init(k, this);
  }
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                   public debugging   methods                   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

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

void MEDDLY::expert_forest::validateIncounts(bool exact)
{
#ifdef VALIDATE_INCOUNTS
  static int idnum = 0;
  idnum++;

  // Inspect every active node's down pointers to determine
  // the incoming count for every active node.
  
  int sz = getLastNode() + 1;
  if (sz > in_val_size) {
    in_validate = (int*) realloc(in_validate, a_size * sizeof(int));
    in_val_size = a_size;
  }
  MEDDLY_DCASSERT(sz <= in_val_size);
  memset(in_validate, 0, sizeof(int) * sz);
  node_reader P;
  for (int i = 1; i < sz; ++i) {
    MEDDLY_DCASSERT(!isTerminalNode(i));
    if (!isActiveNode(i)) continue;
    initNodeReader(P, i, false);

    // add to reference counts
    for (int z=0; z<P.getNNZs(); z++) {
      if (isTerminalNode(P.d(z))) continue;
      MEDDLY_CHECK_RANGE(0, P.d(z), sz);
      in_validate[P.d(z)]++;
    }
  } // for i

  // Add counts for registered dd_edges
  nodecounter foo(this, in_validate);
  visitRegisteredEdges(foo);
  
  // Validate the incoming count stored with each active node using the
  // in_count array computed above
  for (int i = 1; i < sz; ++i) {
    MEDDLY_DCASSERT(!isTerminalNode(i));
    if (!isActiveNode(i)) continue;
    bool fail = exact
      ?  in_validate[i] != readInCount(i)
      :  in_validate[i] >  readInCount(i);
    if (fail) {
      printf("Validation #%d failed\n", idnum);
      printf("For node %d\n    we got %d\n    node says %d\n",
        i, in_validate[i], readInCount(i));
      dump(stdout);
      throw error(error::MISCELLANEOUS);
    }
    // Note - might not be exactly equal
    // because there could be dd_edges that refer to nodes
    // and we didn't count them.
  }

#ifdef TRACK_DELETIONS
  printf("Incounts validated #%d\n", idnum);
#endif
#endif
}


// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                      public handy methods                      '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

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

  node_reader M;

  // Breadth-first search
  for (int mexpl=0; mexpl<mlen; mexpl++) {
    // explore node marked[mexpl]
    initNodeReader(M, marked[mexpl], false);
    for (int i=0; i<M.getNNZs(); i++) {
      if (isTerminalNode(M.d(i))) continue;
      MEDDLY_CHECK_RANGE(0, M.d(i)-1, a_last);
      if (inList[M.d(i)]) continue;
      // add dn to list
      if (mlen+1 >= msize) { 
          // expand.  Note we're leaving an extra slot
          // at the end, for the terminal 0.
          msize += 1024;
          marked = (int*) realloc(marked, msize*sizeof(int));
          if (0==marked) throw error(error::INSUFFICIENT_MEMORY);
      }
      inList[M.d(i)] = true;
      marked[mlen] = M.d(i);
      mlen++;
    } // for i
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

int MEDDLY::expert_forest::getNodeCount(int p) const
{
  int* list = markNodesInSubgraph(p, true);
  if (0==list) return 0;
  int i;
  for (i=0; list[i]; i++) { }
  free(list);
  return i;
}

int MEDDLY::expert_forest::getEdgeCount(int p, bool countZeroes) const
{
  int* list = markNodesInSubgraph(p, true);
  if (0==list) return 0;
  int ec=0;
  node_reader M;
  for (int i=0; list[i]; i++) {
    initNodeReader(M, list[i], countZeroes);
    ec += countZeroes ? M.getSize() : M.getNNZs();
  }
  free(list);
  return ec;
}

void MEDDLY::expert_forest::showNode(FILE* s, int p, int verbose) const
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
  const nodeData& node = getNode(p);
  const level_data &ld = levels[node.level];
  if (verbose) {
    // node: was already written.
    const variable* v = getDomain()->getVar(ABS(node.level));
    if (v->getName()) {
      fprintf(s, " level: %s", v->getName());
    } else {
      fprintf(s, " level: %d", ABS(node.level));
    }
    if (getNodeLevel(p) < 0)
      fprintf(s, "'");
    else
      fprintf(s, " ");
    fprintf(s, " in: %d", ld.countOf(node.offset));
    fprintf(s, " cc: %d", node.cache_count);
  } else {
    fprintf(s, "node: %d", p);
  }
  node_reader *R = node_reader::useReader();
  initNodeReader(*R, p, ld.isFull(node.offset));
  if (R->isFull()) {
    // Full node
    int size = ld.fullSizeOf(node.offset);
    if (verbose) fprintf(s, " size: %d", size);
    fprintf(s, " down: [");
    for (int i=0; i<size; i++) {
      if (i) fprintf(s, "|"); 
      if (R->hasEdges()) {
        fprintf(s, "<");
        showEdgeValue(s, R->rawEdges(), i);
        fprintf(s, ", ");
      } 
      if (isTerminalNode(R->d(i))) {
        showTerminal(s, R->d(i));
      } else {
        fprintf(s, "%d", R->d(i));
      }
      if (R->hasEdges()) fprintf(s, ">");
    } // for i
    fprintf(s, "]");
  } else {
    // Sparse node
    int nnz = ld.sparseSizeOf(node.offset);
    if (verbose) fprintf(s, " nnz : %d", nnz);
    fprintf(s, " down: (");
    for (int z=0; z<nnz; z++) {
      if (z) fprintf(s, ", ");
      fprintf(s, "%d:", R->i(z));
      if (R->hasEdges()) {
        fprintf(s, "<");
        showEdgeValue(s, R->rawEdges(), z);
        fprintf(s, ", ");
      } 
      if (isTerminalNode(R->d(z))) {
        showTerminal(s, R->d(z));
      } else {
        fprintf(s, "%d", R->d(z));
      }
      if (R->hasEdges()) fprintf(s, ">");
    } // for z
    fprintf(s, ")");
  }
  node_reader::recycle(R);

  // show extra header stuff
  if (ld.unhashedHeader) {
    showUnhashedHeader(s, ld.unhashedHeaderOf(node.offset));
  }
  if (ld.hashedHeader) {
    showHashedHeader(s, ld.hashedHeaderOf(node.offset));
  }
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
::reportMemoryUsage(FILE * s, const char* pad, int verb) const
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

unsigned MEDDLY::expert_forest::hashNode(int p) const 
{
  hash_stream s;
  const nodeData& node = getNode(p);
  MEDDLY_DCASSERT(node.level);
  const level_data &ld = levels[node.level];
  s.start(node.level);

  for (int e=0; e<ld.hashedHeader; e++) {
    s.push(ld.hashedHeaderOf(node.offset)[e]);
  }
  
  if (ld.isSparse(node.offset)) {
    int nnzs = ld.sparseSizeOf(node.offset);
    const int* down = ld.sparseDownOf(node.offset);
    const int* indexes = ld.sparseIndexesOf(node.offset);
    if (areEdgeValuesHashed(node.level)) {
      int edge_bytes = edgeSize(node.level) * sizeof(int);
      const char* edge = (const char*) ld.sparseEdgeOf(node.offset);
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(down[z]);
        int* ep = (int*) (edge + z * edge_bytes);
        s.push(indexes[z], down[z], *ep);
      }
    } else {
      for (int z=0; z<nnzs; z++) {
        MEDDLY_DCASSERT(down[z]);
        s.push(indexes[z], down[z]);
      }
    }
  } else {
    int size = ld.fullSizeOf(node.offset);
    const int* down = ld.fullDownOf(node.offset);
    if (areEdgeValuesHashed(node.level)) {
      int edge_bytes = edgeSize(node.level) * sizeof(int);
      const char* edge = (const char*) ld.fullEdgeOf(node.offset);
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        int* ep = (int*) (edge + i * edge_bytes);
        s.push(i, down[i], *ep);
      }
    } else {
      for (int i=0; i<size; i++) {
        if (0==down[i]) continue;
        s.push(i, down[i]);
      }
    }
  }
  return s.finish();
}

int MEDDLY::expert_forest::getSingletonIndex(int n, int &down) const
{
  const nodeData& node = getNode(n);
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


int MEDDLY::expert_forest::getDownPtr(int p, int i) const 
{
  MEDDLY_DCASSERT(i>=0);
  const nodeData& node = getNode(p);
  MEDDLY_DCASSERT(node.level);
  const level_data &ld = levels[node.level];
  if (ld.isFull(node.offset)) {
    // full node - super easy
    int fs = ld.fullSizeOf(node.offset);
    if (i>=fs) return 0;
    const int* dn = ld.fullDownOf(node.offset);
    return dn[i];
  }
  // Node must be sparse; do a binary search
  int low = 0;  // smallest where i might be
  int high = ld.sparseSizeOf(node.offset);  // smallest where i isn't
  const int* ix = ld.sparseIndexesOf(node.offset);
  while (low < high) {
    int mid = (low+high)/2;
    if (ix[mid] == i) {
      return ld.sparseDownOf(node.offset)[mid];
    }
    if (ix[mid] < i)  low = mid+1;
    else              high = mid;
  }
  return 0;
}


void MEDDLY::expert_forest::getDownPtr(int p, int i, int& ev, int& dn) const 
{
  MEDDLY_DCASSERT(i>=0);
  const nodeData& node = getNode(p);
  MEDDLY_DCASSERT(node.level);
  const level_data &ld = levels[node.level];
  if (ld.isFull(node.offset)) {
    // full node - super easy
    int fs = ld.fullSizeOf(node.offset);
    if (i<fs) {
      dn = (ld.fullDownOf(node.offset))[i];
      ev = ((const int*)ld.fullEdgeOf(node.offset))[i];
    } else {
      dn = 0;
      ev = 0;
    }
    return;
  }
  // Node must be sparse; do a binary search
  int low = 0;  // smallest where i might be
  int high = ld.sparseSizeOf(node.offset);  // smallest where i isn't
  const int* ix = ld.sparseIndexesOf(node.offset);
  while (low < high) {
    int mid = (low+high)/2;
    if (ix[mid] == i) {
      dn = ld.sparseDownOf(node.offset)[mid];
      ev = ((const int*)ld.sparseEdgeOf(node.offset))[mid];
      return;
    }
    if (ix[mid] < i)  low = mid+1;
    else              high = mid;
  }
  dn = 0;
  ev = 0;
}


void MEDDLY::expert_forest::getDownPtr(int p, int i, float& ev, int& dn) const
{
  MEDDLY_DCASSERT(i>=0);
  const nodeData& node = getNode(p);
  MEDDLY_DCASSERT(node.level);
  const level_data &ld = levels[node.level];
  if (ld.isFull(node.offset)) {
    // full node - super easy
    int fs = ld.fullSizeOf(node.offset);
    if (i<fs) {
      dn = (ld.fullDownOf(node.offset))[i];
      ev = ((const float*)ld.fullEdgeOf(node.offset))[i];
    } else {
      dn = 0;
      ev = 0;
    }
    return;
  }
  // Node must be sparse; do a binary search
  int low = 0;  // smallest where i might be
  int high = ld.sparseSizeOf(node.offset);  // smallest where i isn't
  const int* ix = ld.sparseIndexesOf(node.offset);
  while (low < high) {
    int mid = (low+high)/2;
    if (ix[mid] == i) {
      dn = ld.sparseDownOf(node.offset)[mid];
      ev = ((const float*)ld.sparseEdgeOf(node.offset))[mid];
      return;
    }
    if (ix[mid] < i)  low = mid+1;
    else              high = mid;
  }
  dn = 0;
  ev = 0;
}



// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                   methods for  reading nodes                   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


void MEDDLY::expert_forest
::initNodeReader(node_reader &nr, int node, bool full) const
{
  const nodeData &n = getNode(node);
  level_data& ld = levels[n.level];
  nr.resize(this, n.level, getLevelSize(n.level), full);
  if (ld.isFull(n.offset)) {
    int i;
    int stop = ld.fullSizeOf(n.offset);
    const int* dn = ld.fullDownOf(n.offset);
    if (full) {
      memcpy(nr.down, dn, stop * sizeof(int));
      int* nrdext = nr.down + stop;
      memset(nrdext, 0, (nr.size-stop) * sizeof(int));
      if (nr.edge_bytes) {
        const void* ev = ld.fullEdgeOf(n.offset);
        memcpy(nr.edge, ev, stop * nr.edge_bytes);
        void* evext = (char*)nr.edge + (stop * nr.edge_bytes);
        memset(evext, 0, (nr.size-stop) * nr.edge_bytes);
      }
    } else {
      int& z = nr.nnzs;
      z = 0;
      if (nr.edge_bytes) {
        void* nev = nr.edge;
        for (i=0; i<stop; i++) {
          if (0==dn[i]) continue;
          nr.down[z] = dn[i];
          nr.index[z] = i;
          const void* ev = (char*)ld.fullEdgeOf(n.offset) + i * nr.edge_bytes;
          memcpy(nev, ev, nr.edge_bytes);
          nev = (char*)nev + nr.edge_bytes;
          z++;
        } // for i
      } else {
        for (i=0; i<stop; i++) if (dn[i]) {
          nr.down[z] = dn[i];
          nr.index[z] = i;
          z++;
        } // for i
      } // if ev
    }
  } else {
    int i = 0;
    int nnz = ld.sparseSizeOf(n.offset);
    const int* dn = ld.sparseDownOf(n.offset);
    const int* ix = ld.sparseIndexesOf(n.offset);
    if (full) {
      if (nr.edge_bytes) {
        const void* ev = ld.sparseEdgeOf(n.offset);
        memset(nr.down, 0, nr.size * sizeof(int));
        memset(nr.edge, 0, nr.size * nr.edge_bytes);
        for (int z=0; z<nnz; z++) {
          nr.down[ix[z]] = dn[z];
          int off = ix[z] * nr.edge_bytes;
          memcpy((char*)nr.edge + off, (char*)ev + off, nr.edge_bytes);
          i++;
        }
      } else {
        for (int z=0; z<nnz; z++) {
          for (; i<ix[z]; i++) nr.down[i] = 0;
          nr.down[i] = dn[z];
          i++;
        }
        for (; i<nr.size; i++) nr.down[i] = 0;
      } // if ev
    } else {
      nr.nnzs = nnz;
      memcpy(nr.down, dn, nnz * sizeof(int));
      memcpy(nr.index, ix, nnz * sizeof(int));
      if (nr.edge_bytes) {
        memcpy(nr.edge, ld.sparseEdgeOf(n.offset), 
          nnz * nr.edge_bytes);
      }
    }
  }
}

void MEDDLY::expert_forest
::initRedundantReader(node_reader &nr, int k, int node, bool full) const
{
  MEDDLY_DCASSERT(0==levels[k].edgeSize);
  int nsize = getLevelSize(k);
  nr.resize(this, k, nsize, full);
  for (int i=0; i<nsize; i++) 
    nr.down[i] = node;
  if (!full) {
    for (int i=0; i<nsize; i++) nr.index[i] = i;
    nr.nnzs = nsize;
  }
}

void MEDDLY::expert_forest
::initRedundantReader(node_reader &nr, int k, int ev, int np, bool full) const
{
  MEDDLY_DCASSERT(1==levels[k].edgeSize);
  int nsize = getLevelSize(k);
  nr.resize(this, k, nsize, full);
  for (int i=0; i<nsize; i++) {
    nr.down[i] = np;
    ((int*)nr.edge)[i] = ev;
  }
  if (!full) {
    for (int i=0; i<nsize; i++) nr.index[i] = i;
    nr.nnzs = nsize;
  }
}

void MEDDLY::expert_forest
::initIdentityReader(node_reader &nr, int k, int i, int node, bool full) const
{
  MEDDLY_DCASSERT(0==levels[k].edgeSize);
  int nsize = getLevelSize(k);
  if (full) {
    nr.resize(this, k, nsize, full);
    memset(nr.down, 0, nsize * sizeof(int));
    nr.down[i] = node;
  } else {
    nr.resize(this, k, 1, full);
    nr.nnzs = 1;
    nr.down[0] = node;
    nr.index[0] = i;
  }
}


void MEDDLY::expert_forest
::initIdentityReader(node_reader &nr, int k, int i, int ev, int node, 
  bool full) const
{
  MEDDLY_DCASSERT(1==levels[k].edgeSize);
  int nsize = getLevelSize(k);
  if (full) {
    nr.resize(this, k, nsize, full);
    memset(nr.down, 0, nsize * sizeof(int));
    memset(nr.edge, 0, nsize * sizeof(int));
    nr.down[i] = node;
    ((int*)nr.edge)[i] = ev;
  } else {
    nr.resize(this, k, 1, full);
    nr.nnzs = 1;
    nr.down[0] = node;
    ((int*)nr.edge)[0] = ev;
    nr.index[0] = i;
  }
}



// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '   public virtual methods provided here,  no need to override   '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


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


// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '    virtual methods to be overridden by some derived classes    '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

/*
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

int MEDDLY::expert_forest
::createReducedHelper(int in, const node_builder &nb, bool &u)
{
  throw error(error::TYPE_MISMATCH);
}

bool MEDDLY::expert_forest::areDuplicates(int, const node_builder&) const
{
  throw error(error::NOT_IMPLEMENTED);
}

bool MEDDLY::expert_forest::areDuplicates(int, const node_reader&) const
{
  throw error(error::NOT_IMPLEMENTED);
}
*/

void MEDDLY::expert_forest::normalize(node_builder &nb, int& ev) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::normalize(node_builder &nb, float& ev) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::showTerminal(FILE* s, int tnode) const
{
  throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::expert_forest::showEdgeValue(FILE* s, const void* edge, int i) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::showHashedHeader(FILE* s, const int* hh) const
{
  throw error(error::TYPE_MISMATCH);
}

void MEDDLY::expert_forest::showUnhashedHeader(FILE* s, const int* uh) const
{
  throw error(error::TYPE_MISMATCH);
}

// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
// '                                                                '
// '                        private  methods                        '
// '                                                                '
// ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

void MEDDLY::expert_forest::handleNewOrphanNode(int p)
{
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(readInCount(p) == 0);

  // insted of orphan_nodes++ here; do it only when the orphan is not going
  // to get deleted or converted into a zombie

  // Two possible scenarios:
  // (1) a reduced node, or
  // (2) a temporary node ready to be deleted.
  // MEDDLY_DCASSERT(isReducedNode(p) || getCacheCount(p) == 0);

  if (getCacheCount(p) == 0) {
    // delete node
    // this should take care of the temporary nodes also
#ifdef TRACK_DELETIONS
    printf("Deleting node %d\n\t", p);
    showNode(stdout, p);
    printf("\n");
    fflush(stdout);
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

void MEDDLY::expert_forest::deleteOrphanNode(int p) 
{
  MEDDLY_DCASSERT(!isPessimistic());
  MEDDLY_DCASSERT(getCacheCount(p) == 0 && readInCount(p) == 0);
#ifdef TRACK_DELETIONS
  printf("Deleting node %d\n\t", p);
  showNode(stdout, p);
  printf("\n");
  fflush(stdout);
#endif
  stats.orphan_nodes--;
  deleteNode(p);
}

void MEDDLY::expert_forest::deleteNode(int p)
{
  MEDDLY_DCASSERT(!isTerminalNode(p));
  MEDDLY_DCASSERT(readInCount(p) == 0);
  MEDDLY_DCASSERT(isActiveNode(p));

#ifdef VALIADE_INCOUNTS_ON_DELETE
  validateIncounts(false);
#endif

  int* foo = getNodeAddress(p);
  int k = getNodeLevel(p);

  // remove from unique table (only applicable to reduced nodes)
  // if (isReducedNode(p)) {
    unsigned h = hashNode(p);
#ifdef DEVELOPMENT_CODE
    // node_finder key(this, p);
    node_reader key;
    initNodeReader(key, p, false);
    key.setHash(hashNode(p));
    if (unique->find(key) != p) {
      fprintf(stderr, "Error in deleteNode\nFind: %d\np: %d\n",
        unique->find(key), p);
      dumpInternal(stdout);
      MEDDLY_DCASSERT(false);
    }
    int x = unique->remove(h, p);
    MEDDLY_DCASSERT(p == x);
#else
    unique->remove(h, p);
#endif

#ifdef TRACK_DELETIONS
    printf("%s: p = %d, unique->remove(p) = %d\n", __func__, p, x);
    fflush(stdout);
#endif

    MEDDLY_DCASSERT(address[p].cache_count == 0);
  // }
  /*
  else {
    // Temporary node
    // TODO:
    // clear cache of corresponding temporary node?
    decrTempNodeCount(k);
  }
  */

  // unlink children
  const int nDptrs = ABS(foo[2]);
  int* downptr = foo + 3 + (foo[2] < 0? nDptrs: 0);
  int* stop = downptr + nDptrs;
#ifdef VALIDATE_INCOUNTS
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

  // recycle the index
  freeActiveNode(p);

  if (levels[k].compactLevel) levels[k].compact(address);

#ifdef VALIDATE_INCOUNTS_ON_DELETE
  validateIncounts(false);
#endif

}

void MEDDLY::expert_forest::zombifyNode(int p)
{
  MEDDLY_DCASSERT(isActiveNode(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  // MEDDLY_DCASSERT(isReducedNode(p));
  MEDDLY_DCASSERT(getCacheCount(p) > 0);  // otherwise this node should be deleted
  MEDDLY_DCASSERT(readInCount(p) == 0);
  MEDDLY_DCASSERT(address[p].cache_count > 0);

  stats.zombie_nodes++;
  levels[getNodeLevel(p)].zombie_nodes++;
  stats.decActive(1);

  // mark node as zombie
  address[p].cache_count = -address[p].cache_count;

  unsigned h = hashNode(p);
#ifdef DEVELOPMENT_CODE 
  // node_finder key(this, p);
  node_reader key;
  initNodeReader(key, p, false);
  key.setHash(hashNode(p));
  if (unique->find(key) != p) {
    fprintf(stderr, "Fail: can't find reduced node %d; got %d\n", p, unique->find(key));
    dumpInternal(stderr);
    MEDDLY_DCASSERT(false);
  }
  int x = unique->remove(h, p);
  MEDDLY_DCASSERT(x==p);
#else
  unique->remove(h, p);
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
#ifdef VALIDATE_INCOUNTS
      int temp = *downptr;
      *downptr = 0;
      unlinkNode(temp);
#else
      unlinkNode(*downptr);
#endif
    }
  } else {
    // Full encoding
    int* downptr = foo + 3;
    int* stop = downptr + foo[2];
    for (; downptr < stop; ++downptr) {
#ifdef VALIDATE_INCOUNTS
      int temp = *downptr;
      *downptr = 0;
      unlinkNode(temp);
#else
      unlinkNode(*downptr);
#endif
    }
  }
  // Recycle node memory
  levels[node_level].recycleNode(node_offset);
}

int MEDDLY::expert_forest
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
    if (nb.hasEdges()) {
      nb.copyIntoSparse(down, index, ld.sparseEdgeOf(address[p].offset), nnz);
    } else {
      nb.copyIntoSparse(down, index, nnz);
    }
  } else {
    // full node wins
    address[p].offset = ld.allocNode(truncsize, p, false);
    MEDDLY_DCASSERT(1==ld.countOf(address[p].offset));
    int* down = ld.fullDownOf(address[p].offset);
    if (nb.hasEdges()) {
      nb.copyIntoFull(down, ld.fullEdgeOf(address[p].offset), truncsize);
    } else {
      nb.copyIntoFull(down, truncsize);
    }
  } // if 

  // add to UT 
  unique->add(nb.hash(), p);
  
#ifdef DEVELOPMENT_CODE
  // node_finder key(this, p);
  node_reader key;
  initNodeReader(key, p, false);
  key.setHash(hashNode(p));
  MEDDLY_DCASSERT(key.hash() == nb.hash());
  MEDDLY_DCASSERT(unique->find(key) == p);
#endif
#ifdef DEBUG_CREATE_REDUCED
  printf("Created node %d\n", p);
  dump(stdout);
#endif
  return p;
}

void MEDDLY::expert_forest::validateDownPointers(const node_builder &nb) const
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

//
// ------------------------------------------------------------------
//  Still need to be organized:
// ------------------------------------------------------------------
//

bool MEDDLY::expert_forest::isStale(int h) const {
  return
    isMarkedForDeletion() || (
      isTerminalNode(h)
      ? terminalNodesAreStale
      : isPessimistic()
        ? isZombieNode(h)
        : (readInCount(h) == 0)
    );
}

/*
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
*/


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
  address = (nodeData*) realloc(address, (a_size+delta) * sizeof(nodeData));
  if (0==address) {
    throw error(error::INSUFFICIENT_MEMORY);
  }
  stats.incMemAlloc(delta * sizeof(nodeData));
  memset(address + a_size, 0, delta * sizeof(nodeData));
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
  address = (nodeData*) realloc(address, new_size * sizeof(nodeData));
  if (0==address) {
    throw error(error::INSUFFICIENT_MEMORY);
  }
  stats.decMemAlloc(delta * sizeof(nodeData));
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

#ifdef USE_OLD_TEMPNODES

int MEDDLY::expert_forest::createTempNode(int lh, int size, bool clear)
{
  throw error(error::NOT_IMPLEMENTED);
}

#endif


#ifdef USE_OLD_EVMDDS
bool MEDDLY::expert_forest
::getDownPtrsAndEdgeValues(int, std::vector<int>&, std::vector<int>&) const
{
  throw error(error::TYPE_MISMATCH);
}

bool MEDDLY::expert_forest
::getDownPtrsAndEdgeValues(int, std::vector<int>&, std::vector<float>&) const
{
  throw error(error::TYPE_MISMATCH);
}

int MEDDLY::expert_forest
::createTempNode(int lh, std::vector<int>&, std::vector<int>&)
{
  throw error(error::TYPE_MISMATCH);
}

int MEDDLY::expert_forest
::createTempNode(int lh, std::vector<int>&, std::vector<float>&)
{
  throw error(error::TYPE_MISMATCH);
}
#endif


