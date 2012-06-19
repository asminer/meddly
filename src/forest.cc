
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
#include <set>
#include <queue>
#include <vector>

// #define DEBUG_CLEANUP
// #define MERGE_RIGHT
// #define MERGE_LEFT
// #define ENABLE_BREAKING_UP_HOLES

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
// *               expert_forest::level_data  methods               *
// *                                                                *
// ******************************************************************

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

void MEDDLY::expert_forest::level_data::clear(expert_forest* p)
{
  parent = p;
  data = 0;
  size = 0;
  last = 0;
  compactLevel = false;
  holes_top = 0;
  holes_bottom = 0;
  hole_slots = 0;
  max_hole_chain = 0;
  zombie_nodes = 0;
  temp_nodes = 0;
  num_compactions = 0;
  levelNode = 0;
}

int MEDDLY::expert_forest::level_data::getHole(int slots, bool search_holes)
{
  MEDDLY_DCASSERT(parent);
  parent->stats.incMemUsed(slots * sizeof(int));

  if (search_holes && parent->areHolesRecycled()) {
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
    const int min_node_size = 5;  // TBD!!!
    /*
      isEVPlus()
        ? 7
        : isEVTimes()
          ? 6
          : 5;
          */

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
  }

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
  cout << "Calling makeHole(" addr << ", " << slots << ")\n";
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
    cout << "Made Last Hole " << addr << "\n";
    dumpInternal(stdout);
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
  cout << "Made Last Hole " << addr << "\n";
  dumpInternal(stdout);
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
  cout << __func__ << "(" << k << ", " << p_offset << ")\n";
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

void MEDDLY::expert_forest::level_data::compact(mdd_node_data* address)
{
  if (0==data) {
    compactLevel = false;
    return;
  }
  if (0 == hole_slots ||  // Already compact
      !needsCompaction()) {  // Level is compact enough!
    compactLevel = false;
#if 0
    printf("%s: level %d... compact enough\n", __func__, k);
#endif
    return;
  }

  if (0 < temp_nodes) return;   // Temp nodes; do not compact
#if 0
  printf("%s: level %d\n", __func__, k);
#endif

#if 0
  printf("Before compaction:\n");
  dumpInternalLevel(stdout, k);
  printf("\n");
#endif

#ifdef DEBUG_SLOW
  fprintf(stderr, "Compacting forest level %d\n", k);
#endif

  // alternate algorithm -- since we now have the node ids in the node data
  int *node_ptr = data + 1;  // since we leave [0] empty
  int *end_ptr = data + last + 1;
  int *curr_ptr = node_ptr;
  int node_size = 0;
  int curr_node = 0;

  int sparseMultiplier = parent->isMultiTerminal() ? -2 : -3;
  int fullMultiplier = parent->isMultiTerminal() ? 1 : 2;

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

      node_size = *(node_ptr + 2);  // [2] = size
      MEDDLY_DCASSERT (node_size != 0);      // assuming zombies have been deleted

      node_size = parent->getDataHeaderSize() +
        (node_size * (node_size < 0? sparseMultiplier: fullMultiplier));

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
  a_size = 1024;
  address = (mdd_node_data *) malloc(a_size * sizeof(mdd_node_data));
  if (0 == address) throw error(error::INSUFFICIENT_MEMORY);
  stats.incMemAlloc(a_size * sizeof(mdd_node_data));
  memset(address, 0, a_size * sizeof(mdd_node_data));
  a_last = peak_nodes = a_unused = 0;
  

  //
  // Initialize level array
  //
  int N = getNumVariables();
  if (rel) {
    raw_levels = new level_data[2*N+1];
    levels = raw_levels + N;
    stats.incMemAlloc(2*N+1 * sizeof(level_data));
    for (int i=2*N; i>=0; i--) raw_levels[i].clear(this);
  } else {
    raw_levels = new level_data[N+1];
    levels = raw_levels;
    stats.incMemAlloc(N+1 * sizeof(level_data));
    for (int i=N; i>=0; i--) raw_levels[i].clear(this);
  }
}


MEDDLY::expert_forest::~expert_forest() 
{
  // Address array
  free(address);

  // Level array
  delete[] raw_levels;
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


