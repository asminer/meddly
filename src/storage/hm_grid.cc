
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

#include "hm_grid.h"


#define MERGE_AND_SPLIT_HOLES
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE


// ******************************************************************
// *                                                                *
// *                                                                *
// *                        hm_grid  methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::hm_grid::hm_grid(node_storage* p) : holeman(5, p)
{
  data = 0;
  size = 0;
  max_request = 0;
  large_holes = 0;
  holes_top = 0;
  holes_bottom = 0;
  updateData(data);
}

MEDDLY::hm_grid::~hm_grid()
{
  decMemAlloc(size*sizeof(node_handle));
  free(data);
}

MEDDLY::node_address MEDDLY::hm_grid::requestChunk(int slots)
{
  if (slots > max_request) {
    max_request = slots;
    // Traverse the large hole list, and re-insert
    // all the holes, in case some are now smaller
    // than max_request.
    node_handle curr = large_holes;
    large_holes = 0;
    for (; curr; ) {
      node_handle next = Next(curr);
      gridInsert(curr);
      curr = next;
    }
  }

  // First, try for a hole exactly of this size
  // by traversing the index nodes in the hole grid
  node_handle chain = 0;
  node_handle curr = holes_bottom;
  while (curr) {
    if (slots == -(data[curr])) break;
    if (slots < -(data[curr])) {
      // no exact match possible
      curr = 0;
      break;
    }
    // move up the hole grid
    curr = Up(curr);
    chain++;
  }

  if (curr) {
    // perfect fit
    hole_slots -= slots;
    // try to not remove the "index" node
    node_handle next = Next(curr);
    if (next) {
      midRemove(next);
#ifdef MEMORY_TRACE
      printf("Removed Non-Index hole %d\n", next);
#ifdef DEEP_MEMORY_TRACE
      dumpInternal(stdout);
#else 
      dumpInternal(stdout, next);
#endif
#endif
      return next;
    }
    indexRemove(curr);
#ifdef MEMORY_TRACE
    printf("Removed Index hole %d\n", curr);
#ifdef DEEP_MEMORY_TRACE
    dumpInternal(stdout);
#else 
    dumpInternal(stdout, curr);
#endif
#endif
    return curr;
  }

#ifdef MERGE_AND_SPLIT_HOLES
  // No hole with exact size.
  // Take the first large hole.
  if (large_holes) {
    curr = large_holes;
    large_holes = Next(curr);
    if (large_holes) Prev(large_holes) = 0;

    // Sanity check: the hole is large enough
    MEDDLY_DCASSERT(slots < -data[curr]);
    
    // Is there space for a leftover hole?
    const int min_node_size = smallestChunk();
    if (slots + min_node_size >= -data[curr]) {
      // leftover space is too small to be useful,
      // just send the whole thing.
#ifdef MEMORY_TRACE
      printf("Removed entire large hole %d\n", curr);
#ifdef DEEP_MEMORY_TRACE
      dumpInternal(stdout);
#else 
      dumpInternal(stdout, curr);
#endif
#endif
      hole_slots += data[curr]; // SUBTRACTS the number of slots
      return curr;
    }

    // This is a large hole, save the leftovers
    node_handle newhole = curr + slots;
    node_handle newsize = -(data[curr]) - slots;
    data[newhole] = -newsize;
    data[newhole + newsize - 1] = -newsize;
    gridInsert(newhole);

#ifdef MEMORY_TRACE
    printf("Removed part of hole %d\n", curr);
#ifdef DEEP_MEMORY_TRACE
    dumpInternal(stdout);
#else 
    dumpInternal(stdout, curr);
#endif
#endif
    data[curr] = -slots;
    hole_slots -= slots;
    return curr;
  }
#endif

  // 
  // Still here?  We couldn't recycle a node.
  // 

  //
  // First -- try to compact if we need to expand
  //
  if (getForest()->getPolicies().compactBeforeExpand) {
    if (last_slot + slots >= size) getParent()->collectGarbage(false);
  }

  //
  // Do we need to expand?
  //
  if (last_slot + slots >= size) {
    // new size is 50% more than previous 
    node_handle want_size = last_slot + slots;
    node_handle new_size = MAX(size, want_size) * 1.5;  // TBD: WTF?

    resize(new_size);
  }

  // 
  // Grab node from the end
  //
  node_handle h = last_slot + 1;
  last_slot += slots;
  data[h] = -slots;
  return h;
}

// ******************************************************************

void MEDDLY::hm_grid::recycleChunk(node_address addr, int slots)
{
#ifdef MEMORY_TRACE
  printf("Calling recycleChunk(%d, %d)\n", addr, slots);
#endif

  decMemUsed(slots * sizeof(node_handle));

  hole_slots += slots;
  data[addr] = data[addr+slots-1] = -slots;

  if (!getForest()->getPolicies().recycleNodeStorageHoles) return;

  // Check for a hole to the left
#ifdef MERGE_AND_SPLIT_HOLES
  if (data[addr-1] < 0) {
    // Merge!
#ifdef MEMORY_TRACE
    printf("Left merging\n");
#endif
    node_handle lefthole = addr + data[addr-1];
    MEDDLY_DCASSERT(data[lefthole] == data[addr-1]);
    if (non_index_hole == Up(lefthole)) midRemove(lefthole);
    else indexRemove(lefthole);
    slots += (-data[lefthole]);
    addr = lefthole;
    data[addr] = data[addr+slots-1] = -slots;
  }
#endif

  // if addr is the last hole, absorb into free part of array
  MEDDLY_DCASSERT(addr + slots - 1 <= last_slot);
  if (addr+slots-1 == last_slot) {
    last_slot -= slots;
    hole_slots -= slots;
    if (size > min_size && (last_slot + 1) < size/2) {
      node_handle new_size = size/2;
      while (new_size > (last_slot + 1) * 2) new_size /= 2;
      if (new_size < min_size) new_size = min_size;
      resize(new_size);
    }
#ifdef MEMORY_TRACE
    printf("Made Last Hole %d, last %d\n", addr, last);
#ifdef DEEP_MEMORY_TRACE
    dumpInternal(stdout);
#else 
    dumpInternal(stdout, 0);
#endif
#endif
    return;
  }

#ifdef MERGE_AND_SPLIT_HOLES
  // Check for a hole to the right
  if (data[addr+slots]<0) {
    // Merge!
#ifdef MEMORY_TRACE
    printf("Right merging\n");
#endif
    node_handle righthole = addr+slots;
    if (non_index_hole == Up(righthole)) midRemove(righthole);
    else indexRemove(righthole);
    slots += (-data[righthole]);
    data[addr] = data[addr+slots-1] = -slots;
  }
#endif

  // Add hole to grid
  gridInsert(addr); 

#ifdef MEMORY_TRACE
  printf("Made Hole %d\n", addr);
#ifdef DEEP_MEMORY_TRACE
  dumpInternal(stdout);
#else
  dumpInternal(stdout, addr);
#endif
#endif
}

// ******************************************************************

MEDDLY::node_address MEDDLY::hm_grid::chunkAfterHole(node_address addr) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(data[addr]<0);
  MEDDLY_DCASSERT(data[addr-data[addr]-1] == data[addr]);
  return addr - data[addr];
}

// ******************************************************************

void MEDDLY::hm_grid::dumpInternalInfo(FILE* s) const
{
  fprintf(s, "Last slot used: %ld\n", long(last_slot));
  fprintf(s, "large_holes: %ld\n", long(large_holes));
  fprintf(s, "Grid: top = %ld bottom = %ld\n", long(holes_top), long(holes_bottom));
}

// ******************************************************************

void MEDDLY::hm_grid::dumpHole(FILE* s, node_address a) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_CHECK_RANGE(1, a, last_slot);
  fprintf(s, "[%ld, ", long(data[a]));
  if (data[a+1]<0)  {
    fprintf(s, "non-index, p: %ld, ", long(data[a+2]));
  } else {
    fprintf(s, "u: %ld, d: %ld, ", long(data[a+1]), long(data[a+2])); 
  }
  long aN = chunkAfterHole(a)-1;
  fprintf(s, "n: %ld, ..., %ld]\n", long(data[a+3]), long(data[aN]));
}

// ******************************************************************

void MEDDLY::hm_grid
::reportMemoryUsage(FILE* s, const char* pad, int verb) const
{
  if (verb>7) {
    long holemem = hole_slots * sizeof(node_handle);
    fprintf(s, "%s  Hole Memory Usage:\t%ld\n", pad, holemem);
  }

  // Compute chain length histogram
  std::map<int, int> chainLengths;
  for (int curr = holes_bottom; curr; curr = Up(curr)) {
    int currHoleOffset = curr;
    int count = 0;
    // traverse this chain to get its length
    for (count = 0; currHoleOffset; count++) {
      currHoleOffset = Next(currHoleOffset);
    }
    int currHoleSize = -data[curr];
    chainLengths[currHoleSize] += count;
  }
  // add the large hole list
  int count = 0;
  for (int curr = large_holes; curr; curr = Next(curr)) {
    count++;
  }
  chainLengths[-1] += count;

  // Display the histogram

  fprintf(s, "%s  Hole Chains (size, count):\n", pad);
  for (std::map<int, int>::iterator iter = chainLengths.begin();
      iter != chainLengths.end(); ++iter)
    {
      if (iter->first<0)
        fprintf(s, "%s\tlarge: %d\n", pad, iter->second);
      else
        fprintf(s, "%s\t%5d: %d\n", pad, iter->first, iter->second);
    }
  fprintf(s, "%s  End of Hole Chains\n", pad);
}

// ******************************************************************

void MEDDLY::hm_grid::clearHolesAndShrink(node_address new_last, bool shrink)
{
  last_slot = new_last;

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  hole_slots = 0;
  large_holes = 0;

  if (shrink && size > min_size && last_slot < size/2) {
    node_handle new_size = size/2;
    while (new_size > min_size && new_size > last_slot * 3) { new_size /= 2; }
    resize(new_size);
  }
}

// ******************************************************************
// *                                                                *
// *                        private  helpers                        *
// *                                                                *
// ******************************************************************

void MEDDLY::hm_grid::gridInsert(node_handle p_offset)
{
#ifdef MEMORY_TRACE
  printf("gridInsert(%d)\n", p_offset);
#endif

  // sanity check to make sure that the first and last slots in this hole
  // have the same value, i.e. -(# of slots in the hole)
  MEDDLY_DCASSERT(data[p_offset] == data[p_offset - data[p_offset] - 1]);

  // Check if we belong in the grid, or the large hole list
  if (-data[p_offset] > max_request) {
#ifdef MEMORY_TRACE
    printf("\tAdding to large_holes: %d\n", large_holes);
#endif
    // add to the large hole list
    Up(p_offset) = non_index_hole;
    Next(p_offset) = large_holes;
    Prev(p_offset) = 0;
    if (large_holes) Prev(large_holes) = p_offset;
    large_holes = p_offset;
    return;
  }

  // special case: empty
  if (0 == holes_bottom) {
#ifdef MEMORY_TRACE
    printf("\tAdding to empty grid\n");
#endif
    // index hole
    Up(p_offset) = 0;
    Down(p_offset) = 0;
    Next(p_offset) = 0;
    holes_top = holes_bottom = p_offset;
    return;
  }
  // special case: at top
  if (data[p_offset] < data[holes_top]) {
#ifdef MEMORY_TRACE
    printf("\tAdding new chain at top\n");
#endif
    // index hole
    Up(p_offset) = 0;
    Next(p_offset) = 0;
    Down(p_offset) = holes_top;
    Up(holes_top) = p_offset;
    holes_top = p_offset;
    return;
  }
  node_handle above = holes_bottom;
  node_handle below = 0;
  while (data[p_offset] < data[above]) {
    below = above;
    above = Up(below);
    MEDDLY_DCASSERT(Down(above) == below);
    MEDDLY_DCASSERT(above);  
  }
  if (data[p_offset] == data[above]) {
#ifdef MEMORY_TRACE
    printf("\tAdding to chain\n");
#endif
    // Found, add this to chain
    // making a non-index hole
    node_handle right = Next(above);
    Up(p_offset) = non_index_hole;
    Prev(p_offset) = above;
    Next(p_offset) = right;
    if (right) Prev(right) = p_offset;
    Next(above) = p_offset;
    return; 
  }
#ifdef MEMORY_TRACE
  printf("\tAdding new chain\n");
#endif
  // we should have above < p_offset < below  (remember, -sizes)
  // create an index hole since there were no holes of this size
  Up(p_offset) = above;
  Down(p_offset) = below;
  Next(p_offset) = 0;
  Down(above) = p_offset;
  if (below) {
    Up(below) = p_offset;
  } else {
    MEDDLY_DCASSERT(above == holes_bottom);
    holes_bottom = p_offset;
  }
}

// ******************************************************************

void MEDDLY::hm_grid::midRemove(node_handle p_offset)
{
#ifdef MEMORY_TRACE
  printf("midRemove(%d)\n", p_offset);
#endif

  // Remove a "middle" node, either in the grid
  // or in the large hole list.
  //
  MEDDLY_DCASSERT(isHoleNonIndex(p_offset));

  node_handle left = Prev(p_offset); 
  node_handle right = Next(p_offset);

#ifdef MEMORY_TRACE
  printf("\tIN left: %d  right: %d  large_holes: %d\n", 
    left, right, large_holes
  );
#endif

  if (left) {
    Next(left) = right;
  } else {
    // MUST be head of the large hole list
    MEDDLY_DCASSERT(large_holes == p_offset);
    large_holes = right;
  }
  if (right) Prev(right) = left;

#ifdef MEMORY_TRACE
  printf("\tOUT large_holes: %d\n", large_holes);
#endif

  // Sanity checks
#ifdef DEVELOPMENT_CODE
  if (large_holes) {
    MEDDLY_CHECK_RANGE(1, large_holes, last_slot+1);
    MEDDLY_DCASSERT(data[large_holes] < 0);
  }
#endif
}

// ******************************************************************

void MEDDLY::hm_grid::indexRemove(node_handle p_offset)
{
#ifdef MEMORY_TRACE
  printf("indexRemove(%d)\n", p_offset);
#endif

  MEDDLY_DCASSERT(!isHoleNonIndex(p_offset));
  node_handle above = Up(p_offset); 
  node_handle below = Down(p_offset); 
  node_handle right = Next(p_offset); 

  if (right >= 1) {
    // there are nodes to the right!
    MEDDLY_DCASSERT(Up(right) < 0);
    Up(right) = above;
    Down(right) = below;
    // update the pointers of the holes (index) above and below it
    if (above) {
      Down(above) = right;
    } else {
      holes_top = right;
    }

    if (below) {
      Up(below) = right;
    } else {
      holes_bottom = right;
    }
    
  } else {
    // there are no non-index nodes
    MEDDLY_DCASSERT(right < 1);

    // this was the last node of its size
    // update the pointers of the holes (index) above and below it
    if (above) {
      Down(above) = below;
    } else {
      holes_top = below;
    }

    if (below) {
      Up(below) = above;
    } else {
      holes_bottom = above;
    }
  }
  // Sanity checks
#ifdef DEVELOPMENT_CODE
  if (large_holes) {
    MEDDLY_CHECK_RANGE(1, large_holes, last_slot+1);
    MEDDLY_DCASSERT(data[large_holes] < 0);
  }
#endif
}

// ******************************************************************

void MEDDLY::hm_grid::resize(node_handle new_slots)
{
  node_handle* new_data 
    = (node_handle*) realloc(data, new_slots * sizeof(node_handle));

  if (0 == new_data) throw error(error::INSUFFICIENT_MEMORY);
  if (0 == data) new_data[0] = 0;
  if (new_slots > size) {
    incMemAlloc((new_slots - size) * sizeof(node_handle));
  } else {
    decMemAlloc((size - new_slots) * sizeof(node_handle));
  }
#ifdef MEMORY_TRACE
    printf("Resized data[]. Old size: %d, New size: %d, Last: %d.\n", 
      size, new_slots, last
    );
#endif
  size = new_slots;
  if (data != new_data) {
    // update pointers
    data = new_data;
    updateData(data);
  }
}


