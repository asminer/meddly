
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

MEDDLY::node_handle MEDDLY::hm_grid::verify_hole_slots;

MEDDLY::hm_grid::hm_grid(node_storage* p) : holeman(5, p)
{
  max_request = 0;
  large_holes = 0;
  holes_top = 0;
  holes_bottom = 0;
}

MEDDLY::hm_grid::~hm_grid()
{
}

MEDDLY::node_address MEDDLY::hm_grid::requestChunk(int slots)
{
#ifdef MEMORY_TRACE
  printf("Requesting %d slots\n", slots);
#endif

  if (slots > max_request) {
#ifdef MEMORY_TRACE
    printf("New max request; shifting holes from large list to grid\n");
#endif
    max_request = slots;
    // Traverse the large hole list, and re-insert
    // all the holes, in case some are now smaller
    // than max_request.
    node_handle curr = large_holes;
    large_holes = 0;
    for (; curr; ) {
      node_handle next = Next(curr);
      insertHole(curr);
      curr = next;
    }
  }

  // find us a good hole
  node_handle found = 0;

  // First, try for a hole exactly of this size
  // by traversing the index nodes in the hole grid
  node_handle chain = 0;
  node_handle curr = holes_bottom;
  while (curr) {
    if (slots == -(data[curr])) {
      if (Next(curr)) {
        found = Next(curr);
      } else {
        found = curr;
      }
      break;
    }
    if (slots < -(data[curr])) {
      // no exact match possible
      curr = 0;
      break;
    }
    // move up the hole grid
    curr = Up(curr);
    chain++;
  }

#ifdef MERGE_AND_SPLIT_HOLES
  // If that failed, try the large hole list
  if (!found && large_holes) {
    // should be large enough
    found = large_holes;
  }

  // If that failed, try the largest hole in the grid
  if (!found && holes_top) if (slots < -data[holes_top]) {
    if (Next(holes_top)) {
      found = Next(holes_top);
    } else {
      found = holes_top;
    }
  }
#endif

  if (found) {
    // We have a hole to recycle
    // Sanity check:
    MEDDLY_DCASSERT(slots <= -data[found]);

    // Remove the hole
    removeHole(found);
    useHole(-data[found]);

    // hole might be larger than requested; decide what to do with leftovers.
    // Is there space for a leftover hole?
    const int min_node_size = smallestChunk();
    if (slots + min_node_size >= -data[found]) {
      // leftover space is too small to be useful,
      // just send the whole thing.
      incFragments(-data[found] - slots);
      return found;
    }

    // Save the leftovers - make a new hole!
    node_handle newhole = found + slots;
    node_handle newsize = -(data[found]) - slots;
    data[found] = -slots;
    data[newhole] = -newsize;
    data[newhole + newsize - 1] = -newsize;
    insertHole(newhole);
    newHole(newsize);
    return found;
  }
  
  // 
  // Still here?  We couldn't recycle a node.
  // 
  return allocFromEnd(slots);
}

// ******************************************************************

void MEDDLY::hm_grid::recycleChunk(node_address addr, int slots)
{
#ifdef MEMORY_TRACE
  printf("Calling recycleChunk(%d, %d)\n", addr, slots);
#endif

  decMemUsed(slots * sizeof(node_handle));

  newHole(slots);
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
    removeHole(lefthole);
    slots += (-data[lefthole]);
    addr = lefthole;
    data[addr] = data[addr+slots-1] = -slots;
    useHole(0);
  }
#endif

  // if addr is the last hole, absorb into free part of array
  MEDDLY_DCASSERT(addr + slots - 1 <= lastSlot());
  if (addr+slots-1 == lastSlot()) {
    releaseToEnd(addr, slots);
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
    removeHole(righthole);
    slots += (-data[righthole]);
    data[addr] = data[addr+slots-1] = -slots;
    useHole(0);
  }
#endif

  // Add hole to grid
  insertHole(addr);

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

void MEDDLY::hm_grid::dumpInternalInfo(FILE* s) const
{
  fprintf(s, "Last slot used: %ld\n", long(lastSlot()));
  fprintf(s, "Total hole slots: %ld\n", holeSlots());
  fprintf(s, "large_holes: %ld\n", long(large_holes));
  fprintf(s, "Grid: top = %ld bottom = %ld\n", long(holes_top), long(holes_bottom));
  verify_hole_slots = 0;
}

// ******************************************************************

void MEDDLY::hm_grid::dumpHole(FILE* s, node_address a) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_CHECK_RANGE(1, a, lastSlot());
  fprintf(s, "[%ld, ", long(data[a]));
  if (data[a+1]<0)  {
    fprintf(s, "non-index, p: %ld, ", long(data[a+2]));
  } else {
    fprintf(s, "u: %ld, d: %ld, ", long(data[a+1]), long(data[a+2])); 
  }
  long aN = chunkAfterHole(a)-1;
  fprintf(s, "n: %ld, ..., %ld]\n", long(data[a+3]), long(data[aN]));
  verify_hole_slots += aN+1-a;
}

// ******************************************************************

void MEDDLY::hm_grid::dumpInternalTail(FILE* s) const
{
  if (verify_hole_slots != holeSlots()) {
    fprintf(s, "Counted hole slots: %ld\n", long(verify_hole_slots));
    MEDDLY_DCASSERT(false);
  }
}

// ******************************************************************

void MEDDLY::hm_grid
::reportStats(FILE* s, const char* pad, unsigned flags) const
{
  static unsigned HOLE_MANAGER =
    expert_forest::HOLE_MANAGER_STATS | expert_forest::HOLE_MANAGER_DETAILED;

  if (! (flags & HOLE_MANAGER)) return;

  fprintf(s, "%sStats for grid hole management\n", pad);

  holeman::reportStats(s, pad, flags);

  if (! (flags & expert_forest::HOLE_MANAGER_DETAILED)) return;

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

  fprintf(s, "%s    Hole Chains (size, count):\n", pad);
  for (std::map<int, int>::iterator iter = chainLengths.begin();
      iter != chainLengths.end(); ++iter)
    {
      if (iter->first<0)
        fprintf(s, "%s\tlarge: %d\n", pad, iter->second);
      else
        fprintf(s, "%s\t%5d: %d\n", pad, iter->first, iter->second);
    }
  fprintf(s, "%s    End of Hole Chains\n", pad);
}

// ******************************************************************

void MEDDLY::hm_grid::clearHolesAndShrink(node_address new_last, bool shrink)
{
  holeman::clearHolesAndShrink(new_last, shrink);

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  large_holes = 0;
}

// ******************************************************************
// *                                                                *
// *                        private  helpers                        *
// *                                                                *
// ******************************************************************

void MEDDLY::hm_grid::insertHole(node_handle p_offset)
{
#ifdef MEMORY_TRACE
  printf("insertHole(%ld)\n", long(p_offset));
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
    makeMiddle(p_offset, 0, large_holes);
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
    makeIndex(p_offset, 0, 0, 0);
    holes_top = holes_bottom = p_offset;
    return;
  }

  // special case: at top
  if (data[p_offset] < data[holes_top]) {
#ifdef MEMORY_TRACE
    printf("\tAdding new chain at top\n");
#endif
    // index hole
    makeIndex(p_offset, 0, holes_top, 0);
    Up(holes_top) = p_offset;
    holes_top = p_offset;
    return;
  }

  // find our vertical position in the grid
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
    makeMiddle(p_offset, above, right);
    if (right) Prev(right) = p_offset;
    Next(above) = p_offset;
    return; 
  }
#ifdef MEMORY_TRACE
  printf("\tAdding new chain\n");
#endif
  // we should have above < p_offset < below  (remember, -sizes)
  // create an index hole since there were no holes of this size
  makeIndex(p_offset, above, below, 0);
  Down(above) = p_offset;
  if (below) {
    Up(below) = p_offset;
  } else {
    MEDDLY_DCASSERT(above == holes_bottom);
    holes_bottom = p_offset;
  }
}

// ******************************************************************

void MEDDLY::hm_grid::removeHole(node_handle p_offset)
{
  if (isIndexHole(p_offset)) {
      //
      // Index node
      //
#ifdef MEMORY_TRACE
      printf("indexRemove(%ld)\n", long(p_offset));
#endif
      node_handle above, below, right;
      getIndex(p_offset, above, below, right);

      if (right) {
        // convert right into an index node
        middle2index(right, above, below);

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
        // remove this row completely

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
#ifdef MEMORY_TRACE
      printf("Removed Index hole %ld\n", long(p_offset));
#ifdef DEEP_MEMORY_TRACE
      dumpInternal(stdout);
#else 
      dumpInternal(stdout, p_offset);
#endif
#endif
      //
      // Done with index node removal
      //
  } else {
      //
      // Not an index node
      //
#ifdef MEMORY_TRACE
      printf("midRemove(%ld)\n", long(p_offset));
#endif
      // Remove a "middle" node, either in the grid
      // or in the large hole list.
      //
      node_handle left, right;
      getMiddle(p_offset, left, right);

      if (left) {
        Next(left) = right;
      } else {
        // MUST be head of the large hole list
        MEDDLY_DCASSERT(large_holes == p_offset);
        large_holes = right;
      }
      if (right) Prev(right) = left;

#ifdef MEMORY_TRACE
      printf("Removed Non-Index hole %ld\n", long(p_offset));
#ifdef DEEP_MEMORY_TRACE
      dumpInternal(stdout);
#else 
      dumpInternal(stdout, p_offset);
#endif
#endif
  }
}


