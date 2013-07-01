
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

#include "hm_array.h"


#define MERGE_AND_SPLIT_HOLES
// #define DEBUG_COMPACTION
// #define DEBUG_SLOW
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE


// ******************************************************************
// *                                                                *
// *                                                                *
// *                        hm_array methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_handle MEDDLY::hm_array::verify_hole_slots;

MEDDLY::hm_array::hm_array(node_storage* p) : holeman(4, p)
{
  data = 0;
  size = 0;
  large_holes = 0;
  for (int i=0; i<LARGE_SIZE; i++) {
    small_holes[i] = 0;
  }
  updateData(data);
#ifdef MEASURE_LARGE_HOLE_STATS
  num_large_hole_traversals = 0;
  count_large_hole_visits = 0;
#endif
}

// ******************************************************************

MEDDLY::hm_array::~hm_array()
{
  decMemAlloc(size*sizeof(node_handle));
  free(data);
}

// ******************************************************************

MEDDLY::node_address MEDDLY::hm_array::requestChunk(int slots)
{
#ifdef MEMORY_TRACE
  printf("Requesting %d slots\n", slots);
#endif
  node_handle found = 0;
  // Try for an exact fit first
  if (slots < LARGE_SIZE) {
    if (small_holes[slots]) {
      found = small_holes[slots];
      listRemove(small_holes[slots], found);
    }
  }
  // Try the large hole list and check for a hole large enough
  if (!found) {
#ifdef MEASURE_LARGE_HOLE_STATS
    num_large_hole_traversals++;
#endif
    node_handle curr = large_holes;
    while (curr) {
#ifdef MEASURE_LARGE_HOLE_STATS
      count_large_hole_visits;
#endif
      if (-data[curr] >= slots) {
        // we have a large enough hole, grab it
        found = curr;
        // remove the hole from the list
        listRemove(large_holes, found);
        break;
      }
    }
  }

  if (found) {
    // We found a hole to recycle.
#ifdef MEMORY_TRACE
    printf("Removed hole %d\n", found);
#ifdef DEEP_MEMORY_TRACE
    getParent()->dumpInternal(stdout);
#else 
    dumpHole(stdout, found);
#endif
#endif
    
    // Is there space for a leftover hole?
    const int min_node_size = smallestChunk();
    if (slots + min_node_size >= -data[found]) {
      // leftover space is too small to be useful (or 0),
      // just send the whole thing.

      // Update stats
      hole_slots += data[found]; // SUBTRACTS the number of slots
    } else {
      // There is space for the leftover slots to make a hole
      // This is a large hole, save the leftovers
      node_handle newhole = found + slots;
      node_handle newsize = -(data[found]) - slots;
      data[newhole] = -newsize;
      data[newhole + newsize - 1] = -newsize;
      if (newsize < LARGE_SIZE) {
        listInsert(small_holes[newsize], newhole);
      } else {
        listInsert(large_holes, newhole);
      }
#ifdef MEMORY_TRACE
      printf("Made leftover hole %ld\n", long(newhole));
#ifdef DEEP_MEMORY_TRACE
      getParent()->dumpInternal(stdout);
#else
      dumpHole(stdout, newhole);
#endif
#endif
      // set size of allocated memory to exactly what was requested
      data[found] = -slots;

      // Update stats
      hole_slots -= slots;
    }
    return found;
  }
  
  // 
  // Still here?  We couldn't recycle a node.
  // 

#ifdef MEMORY_TRACE
  printf("No exact or large hole available\n");
#endif

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

void MEDDLY::hm_array::recycleChunk(node_address addr, int slots)
{
#ifdef MEMORY_TRACE
  printf("Calling recycleChunk(%ld, %d)\n", long(addr), slots);
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
    // find the left hole address
    node_handle lefthole = addr + data[addr-1];
    MEDDLY_DCASSERT(data[lefthole] == data[addr-1]);

    // remove the left hole
    if (-data[lefthole] < LARGE_SIZE) {
      listRemove(small_holes[-data[lefthole]], lefthole);
    } else {
      listRemove(large_holes, lefthole);
    }

    // merge with us
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
    printf("Made Last Hole %ld, last %ld\n", long(addr), long(last_slot));
#ifdef DEEP_MEMORY_TRACE
    getParent()->dumpInternal(stdout);
#endif
#endif
    return;
  }

  // Still here? Wasn't the last hole.

#ifdef MERGE_AND_SPLIT_HOLES
  // Check for a hole to the right
  if (data[addr+slots]<0) {
    // Merge!
#ifdef MEMORY_TRACE
    printf("Right merging\n");
#endif
    // find the right hole address
    node_handle righthole = addr+slots;

    // remove the right hole
    if (-data[righthole] < LARGE_SIZE) {
      listRemove(small_holes[-data[righthole]], righthole);
    } else {
      listRemove(large_holes, righthole);
    }
    
    // merge with us
    slots += (-data[righthole]);
    data[addr] = data[addr+slots-1] = -slots;
  }
#endif

  // Add hole to the proper list
  if (-data[addr] < LARGE_SIZE) {
    listInsert(small_holes[-data[addr]], addr);
  } else {
    listInsert(large_holes, addr);
  }

#ifdef MEMORY_TRACE
  printf("Made Hole %ld\n", long(addr));
#ifdef DEEP_MEMORY_TRACE
  getParent()->dumpInternal(stdout);
#else
  dumpHole(stdout, addr);
#endif
#endif
}

// ******************************************************************

MEDDLY::node_address MEDDLY::hm_array::chunkAfterHole(node_address addr) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(data[addr]<0);
  MEDDLY_DCASSERT(data[addr-data[addr]-1] == data[addr]);
  return addr - data[addr];
}

// ******************************************************************

void MEDDLY::hm_array::dumpInternalInfo(FILE* s) const
{
  fprintf(s, "Last slot used: %ld\n", long(last_slot));
  fprintf(s, "Total hole slots: %ld\n", long(hole_slots));
  fprintf(s, "small_holes: (");
  bool printed = false;
  for (int i=0; i<LARGE_SIZE; i++) if (small_holes[i]) {
    if (printed) fprintf(s, ", ");
    fprintf(s, "%d:%d", i, small_holes[i]);
    printed = true;
  }
  fprintf(s, ")\n");
  fprintf(s, "large_holes: %ld\n", long(large_holes));
  verify_hole_slots = 0;
}

// ******************************************************************

void MEDDLY::hm_array::dumpHole(FILE* s, node_address a) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_CHECK_RANGE(1, a, last_slot);
  long aN = chunkAfterHole(a)-1;
  fprintf(s, "[%ld, p: %ld, n: %ld, ..., %ld]\n", 
      long(data[a]), long(data[a+1]), long(data[a+2]), long(data[aN])
  );
  verify_hole_slots += aN+1-a;
}

// ******************************************************************

void MEDDLY::hm_array::dumpInternalTail(FILE* s) const
{
  if (verify_hole_slots != hole_slots) {
    fprintf(s, "Counted hole slots: %ld\n", long(verify_hole_slots));
    MEDDLY_DCASSERT(false);
  }
}

// ******************************************************************

void MEDDLY::hm_array
::reportMemoryUsage(FILE* s, const char* pad, int verb) const
{
  if (verb>7) {
    long holemem = hole_slots * sizeof(node_handle);
    fprintf(s, "%s  Hole Memory Usage:\t%ld\n", pad, holemem);
  }

  fprintf(s, "%sLength of non-empty chains:\n", pad);
  for (int i=0; i<LARGE_SIZE; i++) {
    long L = listLength(small_holes[i]);
    if (L) {
      fprintf(s, "%s%9d: %ld\n", pad, i, L);
    }
  }
  long LL = listLength(large_holes);
  if (LL) fprintf(s, "%s    large: %ld\n", pad, listLength(large_holes));

#ifdef MEASURE_LARGE_HOLE_STATS
  fprintf(s, "%s#traversals large_holes: %ld\n", pad, num_large_hole_traversals);
  if (num_large_hole_traversals) {
    double avg = count_large_hole_visits;
    avg /= num_large_hole_traversals;
    fprintf(s, "%sAvg cost per traversal : %lf\n", pad, avg);
  }
#endif
}

// ******************************************************************

void MEDDLY::hm_array::clearHolesAndShrink(node_address new_last, bool shrink)
{
  last_slot = new_last;

  // set up hole pointers and such
  hole_slots = 0;
  large_holes = 0;
  for (int i=0; i<LARGE_SIZE; i++) small_holes[i] = 0;

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

void MEDDLY::hm_array::resize(node_handle new_slots)
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
    printf("Resized data[]. Old size: %ld, New size: %ld, Last: %ld.\n", 
      long(size), long(new_slots), long(last_slot)
    );
#endif
  size = new_slots;
  if (data != new_data) {
    // update pointers
    data = new_data;
    updateData(data);
  }
}


