// $Id: hm_grid.cc 436 2013-07-01 21:43:08Z asminer $

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



#include "hm_heap.h"

// #define VERIFY_HEAP
// #define MEMORY_TRACE
// #define DEEP_MEMORY_TRACE
// #define INTERNAL_CODE 0x02

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        hm_heap  methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::hm_heap::hm_heap(node_storage* p) : holeman(6, p)
{
  max_request = 0;
  large_holes = 0;
  large_holes_size = 0;
  holes_top = 0;
  holes_bottom = 0;
}

// ******************************************************************

MEDDLY::hm_heap::~hm_heap()
{
}

// ******************************************************************

MEDDLY::node_address
MEDDLY::hm_heap::requestChunk(int req_slots)
{
  int slots = MAX(smallestChunk(), req_slots);
#ifdef MEMORY_TRACE
  printf("Requesting %d slots, adjusted to %d\n", req_slots, slots);
#endif

  if (slots > max_request) {
#ifdef MEMORY_TRACE
    printf("New max request; shifting holes from large list to grid\n");
#endif
    max_request = slots;

    // Convert large hole heap into a list
    heapToList(large_holes);
    node_handle curr = large_holes;
    large_holes = 0;
    large_holes_size = 0;

    // Traverse the large hole list, and re-insert
    // all the holes, in case some are now smaller
    // than max_request.
    for (; curr; ) {
      node_handle next = Parent(curr);
      insertHole(curr);
      curr = next;
    }
  }

  // Find earliest hole that is large enough for us
  node_handle found = large_holes;
  node_handle curr = holes_top;
  while (curr) {
    if (slots > -data[curr]) {
      break;
    }
    if (slots == -data[curr]) {
      found = curr;
      break;
    }
    if (found) {
      if (curr < found) {
        if (slots + smallestChunk() < -data[curr]) {  // won't cause a fragment
          found = curr;
        }
      }
    } else {
      found = curr;
    }
    curr = Down(curr);
  } // while curr


  if (found) {
    // We have a hole to recycle
    // Sanity check:
    MEDDLY_DCASSERT(slots <= -data[found]);

    // Remove the hole
    removeHole(found);

    node_handle newhole = found + slots;
    node_handle newsize = -(data[found]) - slots;
    data[found] = -slots;
    if (newsize > 0) {
      // Save the leftovers - make a new hole!
      data[newhole] = -newsize;
      data[newhole + newsize - 1] = -newsize;
      if (insertHole(newhole)) {
        newHole(newsize);
      } else {
        newUntracked(newsize);
      }
    }
#ifdef MEMORY_TRACE
    printf("\trecycled %d slots, addr %ld\n", -data[found], long(found));
#endif
    return found;
  }
  
  // 
  // Still here?  We couldn't recycle a node.
  // 
  found = allocFromEnd(slots);
#ifdef MEMORY_TRACE
  printf("\tallocated %d slots, addr %ld\n", -data[found], long(found));
#endif
  MEDDLY_DCASSERT(-data[found] >= smallestChunk());
  return found;
}

// ******************************************************************

void MEDDLY::hm_heap::recycleChunk(node_address addr, int slots)
{
#ifdef MEMORY_TRACE
  printf("Calling recycleChunk(%ld, %d)\n\tchunk contains ", addr, slots);
  dumpInternalNode(stdout, addr, 0x03);
#endif

  decMemUsed(slots * sizeof(node_handle));

  data[addr] = data[addr+slots-1] = -slots;

  if (!getForest()->getPolicies().recycleNodeStorageHoles) return;

  // Check for a hole to the left
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
  }

  // if addr is the last hole, absorb into free part of array
  MEDDLY_DCASSERT(addr + slots - 1 <= lastSlot());
  if (addr+slots-1 == lastSlot()) {
    releaseToEnd(addr, slots);
    return;
  }

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
  }

  // Add hole to grid
  if (insertHole(addr)) {
    newHole(slots);
  } else {
    newUntracked(slots);
  }

#ifdef MEMORY_TRACE
  printf("Made Hole %ld\n", addr);
#ifdef DEEP_MEMORY_TRACE
  dumpInternal(stdout, INTERNAL_CODE);
#else
  dumpHole(stdout, addr);
#endif
#endif
}

// ******************************************************************

void MEDDLY::hm_heap::dumpInternalInfo(FILE* s) const
{
  fprintf(s, "Last slot used: %ld\n", long(lastSlot()));
  fprintf(s, "Total hole slots: %ld\n", holeSlots());
  fprintf(s, "large_holes: %ld\n", long(large_holes));
  fprintf(s, "large_holes_size: %ld\n", long(large_holes_size));
  fprintf(s, "Grid: top = %ld bottom = %ld\n", long(holes_top), long(holes_bottom));
}

// ******************************************************************

void MEDDLY::hm_heap::dumpHole(FILE* s, node_address a) const
{
  MEDDLY_DCASSERT(data);
  MEDDLY_CHECK_RANGE(1, a, lastSlot());
  fprintf(s, "[%ld, ", long(data[a]));
  if (isIndexHole(a)) {
    node_handle up, down, root, size;
    getIndex(a, up, down, root, size);
    fprintf(s, "u: %ld, d: %ld, r: %ld, s: %ld", 
      long(up), long(down), long(root), long(size)
    );
  } else {
    node_handle left, right, parent, id;
    getHeapNode(a, id, parent, left, right);
    fprintf(s, "l: %ld, r: %ld, p: %ld, id: %ld",
      long(left), long(right), long(parent), long(id)
    );
  }
  long aN = chunkAfterHole(a)-1;
  fprintf(s, ", ..., %ld]\n", long(data[aN]));
}

// ******************************************************************

void MEDDLY::hm_heap
::reportStats(FILE* s, const char* pad, unsigned flags) const
{
  static unsigned HOLE_MANAGER =
    expert_forest::HOLE_MANAGER_STATS | expert_forest::HOLE_MANAGER_DETAILED;

  if (! (flags & HOLE_MANAGER)) return;

  fprintf(s, "%sStats for grid hole management\n", pad);

  holeman::reportStats(s, pad, flags);

  if (! (flags & expert_forest::HOLE_MANAGER_DETAILED)) return;

  // Any detailed info?
}

// ******************************************************************

void MEDDLY::hm_heap::clearHolesAndShrink(node_address new_last, bool shrink)
{
  holeman::clearHolesAndShrink(new_last, shrink);

  // set up hole pointers and such
  holes_top = holes_bottom = 0;
  large_holes = 0;
  large_holes_size = 0;
}

// ******************************************************************
// *                                                                *
// *                        private  helpers                        *
// *                                                                *
// ******************************************************************

bool MEDDLY::hm_heap::insertHole(node_handle p_offset)
{
#ifdef MEMORY_TRACE
  printf("insertHole(%ld): ", long(p_offset));
  dumpInternalNode(stdout, p_offset, 0x03);
#endif

  // sanity check to make sure that the first and last slots in this hole
  // have the same value, i.e. -(# of slots in the hole)
  MEDDLY_DCASSERT(data[p_offset] == data[p_offset - data[p_offset] - 1]);

  // If the hole is too small, don't bother to track it
  if (-data[p_offset] < smallestChunk()) {
#ifdef MEMORY_TRACE
    printf("\thole size %d, too small to track\n", -data[p_offset]);
#endif
    return false; 
    // This memory can still be reclaimed 
    // when it is merged with neighboring holes :^)
  }

  // Check if we belong in the grid, or the large hole heap
  if (-data[p_offset] > max_request) {
#ifdef MEMORY_TRACE
    printf("\tAdding to large_holes: %d\n", large_holes);
#endif

    addToHeap(large_holes, large_holes_size, p_offset);
    return true;
  }

  // special case: empty
  if (0 == holes_bottom) {
#ifdef MEMORY_TRACE
    printf("\tAdding to empty grid\n");
#endif
    // index hole
    makeIndex(p_offset, 0, 0, 0, 0);
    return true;
  }

  // special case: at top
  if (data[p_offset] < data[holes_top]) {
#ifdef MEMORY_TRACE
    printf("\tAdding new chain at top\n");
#endif
    // index hole
    makeIndex(p_offset, 0, holes_top, 0, 0);
    return true;
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

    if (p_offset < above) {
      // This node is smaller than the index node.
      // Make this the new index node, and add the current index node
      // to the heap.

      node_handle root = Root(above);
      node_handle size = Size(above);
      node_handle upptr = Up(above);
      addToHeap(root, size, above);
      makeIndex(p_offset, upptr, below, root, size);

    } else {
      // This node is larger than the index node.
      // Add it to the heap.

      addToHeap(Root(above), Size(above), p_offset); 
    }
    return true;
  }
#ifdef MEMORY_TRACE
  printf("\tAdding new chain\n");
#endif
  // we should have above < p_offset < below  (remember, -sizes)
  // create an index hole since there were no holes of this size
  makeIndex(p_offset, above, below, 0, 0);
  return true;
}

// ******************************************************************

void MEDDLY::hm_heap::removeHole(node_handle p_offset)
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(data[p_offset] < 0);

  if (-data[p_offset] < smallestChunk()) {
    useUntracked(-data[p_offset]);
    return; // not tracked
  }

  useHole(-data[p_offset]);

#ifdef MEMORY_TRACE
  printf("removeHole(%ld): ", long(p_offset));
  dumpHole(stdout, p_offset);
#endif
  if (isIndexHole(p_offset)) {
    //
    // Index hole.
    // Positive: no need to search for our location.
    // Negative: need to build another index node
    //
#ifdef MEMORY_TRACE
    printf("\tremoving index hole\n");
#endif
    node_handle above, below, root, size;
    getIndex(p_offset, above, below, root, size);

    if (root) {
      // Heap is not empty.
      // Remove smallest from the heap,
      // and make that the new index node.

      node_handle newindex = root;
      removeFromHeap(root, size, root);
      makeIndex(newindex, above, below, root, size);

#ifdef MEMORY_TRACE
      printf("New index %ld: ", long(newindex));
      dumpHole(stdout, newindex);
#endif

    } else {
      // Heap is empty.
      // No more holes of this size; remove the row completely

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
#endif
    //
    // Done removing index hole.
    //
  } else {
    //
    // Hole is in some heap.
    // Positive: no need to update index node
    // Negative: need to determine the index node that holds the 
    //           heap we're in (unless we're in the large_holes heap)
    //

    if (-data[p_offset] > max_request) {
      // Easy case!
#ifdef MEMORY_TRACE
      printf("removing large hole\n");
#endif
      MEDDLY_DCASSERT(large_holes);
      removeFromHeap(large_holes, large_holes_size, p_offset);

    } else {
      // Hard case
#ifdef MEMORY_TRACE
      printf("removing from heap\n");
#endif
      // Find the correct heap
      node_handle curr;
      for (curr = holes_bottom; curr; curr = Up(curr)) {
        if (data[p_offset] == data[curr]) break;
      }
      MEDDLY_DCASSERT(curr);

      removeFromHeap(Root(curr), Size(curr), p_offset);

    }

#ifdef MEMORY_TRACE
    printf("Removed Heap hole %ld\n", long(p_offset));
#endif
    //
    // Done removing heap hole.
    //
  }

#ifdef DEEP_MEMORY_TRACE
  dumpInternal(stdout, INTERNAL_CODE);
#endif
}

// ******************************************************************

void MEDDLY::hm_heap
::addToHeap(node_handle& root, node_handle& size, node_handle p)
{
#ifdef MEMORY_TRACE
  printf("    addToHeap(%ld, %ld, %ld)\n", long(root), long(size), long(p));
#endif
  MEDDLY_DCASSERT(p);
  MEDDLY_DCASSERT(size<=0);
  // size remains negative
  size--;
  node_handle IDp = -size;
  if (1==IDp) {
    // empty heap
    makeHeapNode(p, IDp, 0, 0);
    rawParent(p) = 0;
    root = p;
    return;
  }
  node_handle pp = takePathToID(root, IDp);
  MEDDLY_DCASSERT(pp);
  MEDDLY_DCASSERT(ID(pp) == IDp/2);
  makeHeapNode(p, IDp, 0, 0);
  upHeap(root, pp, p);
#ifdef DEVELOPMENT_CODE
  verifyHeap(0, root);
#endif
}

// ******************************************************************

void MEDDLY::hm_heap
::removeFromHeap(node_handle& root, node_handle& size, node_handle p)
{
#ifdef MEMORY_TRACE
  printf("    removeFromHeap(%ld, %ld, %ld)\n", long(root), long(size), long(p));
  verifyHeap(0, root);
#endif
  MEDDLY_DCASSERT(p);
  MEDDLY_DCASSERT(p == takePathToID(root, ID(p)));

  // use rightmost leaf node as the replacement node
  node_handle replace = takePathToID(root, -size);
  size++;

  // disconnect rightmost leaf node
  node_handle rmp = Parent(replace);
  if (rmp) {
    if (replace == Right(rmp)) {
      setRight(rmp, 0);
    } else {
      MEDDLY_DCASSERT(replace = Left(rmp));
      setLeft(rmp, 0);
    }
  }

  // Lucky special case:
  // if we're deleting the rightmost leaf node, we're done!
  if (replace == p) {
    if (0==size) root = 0;  // empty heap now
    return;
  }
  
  // repair the heap downward:
  // where p was, put either left child, right child, or replacement
  node_handle* pptr = holeOf(p);
  node_handle newp = downHeap(pptr[ID_index], 
    pptr[left_index], pptr[right_index], replace);

  // repair the heap upward
  upHeap(root, pptr[parent_index], newp);

#ifdef DEVELOPMENT_CODE
  verifyHeap(0, root);
#endif
}

// ******************************************************************

void MEDDLY::hm_heap::upHeap(node_handle& root, node_handle pp, node_handle p)
{
  // p is in position; now perform an "upheap" operation
  while (pp) {
    bool right_child = ID(p) % 2;

    if (right_child) {
      setRight(pp, p);
    } else {
      setLeft(pp, p);
    }

    if (pp < p) break;  

    // sanity checks
    MEDDLY_DCASSERT(ID(pp) > 0);
    MEDDLY_DCASSERT(ID(p) > 0);
    MEDDLY_DCASSERT(ID(p) / 2 == ID(pp));
    MEDDLY_DCASSERT(Parent(p) == pp);

    // rotate
    node_handle newpp = Parent(pp);
    if (right_child) {
      MEDDLY_DCASSERT(Right(pp) == p);
      node_handle L = Left(pp);
      setLeft(pp, Left(p));
      setRight(pp, Right(p));
      setLeft(p, L);
      setRight(p, pp);
    } else {
      MEDDLY_DCASSERT(Left(pp) == p);
      node_handle R = Right(pp);
      setLeft(pp, Left(p));
      setRight(pp, Right(p));
      setLeft(p, pp);
      setRight(p, R);
    }
    rawID(pp) = ID(p);
    rawID(p) /= 2;

    pp = newpp;
  }

  if (0==pp) {
    root = p;
    rawParent(p) = 0;
  }
}

// ******************************************************************

MEDDLY::node_handle MEDDLY::hm_heap::downHeap(
  node_handle ID, node_handle left, node_handle right, node_handle replace)
{
  MEDDLY_DCASSERT(replace);
  node_handle kl = left ? left : replace+1;
  node_handle kr = right ? right : replace+1;

  if (kl < kr) {
    if (replace < kl) {
        // replace is smallest
        makeHeapNode(replace, ID, left, right);
        return replace;
    } else {
        // kl is smallest
        node_handle newleft = 
          downHeap(2*ID, Left(left), Right(left), replace);
        makeHeapNode(left, ID, newleft, right);
        return left;
      }
  } else {
      if (replace < kr) {
        // replace is smallest
        makeHeapNode(replace, ID, left, right);
        return replace;
      } else {
        // kr is smallest
        node_handle newright = 
          downHeap(2*ID+1, Left(right), Right(right), replace);
        makeHeapNode(right, ID, left, newright);
        return right;
      }
  }
}

// ******************************************************************

void MEDDLY::hm_heap::heapToList(node_handle root)
{
#ifdef MEMORY_TRACE
  printf("Converting heap into list...\n");
#endif
  node_handle front, back;
  back = root;
  if (back) {
    rawParent(back) = 0;
  }
  for (front = root; front; front = Parent(front)) {
    if (Left(front)) {
      back = ( rawParent(back) = Left(front) );
      if (Right(front)) {
        back = ( rawParent(back) = Right(front) );
      } // if right child
      rawParent(back) = 0;
    } // if left child
  }
#ifdef MEMORY_TRACE
  printf("Done.  List:  ");
  for (node_handle ptr = root; ptr; ptr = Parent(ptr)) {
    if (ptr != root) printf(", ");
    printf("%ld", long(ptr));
  }
  printf("\n");
#endif
}

// ******************************************************************

#ifdef DEVELOPMENT_CODE
void MEDDLY::hm_heap::verifyHeap(node_handle parent, node_handle node) const
{
#ifdef VERIFY_HEAP
  if (0==node) return;
  MEDDLY_DCASSERT(parent < node);
  MEDDLY_DCASSERT(Parent(node) == parent);
  if (parent) {
    MEDDLY_DCASSERT(ID(parent) == ID(node)/2);
  } else {
    MEDDLY_DCASSERT(1 == ID(node));
  }
  verifyHeap(node, Left(node));
  verifyHeap(node, Right(node));
#endif
}
#endif
