
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


#ifndef HM_HEAP_H
#define HM_HEAP_H

#include "holeman.h"

namespace MEDDLY {
  class hm_heap;
};


/** Grid-based hole management with heaps.

    We have a doubly-linked list of heaps of holes.
    Each heap contains holes of the same size.
    Holes are ordered by address in the heaps,
    and the policy is that we use the smallest address hole 
    that is large enough to satisfy a request.
    Note that the index node has lower address than
    all index nodes in the heap it points to.

    Index node representation:
    ---------------------------------------
    [0] -size (number of slots in hole)     
    [1] up (to larger holes)
    [2] down (to smaller holes)
    [3] root (pointer to heap)
    [4] -heapsize
    :
    :
    [size-1] -size
    ---------------------------------------


    Heap node representation
    ---------------------------------------
    [0] -size (number of slots in hole)     
    [1] left child
    [2] right child
    [3] parent (because we can)
    [4] ID (gives position in the heap)
    :
    :
    [size-1] -size
    ---------------------------------------

  
    Heaps:
    Because we cannot have an array implementation of each heap,
    we use an ID to specify where each node is in the heap, where
    the root node has ID 1, and a generic node with ID i has a parent
    with ID i/2 and children with IDs 2i and 2i+1:

                            1
                  2                   3 
             4         5         6         7
          08   09   10   11   12   13   14   15
*/
class MEDDLY::hm_heap : public holeman {
  public:
    hm_heap(node_storage* p);
    virtual ~hm_heap();

  public:
    // virtual node_address requestChunk(int slots);
    // virtual void recycleChunk(node_address addr, int slots);
    // virtual void dumpInternalInfo(FILE* s) const;
    // virtual void dumpHole(FILE* s, node_address a) const;
    // virtual void dumpInternalTail(FILE* s) const;
    // virtual void reportStats(FILE* s, const char* pad, unsigned flags) const;
    // virtual void clearHolesAndShrink(node_address new_last, bool shrink);

  private:
    // add a hole to the grid or large hole list, as appropriate
    void insertHole(node_handle p_offset);

    // remove a hole from the grid
    void removeHole(node_handle p_offset);

  private:
    static const int left_index = 1;
    static const int right_index = 2;
    static const int parent_index = 3;
    static const int ID_index = 4;

    static const int up_index = 1;
    static const int down_index = 2;
    static const int root_index = 3;
    static const int size_index = 4;

  private:
      inline node_handle* holeOf(node_handle addr) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_CHECK_RANGE(1, addr, lastSlot()+1);
        MEDDLY_DCASSERT(data[addr] < 0);  // it's a hole
        return data + addr;
      }
      inline bool isIndexHandle(node_handle addr) const {
        return holeOf(addr)[size_index] <= 0;
      }
      inline bool isNotIndexHandle(node_handle addr) const {
        return holeOf(addr)[ID_index] > 0;
      }

      inline node_handle& rawLeft(node_handle addr) const {
        MEDDLY_DCASSERT(isNotIndexHandle(addr));
        return holeOf(addr)[left_index];
      }
      inline node_handle& rawRight(node_handle addr) const {
        MEDDLY_DCASSERT(isNotIndexHandle(addr));
        return holeOf(addr)[right_index];
      }
      inline node_handle& rawParent(node_handle addr) const {
        MEDDLY_DCASSERT(isNotIndexHandle(addr));
        return holeOf(addr)[parent_index];
      }
      inline node_handle& rawID(node_handle addr) const {
        MEDDLY_DCASSERT(isNotIndexHandle(addr));
        return holeOf(addr)[ID_index];
      }

      inline node_handle Left(node_handle addr) const {
        return rawLeft(addr); 
      }
      inline node_handle& Left(node_handle addr) {
        return rawLeft(addr); 
      }
      inline node_handle Right(node_handle addr) const {
        return rawRight(addr); 
      }
      inline node_handle& Right(node_handle addr) {
        return rawRight(addr); 
      }
      inline node_handle Parent(node_handle addr) const {
        return rawParent(addr); 
      }
      inline node_handle ID(node_handle addr) const {
        return rawID(addr); 
      }

      inline node_handle& rawUp(node_handle addr) const {
        MEDDLY_DCASSERT(isIndexHandle(addr));
        return holeOf(addr)[up_index];
      }
      inline node_handle& rawDown(node_handle addr) const {
        MEDDLY_DCASSERT(isIndexHandle(addr));
        return holeOf(addr)[down_index];
      }
      inline node_handle& rawRoot(node_handle addr) const {
        MEDDLY_DCASSERT(isIndexHandle(addr));
        return holeOf(addr)[root_index];
      }
      inline node_handle& rawSize(node_handle addr) const {
        MEDDLY_DCASSERT(isIndexHandle(addr));
        return holeOf(addr)[size_index];
      }

      inline node_handle Up(node_handle addr) const {
        return rawUp(addr); 
      }
      inline node_handle Down(node_handle addr) const {
        return rawDown(addr); 
      }
      inline node_handle Root(node_handle addr) const {
        return rawRoot(addr); 
      }
      inline node_handle Size(node_handle addr) const {
        return rawSize(addr); 
      }

      // set all heap pointers at once
      inline void fillHeapNode(node_handle* node, node_handle ID,
        node_handle parent, node_handle left, node_handle right) 
      {
        node[left_index] = left;
        node[right_index] = right;
        node[parent_index] = parent;
        node[ID_index] = ID;
        MEDDLY_DCASSERT(node[ID_index] > 0);
      }
      // set all heap pointers at once
      inline void makeHeapNode(node_handle addr, node_handle ID, 
        node_handle parent, node_handle left, node_handle right) 
      {
        fillHeapNode(holeOf(addr), ID, parent, left, right);
      }

  private:

      // find the last non-null node on the path to the given ID.
      inline node_handle takePathToID(node_handle root, unsigned long ID) const
      {
        // ID bit pattern gives the path :^)
        // determine which bit to start from...
        unsigned long two2b = 0x1;
        int bit = 0;
        for (; two2b <= ID; two2b <<=1) { bit++; }
        two2b >= 1;
        // now, traverse
        for (node_handle n=root; n; ) {
          two2b >>= 1;
          if (0==--bit) return n;   // we're at the right position
          node_handle child = (ID & two2b) ? Right(n) : Left(n);
          if (!child) return n;   // return before we become null
          n = child;
        }
      }

      // add node (already with the proper ID) to the given heap
      void addToHeap(node_handle& root, node_handle& size, node_handle p);

      // remove a node known to be in this heap
      void removeFromHeap(node_handle& root, node_handle& size, node_handle p);

      // repair heap downwards
      node_handle downHeap(node_handle parent, node_handle ID, 
        node_handle left, node_handle right, node_handle replace);

};

#endif

