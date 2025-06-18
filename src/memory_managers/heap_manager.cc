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

// #define TRACE_REALLOCS
// #define MEMORY_TRACE
// #define MEMORY_TRACE_DETAILS
// #define DEBUG_HEAP

#include "hole_base.h"
#include "heap_manager.h"
#include "../operators.h"

namespace MEDDLY {

  // ******************************************************************
  // *                                                                *
  // *                 heap_manager  (template) class                 *
  // *                                                                *
  // ******************************************************************

  /**
    Grid structure for hole management.
    INT should be a signed integer type.

    Small hole (fewer than 5 slots):
    -------------------------------------------------------------
    [0] size with MSB set
    :
    : unused
    :
    [size-1] size with MSB set


    Heap node hole (5 or more slots):
    -------------------------------------------------------------
    [0] size with MSB set
    [1] parent, >= 0: pointer to parent in heap
    [2] left, >= 0: pointer to left child in heap
    [3] right, >= 0: pointer to right child in heap
    :
    : unused
    :
    [size-1] size with MSB set


  */
  template <class INT>
  class heap_manager : public hole_manager<INT> {

    public:
      heap_manager(const char* n, memstats &stats);
      virtual ~heap_manager();

      virtual node_address requestChunk(size_t &numSlots);
      virtual void recycleChunk(node_address h, size_t numSlots);

      virtual void reportStats(output &s, const char* pad, bool human, bool details) const;
      virtual void dumpInternal(output &s) const;
      virtual void dumpInternalUnused(output &s, node_address addr) const;

    private:
      /**
          Find and return the hole at position id in the heap.
          The root node is position 1.
          If a node has position n, it's parent has position n/2,
          its left child has position 2n,
          and its right child has position 2n+1.
            @param  id    Position number to find
            @param  Address of hole at that position, or 0 if none.
      */
      node_address findNodeAtPosition(unsigned long id) const;


      /**
          Remove the "last" node in the heap.
          That's the node with the highest position,
          and the rightmost node at the bottom level.
      */
      node_address removeLastHeapNode();

      /**
          Remove the node corresponding to the given address.
      */
      void removeHeapNode(node_address n);

      /**
          Up heap operation starting at node n.
          Move n up in the heap, while its size
          is larger than its parent's size.
      */
      void upHeap(node_address n);

      /**
          Down heap operation starting at node n.
          Move n down in the heap, while its size
          is smaller than at least one of its children.
      */
      void downHeap(node_address n);

    private:
      inline INT Parent(node_address h) const {
        return hole_manager<INT>::readSlot(h, 1);
      }

      inline INT Left(node_address h) const {
        return hole_manager<INT>::readSlot(h, 2);
      }

      inline INT Right(node_address h) const {
        return hole_manager<INT>::readSlot(h, 3);
      }

      inline void zeroPointers(node_address h) {
        hole_manager<INT>::refSlot(h, 1) = 0; // parent
        hole_manager<INT>::refSlot(h, 2) = 0; // left
        hole_manager<INT>::refSlot(h, 3) = 0; // right
      }

      inline void makeRoot(node_address h) {
        hole_manager<INT>::refSlot(h, 1) = 0;
        heap_root = h;
      }

      inline void setLeft(node_address h, node_address L) {
        hole_manager<INT>::refSlot(h, 2) = L;
        if (L) {
          // Set L's parent to h
          hole_manager<INT>::refSlot(L, 1) = h;
        }
      }

      inline void setRight(node_address h, node_address R) {
        hole_manager<INT>::refSlot(h, 3) = R;
        if (R) {
          // Set R's parent to h
          hole_manager<INT>::refSlot(R, 1) = h;
        }
      }

      inline static INT smallestChunk() {
        return 5;
      }

      inline bool isHole(node_address h) const {
        return hole_manager<INT>::isHole(h);
      }
      inline INT getHoleSize(node_address h) const {
        return hole_manager<INT>::getHoleSize(h);
      }
      inline void setHoleSize(node_address h, INT hs) {
        hole_manager<INT>::setHoleSize(h, hs);
      }

      inline bool isSmallHole(node_address addr) const {
        return getHoleSize(addr) < smallestChunk();
      }

      inline void incHeapNodes() {
        num_heap_nodes++;
        if (num_heap_nodes > max_heap_nodes) {
          max_heap_nodes = num_heap_nodes;
        }
      }
      inline void decHeapNodes() {
        num_heap_nodes--;
      }
      inline void incHeapSlots(long a) {
        MEDDLY_DCASSERT(a>0);
        num_heap_slots += a;
        if (num_heap_slots > max_heap_slots) {
          max_heap_slots = num_heap_slots;
        }
      }
      inline void decHeapSlots(long a) {
        MEDDLY_DCASSERT(a>0);
        num_heap_slots -= a;
      }
      inline void incSmallSlots(long a) {
        if (0==a) return;
        MEDDLY_DCASSERT(a>0);
        num_small_slots += a;
        if (num_small_slots > max_small_slots) {
          max_small_slots = num_small_slots;
        }
        num_small_holes++;
        if (num_small_holes > max_small_holes) {
          max_small_holes = num_small_holes;
        }
      }
      inline void decSmallSlots(long a) {
        MEDDLY_DCASSERT(a>0);
        num_small_slots -= a;
        num_small_holes--;
      }

    private:
      /// Current number of nodes in the heap.
      long num_heap_nodes;

      /// Pointer to root node of heap
      node_address heap_root;

      /// Current hole (not in heap) we're trying to use
      node_address current_hole;

    private:    // stats
      long max_heap_nodes;

      long num_small_holes;
      long max_small_holes;

      long num_small_slots;
      long max_small_slots;
      long num_heap_slots;
      long max_heap_slots;

  }; // class original_grid

  // ******************************************************************
  // *                                                                *
  // *                      heap_manager methods                      *
  // *                                                                *
  // ******************************************************************


  template <class INT>
  heap_manager<INT>::heap_manager(const char* n, memstats &stats)
  : hole_manager<INT>(n, stats)
  {
    num_heap_nodes = 0;
    heap_root = 0;
    current_hole = 0;

    max_heap_nodes = 0;

    num_small_holes = 0;
    max_small_holes = 0;

    num_small_slots = 0;
    max_small_slots = 0;
    num_heap_slots = 0;
    max_heap_slots = 0;
  }

  template <class INT>
  heap_manager<INT>::~heap_manager()
  {
  }


  // ******************************************************************


  template <class INT>
  MEDDLY::node_address heap_manager<INT>::requestChunk(size_t &numSlots)
  {
#ifdef MEMORY_TRACE
    printf("requestChunk(%lu)\n", numSlots);
#endif

    //
    // If current hole is too small, see if it makes sense
    // to grab from the heap
    //
    if (0==current_hole || (size_t(getHoleSize(current_hole)) < numSlots)) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tcurrent_hole is too small, trying heap\n");
#endif
      if (heap_root && (size_t(getHoleSize(heap_root)) >= numSlots)) {
#ifdef MEMORY_TRACE_DETAILS
        printf("\treplacing current_hole with heap root\n");
#endif

        if (1==num_heap_nodes) {
          //
          // Special case - heap contains only the root
          //
          decHeapNodes();
          decHeapSlots(getHoleSize(heap_root));
          node_address tmp = heap_root;
          heap_root = 0;

          if (current_hole) {
            //
            // Make a heap out of current hole
            //
            MEDDLY_DCASSERT(getHoleSize(current_hole) >= smallestChunk());
            zeroPointers(current_hole);
            incHeapNodes();
            incHeapSlots(getHoleSize(current_hole));
          }
          heap_root = current_hole;
          current_hole = tmp;
        } else {
          //
          // Remove root from heap and simultaneously add
          // current_hole to heap.
          // We do this by making current_hole the root and
          // then doing a "downHeap".
          // Unless the current_hole is null, in which case
          // we replace from the end.
          //
          decHeapSlots(getHoleSize(heap_root));
          node_address replace;
          if (0==current_hole) {
            replace = removeLastHeapNode();
          } else {
            replace = current_hole;
          }
          current_hole = heap_root;

          //
          // We have a replacement node, add it to the heap
          //
          incHeapSlots(getHoleSize(replace));

          //
          // Update children
          //
          setLeft(replace, Left(heap_root));
          setRight(replace, Right(heap_root));
          makeRoot(replace);

          //
          // Put new root node in proper position
          //
          downHeap(heap_root);
        }
#ifdef DEBUG_HEAP
        FILE_output out(stdout);
        dumpInternal(out);
#endif
      }
    }

    //
    // Check current hole: if large enough, pull from it
    //
    if (current_hole && (size_t(getHoleSize(current_hole)) >= numSlots)) {

      node_address h = current_hole;

      //
      // Deal with any leftover slots in the hole
      //
      size_t leftover_slots = size_t(getHoleSize(current_hole)) - numSlots;
#ifdef MEMORY_TRACE_DETAILS
      printf("\t %lu remaining slots in current\n", leftover_slots);
#endif
      if (leftover_slots > 0) {
        //
        // Note - we even recycle holes of size 1,
        // because they can be merged to the left or right,
        // depending on who is recycled first
        //
        current_hole += numSlots;
        setHoleSize(current_hole, leftover_slots);
#ifdef MEMORY_TRACE_DETAILS
        printf("\tcurrent hole %ld has size %lu\n", current_hole, leftover_slots);
#endif
      }
      if (leftover_slots < size_t(smallestChunk())) {
        // cannot track this hole, so leave it
        incSmallSlots(leftover_slots);
        current_hole = 0;
#ifdef MEMORY_TRACE_DETAILS
        printf("\tcurrent hole too small to track\n");
#endif
      }

#ifdef MEMORY_TRACE_DETAILS
      printf("requestChunk(%lu) grabbed from current, returned hole %ld\n", numSlots, h);
#endif
      memory_manager::incMemUsed(getHoleSize(h) * sizeof(INT));
      return h;
    }

    //
    // Still here?  We couldn't recycle a chunk.
    //
#ifdef MEMORY_TRACE_DETAILS
    printf("\tNo recycleable chunks large enough, grabbing from end\n");
#endif
    node_address h = hole_manager<INT>::allocateFromArray(numSlots);
    if (0==h) {
      numSlots = 0;
      return 0;
    }
#ifdef MEMORY_TRACE_DETAILS
    printf("requestChunk(%lu) returned %ld\n", numSlots, h);
#endif
    memory_manager::incMemUsed(numSlots * sizeof(INT));
    return h;
  }


  // ******************************************************************


  template <class INT>
  void heap_manager<INT>::recycleChunk(node_address h, size_t numSlots)
  {
#ifdef MEMORY_TRACE
    printf("recycling chunk %lu size %lu\n", h, numSlots);
#endif
    memory_manager::decMemUsed(numSlots * sizeof(INT));

    setHoleSize(h, numSlots);

    //
    // Check to the left for another hole
    //
    if (isHole(h-1)) {
      MEDDLY_DCASSERT(node_address(getHoleSize(h-1)) < h);
      node_address hleft = h - getHoleSize(h-1);
#ifdef MEMORY_TRACE_DETAILS
      printf("\tMerging to the left, holes %lu and %lu\n", hleft, h);
#endif
      if (current_hole != hleft) {
        if (isSmallHole(hleft)) {
          decSmallSlots(getHoleSize(hleft));
        } else {
          // hleft must be in the heap; remove it
          removeHeapNode(hleft);
        }
      }
      numSlots += getHoleSize(hleft);
      h = hleft;
      setHoleSize(h, numSlots);
    }

    //
    // Can we absorb this hole at the end?
    //
    if (hole_manager<INT>::recycleHoleInArray(h, numSlots)) {
      if (h==current_hole) {
        current_hole = 0;
      }
      return;
    }

    //
    // Check to the right for another hole
    //
    if (isHole(h + numSlots)) {
      node_address hright = h + numSlots;
#ifdef MEMORY_TRACE_DETAILS
      printf("\tMerging to the right, holes %lu and %lu\n", h, hright);
#endif
      if (current_hole != hright) {
        if (isSmallHole(hright)) {
          decSmallSlots(getHoleSize(hright));
        } else {
          // hright must be in the heap; remove it
          removeHeapNode(hright);
        }
      }
      numSlots += getHoleSize(hright);
      setHoleSize(h, numSlots);
      if (hright == current_hole) {
        current_hole = h;
      }
    }

    //
    // Hole is ready
    //
    if (h == current_hole) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tRecycled chunk merged with current hole: %ld\n", current_hole);
#endif
      return;
    }
    if (isSmallHole(h)) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tRecycled chunk is too small for the heap\n");
#endif
      incSmallSlots(getHoleSize(h));
      return;
    }

    //
    // Add to heap
    //
    incHeapNodes();
    incHeapSlots(numSlots);
    zeroPointers(h);
    if (0==heap_root) {
      // special case - start a new heap
#ifdef MEMORY_TRACE_DETAILS
      printf("\tCreating heap from hole\n");
#endif
      MEDDLY_DCASSERT(1==num_heap_nodes);
      makeRoot(h);
    } else {
      // find parent of where we will go
#ifdef MEMORY_TRACE_DETAILS
      printf("\tAdding hole to heap\n");
#endif
      node_address p = findNodeAtPosition(num_heap_nodes/2);
      MEDDLY_DCASSERT(p);
      MEDDLY_DCASSERT(0==Right(p));
      if (0==Left(p)) {
        setLeft(p, h);
      } else {
        setRight(p, h);
      }
      upHeap(h);
    }

#ifdef DEBUG_HEAP
    FILE_output out(stdout);
    dumpInternal(out);
#endif
  }


  // ******************************************************************


  template <class INT>
  void heap_manager<INT>::reportStats(output &s, const char* pad, bool human, bool details) const
  {
    s << pad << "Report for heap memory manager:\n";
    s << pad << "  Current #holes: ";
    s << num_small_holes+num_heap_nodes << "\n";
    s << pad << "  Current bytes in holes: ";
    s.put_mem( (num_small_slots+num_heap_slots) * sizeof(INT), human );
    s << "\n";
    if (details) {
      s << pad << "  Hole counts:\n";
      s << pad << "      " << num_small_holes << " current small (untracked)\n";
      s << pad << "      " << num_heap_nodes << " current in heap\n";
      s << pad << "      " << max_small_holes << " maximum small (untracked)\n";
      s << pad << "      " << max_heap_nodes << " maximum in heap\n";
      s << pad << "  Hole bytes:\n";
      s << pad << "      ";
      s.put_mem( num_small_slots * sizeof(INT), human );
      s << " current small (untracked)\n";
      s << pad << "      ";
      s.put_mem( num_heap_slots * sizeof(INT), human );
      s << " current in Heap\n";
      s << pad << "      ";
      s.put_mem( max_small_slots * sizeof(INT), human );
      s << " maximum small (untracked)\n";
      s << pad << "      ";
      s.put_mem( max_heap_slots * sizeof(INT), human );
      s << " maximum in Heap\n";
    }
  }


  // ******************************************************************


  template <class INT>
  void heap_manager<INT>::dumpInternal(output &s) const
  {
    s << "Internal storage for heap:\n";
    hole_manager<INT>::showInternal(s);
    if (0==heap_root) {
      s << "  Empty heap\n";
    } else {
      s << "  Root is " << heap_root << "\n";
      for (long id=1; id <= num_heap_nodes; id++) {
        node_address addr = findNodeAtPosition(id);
        s << "    Node at position " << long(id) << ": " << addr << "\n";
        s << "        Size  : " << getHoleSize(addr) << "\n";
        s << "        Parent: " << Parent(addr) << "\n";
        s << "        Left  : " << Left(addr) << "\n";
        s << "        Right : " << Right(addr) << "\n";
        s << "\n";
      }
    }
    s << "  Current hole: " << current_hole;
    if (current_hole) s << " size " << getHoleSize(current_hole);
    s << "\n";
    s << "  num_small_holes: " << num_small_holes << "\n";
    s << "  max_small_holes: " << max_small_holes << "\n";
    s << "  num_small_slots: " << num_small_slots << "\n";
    s << "  max_small_slots: " << max_small_slots << "\n";
    s << "  num_heap_nodes: " << num_heap_nodes << "\n";
    s << "  max_heap_nodes: " << max_heap_nodes << "\n";
    s << "  num_heap_slots: " << num_heap_slots << "\n";
    s << "  max_heap_slots: " << max_heap_slots << "\n";
  }


  // ******************************************************************


  template <class INT>
  void heap_manager<INT>::dumpInternalUnused(output &s, node_address addr) const
  {
    hole_manager<INT>::showInternalAddr(s, addr, 3);
  }


  // ******************************************************************

  template <class INT>
  MEDDLY::node_address heap_manager<INT>::findNodeAtPosition(unsigned long id) const
  {
    if (0==heap_root) return 0;
    //
    // The bit pattern of id gives the path:
    // Ignore the MSB, then read the bit pattern from most
    // to least significant, 0 goes left, 1 goes right.
    //

    //
    // Determine which bit to start from
    //
    unsigned long two2b = 0x1;
    int bit = 0;
    for (; two2b <= id; two2b <<=1) { bit++; }
    two2b >>= 1;
    //
    // Traverse
    //
    node_address n = heap_root;
    for (;;) {
      if (0==--bit) return n;
      two2b >>= 1;
      node_address child = (id & two2b) ? Right(n) : Left(n);
      MEDDLY_DCASSERT(child);
      n = child;
    } // for
    return 0;   // We never get here; keep compilers happy
  }

  // ******************************************************************

  template <class INT>
  MEDDLY::node_address heap_manager<INT>::removeLastHeapNode()
  {
    MEDDLY_DCASSERT(heap_root);
    MEDDLY_DCASSERT(num_heap_nodes);
    node_address n = findNodeAtPosition(num_heap_nodes);
    decHeapNodes();
    decHeapSlots(getHoleSize(n));
    MEDDLY_DCASSERT(n);
    MEDDLY_DCASSERT(Left(n)==0);
    MEDDLY_DCASSERT(Right(n)==0);
    node_address p = Parent(n);
    if (0==p) {
      // this is the root
      MEDDLY_DCASSERT(n==heap_root);
      heap_root = 0;
      return n;
    }
    // Not the root, update parent's pointers
    if (Right(p)==0) {
      // we must be the left child
      MEDDLY_DCASSERT(node_address(Left(p))==n);
      setLeft(p, 0);
    } else {
      // we must be the right child
      MEDDLY_DCASSERT(node_address(Right(p))==n);
      setRight(p, 0);
    }
    return n;
  }

  // ******************************************************************

  template <class INT>
  void heap_manager<INT>::removeHeapNode(node_address target)
  {
    MEDDLY_DCASSERT(target);
    MEDDLY_DCASSERT(heap_root);

    node_address replace = removeLastHeapNode();
    MEDDLY_DCASSERT(replace);

    if (replace == target) return;

    //
    // Jam the removed last node in place of target.
    //
    incHeapSlots(getHoleSize(replace));
    decHeapSlots(getHoleSize(target));

    //
    // Fix parent's pointer
    //
    node_address p = Parent(target);
    if (0==p) {
      MEDDLY_DCASSERT(heap_root == target);
      makeRoot(replace);
    } else {
      if (node_address(Left(p)) == target) {
        setLeft(p, replace);
      } else {
        MEDDLY_DCASSERT(node_address(Right(p)) == target);
        setRight(p, replace);
      }
    }

    //
    // Update replace's pointers
    //
    setLeft(replace, Left(target));
    setRight(replace, Right(target));

    //
    // target is out of the heap.
    // Now, put replacement node in the proper position.
    //
    downHeap(replace);
  }

  // ******************************************************************

  template <class INT>
  void heap_manager<INT>::upHeap(node_address n)
  {
    MEDDLY_DCASSERT(n);
    for (;;) {
      node_address p = Parent(n);
      if (0==p) {
        MEDDLY_DCASSERT(heap_root == n);
        return;
      }
      if (getHoleSize(n) <= getHoleSize(p)) return; // done swapping

      //
      // Rotation.
      // When we're done, n and p swap positions in the tree.
      //

      //
      // Make p's parent point to n instead of p
      //
      node_address gp = Parent(p);
      if (0==gp) {
        MEDDLY_DCASSERT(heap_root == p);
        makeRoot(n);
      } else {
        if (node_address(Left(gp)) == p) {
          setLeft(gp, n);
        } else {
          MEDDLY_DCASSERT(node_address(Right(gp)) == p);
          setRight(gp, n);
        }
      }

      //
      // Fix p and n's children pointers
      //
      if (node_address(Left(p)) == n) {
        //
        // n was the left child of p.
        // Make p the left child of n.
        //
        setLeft(p, Left(n));
        setLeft(n, p);
        node_address tmp = Right(p);
        setRight(p, Right(n));
        setRight(n, tmp);
      } else {
        MEDDLY_DCASSERT(node_address(Right(p)) == n);
        //
        // n was the right child of p.
        // make p the right child of n.
        setRight(p, Right(n));
        setRight(n, p);
        node_address tmp = Left(p);
        setLeft(p, Left(n));
        setLeft(n, tmp);
      }
    } // infinite loop
  }

  // ******************************************************************

  template <class INT>
  void heap_manager<INT>::downHeap(node_address n)
  {
    MEDDLY_DCASSERT(n);

    for (;;) {

      node_address l = Left(n);
      node_address r = Right(n);

      if (0==l) return; // both children are null

      bool promote_left;

      if (getHoleSize(l) <= getHoleSize(n)) {
        // Might be done swapping...
        if (0==r) return;
        if (getHoleSize(r) <= getHoleSize(n)) return;

        //
        // Need to rotate, but we know exactly how
        //
        promote_left = false;
      } else {
        MEDDLY_DCASSERT(getHoleSize(l) > getHoleSize(n));
        if (0==r || (getHoleSize(r) <= getHoleSize(n))) {
          promote_left = true;
        } else {
          MEDDLY_DCASSERT(r);
          MEDDLY_DCASSERT(getHoleSize(r) > getHoleSize(n));
          promote_left = getHoleSize(l) > getHoleSize(r);
        }
      }

      //
      // Need to promote a child, and we even know which one :^)
      //
      node_address p = Parent(n);
      node_address newn;

      if (promote_left) {
        // promote left child
        newn = Left(n);
        setLeft(n, Left(newn));
        setLeft(newn, n);
        node_address tmp = Right(newn);
        setRight(newn, Right(n));
        setRight(n, tmp);
        // done promote left child
      } else {
        // prmote right child
        newn = Right(n);
        setRight(n, Right(newn));
        setRight(newn, n);
        node_address tmp = Left(newn);
        setLeft(newn, Left(n));
        setLeft(n, tmp);
        // done promote right child
      }
      // update parent
      if (p) {
        if (node_address(Left(p)) == n) {
          setLeft(p, newn);
        } else {
          MEDDLY_DCASSERT(node_address(Right(p)) == n);
          setRight(p, newn);
        }
      } else {
        makeRoot(newn);
      }
    } // infinite loop
  }

  // ******************************************************************

}; // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                                                                *
// *                       heap_style methods                       *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::heap_style::heap_style(const char* n)
: memory_manager_style(n)
{
}

MEDDLY::heap_style::~heap_style()
{
}

MEDDLY::memory_manager*
MEDDLY::heap_style::initManager(unsigned char granularity,
  unsigned char minsize, memstats &stats) const
{
  if (sizeof(int) == granularity) {
    return new heap_manager <int>(getName(), stats);
  }

  if (sizeof(long) == granularity) {
    return new heap_manager <long>(getName(), stats);
  }

  if (sizeof(short) == granularity) {
    return new heap_manager <short>(getName(), stats);
  }

  // unsupported granularity

  return 0;
}

