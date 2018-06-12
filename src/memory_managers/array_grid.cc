
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
// #define DEBUG_GRID
// #define REPORT_ALL_MEDIUM_LISTS

#include "hole_base.h"
#include "array_grid.h"

namespace MEDDLY {

  // ******************************************************************
  // *                                                                *
  // *                  array_grid  (template) class                  *
  // *                                                                *
  // ******************************************************************

  /**
    Hole management using an array of lists for smaller holes,
    and a grid structure for larger holes.
    Huge holes (larger than ever requested) go in their own doubly-linked list.
    INT should be a signed integer type.

    Small hole (untracked):
    -------------------------------------------------------------
    Minimum size: 1
    [0] size with MSB set
    :
    : unused
    :
    [size-1] size with MSB set


    Medium hole (in a list with other holes of this size):
    -------------------------------------------------------------
    Minimum size: 4
    [0] size with MSB set
    [1] address of previous
    [2] address of next
    :
    : unused
    :
    [size-1] size with MSB set

    
    Array structure for medium holes:
      medium_list[i] : pointer to front of list with holes of size i.
    Array goes from element 0 (unused) to element LargeNodeSize-1.


    Large index hole (in grid):
    -------------------------------------------------------------
    [0] size with MSB set
    [1] prev pointer
    [2] next pointer
    [3] up, >= 0, goes to larger holes
    [4] down, >= 0, goes to smaller holes
    :
    : unused
    :
    [size-1] size with MSB set


    Large non-index hole (in grid):
    -------------------------------------------------------------
    [0] size with MSB set
    [1] prev pointer
    [2] next pointer
    [3] any negative value (only checked in development code assertions)
    [4] any negative value (only checked in development code assertions)
    :
    : unused
    :
    [size-1] size with MSB set


    The grid structure for large holes:
      (holes_bottom)
      holes_of_size_0 (index) <--> (non-index) <--> ... <--> (non-index) --> NULL
      ^
      |
      v
      holes_of_size_1 (index) <--> (non-index) <--> ... --> NULL
      ^
      |
      v
      :
      :
      (holes_top)
  */
  template <class INT>
  class array_plus_grid : public hole_manager<INT> {

    public:
      array_plus_grid(const char* n, forest::statset &stats);
      virtual ~array_plus_grid();

      virtual node_address requestChunk(size_t &numSlots);
      virtual void recycleChunk(node_address h, size_t numSlots);

      virtual void reportStats(output &s, const char* pad, bool human, bool details) const;
      virtual void dumpInternal(output &s) const;
      virtual void dumpInternalUnused(output &s, node_address addr) const;

    private:
      /**
          Move "grid_current" to the row with the requested size.
          If none, grid_current will be either one row below or one row above
          the requested size.
            @param  size    Size of hole we want
            @param  curr    Pointer to move

            @return   0,  if we found a row with exactly this size;
                          otherwise, no row with this size exists, and we return:
                     -1,  if we point to a row with smaller size (or null);
                     +1,  if we point to a row with larger size (or null).
      */
      int moveCurrentToRow(INT size, INT& curr) const;

      /**
          Stop tracking hole with handle h.
          If it is a small hole, do nothing.
          If it is a medium hole, remove it from its list.
          If it is a large hole, remove it from the grid.
          If it is a huge hole, remove it from the huge hole list.
      */
      void stopTrackingHole(node_address h);

      /**
          Start tracking hole with handle h.
          If it is a small hole, do nothing.
          If it is a medium hole, add it to the list for this hole's size.
          If it is a large hole, add it to the grid.
          If it is a huge hole, add it to the huge hole list.
      */
      void startTrackingHole(node_address h);

    private:
      static const INT MediumHoleSize = 4;
      static const INT LargeHoleSize = 6;   // smallest possible value

      inline bool isHole(node_address h) const {
        return hole_manager<INT>::isHole(h);
      }
      inline INT getHoleSize(node_address h) const {
        // stupid compiler
        return hole_manager<INT>::getHoleSize(h);
      }
      inline void setHoleSize(node_address h, INT hs) {
        hole_manager<INT>::setHoleSize(h, hs);
      }

      inline bool matchingHoleSizes(node_address h) const {
        return hole_manager<INT>::matchingHoleSizes(h);
      }

      inline bool isSmallHole(node_address h) const {
        return getHoleSize(h) < MediumHoleSize;
      }
      inline bool isLargeHole(node_address h) const {
        return getHoleSize(h) >= LargeHoleSize;
      }

      // Pointers common to all types of holes except small

      inline INT Prev(node_address h) const {
        MEDDLY_DCASSERT(!isSmallHole(h));
        return hole_manager<INT>::readSlot(h, 1);
      }
      inline void setPrev(node_address h, INT v) {
        MEDDLY_DCASSERT(!isSmallHole(h));
        hole_manager<INT>::refSlot(h, 1) = v;
      }

      inline INT Next(node_address h) const {
        MEDDLY_DCASSERT(!isSmallHole(h));
        return hole_manager<INT>::readSlot(h, 2);
      }
      inline void setNext(node_address h, INT v) {
        MEDDLY_DCASSERT(!isSmallHole(h));
        hole_manager<INT>::refSlot(h, 2) = v;
      }
      

      // Large hole handy methods
      inline void setNonIndex(node_address h) {
        MEDDLY_DCASSERT(isLargeHole(h));
        hole_manager<INT>::refSlot(h, 3) = -1;
        // hole_manager<INT>::refSlot(h, 4) = -1;
      }
      inline bool isIndexHole(node_address h) {
        MEDDLY_DCASSERT(isLargeHole(h));
        return hole_manager<INT>::readSlot(h, 3) >= 0;
      }

      inline INT Up(node_address h) const {
        MEDDLY_DCASSERT(isLargeHole(h));
        return hole_manager<INT>::readSlot(h, 3);
      }
      inline void setUp(node_address h, INT v) {
        MEDDLY_DCASSERT(isLargeHole(h));
        hole_manager<INT>::refSlot(h, 3) = v;
      }
      
      inline INT Down(node_address h) const {
        MEDDLY_DCASSERT(isLargeHole(h));
        return hole_manager<INT>::readSlot(h, 4);
      }
      inline void setDown(node_address h, INT v) {
        MEDDLY_DCASSERT(isLargeHole(h));
        hole_manager<INT>::refSlot(h, 4) = v;
      }
      

    private:

      // ---------- medium hole stuff ----------

      INT medium_hole_list[LargeHoleSize];

      // ---------- large hole stuff ----------

      /// Pointer to top of holes grid
      INT grid_bottom;
      /// Pointer to bottom of holes grid
      INT grid_top;
      /// Pointer where we left off the last search
      INT grid_current;

      // ---------- huge hole stuff ----------

      /// Largest hole ever requested
      size_t max_request;
      /// Pointer to list of huge holes
      INT huge_holes;


    private:    // stats
      long num_small_holes;
      long num_medium_holes[LargeHoleSize]; // element 0 counts for all sizes.
      long num_grid_holes;
      long num_huge_holes;

      long num_small_slots;
      long num_medium_slots;
      long num_grid_slots;
      long num_huge_slots;

  }; // class original_grid
}; // namespace MEDDLY


// ******************************************************************
// *                                                                *
// *                    array_plus_grid  methods                    *
// *                                                                *
// ******************************************************************

template <class INT>
MEDDLY::array_plus_grid<INT>::array_plus_grid(const char* n, forest::statset &stats)
 : hole_manager<INT>(n, stats)
{
  // small hole stuff
  num_small_holes = 0;
  num_small_slots = 0;

  // medium hole stuff
  for (INT i=0; i<LargeHoleSize; i++) {
    medium_hole_list[i] = 0;
    num_medium_holes[i] = 0;
  }

  // large hole stuff
  grid_bottom = 0;
  grid_top = 0;
  grid_current = 0;
  num_grid_holes = 0;
  num_grid_slots = 0;

  // huge hole stuff
  max_request = 0;
  // max_request = LargeHoleSize-1;  
  huge_holes = 0;
  num_huge_holes = 0;
  num_huge_slots = 0;
}

// ******************************************************************

template <class INT>
MEDDLY::array_plus_grid<INT>::~array_plus_grid()
{
}


// ******************************************************************

template <class INT>
MEDDLY::node_address 
MEDDLY::array_plus_grid<INT>::requestChunk(size_t &numSlots)
{
#ifdef MEMORY_TRACE
  printf("requestChunk(%lu)\n", numSlots);
#endif

  //
  // See if we need to update the huge list
  //
  if (numSlots > max_request) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\tNew max request %lu; shifting holes from huge list to grid\n", numSlots);
#endif
    max_request = numSlots;

    //
    // Traverse the huge hole list, and re-insert everything.
    // Holes that are still huge will go back into the huge hole list,
    // and the rest will go where they belong (grid or medium list).
    //

    INT curr = huge_holes;
    huge_holes = 0;
    for (; curr; ) {
      INT next = Next(curr);
      startTrackingHole(curr);  
      curr = next;
    }
#ifdef DEBUG_GRID
    FILE_output out(stdout);
    dumpInternal(out);
#endif
  }

  //
  // See if there is a hole of this size.
  //
  node_address h = 0;
  if (numSlots < LargeHoleSize) {
    // Medium sized hole - try list
    h = medium_hole_list[numSlots];
  } else {
    //
    // Large or very small, it's the same
    // Check the grid for a hole of exactly this size
    //
    if (0==moveCurrentToRow(numSlots, grid_current)) {
      MEDDLY_DCASSERT(grid_current);
      MEDDLY_DCASSERT(size_t(getHoleSize(grid_current)) == numSlots);
      // It's easier to remove a non-index node, try that first
      h = node_address(Next(grid_current));
      if (0==h) {
        h = grid_current;
      }
    }
  }

  //
  // If nothing found yet, try the huge hole list.
  // Anything there should be large enough.
  //
  if (0==h) {
    h = huge_holes;
  }

  
  //
  // Ok, that's everywhere.  If we have a hole,
  // use it and recycle any leftovers.
  //
  if (h) {
    stopTrackingHole(h);
    MEDDLY_DCASSERT(size_t(getHoleSize(h)) >= numSlots);
    memory_manager::incMemUsed(getHoleSize(h) * sizeof(INT));
    size_t leftover_slots = size_t(getHoleSize(h)) - numSlots;
    if (leftover_slots > 0) {
      hole_manager<INT>::clearHole(h, numSlots); // otherwise we'll just merge them again
      recycleChunk(h+numSlots, leftover_slots);
    }
#ifdef MEMORY_TRACE_DETAILS
    printf("requestChunk(%lu) returned hole %ld\n", numSlots, h);
#endif
    return h;
  }

  //
  // Still here?  We couldn't recycle a chunk.
  //
#ifdef MEMORY_TRACE_DETAILS
  printf("\tNo recycleable chunks, grabbing from end\n");
#endif
  h = hole_manager<INT>::allocateFromArray(numSlots);
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
void MEDDLY::array_plus_grid<INT>
::recycleChunk(node_address h, size_t numSlots)
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
    stopTrackingHole(hleft);
    numSlots += getHoleSize(hleft); 
    h = hleft;
    setHoleSize(h, numSlots);
#ifdef DEBUG_GRID
    FILE_output out(stdout);
    dumpInternal(out);
#endif
  }

  //
  // Can we absorb this hole at the end?
  //
  if (hole_manager<INT>::recycleHoleInArray(h, numSlots)) {
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
    stopTrackingHole(hright);
    numSlots += getHoleSize(hright);
    setHoleSize(h, numSlots);
#ifdef DEBUG_GRID
    FILE_output out(stdout);
    dumpInternal(out);
#endif
  }

  //
  // Hole is ready
  //
  startTrackingHole(h);
#ifdef DEBUG_GRID
  FILE_output out(stdout);
  dumpInternal(out);
#endif
  }


// ******************************************************************

template <class INT>
void MEDDLY::array_plus_grid<INT>
::reportStats(output &s, const char* pad, bool human, bool details) const
{
  s << pad << "Report for array_plus_grid memory manager:\n";
  s << pad << "  largest request: " << long(max_request) << " slots\n";
  long total_holes = num_small_holes + num_grid_holes 
    + num_huge_holes + num_medium_holes[0];
  s << pad << "  Current #holes: " << total_holes << "\n";
  if (details) {
    s << pad << "      " << num_small_holes << " small (untracked)\n";
#ifdef REPORT_ALL_MEDIUM_LISTS
    for (int i=1; i<LargeHoleSize; i++) {
      if (num_medium_holes[i]) {
        s << pad << "      " << num_medium_holes[i] << " in medium " << i << " list\n";
      }
    }
#else
    s << pad << "      " << num_medium_holes[0] << " in medium lists\n";
#endif
    s << pad << "      " << num_grid_holes << " in grid\n";
    s << pad << "      " << num_huge_holes << " in huge hole list\n";
  }
  long total_slots = num_small_slots + num_grid_slots
    + num_medium_slots + num_huge_slots;
  s << pad << "  Current bytes in holes: ";
  s.put_mem(total_slots * sizeof(INT), human);
  s << "\n";
  if (details) {
    s << pad << "      ";
    s.put_mem( num_small_slots * sizeof(INT), human );
    s << " small (untracked)\n";
#ifdef REPORT_ALL_MEDIUM_LISTS
    for (int i=1; i<LargeHoleSize; i++) {
      if (num_medium_holes[i]) {
        s << pad << "      ";
        s.put_mem( num_medium_holes[i] * i * sizeof(INT), human );
        s << " in medium " << i << " list\n";
      }
    }
#else
    s << pad << "      ";
    s.put_mem( num_medium_slots * sizeof(INT), human );
    s << " in medium lists\n";
#endif
    s << pad << "      ";
    s.put_mem( num_grid_slots * sizeof(INT), human );
    s << " in grid\n";
    s << pad << "      ";
    s.put_mem( num_huge_slots * sizeof(INT), human );
    s << " in huge hole list\n";
  }
}

// ******************************************************************

template <class INT>
void MEDDLY::array_plus_grid<INT>
::dumpInternal(output &s) const
{
  s << "Internal storage for array_plus_grid:\n";
  hole_manager<INT>::showInternal(s);

  s << "\n  Non-empty medium hole lists:\n";
  for (INT i=0; i<LargeHoleSize; i++) {
    if (0==medium_hole_list[i]) continue;
    s << "    [size " << i << "]";
    INT prev = 0;
    for (INT curr = medium_hole_list[i]; curr; curr=Next(curr)) {
      s << "  <-->  ";
      s << "#" << curr; 
      if (getHoleSize(curr) != i) {
        s << " (size = " << getHoleSize(curr) << " != " << i << " !!!) ";
      }
      if (prev != Prev(curr)) {
        s << " (prev = " << Prev(curr) << " != " << prev << "!!!) ";
      }
      prev = curr;
    }
    s << "\n";
  }

  s << "\n  Large hole grid:\n";
  if (0==grid_bottom) {
    s << "  Empty Grid\n";
  } else {
    s << "  Grid:\n";
    s << "  bottom\n";
    s << "    |\n";
    s << "    v\n";
    
    INT down = 0;
    for (INT index=grid_bottom; index; index=Up(index)) {
      if (grid_current == index) {
        s << "curr-> ";
      }
      s << "  #" << index << "  (size ";
      s << getHoleSize(index) << ")";

      if (down != Down(index)) {
        s << " (down = " << Down(index) << " != " << down << "!!!) ";
      }
      down = index;

      INT prev = index;
      for (INT ch=Next(index); ch; ch=Next(ch)) {
        s << "  <-->  #" << ch;
        if (getHoleSize(ch) != getHoleSize(index)) {
          s << "  (size " << getHoleSize(ch) << " != " << getHoleSize(index) << "!!!)";
        }
        if (prev != Prev(ch)) {
          s << " (prev = " << Prev(ch) << " != " << prev << "!!!) ";
        }
        prev = ch;
      }

      s << "\n";
      if (Up(index)) {
        s << "    ^\n";
        s << "    |\n";
        s << "    v\n";
      }

    } // for index

    s << "    ^\n";
    s << "    |\n";
    s << "  top\n";
    if (grid_top != down) {
      s << "  (top = " << grid_top << " != " << down << "!!!)\n";
    }
  } // non-empty grid

  s << "\n  Huge hole list:\n";
  if (0==huge_holes) {
    s << "empty\n";
  } else {
    INT prev = 0;
    for (INT index = huge_holes; index; index=Next(index)) {
      if (index != huge_holes) s << "  <-->  ";
      s << "#" << index << " (size ";
      s << getHoleSize(index) << ")";

      if (prev != Prev(index)) {
        s << " (prev = " << Prev(index) << " != " << prev << "!!!) ";
      }
      prev = index;
    }
    s << "\n";
  }

  s << "\n  max request: " << long(max_request) << "\n";
}


// ******************************************************************

template <class INT>
void MEDDLY::array_plus_grid<INT>
::dumpInternalUnused(output &s, node_address addr) const
{
  hole_manager<INT>::showInternalAddr(s, addr, 3);
}

// ******************************************************************

template <class INT>
int MEDDLY::array_plus_grid<INT>
::moveCurrentToRow(INT size, INT &current) const
{
  if (0==current) {
    current = grid_bottom;
  }
  if (0==current) return -1;

  if (getHoleSize(current) == size) {
    return 0;
  }

  if (getHoleSize(current) < size) {
    // go to larger holes until not less than size
    for (;;) {
      INT up = Up(current);
      if (0==up) return -1; // we're still smaller
      current = up;
      INT delta = getHoleSize(current) - size;
      if (delta >= 0) return delta;
    }
  }

  MEDDLY_DCASSERT(getHoleSize(current) > size);
  // go to smaller holes until not greater than size
  for (;;) {
    INT down = Down(current);
    if (0==down) return +1; // we're still larger
    current = down;
    INT delta = getHoleSize(current) - size;
    if (delta <= 0) return delta;
  }
  
  // should never get here.
  MEDDLY_DCASSERT(0);
  return -1;
}


// ******************************************************************


template <class INT>
void MEDDLY::array_plus_grid<INT>
::stopTrackingHole(node_address h) 
{
#ifdef MEMORY_TRACE_DETAILS
  printf("stopTrackingHole(%lu)\n", h);
#endif

  MEDDLY_DCASSERT(getHoleSize(h)>0);
  MEDDLY_DCASSERT(matchingHoleSizes(h));

  //
  // Is this a small hole?  If so, we don't track them,
  // so there is nothing to do except update stats.
  //
  if (isSmallHole(h)) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\thole size %ld, too small to track\n", long(getHoleSize(h)));
#endif
    // update stats
    num_small_holes--;
    num_small_slots -= getHoleSize(h);
    return;
  }

  //
  // Is this a medium hole?  If so, remove it from
  // the list containing it and update stats.
  //
  if (!isLargeHole(h)) {
    INT size = getHoleSize(h);
#ifdef MEMORY_TRACE_DETAILS
    printf("\tremoving from medium hole list [size=%ld]\n", long(size));
#endif
    INT left = Prev(h);
    INT right = Next(h);

    if (left) {
      MEDDLY_DCASSERT(h == node_address(Next(left)));
      setNext(left, right);
    } else {
      MEDDLY_DCASSERT(node_address(medium_hole_list[size]) == h);
      medium_hole_list[size] = right;
    }

    if (right) {
      MEDDLY_DCASSERT(h == node_address(Prev(right)));
      setPrev(right, left);
    }

    // update stats
    num_medium_holes[0]--;
    num_medium_holes[size]--;
    num_medium_slots -= size;
    return;
  }

  //
  // Is this a huge hole?  If so, remove from the huge hole list
  //
  if (size_t(getHoleSize(h)) > max_request) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\tremoving from huge hole list %ld\n", long(huge_holes));
#endif
    INT left = Prev(h);
    INT right = Next(h);

    if (left) {
      MEDDLY_DCASSERT(h == node_address(Next(left)));
      setNext(left, right);
    } else {
      MEDDLY_DCASSERT(node_address(huge_holes) == h);
      huge_holes = right;
    }

    if (right) {
      MEDDLY_DCASSERT(h == node_address(Prev(right)));
      setPrev(right, left);
    }

    // update stats
    num_huge_holes--;
    num_huge_slots -= getHoleSize(h);
    return;
  }

  //
  // Must be a large hole, in the grid.
  // Update stats, then split out the various cases.
  //

  // update stats
  num_grid_holes--;
  num_grid_slots -= getHoleSize(h);

  // 
  // Easier case - not an index hole
  //
  if (!isIndexHole(h)) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\tremoving non-index hole from grid\n");
#endif
    INT left = Prev(h);
    INT right = Next(h);

    MEDDLY_DCASSERT(left);

    MEDDLY_DCASSERT(h == node_address(Next(left)));
    setNext(left, right);

    if (right) {
      MEDDLY_DCASSERT(h == node_address(Prev(right)));
      setPrev(right, left);
    }
    return;
  }


#ifdef MEMORY_TRACE_DETAILS
  printf("\tremoving index hole from grid\n");
#endif

  //
  // We're removing an index node.
  // Two main cases to consider:
  //   (1) index node with empty chain: 
  //          the entire row is killed
  //   (2) index node with non-empty chain: 
  //          the front of chain becomes a new index node
  //

  if (0==Next(h)) {
    //
    // Empty chain.  Remove this "row" completely from the grid.
    //
#ifdef MEMORY_TRACE_DETAILS
    printf("\tremoving entire row\n");
#endif

    INT above = Up(h);
    INT below = Down(h);
    if (node_address(grid_current) == h) {
      if (above)  grid_current = above;
      else        grid_current = below;
    }

    if (below) {
      MEDDLY_DCASSERT(node_address(Up(below)) == h);
      setUp(below, above);
    } else {
      MEDDLY_DCASSERT(node_address(grid_bottom) == h);
      grid_bottom = above;
    }

    if (above) {
      MEDDLY_DCASSERT(node_address(Down(above)) == h);
      setDown(above, below);
    } else {
      MEDDLY_DCASSERT(node_address(grid_top) == h);
      grid_top = below;
    }

    return; 
  }

  // 
  // Non-empty chain.
  // Shift the front over to become the new index node
  //
#ifdef MEMORY_TRACE_DETAILS
  printf("\tshifting row\n");
#endif

  INT above = Up(h);
  INT below = Down(h);
  INT next = Next(h);
  MEDDLY_DCASSERT(next);

  if (node_address(grid_current) == h) {
    grid_current = next;
  }
  setPrev(next, 0);
  setUp(next, above);
  setDown(next, below);

  if (below) {
    MEDDLY_DCASSERT(node_address(Up(below)) == h);
    setUp(below, next);
  } else {
    MEDDLY_DCASSERT(node_address(grid_bottom) == h);
    grid_bottom = next;
  }

  if (above) {
    MEDDLY_DCASSERT(node_address(Down(above)) == h);
    setDown(above, next);
  } else {
    MEDDLY_DCASSERT(node_address(grid_top) == h);
    grid_top = next;
  }
}


// ******************************************************************

template <class INT>
void MEDDLY::array_plus_grid<INT>
::startTrackingHole(node_address h) 
{
#ifdef MEMORY_TRACE_DETAILS
  printf("startTrackingHole(%lu) size %ld\n", h, long(getHoleSize(h)));
#endif

  MEDDLY_DCASSERT(getHoleSize(h)>0);
  MEDDLY_DCASSERT(matchingHoleSizes(h));
    
  //
  // Check if the hole is too small to track.
  // If so, just leave it, and hope it gets merged later.
  //
  if (isSmallHole(h)) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\thole size %ld, too small to track\n", long(getHoleSize(h)));
#endif
    // update stats
    num_small_holes++;
    num_small_slots += getHoleSize(h);
    return;
  }

  //
  // Is this a medium hole?  If so, add it to the appropriate
  // list and update stats.
  //
  if (!isLargeHole(h)) {
    INT size = getHoleSize(h);
#ifdef MEMORY_TRACE_DETAILS
    printf("\tadding medium hole to list [size=%ld]\n", long(size));
#endif
    setNext(h, medium_hole_list[size]);
    setPrev(h, 0);

    if (medium_hole_list[size]) {
      setPrev(medium_hole_list[size], h);
    }
    medium_hole_list[size] = h;

    // update stats
    num_medium_holes[0]++;
    num_medium_holes[size]++;
    num_medium_slots += size;
    return;
  }

  //
  // Is this a huge hole?  If so, add it to the huge hole list
  //
  if (size_t(getHoleSize(h)) > max_request) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\tadding to huge hole list %ld\n", long(huge_holes));
#endif
    setNonIndex(h);
    setPrev(h, 0);
    setNext(h, huge_holes);
    if (huge_holes) {
      setPrev(huge_holes, h);
    }
    huge_holes = h;

    // update stats
    num_huge_holes++;
    num_huge_slots += getHoleSize(h);
    return;
  }

  //
  // Must be a large hole; add it to the grid
  //
  MEDDLY_DCASSERT(isLargeHole(h));
  // update stats
  num_grid_holes++;
  num_grid_slots += getHoleSize(h);

  //
  // Special case: empty grid
  //
  if (0 == grid_bottom) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\tadding to empty grid\n");
#endif
    setPrev(h, 0);
    setNext(h, 0);
    setUp(h, 0);
    setDown(h, 0);
    grid_bottom = grid_top = h;
    return;
  }

  //
  // Special case: at top
  // We check this for speed, because the general case
  // is to search from the bottom up.
  //
  if (getHoleSize(h) > getHoleSize(grid_top)) {
#ifdef MEMORY_TRACE_DETAILS
    printf("\tadding new chain at top\n");
#endif
    setPrev(h, 0);
    setNext(h, 0);
    setUp(h, 0);  // nothing larger
    setDown(h, grid_top);
    if (grid_top) setUp(grid_top, h);
    grid_top = h;
    return;
  }

  //
  // Add to grid, general case
  //
  INT curr = grid_bottom;
  int status = moveCurrentToRow(getHoleSize(h), curr);
  MEDDLY_DCASSERT(curr);
  if (0==status) {
    //
    // There's a chain of our size, add to it
    //
#ifdef MEMORY_TRACE_DETAILS
    printf("\tadding to chain\n");
#endif
    INT right = Next(curr);
    setNonIndex(h);
    setPrev(h, curr);
    setNext(h, right);
    if (right) {
      setPrev(right, h);
    }
    setNext(curr, h);
    return;
  }

  //
  // Nothing in our size, determine above and below rows
  //
  INT above, below;
  if (status<0) {
    below = curr;
    above = Up(curr); 
  } else {
    above = curr;
    below = Down(curr);
  }

  //
  // Make a new chain, here
  //
#ifdef MEMORY_TRACE_DETAILS
  printf("\tadding new chain\n");
#endif
  setPrev(h, 0);
  setNext(h, 0);
  setUp(h, above);
  setDown(h, below);
  if (above) {
    setDown(above, h);
  } else {
    grid_top = h;
  }
  if (below) {
    setUp(below, h);
  } else {
    grid_bottom = h;
  }
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    array_grid_style methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::array_grid_style::array_grid_style(const char* n)
: memory_manager_style(n)
{
}

MEDDLY::array_grid_style::~array_grid_style()
{
}

MEDDLY::memory_manager*
MEDDLY::array_grid_style::initManager(unsigned char granularity, 
  unsigned char minsize, forest::statset &stats) const
{
  if (sizeof(int) == granularity) {
    return new array_plus_grid <int>(getName(), stats);
  }

  if (sizeof(long) == granularity) {
    return new array_plus_grid <long>(getName(), stats);
  }

  if (sizeof(short) == granularity) {
    return new array_plus_grid <short>(getName(), stats);
  }

  // unsupported granularity

  return 0;
}

