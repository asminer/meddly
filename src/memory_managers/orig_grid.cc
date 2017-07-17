
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
// #define DONT_RECYCLE_FROM_GRID
// #define DONT_RECYCLE_FROM_LARGE

#include "hole_base.h"
#include "orig_grid.h"

namespace MEDDLY {

  // ******************************************************************
  // *                                                                *
  // *                 original_grid (template) class                 *
  // *                                                                *
  // ******************************************************************

  /**
    Grid structure for hole management.
    INT should be a signed integer type.

    Index hole:
    -------------------------------------------------------------
    [0] size with MSB set
    [1] up, >= 0, goes to larger holes
    [2] down, >= 0, goes to smaller holes
    [3] next pointer (list of holes of the same size)
    :
    : unused
    :
    [size-1] size with MSB set


    Non-index hole:
    -------------------------------------------------------------
    [0] size with MSB set
    [1] magic value < 0 to indicate non-index hole
    [2] prev pointer
    [3] next pointer
    :
    : unused
    :
    [size-1] size with MSB set


    The grid structure:
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
  class original_grid : public hole_manager<INT> {

    public:
      original_grid(const char* n, forest::statset &stats);
      virtual ~original_grid();

      virtual node_address requestChunk(size_t &numSlots);
      virtual void recycleChunk(node_address h, size_t numSlots);

      virtual void reportStats(output &s, const char* pad, bool human, bool details) const;
      virtual void dumpInternal(output &s) const;
      virtual void dumpInternalUnused(output &s, node_address addr) const;

    private:
      /**
          Remove a hole, at handle h, from the grid.
      */
      void removeFromGrid(node_address h);

      /**
          Add a hole, at handle h, to the grid.
      */
      void addToGrid(node_address h);

    private:
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

      inline bool isIndexHole(node_address h) const {
        return hole_manager<INT>::readSlot(h, 1) >= 0;
      }

      inline void setNonIndexHole(node_address h) {
        hole_manager<INT>::refSlot(h, 1) = -1;
      }

      inline void createIndexHoleUp(node_address h, INT up) {
        hole_manager<INT>::refSlot(h, 1) = up;
      }

      inline INT Up(node_address h) const {
        MEDDLY_DCASSERT(isIndexHole(h));
        return hole_manager<INT>::readSlot(h, 1);
      }
      inline INT& Up(node_address h) {
        MEDDLY_DCASSERT(isIndexHole(h));
        return hole_manager<INT>::refSlot(h, 1);
      }
      
      inline INT Down(node_address h) const {
        MEDDLY_DCASSERT(isIndexHole(h));
        return hole_manager<INT>::readSlot(h, 2);
      }
      inline INT& Down(node_address h) {
        MEDDLY_DCASSERT(isIndexHole(h));
        return hole_manager<INT>::refSlot(h, 2);
      }
      

      inline INT Prev(node_address h) const {
        MEDDLY_DCASSERT(!isIndexHole(h));
        return hole_manager<INT>::readSlot(h, 2);
      }
      inline INT& Prev(node_address h) {
        MEDDLY_DCASSERT(!isIndexHole(h));
        return hole_manager<INT>::refSlot(h, 2);
      }

      inline INT Next(node_address h) const {
        return hole_manager<INT>::readSlot(h, 3);
      }
      inline INT& Next(node_address h) {
        return hole_manager<INT>::refSlot(h, 3);
      }
      
      inline static INT smallestChunk() {
        return 5;
      }

    private:
      /// Largest hole ever requested
      size_t max_request;

      /// Pointer to list of large holes
      INT large_holes;

      /// Pointer to top of holes grid
      INT holes_bottom;

      /// Pointer to bottom of holes grid
      INT holes_top;

      /// Pointer where we left off the last search
      INT holes_current;

    private:    // stats
      long num_small_holes;
      long num_grid_holes;
      long num_large_holes;

      long num_small_slots;
      long num_grid_slots;
      long num_large_slots;

  }; // class original_grid

  // ******************************************************************
  // *                                                                *
  // *                     original_grid  methods                     *
  // *                                                                *
  // ******************************************************************

  template <class INT>
  original_grid<INT>::original_grid(const char* n, forest::statset &stats)
  : hole_manager<INT>(n, stats)
  {
    max_request = 0;
    large_holes = 0;
    holes_bottom = 0;
    holes_top = 0;
    holes_current = 0;

    num_small_holes = 0;
    num_grid_holes = 0;
    num_large_holes = 0;

    num_small_slots = 0;
    num_grid_slots = 0;
    num_large_slots = 0;
  }

  template <class INT>
  original_grid<INT>::~original_grid()
  {
  }


  // ******************************************************************


  template <class INT>
  MEDDLY::node_address original_grid<INT>::requestChunk(size_t &numSlots)
  {
#ifdef MEMORY_TRACE
    printf("requestChunk(%lu)\n", numSlots);
#endif

    if (numSlots > max_request) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tNew max request %lu; shifting holes from large list to grid\n", numSlots);
#endif
      max_request = numSlots;

      //
      // Traverse the large hole list, and re-insert everything.
      // Holes that are still large will go back into the large hole list,
      // and the rest will go into the grid where they belong
      //

      INT curr = large_holes;
      large_holes = 0;
      for (; curr; ) {
        INT next = Next(curr);
        addToGrid(curr);  
        curr = next;
      }
#ifdef DEBUG_GRID
      FILE_output out(stdout);
      dumpInternal(out);
#endif
    }

    //
    // Check a few places for holes.
    // Update h with the hole we're going to use.
    //
    node_address h = 0;

#ifndef DONT_RECYCLE_FROM_GRID
    //
    // Search for an existing chunk with exactly this size
    // 
    if (0==holes_current) {
      // Start at the bottom
      holes_current = holes_bottom;
    }

    //
    // Don't bother searching if grid is empty
    //
    if (holes_current) {

      while (getHoleSize(holes_current) > numSlots) {
        INT next = Down(holes_current);
        if (0==next) break;
        holes_current = next;
      }

      while (getHoleSize(holes_current) < numSlots) {
        INT next = Up(holes_current);
        if (0==next) break;
        holes_current = next;
      }

      //
      // holes_current either points to a row of the requested size,
      // or no such row exists.
      //
      if (getHoleSize(holes_current) >= numSlots) {
        //
        // We have a hole we can use.
        //

#ifdef MEMORY_TRACE_DETAILS
        printf("\tUsing hole of ");
        if (getHoleSize(holes_current) == numSlots) printf("requested ");
        printf("size %ld\n", long(getHoleSize(holes_current)));
#endif
       
        //
        // Remove the next hole because it's not an index hole
        // (that's easier), if we can
        //
        h = Next(holes_current);
        if (0==h) h = holes_current;
      } // if hole is large enough
    } // if holes_current
#endif  // #ifndef DONT_RECYCLE_FROM_GRID


#ifndef DONT_RECYCLE_FROM_LARGE
    if (0==h && large_holes) {
      //
      // Couldn't recycle from grid.  Try large hole list.
      //

      MEDDLY_DCASSERT(getHoleSize(large_holes) >= numSlots);

#ifdef MEMORY_TRACE_DETAILS
      printf("\tUsing large hole of size %ld\n", long(getHoleSize(large_holes)));
#endif

      h = large_holes;
    }
#endif // #ifndef DONT_RECYCLE_FROM_LARGE


    if (h) {
      //
      // We found a hole to use.  Remove it from the grid.
      //
      removeFromGrid(h);
      memory_manager::incMemUsed(getHoleSize(h) * sizeof(INT));

      //
      // Deal with leftover slots, if any
      //
      MEDDLY_DCASSERT(getHoleSize(h) >= numSlots);
      size_t leftover_slots = getHoleSize(h) - numSlots;
      if (leftover_slots > 0) {
        // 
        // Note - we even recycle holes of size 1,
        // because they can be merged to the left or right,
        // depending on who is recycled first
        //
        hole_manager<INT>::clearHole(h, numSlots); // otherwise we'll just merge them again
        recycleChunk(h+numSlots, leftover_slots);
      }

      // TBD - update stats
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
  void original_grid<INT>::recycleChunk(node_address h, size_t numSlots)
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
      MEDDLY_DCASSERT(getHoleSize(h-1) < h);
      node_address hleft = h - getHoleSize(h-1);
#ifdef MEMORY_TRACE_DETAILS
      printf("\tMerging to the left, holes %lu and %lu\n", hleft, h);
#endif
      removeFromGrid(hleft);
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
      removeFromGrid(hright);
      numSlots += getHoleSize(hright);
      setHoleSize(h, numSlots);
#ifdef DEBUG_GRID
      FILE_output out(stdout);
      dumpInternal(out);
#endif
    }

    //
    // Hole is ready, add to grid
    //
    addToGrid(h); 
#ifdef DEBUG_GRID
    FILE_output out(stdout);
    dumpInternal(out);
#endif
  }


  // ******************************************************************


  template <class INT>
  void original_grid<INT>::reportStats(output &s, const char* pad, bool human, bool details) const
  {
    s << pad << "Report for original_grid memory manager:\n";
    s << pad << "  Current #holes: ";
    s << num_small_holes+num_grid_holes+num_large_holes << "\n";
    if (details) {
      s << pad << "      " << num_small_holes << " small (untracked)\n";
      s << pad << "      " << num_grid_holes << " in grid\n";
      s << pad << "      " << num_large_holes << " in large hole list\n";
    }
    s << pad << "  Current bytes in holes: ";
    s.put_mem( (num_small_slots+num_grid_slots+num_large_slots) * sizeof(INT), human );
    s << "\n";
    if (details) {
      s << pad << "      ";
      s.put_mem( num_small_slots * sizeof(INT), human );
      s << " small (untracked)\n";
      s << pad << "      ";
      s.put_mem( num_grid_slots * sizeof(INT), human );
      s << " in grid\n";
      s << pad << "      ";
      s.put_mem( num_large_slots * sizeof(INT), human );
      s << " in large hole list\n";
    }
  }


  // ******************************************************************


  template <class INT>
  void original_grid<INT>::dumpInternal(output &s) const
  {
    s << "Internal storage for original_grid:\n";
    hole_manager<INT>::showInternal(s);

    if (0==holes_bottom) {
      s << "  Empty Grid\n";
    } else {

      s << "  Grid:\n";
      s << "  bottom\n";
      s << "    |\n";
      s << "    v\n";
    
      INT down = 0;
      for (INT index=holes_bottom; index; index=Up(index)) {
        if (holes_current == index) {
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
      if (holes_top != down) {
        s << "  (top = " << holes_top << " != " << down << "!!!)\n";
      }
    }

    s << "\n";
    s << "  max request: " << long(max_request) << "\n";
    s << "  large hole list:\n    ";
    
    if (0==large_holes) {
      s << "empty\n";
    } else {

      INT prev = 0;
      for (INT index = large_holes; index; index=Next(index)) {

          if (index != large_holes) s << "  <-->  ";
          s << "#" << index << " (size ";
          s << getHoleSize(index) << ")";

          if (prev != Prev(index)) {
            s << " (prev = " << Prev(index) << " != " << prev << "!!!) ";
          }
          prev = index;
      }
      s << "\n";
    }

  }


  // ******************************************************************


  template <class INT>
  void original_grid<INT>::dumpInternalUnused(output &s, node_address addr) const
  {
    hole_manager<INT>::showInternalAddr(s, addr, 3);
  }


  // ******************************************************************


  template <class INT>
  void original_grid<INT>::removeFromGrid(node_address h) 
  {
#ifdef MEMORY_TRACE_DETAILS
    printf("removeFromGrid(%lu)\n", h);
#endif

    MEDDLY_DCASSERT(getHoleSize(h)>0);
    MEDDLY_DCASSERT(matchingHoleSizes(h));

    //
    // Is this a small hole?  If so, it's not in the grid
    //
    if (getHoleSize(h) < smallestChunk()) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\thole size %ld, too small to track\n", long(getHoleSize(h)));
#endif
      // update stats
      num_small_holes--;
      num_small_slots -= getHoleSize(h);
      return;
    }

    //
    // Is this a large hole?  If so, remove from large hole list
    //
    if (getHoleSize(h) > max_request) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tremoving from large hole list %ld\n", long(large_holes));
#endif
      INT left = Prev(h);
      INT right = Next(h);

      if (left) {
        Next(left) = right;
      } else {
        MEDDLY_DCASSERT(large_holes == h);
        large_holes = right;
      }

      if (right) {
        Prev(right) = left;
      }

      // update stats
      num_large_holes--;
      num_large_slots -= getHoleSize(h);
      return;
    }

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

      Next(left) = right;
      if (right) {
        Prev(right) = left;
      }

      return;
    }

#ifdef MEMORY_TRACE_DETAILS
      printf("\tremoving index hole from grid\n");
#endif

    //
    // We're removing an index node.
    // Two main cases to consider:
    //   (1) index node with empty chain
    //   (2) index node with non-empty chain
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

      if (below) {
        Up(below) = above;
      } else {
        holes_bottom = above;
      }

      if (above) {
        Down(above) = below;
      } else {
        holes_top = below;
      }

      if (holes_current == h) {
        holes_current = below;
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

    if (above) {
      Down(above) = next;
    } else {
      holes_top = next;
    }
    createIndexHoleUp(next, above);

    if (below) {
      Up(below) = next;
    } else {
      holes_bottom = next;
    }
    Down(next) = below;

    if (holes_current == h) {
      holes_current = next;
    }
  }


  // ******************************************************************


  template <class INT>
  void original_grid<INT>::addToGrid(node_address h) 
  {
#ifdef MEMORY_TRACE_DETAILS
    printf("addToGrid(%lu) size %ld\n", h, long(getHoleSize(h)));
#endif

    MEDDLY_DCASSERT(getHoleSize(h)>0);
    MEDDLY_DCASSERT(matchingHoleSizes(h));
    
    //
    // Check if the hole is too small to track.
    // If so, just leave it, and hope it gets merged later.
    //
    if (getHoleSize(h) < smallestChunk()) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\thole size %ld, too small to track\n", long(getHoleSize(h)));
#endif
      // update stats
      num_small_holes++;
      num_small_slots += getHoleSize(h);
      return;
    }

    //
    // If the hole is large enough, just add it to the large hole list
    //
    if (getHoleSize(h) > max_request) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tadding to large hole list %ld\n", long(large_holes));
#endif
      setNonIndexHole(h);
      Prev(h) = 0;
      Next(h) = large_holes;
      if (large_holes) Prev(large_holes) = h;
      large_holes = h;
      // update stats
      num_large_holes++;
      num_large_slots += getHoleSize(h);
      return;
    }

    // update stats
    num_grid_holes++;
    num_grid_slots += getHoleSize(h);
    
    //
    // Special case: empty grid
    //
    if (0 == holes_bottom) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tadding to empty grid\n");
#endif
      createIndexHoleUp(h, 0);
      Down(h) = 0;
      Next(h) = 0;
      holes_bottom = holes_top = h;
      return;
      // TBD - update stats
    }

    //
    // Special case: at top
    //
    if (getHoleSize(h) > getHoleSize(holes_top)) {
#ifdef MEMORY_TRACE_DETAILS
      printf("\tadding new chain at top\n");
#endif
      createIndexHoleUp(h, 0);
      Down(h) = holes_top;
      Next(h) = 0;
      if (holes_top) Up(holes_top) = h;
      holes_top = h;
      return;
      // TBD - update stats
    }

    //
    // Add to the grid now.
    //

    //
    // First step: find our vertical position in the grid.
    // 
    INT above = holes_bottom;
    INT below = 0;
    while (getHoleSize(h) > getHoleSize(above)) {
      below = above;
      above = Up(below);
      MEDDLY_DCASSERT(Down(above) == below);
      MEDDLY_DCASSERT(above);
    }
    //
    // Found where we belong, check if it's the exact size
    //
    if (getHoleSize(h) == getHoleSize(above)) {
      //
      // There's a chain of our size, add to it
      //
#ifdef MEMORY_TRACE_DETAILS
      printf("\tadding to chain\n");
#endif
      INT right = Next(above);
      setNonIndexHole(h);
      Prev(h) = above;
      Next(h) = right;
      if (right) Prev(right) = h;
      Next(above) = h;
      return;
      // TBD - update stats
    }

    //
    // Make a new chain, here
    //
#ifdef MEMORY_TRACE_DETAILS
    printf("\tadding new chain\n");
#endif
    createIndexHoleUp(h, above);
    Down(h) = below;
    Next(h) = 0;
    if (above) {
      Down(above) = h;
    } else {
      holes_top = h;
    }
    if (below) {
      Up(below) = h;
    } else {
      holes_bottom = h;
    }

    // TBD - update stats
  }

}; // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    orig_grid_style  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::orig_grid_style::orig_grid_style(const char* n)
: memory_manager_style(n)
{
}

MEDDLY::orig_grid_style::~orig_grid_style()
{
}

MEDDLY::memory_manager*
MEDDLY::orig_grid_style::initManager(unsigned char granularity, 
  unsigned char minsize, forest::statset &stats) const
{
  if (sizeof(int) == granularity) {
    return new original_grid <int>(getName(), stats);
  }

  if (sizeof(long) == granularity) {
    return new original_grid <long>(getName(), stats);
  }

  if (sizeof(short) == granularity) {
    return new original_grid <short>(getName(), stats);
  }

  // unsupported granularity

  return 0;
}

