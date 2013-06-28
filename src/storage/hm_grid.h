
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


#ifndef HM_GRID_H
#define HM_GRID_H

#include "holeman.h"

namespace MEDDLY {
  class hm_grid;
};

/** Grid-based hole management.
    Scheme from early versions of this library.

    When nodes are deleted, the memory slots are marked as a "hole",
    using the following format.

          slot[0] : -numslots, the number of slots in the hole
            .
            .
            .
          slot[L] : -numslots, with L = numslots-1
        
    The first few slots of the hole are used for a hole management
    data structure, described below.
    Note that a hole is guaranteed to be at least 5 slots long
    (assuming a node of size 0, with no extra header info, is impossible).


    Hole management.
    ==============================================================
    There are two kinds of holes depending on their location in the grid:
    Index Holes and Non-index Holes.
    Rows in the grid correspond to holes of the same size.
    The left-most column of the grid is a (vertical) list,
    and these are the index holes.
    Nodes in the middle are not connected vertically,
    and these are the non-index holes.

    The hole grid structure:
    ------------------------
    (holes_bottom)
    holes_of_size_0 (index) -- (non_index) -- (non_index) -- NULL
    |
    holes_of_size_1 -- ""
    |
    :
    :
    (holes_top)

    TBD:
    Note that the grid only stores holes up to the largest hole
    requested.  Larger holes are stored in the "large holes list".

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
class MEDDLY::hm_grid : public holeman {
  public:
    hm_grid(node_storage* p);
    virtual ~hm_grid();

  public:
    virtual node_address requestChunk(int slots);
    virtual void recycleChunk(node_address addr, int slots);
    virtual node_address chunkAfterHole(node_address addr) const;
    virtual void dumpInternalInfo(FILE* s) const;
    virtual void dumpHole(FILE* s, node_address a) const;
    virtual void reportMemoryUsage(FILE* s, const char* pad, int vL) const;
    virtual void clearHolesAndShrink(node_address new_last, bool shrink);

  private:
      inline node_handle* holeOf(node_handle addr) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_CHECK_RANGE(1, addr, last_slot+1);
        MEDDLY_DCASSERT(data[addr] < 0);  // it's a hole
        return data + addr;
      }
      inline node_handle& h_up(node_handle off) const {
        return holeOf(off)[hole_up_index];
      }
      inline node_handle& h_down(node_handle off) const {
        return holeOf(off)[hole_down_index];
      }
      inline node_handle& h_prev(node_handle off) const {
        return holeOf(off)[hole_prev_index];
      }
      inline node_handle& h_next(node_handle off) const {
        return holeOf(off)[hole_next_index];
      }

      inline node_handle& Up(node_handle off)       { return h_up(off); }
      inline node_handle  Up(node_handle off) const { return h_up(off); }

      inline node_handle& Down(node_handle off)       { return h_down(off); }
      inline node_handle  Down(node_handle off) const { return h_down(off); }

      inline node_handle& Prev(node_handle off)       { return h_prev(off); }
      inline node_handle  Prev(node_handle off) const { return h_prev(off); }

      inline node_handle& Next(node_handle off)       { return h_next(off); }
      inline node_handle  Next(node_handle off) const { return h_next(off); }

      inline bool isHoleNonIndex(node_handle p_offset) const {
          return (non_index_hole == h_up(p_offset));
      }

      // add a hole to the hole grid
      void gridInsert(node_handle p_offset);

      // remove a non-index hole from the hole grid
      void midRemove(node_handle p_offset);

      // remove an index hole from the hole grid
      void indexRemove(node_handle p_offset);

      // resize the data array.
      void resize(node_handle new_slots);


  private:
      static const int min_size = 1024;
      static const int non_index_hole = -2;   // any negative will do
      // hole indexes
      static const int hole_up_index = 1;
      static const int hole_down_index = hole_up_index+1;
      static const int hole_prev_index = hole_down_index;
      static const int hole_next_index = hole_prev_index+1;

  private:
      /// data array; we manage this
      node_handle* data;
      /// Size of data array.
      node_handle size;

      /// Largest hole ever requested
      node_handle max_request;
      /// List of large holes
      node_handle large_holes;
      /// Pointer to top of holes grid
      node_handle holes_top;
      /// Pointer to bottom of holes grid
      node_handle holes_bottom;
};

#endif

