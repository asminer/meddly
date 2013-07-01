
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


#ifndef HM_ARRAY_H
#define HM_ARRAY_H

#include "holeman.h"

namespace MEDDLY {
  class hm_array;
};

#define MEASURE_LARGE_HOLE_STATS

/** Array-based hole management.

    Small holes are stored by size in an array of lists.
    Large holes are stored in a single list.
    Lists are doubly-linked so we may remove an arbitrary node.

    Holes are represented as follows:
    ---------------------------------------
    [0] -size (number of slots in hole)     
    [1] prev pointer
    [2] next pointer
    :
    :
    [size-1] -size

*/
class MEDDLY::hm_array : public holeman {
  public:
    hm_array(node_storage* p);
    virtual ~hm_array();

  public:
    virtual node_address requestChunk(int slots);
    virtual void recycleChunk(node_address addr, int slots);
    virtual node_address chunkAfterHole(node_address addr) const;
    virtual void dumpInternalInfo(FILE* s) const;
    virtual void dumpHole(FILE* s, node_address a) const;
    virtual void dumpInternalTail(FILE* s) const;
    virtual void reportMemoryUsage(FILE* s, const char* pad, int vL) const;
    virtual void clearHolesAndShrink(node_address new_last, bool shrink);

  private:
      inline node_handle* holeOf(node_handle addr) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_CHECK_RANGE(1, addr, last_slot+1);
        MEDDLY_DCASSERT(data[addr] < 0);  // it's a hole
        return data + addr;
      }
      inline node_handle& h_prev(node_handle off) const {
        return holeOf(off)[hole_prev_index];
      }
      inline node_handle& h_next(node_handle off) const {
        return holeOf(off)[hole_next_index];
      }

      inline node_handle& Prev(node_handle off)       { return h_prev(off); }
      inline node_handle  Prev(node_handle off) const { return h_prev(off); }

      inline node_handle& Next(node_handle off)       { return h_next(off); }
      inline node_handle  Next(node_handle off) const { return h_next(off); }

      // resize the data array.
      void resize(node_handle new_slots);

      // Insert a node
      inline void listInsert(node_handle& list, node_handle node) {
        Next(node) = list;
        Prev(node) = 0;
        if (list) {
          Prev(list) = node;
        } 
        list = node; 
      }

      // Remove a node
      inline void listRemove(node_handle& list, node_handle node) {
        if (Prev(node)) {
          Next(Prev(node)) = Next(node);
        } else {
          list = Next(node);
        }
        if (Next(node)) {
          Prev(Next(node)) = Prev(node);
        }
      }

      // Determine list length
      inline long listLength(node_handle list) const {
        long count = 0;
        for (; list; list = Next(list)) count++;
        return count;
      }

  private:
      static const int min_size = 1024;
      // hole indexes
      static const int hole_prev_index = 1;
      static const int hole_next_index = 2;

      // when does a hole become "large"
      static const int LARGE_SIZE = 128;
      // static const int LARGE_SIZE = 1;  // extreme - everything is "large"

  private:
      /// data array; we manage this
      node_handle* data;
      /// Size of data array.
      node_handle size;

      /// List of large holes
      node_handle large_holes;
      /// Lists of holes by size
      node_handle small_holes[LARGE_SIZE];

#ifdef MEASURE_LARGE_HOLE_STATS
      long num_large_hole_traversals;
      long count_large_hole_visits;
#endif

      static node_handle verify_hole_slots;
};

#endif

