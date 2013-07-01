
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


#ifndef COMPACT_STORAGE_H
#define COMPACT_STORAGE_H

// TODO: Testing

#include "../defines.h"
#include "../hash_stream.h"


// ******************************************************************
// *                                                                *
// *                                                                *
// *                    new node_compacted class                    *
// *                                                                *
// *                                                                *
// ******************************************************************

/** New design for node storage in a forest.

    Limits: node size fits in an "int".

    Guarantees: nodes are "node_handle aligned" except for
                portions that are "compacted".

    ------------------------------------------------------------
    Format, "outer", viewed as an array of node_handles:
    ------------------------------------------------------------

        slot[0]:    incoming count, will be >= 0
        slot[1]:    next pointer in unique table, or special value

          ... memory chunk, described below, with padding
          at the end as necessary since the chunk is spread
          over an array of node_handles ...

        slot[N]:    tail (node handle number), will be >= 0

        
        Memory chunk:

          int   size  :   node size.  >=0 for truncated full, <0 for sparse
          byte  style :   how to compact the remaining arrays
                  3 bits - (size of pointers in bytes) - 1
                      000: 1 byte pointers
                      001: 2 byte pointers
                      ...
                      111: 8 byte pointers (longs on 64-bit arch.)

                  2 bits - (size of indexes in bytes) - 1
                      00: 1 byte indexes
                      ..
                      11: 4 byte indexes (ints; max node size anyway)

                  3 bits - unused, reserve for edge value compaction?

          unhashed header info: raw bytes
          hashed header info: raw bytes

          compacted downward pointer array (size * #bytes per pointer)
          compacted index array if sparse (size * #bytes per index)
          compacted edge value array (size * #bytes per edge)

          ignored bytes until we reach node_handle boundary


    ------------------------------------------------------------
    Format of "holes", viewed as an array of node_handles
    ------------------------------------------------------------

        slot[0]:    -size of hole,   <0

        slot[1]:  }
        slot[2]:  } used for hole management
        slot[3]:  }

        slot[N]:    -size of hole,   <0

    ------------------------------------------------------------
    
    Minimum chunk size is 5 slots.

    Alternative to the grid for holes: array of lists/heaps

      small   [4..255]:  holes of exactly this size
      medium  [4..255]:  slot i is for 64*i <= hole < 64*(i+1)
      large   [4..255]:  slot i is for 4096*i <= hole < 4096*(i+1)

      hugelist: holes of size 4096*256 = 1048576 and larger,
                not organized by size

*/

// class goes here

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    old node_compacted class                    *
// *                                                                *
// *                                                                *
// ******************************************************************

#if 0

/** Compacted node storage in a forest.
    TBD!  Experimental!
    Note this allows for 64-bit pointers and node counts.
    Implemented in node_wrappers.cc
*/
class MEDDLY::node_compacted {
      static const int min_size = 1024;

      /// Special values
      static const long non_index_hole = -2;
      static const long temp_node_value = -5;

      // header indexes
      static const int count_index = 0;
      static const int next_index = 1;    
      static const int size_index = 2;

      // Counts for extra slots
      static const int commonHeaderBytes = 2*sizeof(long) + sizeof(int) + 1;
      static const int commonTailBytes = sizeof(long);

      // Parent forest.
      expert_forest* parent;

      /// data array
      unsigned char* data;
      /// shifted data array, for node size
      unsigned char* data_size;
      /// shifted data array, for storage info
      unsigned char* data_info;
      /// shifted data array, for "node start"
      unsigned char* data_down;
      /// Size of data array.
      long size;
      /// Last used data slot.  Also total number of bytes "allocated"
      long last;

  // Holes grid info
  private:
      /// List of large holes
      long large_holes;
      /// Pointer to top of holes grid
      long holes_top;
      /// Pointer to bottom of holes grid
      long holes_bottom;
      /// Total bytes in holes
      long hole_slots;
      /// Largest hole ever requested
      int max_request;

  // header sizes; vary by forest.
  public:
      /// Size of each outgoing edge's value (can be 0).
      char evbytes;
      /// Size of extra unhashed data (typically 0).
      char uhbytes;
      /// Size of extra hashed data (typically 0).
      char hhbytes;

  // --------------------------------------------------------
  // |  Public interface.
  public:
      node_compacted();
      ~node_compacted();

      /// Bind this to a particular forest.
      void init(expert_forest* p);

      /// Compact this level.  (Rearrange, to remove all holes.)
      void compact(bool shrink);



      /// For debugging: dump everything.
      void dumpInternal(FILE* s) const;

      /// For debugging: dump entry starting at given slot.
      void dumpInternal(FILE* s, long a) const;

  // --------------------------------------------------------
  // |  inlined helpers.
  private:
      // get downward pointer bytes from a status value.
      static inline int dpbytesOf(unsigned char d) {
        return (d & 0xf0) >> 4;
      }
      // get index bytes from a status value.
      static inline int ixbytesOf(unsigned char d) {
        return (d & 0x0f);
      }
      // number of bytes used for downward pointers, for node a.
      inline int dpbytes(long a) const {
        MEDDLY_DCASSERT(data_info);
        return dpbytesOf(data_info[a]);
      }
      // number of bytes used for index pointers, for node a.
      // zero means the node is "full".
      inline int ixbytes(long a) const {
        MEDDLY_DCASSERT(data_info);
        return ixbytesOf(data_info[a]);
      }

  // --------------------------------------------------------
  // |  Misc. helpers.
  private:

      // resize the data array.
      void resize(long new_bytes);
};

#endif  // if 0

#endif  // include guard

