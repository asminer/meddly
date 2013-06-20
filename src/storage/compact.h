
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
// *                     node_compacted helpers                     *
// *                                                                *
// *                                                                *
// ******************************************************************

template <int bytes, class INT>
inline void valToData(INT a, unsigned char* b)
{
  //
  // The compiler is smart enough to optimize out these if's
  //
  b[0] = a & 0xff;
  if (bytes > 1)  { a >>= 8;  b[1] = a & 0xff;  }
  if (bytes > 2)  { a >>= 8;  b[2] = a & 0xff;  }
  if (bytes > 3)  { a >>= 8;  b[3] = a & 0xff;  }
  if (bytes > 4)  { a >>= 8;  b[4] = a & 0xff;  }
  if (bytes > 5)  { a >>= 8;  b[5] = a & 0xff;  }
  if (bytes > 6)  { a >>= 8;  b[6] = a & 0xff;  }
  if (bytes > 7)  { a >>= 8;  b[7] = a & 0xff;  }
}

template <int bytes, class INT>
inline void dataToVal(const unsigned char* b, INT &a)
{
  //
  // The compiler is smart enough to optimize out these if's
  //
  if (bytes > 7)  { a |= b[7];  a <<= 8;  }
  if (bytes > 6)  { a |= b[6];  a <<= 8;  }
  if (bytes > 5)  { a |= b[5];  a <<= 8;  }
  if (bytes > 4)  { a |= b[4];  a <<= 8;  }
  if (bytes > 3)  { a |= b[3];  a <<= 8;  }
  if (bytes > 2)  { a |= b[2];  a <<= 8;  }
  if (bytes > 1)  { a |= b[1];  a <<= 8;  }
  a |= b[0];
}

template <int bytes>
inline void longToData(long L, unsigned char* d)
{
  valToData<bytes>(L, d);
}

template <int bytes>
inline void dataToLong(const unsigned char* d, long& L)
{
  // deal with negatives properly
  if (d[bytes-1] & 0x80) {
    L = (~0L) << 8;
  } else {
    L = 0;
  }
  dataToVal<bytes>(d, L);
}

template <int bytes>
inline void intToData(int L, unsigned char* d)
{
  valToData<bytes>(L, d);
}

template <int bytes>
inline void dataToInt(const unsigned char* d, int& L)
{
  // deal with negatives properly
  if (d[bytes-1] & 0x80) {
    L = (~0) << 8;
  } else {
    L = 0;
  }
  dataToVal<bytes>(d, L);
}

template <int bytes>
inline void downToData(long P, unsigned char* d)
{
  // positive P: as usual.
  if (P >= 0) {
    valToData<bytes>(P, d);
    return;
  }

  // negative P: this is a terminal pointer.
  //              next msb set - terminal value is negative.
  //
  //  No conversion necessary because msb propogates when we shift
  static const unsigned long nmsb = (0x40L) << ((sizeof(long)-1)*8);
  if (P & nmsb) {
    valToData<bytes>(P, d);
    return;
  }

  // negative P: this is a terminal pointer.
  //              next msb clr - terminal value is positive.
  //              
  //  The thing to do here is deal with the msb manually:
  //  clear msb, encode, set msb.
  static const unsigned long msboff = ~ ((0x80L) << ((sizeof(long)-1)*8));
  valToData<bytes>(P & msboff, d);
  d[bytes-1] |= 0x80;
}

template <int bytes>
inline void dataToDown(const unsigned char* d, long& P)
{
  // Is this a terminal value?
  if (d[bytes-1] & 0x80) {
    // YES.
    // Is this a negative terminal value?
    if (d[bytes-1] & 0x40) {
      // YES.
      // Easy case: same as ordinary negatives.
      P = (~0L) << 8;
      dataToVal<bytes>(d, P);
      return;
    }
    // NO.
    // Positive terminal value.
    P = 0;
    dataToVal<bytes>(d, P);
    if (bytes != 8) {
      // Move MSB
      static const unsigned long bmsboff = ~ ((0x80L) << ((bytes-1)*8));
      static const unsigned long msbon = (0x80L) << ((sizeof(long)-1)*8);
      P = (P & bmsboff) | msbon;
    }
    return;
  }

  // non-terminal value: as usual
  P = 0;
  dataToVal<bytes>(d, P);
}

template <class INT>
inline void bytesRequiredForVal(INT a, int& bytes)
{
  if (sizeof(INT) == bytes) return;
  long la = (a<0) ? -(a<<1) : a<<1; 
  static unsigned long msbyte = (0xffL) << (sizeof(INT)-1)*8;
  unsigned long mask = msbyte;
  for (int b=sizeof(INT); b>bytes; b--) {
    if (la & mask) {
      bytes = b;
      return;
    }
    mask = (mask >> 8) & ~msbyte;
  }
}

inline void bytesRequiredForDown(long a, int& bytes)
{
  if (a<0) {
    // terminal value
    a <<= 1;
    if (a<0) {
      a <<= 1;
      a = -a;
    } else {
      a <<= 1;
    }
  } else {
    a <<= 1;
  }
  static unsigned long msbyte = (0xffL) << (sizeof(long)-1)*8;
  unsigned long mask = msbyte;
  for (int b=sizeof(long); b>bytes; b--) {
    if (a & mask) {
      bytes = b;
      return;
    }
    mask = (mask >> 8) & ~msbyte;
  }
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_compacted class                      *
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

