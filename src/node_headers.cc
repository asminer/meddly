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


#include "defines.h"
#include "node_headers.h"
#include "io.h"
#include "forest.h"
#include "storage/bytepack.h"
#include <climits>

// #define DEBUG_COUNTER_RESIZE

// #define DEBUG_HANDLE_MGT

// #define DEBUG_HANDLE_FREELIST

// #define ALWAYS_CHECK_FREELISTS

// #define DEBUG_ADDRESS_RESIZE

// #define DEBUG_SWEEP
// #define DEBUG_SWEEP_DETAIL
// #define DEBUG_SWEEP_DETAIL_CT

// ******************************************************************
// *                                                                *
// *                        Helper functions                        *
// *                                                                *
// ******************************************************************


// ******************************************************************

#ifdef DEBUG_HANDLE_FREELIST
void print_sequence(long a)
{
  static bool printed;
  static long first = 0;
  static long last = -1;
  if (a<=0) {
    if (first>last) return;
    if (printed) printf(", ");
    printed = false;
    if (first < last) {
      if (first+1<last) {
        printf("%ld ... %ld", first, last);
      } else {
        printf("%ld, %ld", first, last);
      }
    } else {
      printf("%ld", first);
    }
    first = 0;
    last = -1;
    return;
  }
  // a > 0
  if (0==first) {
    first = last = a;
    return;
  }
  if (last+1 == a) {
    last++;
    return;
  }
  // break in the sequence, we need to print
  if (printed) printf(", "); else printed = true;
  if (first < last) {
    if (first+1<last) {
      printf("%ld ... %ld", first, last);
    } else {
      printf("%ld, %ld", first, last);
    }
  } else {
    printf("%ld", first);
  }
  first = last = a;
}

// ******************************************************************

inline void dump_handle_info(const MEDDLY::node_headers &NH, long size)
{
  printf("Used handles:  ");
  print_sequence(0);
  for (long i=1; i<size; i++) if (!NH.isDeleted(i)) {
    print_sequence(i);
  }
  print_sequence(0);
  printf("\nFree handles:  ");
  for (long i=1; i<size; i++) if (NH.isDeleted(i)) {
    print_sequence(i);
  }
  print_sequence(0);
  printf("\n");
}
#endif

// ******************************************************************
// *                                                                *
// *               node_headers::level_array  methods               *
// *                                                                *
// ******************************************************************

MEDDLY::node_headers::level_array::level_array(node_headers &p, int max_level)
 : parent(p)
{
  data8 = 0;
  data16 = 0;
  data32 = 0;
  size = 0;

  if (max_level < 128) {
    bytes = 1;
    MEDDLY_DCASSERT(1 == sizeof(char));
  } else if (max_level < 32768) {
    bytes = 2;
    MEDDLY_DCASSERT(2 == sizeof(short));
  } else {
    bytes = 4;
    MEDDLY_DCASSERT(4 == sizeof(int));
  }

  parent.changeHeaderSize(0, bytes*8);
}

MEDDLY::node_headers::level_array::~level_array()
{
  free(data8);
  free(data16);
  free(data32);
  parent.changeHeaderSize(bytes*8, 0);
}

void MEDDLY::node_headers::level_array::expand(size_t ns)
{
  if (ns <= size) return;

  char* d8 = 0;
  short* d16 = 0;
  int* d32 = 0;

  switch (bytes) {
    case 1:   // char array
              d8 = (char*) realloc(data8, ns * bytes);
              if (0==d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              memset(d8 + size, 0, (ns-size) * bytes );
              data8 = d8;
              break;

    case 2:   // short array
              d16 = (short*) realloc(data16, ns * bytes);
              if (0==d16) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              memset(d16 + size, 0, (ns-size) * bytes );
              data16 = d16;
              break;

    case 4:   // int array
              d32 = (int*) realloc(data32, ns * bytes);
              if (0==d32) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              memset(d32 + size, 0, (ns-size) * bytes );
              data32 = d32;
              break;

    default:
              MEDDLY_DCASSERT(0);
  }

  size = ns;
#ifdef DEBUG_LEVEL_RESIZE
  printf("Enlarged level array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::level_array::shrink(size_t ns)
{
  if (ns >= size) return;

  char* d8 = 0;
  short* d16 = 0;
  int* d32 = 0;

  switch (bytes) {
    case 1:   // char array
              d8 = (char*) realloc(data8, ns * bytes);
              if (0==d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              data8 = d8;
              break;

    case 2:   // short array
              d16 = (short*) realloc(data16, ns * bytes);
              if (0==d16) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              data16 = d16;
              break;

    case 4:   // int array
              d32 = (int*) realloc(data32, ns * bytes);
              if (0==d32) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              data32 = d32;
              break;

    default:
              MEDDLY_DCASSERT(0);
  }

  size = ns;
#ifdef DEBUG_LEVEL_RESIZE
  printf("Reduced level array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::level_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put(long(get(first)), width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put(long(get(i)), width);
  }
  s << "]";
}

// ******************************************************************
// *                                                                *
// *              node_headers::counter_array  methods              *
// *                                                                *
// ******************************************************************

MEDDLY::node_headers::counter_array::counter_array(node_headers &p) : parent(p)
{
  data8 = 0;
  data16 = 0;
  data32 = 0;
  size = 0;
  counts_09bit = 0;
  counts_17bit = 0;
  bytes = sizeof(unsigned char);
  parent.changeHeaderSize(0, bytes*8);
}

MEDDLY::node_headers::counter_array::~counter_array()
{
  free(data8);
  free(data16);
  free(data32);
  parent.changeHeaderSize(bytes*8, 0);
}

void MEDDLY::node_headers::counter_array::expand(size_t ns)
{
  if (ns <= size) return;

  unsigned char* d8 = 0;
  unsigned short* d16 = 0;
  unsigned int* d32 = 0;

  switch (bytes) {
    case 1:   // char array
              MEDDLY_DCASSERT(0==data16);
              MEDDLY_DCASSERT(0==data32);
              d8 = (unsigned char*) realloc(data8, ns * bytes);
              if (0==d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              memset(d8 + size, 0, (ns-size) * bytes );
              data8 = d8;
              break;

    case 2:   // short array
              if (counts_09bit) {
                d16 = (unsigned short*) realloc(data16, ns * bytes);
                if (0==d16) {
                  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                memset(d16 + size, 0, (ns-size) * bytes );
                data16 = d16;
                MEDDLY_DCASSERT(0==data8);
                MEDDLY_DCASSERT(0==data32);
              } else {
                shrink16to8(ns);
                memset(data8 + size, 0, (ns-size) * bytes );
              }
              break;

    case 4:   // int array
              if (counts_17bit) {
                d32 = (unsigned int*) realloc(data32, ns * bytes);
                if (0==d32) {
                  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                memset(d32 + size, 0, (ns-size) * bytes );
                data32 = d32;
                MEDDLY_DCASSERT(0==data8);
                MEDDLY_DCASSERT(0==data16);
              } else if (counts_09bit) {
                shrink32to16(ns);
                memset(data16 + size, 0, (ns-size) * bytes );
              } else {
                shrink32to8(ns);
                memset(data8 + size, 0, (ns-size) * bytes );
              }
              break;

    default:
              MEDDLY_DCASSERT(0);
  }; // switch bytes

  size = ns;
#ifdef DEBUG_COUNTER_RESIZE
  printf("Enlarged counter array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::counter_array::shrink(size_t ns)
{
  if (ns >= size) return;

  unsigned char* d8 = 0;
  unsigned short* d16 = 0;
  unsigned int* d32 = 0;

  switch (bytes) {
    case 1:   // char array
              MEDDLY_DCASSERT(0==counts_09bit);
              MEDDLY_DCASSERT(0==counts_17bit);
              d8 = (unsigned char*) realloc(data8, ns * bytes);
              if (0==d8) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
              }
              data8 = d8;
              break;

    case 2:   // short array
              MEDDLY_DCASSERT(0==counts_17bit);
              if (counts_09bit) {
                d16 = (unsigned short*) realloc(data16, ns * bytes);
                if (0==d16) {
                  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                data16 = d16;
              } else {
                shrink16to8(ns);
              }
              break;

    case 4:   // int array
              if (counts_17bit) {
                d32 = (unsigned int*) realloc(data32, ns * bytes);
                if (0==d32) {
                  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
                }
                data32 = d32;
              } else if (counts_09bit) {
                shrink32to16(ns);
              } else {
                shrink32to8(ns);
              }
              break;

    default:
              MEDDLY_DCASSERT(0);
  }; // switch bytes

  size = ns;
#ifdef DEBUG_COUNTER_RESIZE
  printf("Reduced counter array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::counter_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put((unsigned long) get(first), width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put((unsigned long) get(i), width);
  }
  s << "]";
}

void MEDDLY::node_headers::counter_array::expand8to16(size_t j)
{
#ifdef DEBUG_COUNTER_RESIZE
  printf("Expanding counter array from 8 to 16 bits (for element %lu)\n", j);
#endif
  MEDDLY_DCASSERT(1==bytes);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(0==data32);
  MEDDLY_DCASSERT(0==counts_09bit);
  MEDDLY_DCASSERT(0==counts_17bit);

  counts_09bit = 1;
  bytes = sizeof(unsigned short);
  data16 = (unsigned short*) malloc(size * bytes);
  if (0==data16) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  for (size_t i=0; i<size; i++) {
    data16[i] = data8[i];
  }
  data16[j] = 256;
  free(data8);
  data8 = 0;
  parent.changeHeaderSize(8, 16);
}

void MEDDLY::node_headers::counter_array::expand16to32(size_t j)
{
#ifdef DEBUG_COUNTER_RESIZE
  printf("Expanding counter array from 16 to 32 bits (for element %lu)\n", j);
#endif
  MEDDLY_DCASSERT(2==bytes);
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data32);
  MEDDLY_DCASSERT(0==counts_17bit);

  counts_17bit = 1;
  bytes = sizeof(unsigned int);
  data32 = (unsigned int*) malloc(size * bytes);
  if (0==data32) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  for (size_t i=0; i<size; i++) {
    data32[i] = data16[i];
  }
  data32[j] = 65536;
  free(data16);
  data16 = 0;
  parent.changeHeaderSize(16, 32);
}

void MEDDLY::node_headers::counter_array::shrink16to8(size_t ns)
{
#ifdef DEBUG_COUNTER_RESIZE
  printf("Shrinking counter array from 16 to 8 bits\n");
#endif
  MEDDLY_DCASSERT(2==bytes);
  MEDDLY_DCASSERT(data16);
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data32);
  MEDDLY_DCASSERT(0==counts_09bit);
  MEDDLY_DCASSERT(0==counts_17bit);

  parent.changeHeaderSize(16, 8);
  bytes = sizeof(unsigned char);
  data8 = (unsigned char*) malloc(ns * bytes);
  if (0==data8) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size_t stop = (ns < size) ? ns : size;
  for (size_t i=0; i<stop; i++) {
    data8[i] = data16[i];
  }
  free(data16);
  data16 = 0;
}

void MEDDLY::node_headers::counter_array::shrink32to16(size_t ns)
{
#ifdef DEBUG_COUNTER_RESIZE
  printf("Shrinking counter array from 32 to 16 bits\n");
#endif
  MEDDLY_DCASSERT(4==bytes);
  MEDDLY_DCASSERT(data32);
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(0==counts_17bit);

  parent.changeHeaderSize(32, 16);
  bytes = sizeof(unsigned short);
  data16 = (unsigned short*) malloc(ns * bytes);
  if (0==data16) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size_t stop = (ns < size) ? ns : size;
  for (size_t i=0; i<stop; i++) {
    data16[i] = data32[i];
  }
  free(data32);
  data32 = 0;
}

void MEDDLY::node_headers::counter_array::shrink32to8(size_t ns)
{
#ifdef DEBUG_COUNTER_RESIZE
  printf("Shrinking counter array from 32 to 8 bits\n");
#endif
  MEDDLY_DCASSERT(4==bytes);
  MEDDLY_DCASSERT(data32);
  MEDDLY_DCASSERT(0==data8);
  MEDDLY_DCASSERT(0==data16);
  MEDDLY_DCASSERT(0==counts_09bit);
  MEDDLY_DCASSERT(0==counts_17bit);

  parent.changeHeaderSize(32, 8);
  bytes = sizeof(unsigned char);
  data8 = (unsigned char*) malloc(ns * bytes);
  if (0==data8) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size_t stop = (ns < size) ? ns : size;
  for (size_t i=0; i<stop; i++) {
    data8[i] = data32[i];
  }
  free(data32);
  data32 = 0;
}

// ******************************************************************
// *                                                                *
// *              node_headers::address_array  methods              *
// *                                                                *
// ******************************************************************

MEDDLY::node_headers::address_array::address_array(node_headers &p) : parent(p)
{
  data32 = 0;
  num_large_elements = 0;
  bytes = 4;
  data64 = 0;
  size = 0;
  parent.changeHeaderSize(0, bytes*8);
}

MEDDLY::node_headers::address_array::~address_array()
{
  free(data32);
  free(data64);
  parent.changeHeaderSize(bytes*8, 0);
}

void MEDDLY::node_headers::address_array::expand(size_t ns)
{
  if (ns <= size) return;

  if (4==bytes) {
    MEDDLY_DCASSERT(0==data64);
    unsigned int* d = (unsigned int*) realloc(data32, ns * bytes);
    if (0==d) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    memset(d + size, 0, (ns-size) * bytes);
    data32 = d;
  } else {
    MEDDLY_DCASSERT(0==data32);
    MEDDLY_DCASSERT(8==bytes);
    if (num_large_elements) {
      // Stay at 64 bits
      unsigned long* d = (unsigned long*) realloc(data64, ns * bytes);
      if (0==d) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }
      memset(d + size, 0, (ns-size) * sizeof(unsigned long) );
      data64 = d;
    } else {
      shrink64to32(ns);
    }
  }

  size = ns;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::address_array::shrink(size_t ns)
{
  if (ns >= size) return;

  if (4==bytes) {
    MEDDLY_DCASSERT(data32);
    MEDDLY_DCASSERT(0==data64);
    unsigned int* d = (unsigned int*) realloc(data32, ns * bytes);
    if (0==d) {
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    data32 = d;
  } else {
    MEDDLY_DCASSERT(8==bytes);
    MEDDLY_DCASSERT(0==data32);
    if (num_large_elements) {
      // Stay at 64 bits
      unsigned long* d = (unsigned long*) realloc(data64, ns * sizeof(unsigned long));
      if (0==d) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
      }
      data64 = d;
    } else {
      shrink64to32(ns);
    }
  }

  size = ns;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced address array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::address_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put(get(first), width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put(get(i), width);
  }
  s << "]";
}

void MEDDLY::node_headers::address_array
  ::expand32to64()
{
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Expanding address array from 32 to 64 bits\n");
#endif
  MEDDLY_DCASSERT(4==bytes);
  MEDDLY_DCASSERT(0==data64);
  MEDDLY_DCASSERT(0==num_large_elements);

  bytes = 8;
  MEDDLY_DCASSERT(sizeof(unsigned long) == bytes);
  data64 = (unsigned long*) malloc(size * bytes);
  if (0==data64) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  for (size_t i=0; i<size; i++) {
    data64[i] = data32[i];
  }
  free(data32);
  data32 = 0;
  num_large_elements = 1;
  parent.changeHeaderSize(32, 64);
}

void MEDDLY::node_headers::address_array::shrink64to32(size_t ns)
{
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrinking address array from 64 to 32 bits\n");
#endif
  MEDDLY_DCASSERT(8==bytes);
  MEDDLY_DCASSERT(data64);
  MEDDLY_DCASSERT(0==data32);
  MEDDLY_DCASSERT(0==num_large_elements);

  bytes = 4;
  MEDDLY_DCASSERT(sizeof(unsigned int) == bytes);
  data32 = (unsigned int*) malloc(ns * bytes);
  if (0==data32) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size_t stop = (ns < size) ? ns : size;
  for (size_t i=0; i<stop; i++) {
    data32[i] = data64[i];
  }
  free(data64);
  data64 = 0;
  parent.changeHeaderSize(64, 32);
}

// ******************************************************************
// *                                                                *
// *                node_headers::bitvector  methods                *
// *                                                                *
// ******************************************************************

MEDDLY::node_headers::bitvector::bitvector(node_headers &p) : parent(p)
{
  data = 0;
  size = 0;
  parent.changeHeaderSize(0, sizeof(bool)*8);
}

MEDDLY::node_headers::bitvector::~bitvector()
{
  free(data);
  parent.changeHeaderSize(sizeof(bool)*8, 0);
}

void MEDDLY::node_headers::bitvector::expand(size_t ns)
{
  if (ns <= size) return;

  bool* d = (bool*) realloc(data, ns * sizeof(bool));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(bool) );
  size = ns;
  data = d;
#ifdef DEBUG_BITVECTOR_RESIZE
  printf("Enlarged bitvector, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::bitvector::shrink(size_t ns)
{
  if (ns >= size) return;
  bool* d = (bool*) realloc(data, ns * sizeof(bool));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size = ns;
  data = d;
#ifdef DEBUG_BITVECTOR_RESIZE
  printf("Reduced bitvector, new size %lu\n", size);
#endif
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_headers methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

const size_t START_SIZE = 512;
// const size_t MAX_ADD = 65536;
const size_t MAX_ADD = 16777216;

inline size_t next_size(size_t s)
{
  if (0==s) return START_SIZE;
  if (s < MAX_ADD) return s*2;
  return s+MAX_ADD;
}

inline size_t prev_size(size_t s)
{
  if (s <= MAX_ADD) return s/2;
  return s-MAX_ADD;
}

inline size_t next_check(size_t s)
{
  if (0==s) return START_SIZE/4;
  if (s < MAX_ADD) return s*2;
  return s+MAX_ADD;
}

inline size_t prev_check(size_t s)
{
  if (s <= MAX_ADD) return s/2;
  return s-MAX_ADD;
}

MEDDLY::node_headers::node_headers(expert_forest &P)
  : parent(P)
{
  h_bits = 0;
  a_last = 0;
  a_size = 0;
  a_freed = 0;
  a_next_shrink = 0;
  for (int i=0; i<8; i++) a_unused[i] = 0;
  a_lowest_index = 8;
  a_sweep = SIZE_MAX;
  addresses = new address_array(*this);
  levels = new level_array(*this, parent.getNumVariables());
  if (parent.getPolicies().useReferenceCounts) {
    cache_counts = new counter_array(*this);
    is_in_cache = 0;
    incoming_counts = new counter_array(*this);
    is_reachable = 0;
  } else {
    cache_counts = 0;
    is_in_cache = new bitvector(*this);
    incoming_counts = 0;
    is_reachable = new bitvector(*this);
  }
  implicit_bits = 0;
}

// ******************************************************************

MEDDLY::node_headers::~node_headers()
{
  delete addresses;
  delete levels;
  delete cache_counts;
  delete is_in_cache;
  delete incoming_counts;
  delete is_reachable;
  delete implicit_bits;
}

// ******************************************************************

void MEDDLY::node_headers
::reportStats(output &s, const char* pad, unsigned flags) const
{
  if (flags & (expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED) ) {
    s << pad << "Node headers: array based\n";
  }
  if (flags & expert_forest::STORAGE_STATS) {
    s << pad << "  Node header: " << int(h_bits) << " bits\n";
  }

  if (flags & expert_forest::STORAGE_DETAILED) {
    if (addresses)        s << pad << "    address   : " << addresses->entry_bits() << " bits\n";
    if (levels)           s << pad << "    level     : " << levels->entry_bits() << " bits\n";
    if (cache_counts)     s << pad << "    #caches   : " << cache_counts->entry_bits() << " bits\n";
    if (is_in_cache)      s << pad << "    in_cache? : " << is_in_cache->entry_bits() << " bits\n";
    if (incoming_counts)  s << pad << "    #incoming : " << incoming_counts->entry_bits() << " bits\n";
    if (is_reachable)     s << pad << "    reachable?: " << is_reachable->entry_bits() << " bits\n";
    if (implicit_bits)    s << pad << "    implicit? : " << implicit_bits->entry_bits() << " bits\n";
  }

  if (flags & expert_forest::STORAGE_STATS) {
    s << pad << "  " << a_size << " headers allocated\n";
    s << pad << "  " << a_last << " last header\n";
    s << pad << "  " << a_freed << " recycled headers\n";
    s << pad << "  " << a_last - a_freed << " headers in use\n";
  }

}

// ******************************************************************

void MEDDLY::node_headers
::showHeader(output &s, node_handle p) const
{
    MEDDLY_DCASSERT(p>0);

  int k = ABS(getNodeLevel(p));
  const variable* v = parent.getDomain()->getVar(parent.getVarByLevel(k));
  if (v->getName()) {
    s << " level: " << v->getName();
  } else {
    s << " level: " <<  k;
  }
  s.put( (getNodeLevel(p) < 0) ? '\'' : ' ' );

  if (incoming_counts)  s << " in: " << (unsigned long)incoming_counts->get((size_t)p);
  if (cache_counts)     s << " cc: " << (unsigned long)cache_counts->get((size_t)p);
  if (is_reachable)     s << " reach: " << (unsigned long)is_reachable->get((size_t)p);
  if (is_in_cache)      s << " cache: " << (unsigned long)is_in_cache->get((size_t)p);
}

// ******************************************************************

/*
void MEDDLY::node_headers::turnOffCacheCounts()
{
  delete cache_counts;
  cache_counts = 0;
}
*/

// ******************************************************************

/*
void MEDDLY::node_headers::turnOffIncomingCounts()
{
  delete incoming_counts;
  incoming_counts = 0;
}
*/

// ******************************************************************

MEDDLY::node_handle MEDDLY::node_headers::getFreeNodeHandle()
{
  parent.mstats.incMemUsed(h_bits/8);
  node_handle found = 0;
  for (int i=a_lowest_index; i<8; i++) {
    // try the various lists
    while (a_unused[i] > a_last) {
      a_unused[i] = getNextOf(a_unused[i]);
    }
    if (a_unused[i]) {  // get a recycled one from list i
      found = a_unused[i];
      a_unused[i] = getNextOf(a_unused[i]);
      break;
    } else {
      if (i == a_lowest_index) a_lowest_index++;
    }
  }
  //
  // We got something from a free list.
  //
  if (found) {
#ifdef DEBUG_HANDLE_FREELIST
    // address[found].setNotDeleted();
    dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
    printf("Forest %u using recycled handle %d\n", parent.FID(), found);
#endif
    a_freed--;
    return found;
  }

  //
  // Free lists are empty.
  // Use sweep index if it's small enough (this is set
  // if we're using mark and sweep).
  //
  if (a_sweep < a_last) {
    MEDDLY_DCASSERT(addresses);
    MEDDLY_DCASSERT(is_in_cache);
    MEDDLY_DCASSERT(is_reachable);
    for (;;) {
      a_sweep = is_in_cache->firstZero(++a_sweep);
      if (a_sweep > a_last) break;
      if (is_reachable->get(a_sweep)) continue;
      found = a_sweep;
      break;
    }
    if (found) {
      MEDDLY_DCASSERT(0==addresses->get((size_t)found));
    } else {
      // we're done sweeping
      a_sweep = SIZE_MAX;
    }
  }

  if (found) {
#ifdef DEBUG_HANDLE_FREELIST
    // address[found].setNotDeleted();
    dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
    printf("Forest %u using swept handle %d\n", parent.FID(), found);
#endif
    return found;
  }

  if (a_last+1 >= a_size) {
    expandHandleList();
  }
  a_last++;
  MEDDLY_DCASSERT(a_last < a_size);

#ifdef DEBUG_HANDLE_FREELIST
  dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
  printf("Forest %u using end handle %lu\n", parent.FID(), size_t(a_last));
#endif
  return a_last;
}

// ******************************************************************

void MEDDLY::node_headers::recycleNodeHandle(node_handle p)
{
    MEDDLY_DCASSERT(p>0);
    const size_t pp = size_t(p);

  MEDDLY_DCASSERT(0==getNodeCacheCount(p));

#ifdef DEBUG_HANDLE_MGT
  printf("Forest %u recycling handle %d\n", parent.FID(), p);
#endif

  parent.mstats.decMemUsed(h_bits/8);
  a_freed++;

  // Determine which list to add this into,
  // and add it there
  int i = bytesRequiredForDown(p) -1;
  setNextOf(pp, a_unused[i]);
  a_unused[i] = pp;
  a_lowest_index = MIN(a_lowest_index, (char)i);

  // if this was the last node, collapse nodes into the
  // "not yet allocated" pile.  But, we don't remove them
  // from the free list(s); we simply discard any too-large
  // ones when we pull from the free list(s).
  if (pp == a_last) {
    while (a_last && isDeleted(a_last) && (0==getNodeCacheCount(a_last))) {
      a_last--;
      a_freed--;
    }

#ifdef DEBUG_HANDLE_MGT
    printf("Forest %u collapsing handle end from %d to %lu\n", parent.FID(), p, size_t(a_last));
#endif

    if (a_last < a_next_shrink) {
#ifdef DEBUG_HANDLE_MGT
      printf("Forest %u collapsed end less than shrink check %lu\n", parent.FID(), size_t(a_next_shrink));
#endif
      shrinkHandleList();
    }

  }
#ifdef DEBUG_HANDLE_FREELIST
  dump_handle_info(*this, a_last+1);
#endif
#ifdef ALWAYS_CHECK_FREELISTS
  validateFreeLists();
#endif
}

// ******************************************************************

void MEDDLY::node_headers::swapNodes(node_handle p, node_handle q, bool swap_incounts)
{
  MEDDLY_DCASSERT(p!=q);
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(q>0);

  const size_t sp = size_t(p);
  const size_t sq = size_t(q);
  if (addresses)                          addresses->swap(sp, sq);
  if (levels)                             levels->swap(sp, sq);
  if (cache_counts)                       cache_counts->swap(sp, sq);
  if (swap_incounts && incoming_counts)   incoming_counts->swap(sp, sq);
  if (implicit_bits)                      implicit_bits->swap(sp, sq);
}

// ******************************************************************

void MEDDLY::node_headers::sweepAllInCacheBits()
{
  a_sweep = 0;

#ifdef DEBUG_SWEEP
  printf("Forest %u finished cache scan\n", parent.FID());
#endif

#ifdef DEBUG_SWEEP_DETAIL
  printf("Done sweeping; nodes that will be reclaimed:\n\t");

  MEDDLY_DCASSERT(is_in_cache);
  MEDDLY_DCASSERT(is_reachable);
  size_t p = 0;
  for (;;) {
    p = is_in_cache->firstZero(++p);
    if (is_reachable->get(p)) continue;
    printf("%d, ", p);
  } // infinite loop
  printf(".\n");

#ifdef DEBUG_SWEEP_DETAIL_CT
  printf("Compute tables:\n");
  FILE_output s(stdout);
  operation::showAllComputeTables(s, 5);
#endif
#endif

  // TBD - scan backwards and collapse a_last, and then shrink array if needed
}


// ******************************************************************

void MEDDLY::node_headers::dumpInternal(output &s) const
{
  s << "Node headers and management:\n";
  for (int i=0; i<8; i++) {
    s << "    First " << i << "-byte unused node index: " << a_unused[i] << "\n";
  }
  int awidth = digits(a_last);
  s << "    Node# :  ";
  for (unsigned long p=1; p<=a_last; p++) {
    if (p>1) s.put(' ');
    s.put(p, awidth);
  }
  s << "\n    Level  : ";
  if (levels) {
    levels->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n    Offset : ";
  if (addresses) {
    addresses->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n   Incount : ";
  if (incoming_counts) {
    incoming_counts->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n    Cache  : ";
  if (cache_counts) {
    cache_counts->show(s, 1, a_last, awidth);
    s << "\n";
  } else {
    s << "(null)\n";
  }
  s << "\n";
}

// ******************************************************************

void MEDDLY::node_headers::validateFreeLists() const
{
  const int MAX_ERRORS = 25;
  int errors = 0;
  printf("  Validating Free Lists\n");
  bool* inlist = new bool[a_last+1];
  for (size_t i=0; i<=a_last; i++) {
    inlist[i] = false;
  }
  for (unsigned i=0; i<8; i++) {
    for (size_t curr = size_t(a_unused[i]); curr; curr = size_t(getNextOf(curr))) {
      if (curr > a_last) continue;
      inlist[curr] = true;
      if (isDeleted(curr)) continue;
      printf("\tBAD FREE LIST: item %lu in free list %u is not deleted?\n",
        curr, i
      );
      ++errors;
      if (errors >= MAX_ERRORS) {
        printf("\tToo many errors, not printing the rest\n");
        exit(1);
      }
    }
  }
  for (size_t i=1; i<=a_last; i++) {
    if (!isDeleted(i)) continue;
    if (inlist[i]) continue;
    printf("\tMISSING: item %lu is deleted, not in any free list\n", i);
    ++errors;
    if (errors >= MAX_ERRORS) {
      printf("\tToo many errors, not printing the rest\n");
      exit(1);
    }
  }
  printf("  Done validating free lists\n");
  delete[] inlist;
  if (errors>0) exit(1);
}

// ******************************************************************

void MEDDLY::node_headers::changeHeaderSize(unsigned oldbits, unsigned newbits)
{
  if (oldbits < newbits) {
    parent.mstats.incMemUsed((a_last - a_freed)*(newbits-oldbits)/8);
    parent.mstats.incMemAlloc(a_size*(newbits-oldbits)/8);
  } else {
    parent.mstats.decMemUsed((a_last - a_freed)*(oldbits-newbits)/8);
    parent.mstats.decMemAlloc(a_size*(oldbits-newbits)/8);
  }
  h_bits += newbits;
  h_bits -= oldbits;
}

// ******************************************************************

void MEDDLY::node_headers::expandHandleList()
{
  //
  // If we're not using reference counts,
  // determine which nodes are reachable.
  //
  if (!parent.getPolicies().useReferenceCounts) {

#ifdef DEBUG_SWEEP
      printf("Forest %u starting unreachable scan\n", parent.FID());
#endif

      parent.markAllRoots();

#ifdef DEBUG_SWEEP
      printf("Forest %u finished unreachable scan\n", parent.FID());
#endif
#ifdef DEBUG_SWEEP_DETAIL
      printf("Unreachable nodes:\n\t");

      for (node_handle p=1; p<=a_last; p++) {
        MEDDLY_DCASSERT(is_reachable);
        if (is_reachable->get(p)) continue;
        printf("%d, ", p);
      }
      printf(".\n");

#endif  // DEBUG_SWEEP_DETAIL

      //
      // Because there's no way for us to recover
      // an unreachable node, go ahead and delete them now.
      //
      MEDDLY_DCASSERT(is_reachable);
      size_t unrch = 0;
      for (;;) {
          unrch = is_reachable->firstZero(++unrch);
          if (unrch >= a_size) break;
          if (!isDeleted(unrch)) parent.deleteNode(unrch);
      }

  } // no incoming counts

  //
  // Done with GC, now work on expanding
  //

  size_t old_size = a_size;
  do {
    a_next_shrink = next_check(a_next_shrink);
    a_size = next_size(a_size);
  } while (a_last+1 >= a_size);

  parent.mstats.incMemAlloc( ((a_size-old_size)*h_bits)/8 );

  if (addresses)        addresses->expand(a_size);
  if (levels)           levels->expand(a_size);
  if (cache_counts)     cache_counts->expand(a_size);
  if (is_in_cache)      is_in_cache->expand(a_size);
  if (incoming_counts)  incoming_counts->expand(a_size);
  if (is_reachable)     is_reachable->expand(a_size);
  if (implicit_bits)    implicit_bits->expand(a_size);

#ifdef DEBUG_ADDRESS_RESIZE
  printf("Expanded node headers arrays, new size %lu (shrink check %lu)\n", a_size, a_next_shrink);
#endif
}


// ******************************************************************

void MEDDLY::node_headers::shrinkHandleList()
{
#ifdef DEBUG_ADDRESS_RESIZE
  printf("About to shrink node headers arrays (current size %ld)...\n",
    a_size);
  validateFreeLists();
#endif

  // clean out free lists, because we're about
  // to shrink the address list which holds the pointers
  for (int i=0; i<8; i++) {
    //
    // clean list i
    //
    size_t prev = 0;
    size_t curr;
    for (curr = a_unused[i]; curr; curr=getNextOf(curr))
    {
      if (curr > a_last) continue;  // don't add to the list
      if (prev) {
        setNextOf(prev, curr);
      } else {
        a_unused[i] = curr;
      }
      prev = curr;
    }
    if (prev) {
      setNextOf(prev, 0);
    } else {
      a_unused[i] = 0;
    }
  } // for i

  size_t old_size = a_size;
  do {
    a_size = prev_size(a_size);
    a_next_shrink = (a_size > START_SIZE) ? prev_check(a_next_shrink) : 0;
  } while (a_last < a_next_shrink);

  parent.mstats.decMemAlloc( ((old_size-a_size)*h_bits)/8 );

  if (addresses)        addresses->shrink(a_size);
  if (levels)           levels->shrink(a_size);
  if (cache_counts)     cache_counts->shrink(a_size);
  if (is_in_cache)      is_in_cache->shrink(a_size);
  if (incoming_counts)  incoming_counts->shrink(a_size);
  if (is_reachable)     is_reachable->shrink(a_size);
  if (implicit_bits)    implicit_bits->shrink(a_size);

#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank node headers arrays, new size %lu (shrink check %lu)\n", a_size, a_next_shrink);
  printf("Done shrinking node headers arrays (size is now %ld)...\n", a_size);
  validateFreeLists();
#endif
}


// ******************************************************************

void MEDDLY::node_headers::lastUnlink(node_handle p)
{
    //
    // Node p is unreachable.  See if we need to do some cleanup
    //

    if (
        (cache_counts && (0==cache_counts->get(size_t(p))))
        ||
        (is_in_cache && (0==is_in_cache->get(size_t(p))))
    ) {
        //
        // Unreachable and not in any caches.  Delete and recycle handle.
        //
        parent.deleteNode(p);
        recycleNodeHandle(p);
        return;
    }

    if (pessimistic) {
        //
        // Delete; keep handle until caches are cleared.
        //
        parent.deleteNode(p);
    } else {
        //
        // Optimistic.  Keep unreachables around.
        //
#ifdef TRACK_UNREACHABLE_NODES
        parent.stats.unreachable_nodes++;
#endif
    }
}

// ******************************************************************

void MEDDLY::node_headers::reviveNode(node_handle p)
{
    // Reclaim an unreachable
    parent.stats.reclaimed_nodes++;
#ifdef TRACK_UNREACHABLE_NODES
    parent.stats.unreachable_nodes--;
#endif
}

// ******************************************************************

void MEDDLY::node_headers::lastUncache(node_handle p)
{
    //
    // Cache count for p is now zero; might need to clean up
    //

    if (isDeleted(p)) {
        //
        // Must be using pessimistic
        //
#ifdef TRACK_UNREACHABLE_NODES
        parent.stats.unreachable_nodes--;
#endif
        recycleNodeHandle(p);
    } else {
        // we were active.
        // See if we're now completely disconnected
        // and if so, tell parent to recycle node storage
        if (
          (incoming_counts && (0==incoming_counts->get(size_t(p))))
          ||
          (is_reachable && (0==is_reachable->get(size_t(p))))
        ) {
#ifdef TRACK_UNREACHABLE_NODES
            parent.stats.unreachable_nodes--;
#endif
            parent.deleteNode(p);
            recycleNodeHandle(p);
        }
    }
}

// ******************************************************************
