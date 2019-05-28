
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
#include "storage/bytepack.h"

// #define DEBUG_HANDLE_MGT

// #define DEBUG_HANDLE_FREELIST

// #define DEBUG_ADDRESS_RESIZE

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

#ifndef OLD_NODE_HEADERS

MEDDLY::node_headers::level_array::level_array(node_headers &p) : parent(p)
{
  data32 = 0;
  size = 0;
  parent.changeHeaderSize(0, sizeof(int)*8);
}

MEDDLY::node_headers::level_array::~level_array()
{
  free(data32);
  parent.changeHeaderSize(sizeof(int)*8, 0);
}

void MEDDLY::node_headers::level_array::expand(size_t ns)
{
  if (ns <= size) return;

  int* d = (int*) realloc(data32, ns * sizeof(int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(int) );
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged level array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::level_array::shrink(size_t ns)
{
  if (ns >= size) return;
  int* d = (int*) realloc(data32, ns * sizeof(int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced level array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::level_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put(long(data32[first]), width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put(long(data32[i]), width);
  }
  s << "]";
}

#endif // OLD_NODE_HEADERS

// ******************************************************************
// *                                                                *
// *              node_headers::counter_array  methods              *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

MEDDLY::node_headers::counter_array::counter_array(node_headers &p) : parent(p)
{
  data32 = 0;
  size = 0;
  parent.changeHeaderSize(0, sizeof(unsigned int)*8);
}

MEDDLY::node_headers::counter_array::~counter_array()
{
  free(data32);
  parent.changeHeaderSize(sizeof(unsigned int)*8, 0);
}

void MEDDLY::node_headers::counter_array::expand(size_t ns)
{
  if (ns <= size) return;

  unsigned int* d = (unsigned int*) realloc(data32, ns * sizeof(unsigned int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(unsigned int) );
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged counter array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::counter_array::shrink(size_t ns)
{
  if (ns >= size) return;
  unsigned int* d = (unsigned int*) realloc(data32, ns * sizeof(unsigned int));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size = ns;
  data32 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced counter array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::counter_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put((unsigned long) data32[first], width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put((unsigned long) data32[i], width);
  }
  s << "]";
}

#endif

// ******************************************************************
// *                                                                *
// *              node_headers::address_array  methods              *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

MEDDLY::node_headers::address_array::address_array(node_headers &p) : parent(p)
{
  data64 = 0;
  size = 0;
  parent.changeHeaderSize(0, sizeof(unsigned long)*8);
}

MEDDLY::node_headers::address_array::~address_array()
{
  free(data64);
  parent.changeHeaderSize(sizeof(unsigned long)*8, 0);
}

void MEDDLY::node_headers::address_array::expand(size_t ns)
{
  if (ns <= size) return;

  unsigned long* d = (unsigned long*) realloc(data64, ns * sizeof(unsigned long));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(d + size, 0, (ns-size) * sizeof(unsigned long) );
  size = ns;
  data64 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::address_array::shrink(size_t ns)
{
  if (ns >= size) return;
  unsigned long* d = (unsigned long*) realloc(data64, ns * sizeof(unsigned long));
  if (0==d) {
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  size = ns;
  data64 = d;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced address array, new size %lu\n", size);
#endif
}

void MEDDLY::node_headers::address_array
  ::show(output &s, size_t first, size_t last, int width) const
{
  s << "[";
  s.put(data64[first], width);
  for (size_t i=first+1; i<=last; i++) {
    s.put('|');
    s.put(data64[i], width);
  }
  s << "]";
}

#endif

// ******************************************************************
// *                                                                *
// *                node_headers::bitvector  methods                *
// *                                                                *
// ******************************************************************

#ifndef OLD_NODE_HEADERS

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
#ifdef DEBUG_ADDRESS_RESIZE
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
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Reduced bitvector, new size %lu\n", size);
#endif
}

#endif

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_headers methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

const size_t START_SIZE = 512;
const size_t MAX_ADD = 65536;

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
#ifdef OLD_NODE_HEADERS
  //
  // Inltialize address array
  //
  a_size = a_min_size;
  address = (node_header *) malloc(a_size * sizeof(node_header));
  if (0 == address) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  parent.mstats.incMemAlloc(a_size * sizeof(node_header));
  memset(address, 0, a_size * sizeof(node_header));
  a_last = a_next_shrink = a_freed = 0;
  for (int i=0; i<8; i++) a_unused[i] = 0;
  a_lowest_index = 8;

  //
  // hooks for later
  //
  usesCacheCounts = true;
  usesIncomingCounts = true;
#else
  h_bits = 0;
  a_last = 0;
  a_size = 0;
  a_freed = 0;
  a_next_shrink = 0;
  for (int i=0; i<8; i++) a_unused[i] = 0;
  a_lowest_index = 8;
  addresses = new address_array(*this);
  levels = new level_array(*this);
  cache_counts = new counter_array(*this);
  incoming_counts = new counter_array(*this);
  implicit_bits = 0;
#endif
}

// ******************************************************************

MEDDLY::node_headers::~node_headers()
{
#ifdef OLD_NODE_HEADERS
  parent.mstats.decMemAlloc(a_size * sizeof(node_header));

  // Address array
  free(address);
#else
  delete addresses;
  delete levels;
  delete cache_counts;
  delete incoming_counts;
  delete implicit_bits;
#endif
}

// ******************************************************************

void MEDDLY::node_headers
::reportStats(output &s, const char* pad, unsigned flags) const
{
  if (flags & (expert_forest::STORAGE_STATS | expert_forest::STORAGE_DETAILED) ) {
#ifdef OLD_NODE_HEADERS
    s << pad << "Node headers: struct based\n";
#else
    s << pad << "Node headers: array based\n";
#endif
  }
  if (flags & expert_forest::STORAGE_STATS) {
#ifdef OLD_NODE_HEADERS
    unsigned h_bits = sizeof(node_header) * 8;
#endif
    s << pad << "  Node header: " << int(h_bits) << " bits\n";
  }

#ifndef OLD_NODE_HEADERS
  if (flags & expert_forest::STORAGE_DETAILED) {
    if (addresses)        s << pad << "    address   : " << addresses->entry_bits() << " bits\n";
    if (levels)           s << pad << "    level     : " << levels->entry_bits() << " bits\n";
    if (cache_counts)     s << pad << "    #caches   : " << cache_counts->entry_bits() << " bits\n";
    if (incoming_counts)  s << pad << "    #incoming : " << incoming_counts->entry_bits() << " bits\n";
    if (implicit_bits)    s << pad << "    implicit? : " << implicit_bits->entry_bits() << " bits\n";
  }
#endif

  if (flags & expert_forest::STORAGE_STATS) {
    s << pad << "  " << a_size << " headers allocated\n";
    s << pad << "  " << a_last << " last header\n";
    s << pad << "  " << a_freed << " recycled headers\n";
    s << pad << "  " << a_last - a_freed << " headers in use\n";
  }

}

// ******************************************************************

void MEDDLY::node_headers::turnOffCacheCounts()
{
#ifdef OLD_NODE_HEADERS
  usesCacheCounts = false;
  // For now - nothing else to do
#else
  delete cache_counts;
  cache_counts = 0;
#endif
}

// ******************************************************************

void MEDDLY::node_headers::turnOffIncomingCounts()
{
#ifdef OLD_NODE_HEADERS
  usesIncomingCounts = false;
  // For now - nothing else to do
#else
  delete incoming_counts;
  incoming_counts = 0;
#endif
}

// ******************************************************************

MEDDLY::node_handle MEDDLY::node_headers::getFreeNodeHandle() 
{
#ifdef OLD_NODE_HEADERS
  parent.mstats.incMemUsed(sizeof(node_header));
#else
  parent.mstats.incMemUsed(h_bits/8);
#endif
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
  if (found) {  
#ifdef DEBUG_HANDLE_FREELIST
    // address[found].setNotDeleted();
    dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
    fprintf(stderr, "Using recycled handle %d\n", found);
#endif
    a_freed--;
    return found;
  }
  a_last++;
  if (a_last >= a_size) {
    expandHandleList();
  }
  MEDDLY_DCASSERT(a_last < a_size);

#ifdef DEBUG_HANDLE_FREELIST
  dump_handle_info(*this, a_last+1);
#endif
#ifdef DEBUG_HANDLE_MGT
  fprintf(stderr, "Using end handle %d\n", a_last);
#endif
  return a_last;
}

// ******************************************************************

void MEDDLY::node_headers::recycleNodeHandle(node_handle p) 
{
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(0==getNodeCacheCount(p));

#ifdef DEBUG_HANDLE_MGT
  fprintf(stderr, "Recycling handle %d\n", p);
#endif

#ifdef OLD_NODE_HEADERS
  parent.mstats.decMemUsed(sizeof(node_header));
#else
  parent.mstats.decMemUsed(h_bits/8);
#endif
  deactivate(p);
  a_freed++;

  // Determine which list to add this into,
  // and add it there
  int i = bytesRequiredForDown(p) -1;
#ifdef OLD_NODE_HEADERS
  setNextOf(p, a_unused[i]);
  a_unused[i] = p;
#else
  setNextOf(size_t(p), a_unused[i]);
  a_unused[i] = size_t(p);
#endif
  a_lowest_index = MIN(a_lowest_index, (char)i);

  // if this was the last node, collapse nodes into the
  // "not yet allocated" pile.  But, we don't remove them
  // from the free list(s); we simply discard any too-large
  // ones when we pull from the free list(s).
  if (p == a_last) {
    while (a_last && isDeactivated(a_last)) {
      a_last--;
      a_freed--;
    }

#ifdef DEBUG_HANDLE_MGT
    fprintf(stderr, "Collapsing handle end from %d to %lu\n", p, a_last);
#endif

    if (a_last < a_next_shrink) {
#ifdef DEBUG_HANDLE_MGT
      fprintf(stderr, "Collapsed end less than shrink check %lu\n", a_next_shrink);
#endif
      shrinkHandleList();
    }

  }
#ifdef DEBUG_HANDLE_FREELIST
  dump_handle_info(*this, a_last+1);
#endif
}

// ******************************************************************

void MEDDLY::node_headers::swapNodes(node_handle p, node_handle q, bool swap_incounts)
{
  MEDDLY_DCASSERT(p!=q);
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(q>0);
  MEDDLY_DCASSERT(q<=a_last);

#ifdef OLD_NODE_HEADERS
  MEDDLY_DCASSERT(address);
  SWAP(address[p], address[q]);

  if (!swap_incounts) {
    SWAP(address[p].incoming_count, address[q].incoming_count);
  } 
#else
  size_t sp = size_t(p);
  size_t sq = size_t(q);
  if (addresses)                          addresses->swap(sp, sq);
  if (levels)                             levels->swap(sp, sq);
  if (cache_counts)                       cache_counts->swap(sp, sq);
  if (swap_incounts && incoming_counts)   incoming_counts->swap(sp, sq);
  if (implicit_bits)                      implicit_bits->swap(sp, sq);
#endif
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
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put(' ');
    s.put(long(p), awidth);
  }
#ifdef OLD_NODE_HEADERS
  s << "\n    Level  : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].level), awidth);
  }
  s << "]\n";
  s << "\n    Offset : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].offset), awidth);
  }
  s << "]\n";
  s << "\n   Incount : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].incoming_count), awidth);
  }
  s << "]\n";
  s << "\n    Cache  : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].cache_count), awidth);
  }
  s << "]\n\n";
#else
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
#endif
}

// ******************************************************************

#ifndef OLD_NODE_HEADERS

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

#endif

// ******************************************************************

void MEDDLY::node_headers::expandHandleList()
{
#ifdef OLD_NODE_HEADERS
  // increase size by 50%
  int delta = a_size / 2;
  MEDDLY_DCASSERT(delta>=0);
  node_header* new_address = (node_header*) 
    realloc(address, (a_size+delta) * sizeof(node_header));
  if (0==new_address) {
    /*
    fprintf(stderr, "Error in allocating array of size %lu at %s, line %d\n",
        (a_size+delta) * sizeof(node_header), __FILE__, __LINE__);
    */
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  address = new_address;
  parent.mstats.incMemAlloc(delta * sizeof(node_header));
  memset(address + a_size, 0, delta * sizeof(node_header));
  a_size += delta;
  a_next_shrink = a_size / 2;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %d\n", a_size);
#endif

#else // OLD_NODE_HEADERS

  size_t old_size = a_size;
  do {
    a_next_shrink = next_check(a_next_shrink);
    a_size = next_size(a_size);
  } while (a_last >= a_size);

  if (addresses)        addresses->expand(a_size);
  if (levels)           levels->expand(a_size);
  if (cache_counts)     cache_counts->expand(a_size);
  if (incoming_counts)  incoming_counts->expand(a_size);
  if (implicit_bits)    implicit_bits->expand(a_size);

  parent.mstats.incMemAlloc( ((a_size-old_size)*h_bits)/8 );

#ifdef DEBUG_ADDRESS_RESIZE
  printf("Expanded node headers arrays, new size %lu (shrink check %lu)\n", a_size, a_next_shrink);
#endif
#endif
}


// ******************************************************************

void MEDDLY::node_headers::shrinkHandleList()
{
#ifdef OLD_NODE_HEADERS

  // Determine new size
  int new_size = a_min_size;
  while (a_last >= new_size) new_size += new_size/2;
  int delta = a_size - new_size;
  if (0==delta) {
    a_next_shrink = 0;
    return;
  }

  // clean out free lists, because we're about
  // to trash memory beyond new_size.
  for (int i=0; i<8; i++) {
    //
    // clean list i
    //
    node_handle prev = 0;
    node_handle curr;
    for (curr = a_unused[i]; curr; curr=getNextOf(curr)) 
    {
      if (curr > a_last) continue;  // don't add to the list
      if (prev) {
        setNextOf(prev, curr);
      } else {
        a_unused[i] = curr;
      }
    } 
    if (prev) {
      setNextOf(prev, 0);
    } else {
      a_unused[i] = 0;
    }
  } // for i

  // shrink the array
  MEDDLY_DCASSERT(delta>=0);
  MEDDLY_DCASSERT(a_size-delta>=a_min_size);
  node_header* new_address = (node_header*) 
    realloc(address, new_size * sizeof(node_header));
  if (0==new_address) {
    /*
    fprintf(stderr, "Error in allocating array of size %lu at %s, line %d\n",
        new_size*sizeof(node_header), __FILE__, __LINE__);
    */
    throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  address = new_address;
  parent.mstats.decMemAlloc(delta * sizeof(node_header));
  a_size -= delta;
  a_next_shrink = a_size / 2;
  MEDDLY_DCASSERT(a_last < a_size);
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank address array, new size %d\n", a_size);
#endif

#else // OLD_NODE_HEADERS

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

  if (addresses)        addresses->shrink(a_size);
  if (levels)           levels->shrink(a_size);
  if (cache_counts)     cache_counts->shrink(a_size);
  if (incoming_counts)  incoming_counts->shrink(a_size);
  if (implicit_bits)    implicit_bits->shrink(a_size);

  parent.mstats.incMemAlloc( ((old_size-a_size)*h_bits)/8 );

#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank node headers arrays, new size %lu (shrink check %lu)\n", a_size, a_next_shrink);
#endif

#endif
}


// ******************************************************************

