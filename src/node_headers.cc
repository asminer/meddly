
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

inline void dump_handle_info(const MEDDLY::node_header* A, long size)
{
  printf("Used handles:  ");
  print_sequence(0);
  for (long i=1; i<size; i++) if (!A[i].isDeleted()) {
    print_sequence(i);
  }
  print_sequence(0);
  printf("\nFree handles:  ");
  for (long i=1; i<size; i++) if (A[i].isDeleted()) {
    print_sequence(i);
  }
  print_sequence(0);
  printf("\n");
}
#endif


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      node_headers methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::node_headers::node_headers(expert_forest &P)
  : parent(P)
{
  //
  // Inltialize address array
  //
  a_size = a_min_size;
  address = (node_header *) malloc(a_size * sizeof(node_header));
  if (0 == address) throw error(error::INSUFFICIENT_MEMORY);
  parent.stats.incMemAlloc(a_size * sizeof(node_header));
  memset(address, 0, a_size * sizeof(node_header));
  a_last = a_next_shrink = 0;
  for (int i=0; i<8; i++) a_unused[i] = 0;
  a_lowest_index = 8;

  //
  // hooks for later
  //
  usesCacheCounts = true;
  usesIncomingCounts = true;
}

// ******************************************************************

MEDDLY::node_headers::~node_headers()
{
  parent.stats.decMemAlloc(a_size * sizeof(node_header));

  // Address array
  free(address);
}

// ******************************************************************

void MEDDLY::node_headers::turnOffCacheCounts()
{
  usesCacheCounts = false;
  // For now - nothing else to do
}

// ******************************************************************

void MEDDLY::node_headers::turnOffIncomingCounts()
{
  usesIncomingCounts = false;
  // For now - nothing else to do
}

// ******************************************************************

MEDDLY::node_handle MEDDLY::node_headers::getFreeNodeHandle() 
{
  MEDDLY_DCASSERT(address);
  parent.stats.incMemUsed(sizeof(node_header));
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
    address[found].setNotDeleted();
    dump_handle_info(address, a_last+1);
#endif
    return found;
  }
  a_last++;
  if (a_last >= a_size) {
    expandHandleList();
  }
  MEDDLY_DCASSERT(a_last < a_size);
#ifdef DEBUG_HANDLE_FREELIST
  address[a_last].setNotDeleted();
  dump_handle_info(address, a_last+1);
#endif
  return a_last;
}

// ******************************************************************

void MEDDLY::node_headers::recycleNodeHandle(node_handle p) 
{
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(0==address[p].cache_count);
  parent.stats.decMemUsed(sizeof(node_header));
  deactivate(p);

  // Determine which list to add this into,
  // and add it there
  int i = bytesRequiredForDown(p) -1;
  setNextOf(p, a_unused[i]);
  a_unused[i] = p;
  a_lowest_index = MIN(a_lowest_index, (char)i);

  // if this was the last node, collapse nodes into the
  // "not yet allocated" pile.  But, we don't remove them
  // from the free list(s); we simply discard any too-large
  // ones when we pull from the free list(s).
  if (p == a_last) {
    while (a_last && isDeactivated(a_last)) {
      a_last--;
    }

    if (a_last < a_next_shrink) shrinkHandleList();
  }
#ifdef DEBUG_HANDLE_FREELIST
  dump_handle_info(address, a_last+1);
#endif
}

// ******************************************************************

void MEDDLY::node_headers::swapNodes(node_handle p, node_handle q, bool swap_incounts)
{
  MEDDLY_DCASSERT(address);
  MEDDLY_DCASSERT(p>0);
  MEDDLY_DCASSERT(p<=a_last);
  MEDDLY_DCASSERT(q>0);
  MEDDLY_DCASSERT(q<=a_last);

  SWAP(address[p], address[q]);

  if (!swap_incounts) {
    SWAP(address[p].incoming_count, address[q].incoming_count);
  } 
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
  s << "\n    Cache  : [";
  for (node_handle p=1; p<=a_last; p++) {
    if (p>1) s.put('|');
    s.put(long(address[p].cache_count), awidth);
  }
  s << "]\n\n";
}

// ******************************************************************

void MEDDLY::node_headers::expandHandleList()
{
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
    throw error(error::INSUFFICIENT_MEMORY);
  }
  address = new_address;
  parent.stats.incMemAlloc(delta * sizeof(node_header));
  memset(address + a_size, 0, delta * sizeof(node_header));
  a_size += delta;
  a_next_shrink = a_size / 2;
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Enlarged address array, new size %d\n", a_size);
#endif
}

// ******************************************************************

void MEDDLY::node_headers::shrinkHandleList()
{
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
    throw error(error::INSUFFICIENT_MEMORY);
  }
  address = new_address;
  parent.stats.decMemAlloc(delta * sizeof(node_header));
  a_size -= delta;
  a_next_shrink = a_size / 2;
  MEDDLY_DCASSERT(a_last < a_size);
#ifdef DEBUG_ADDRESS_RESIZE
  printf("Shrank address array, new size %d\n", a_size);
#endif
}


// ******************************************************************

