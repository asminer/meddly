
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

#include "freelists.h"
#include "../io.h"
#include "../error.h"

// #define DEBUG_FREELISTS

// #define USE_SLOW_GET_CHUNK_ADDRESS

// ******************************************************************
// *                                                                *
// *                     freelist_manager class                     *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

/**
    Simple hole management based on free lists.
    We don't try to merge holes at all.
    Instead we maintain, for each size, a list of holes.
*/
template <class INT>
class freelist_manager : public memory_manager {
  public:
    freelist_manager(const char* n, memstats &stats);
    virtual ~freelist_manager();

    virtual bool mustRecycleManually() const {
      return false;
    }

    virtual bool firstSlotMustClearMSB() const {
      return false;
    }

    virtual bool lastSlotMustClearMSB() const {
      return false;
    }

    virtual node_address requestChunk(size_t &numSlots);
    virtual void recycleChunk(node_address h, size_t numSlots);
#ifdef USE_SLOW_GET_CHUNK_ADDRESS
    virtual void* slowChunkAddress(node_address h) const;
#endif
    virtual bool isValidHandle(node_address h) const;

    virtual void reportStats(output &s, const char* pad, bool human, bool details) const;
    virtual void dumpInternal(output &s) const;

    virtual node_address getFirstAddress() const;
    virtual bool isAddressInUse(node_address addr) const;
    virtual node_address getNextAddress(node_address addr) const;
    virtual void dumpInternalUnused(output &s, node_address addr) const;

  private:
    INT* entries;
    int entriesSize;
    int entriesAlloc;

    static const int maxEntrySize = 15;
    static const size_t maxEntryBytes = sizeof(int) * maxEntrySize;
    INT* freeList;
}; // class freelist_manager

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                    freelist_manager methods                    *
// *                                                                *
// ******************************************************************

template <class INT>
MEDDLY::freelist_manager<INT>::freelist_manager(const char* n, memstats &stats)
 : memory_manager(n, stats)
{
  freeList = new INT[1+maxEntrySize];
  for (int i=0; i<=maxEntrySize; i++) {
    freeList[i] = 0;
  }
  incMemUsed( (1+maxEntrySize) * sizeof(INT) );
  incMemAlloc( (1+maxEntrySize) * sizeof(INT) );

  entriesAlloc = 1024;
  entriesSize = 1;
  entries = (INT*) malloc(entriesAlloc * sizeof(INT) );
  if (0==entries) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  entries[0] = 0;     // NEVER USED; set here for sanity.

  incMemUsed( entriesSize * sizeof(INT) );
  incMemAlloc( entriesAlloc * sizeof(INT) );

#ifndef USE_SLOW_GET_CHUNK_ADDRESS
  setChunkBase(entries);
  setChunkMultiplier(sizeof(INT));
#endif
}

// ******************************************************************

template <class INT>
MEDDLY::freelist_manager<INT>::~freelist_manager()
{
  delete[] entries;
  delete[] freeList;

  // Clever way to update stats!

  zeroMemUsed();
  zeroMemAlloc();
}

// ******************************************************************

template <class INT>
MEDDLY::node_address MEDDLY::freelist_manager<INT>::requestChunk(size_t &numSlots)
{
  // check free list
  if (numSlots > maxEntrySize) {
    fprintf(stderr, "MEDDLY error: request for compute table entry larger than max size\n");
    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);  // best we can do
  }
  if (numSlots<1) return 0;
  if (freeList[numSlots]) {
    INT h = freeList[numSlots];
    freeList[numSlots] = entries[h];
#ifdef DEBUG_FREELISTS
    fprintf(stderr, "Re-used entry %ld size %ld\n", node_address(h), numSlots);
#endif
    incMemUsed(numSlots * sizeof(INT));
    return h;
  }
  if (entriesSize + numSlots > entriesAlloc) {
    // Expand by a factor of 1.5
    size_t neA = entriesAlloc + (entriesAlloc/2);
    INT* ne = (INT*) realloc(entries, neA * sizeof(INT));
    if (0==ne) {
      fprintf(stderr,
          "Error in allocating array of size %lu at %s, line %d\n",
          neA * sizeof(INT), __FILE__, __LINE__);
      throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    incMemAlloc( (entriesAlloc/2) * sizeof(INT) );

    entries = ne;
    entriesAlloc = neA;

#ifndef USE_SLOW_GET_CHUNK_ADDRESS
    setChunkBase(entries);
#endif
  }
  MEDDLY_DCASSERT(entriesSize + numSlots <= entriesAlloc);
  node_address h = entriesSize;
  entriesSize += numSlots;
#ifdef DEBUG_FREELISTS
  fprintf(stderr, "New entry %ld size %ld\n", h, numSlots);
#endif
  incMemUsed(numSlots * sizeof(INT));
  return h;

}

// ******************************************************************

template <class INT>
void MEDDLY::freelist_manager<INT>::recycleChunk(node_address h, size_t numSlots)
{
#ifdef DEBUG_FREELISTS
  fprintf(stderr, "Recycling entry %ld size %ld\n", h, numSlots);
#endif
  entries[h] = freeList[numSlots];
  freeList[numSlots] = h;
  decMemUsed(numSlots * sizeof(INT));
}

// ******************************************************************

#ifdef USE_SLOW_GET_CHUNK_ADDRESS

template <class INT>
void* MEDDLY::freelist_manager<INT>::slowChunkAddress(node_address h) const
{
  return (void*)(entries + h);
}

#endif

// ******************************************************************

template <class INT>
bool MEDDLY::freelist_manager<INT>::isValidHandle(node_address h) const
{
  return (h >= 1) && (h < entriesSize);
}

// ******************************************************************

template <class INT>
void MEDDLY::freelist_manager<INT>::reportStats(output &s, const char* pad,
  bool human, bool details) const
{
  s << pad << "Report for freelist memory manager:\n";
  if (details) {
    s << pad << "  Current #holes: \n";
  }
  long total = 0;
  long tbytes = 0;
  for (int i=0; i<maxEntrySize; i++) {
    long count=0;
    for (INT ptr = freeList[i]; ptr; ptr = entries[ptr]) {
      count++;
    }
    if (0==count) continue;
    if (details) {
      s << pad << "    size " << i << ": " << count << "\n";
    }
    total += count;
    tbytes += count * i * sizeof(INT);
  }
  if (details) {
    s << pad << "  total " << total << " holes, " << tbytes << " bytes\n";
  } else {
    s << pad << "  Current #holes: " << total << ", " << tbytes << " bytes\n";
  }
}

// ******************************************************************

template <class INT>
void MEDDLY::freelist_manager<INT>::dumpInternal(output &s) const
{
  if (0==entries) {
    s << "Entries: null\n";
  } else {
    s << "Entries: [" << long(entries[0]);
    for (int i=1; i<entriesSize; i++) {
      s << ", " << long(entries[i]);
    }
    s << "]\n";
  }
  if (0==freeList) {
    s << "Free: null\n";
  } else {
    s << "Free: [" << long(freeList[0]);
    for (int i=1; i<=maxEntrySize; i++) {
      s << ", " << long(freeList[i]);
    }
    s << "]\n";
  }
}

// ******************************************************************

template <class INT>
MEDDLY::node_address MEDDLY::freelist_manager<INT>::getFirstAddress() const
{
  return 0;
}

// ******************************************************************

template <class INT>
bool MEDDLY::freelist_manager<INT>::isAddressInUse(node_address addr) const
{
  return false;   // Can't tell, don't even try
}

// ******************************************************************

template <class INT>
MEDDLY::node_address MEDDLY::freelist_manager<INT>::getNextAddress(node_address addr) const
{
  return 0;
}

// ******************************************************************

template <class INT>
void MEDDLY::freelist_manager<INT>::dumpInternalUnused(output &s, node_address addr) const
{
  s << "Free: [ next: " << long(entries[addr]) << " ]\n";
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     freelist_style methods                     *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::freelist_style::freelist_style(const char* n)
: memory_manager_style(n)
{
}

MEDDLY::freelist_style::~freelist_style()
{
}

MEDDLY::memory_manager*
MEDDLY::freelist_style::initManager(unsigned char granularity,
  unsigned char minsize, memstats &stats) const
{
  switch (granularity) {

    case sizeof(short):
      return new freelist_manager<short>(getName(), stats);

    case sizeof(int):
      return new freelist_manager<int>(getName(), stats);

    case sizeof(long):
      return new freelist_manager<long>(getName(), stats);

    default:
      return 0;
  }
}

