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

#ifndef MEDDLY_HOLE_BASE_H
#define MEDDLY_HOLE_BASE_H

#include "../memory.h"
#include "orig_grid.h"
#include "../io.h"

#if 0
#ifdef HAVE_MALLOC_GOOD_SIZE
#include <malloc/malloc.h>
#endif
#endif

// #define USE_SLOW_GET_CHUNK_ADDRESS

namespace MEDDLY {

  // soon:
  // class hole_manager_byte;

  // ******************************************************************
  // *                                                                *
  // *                    holeman (template) class                    *
  // *                                                                *
  // ******************************************************************

  /**
    Abstract template class for hole management.
    INT should be a signed integer type.
    Several of the concrete memoery managers are derived from this class.

    Hole requirements:
      First slot is the hole sizze, with MSB set
      Last slot is the hole size, with MSB set

  */
  template <class INT>
  class hole_manager : public memory_manager {
    public:
      hole_manager(const char* n, memstats &stats);
      virtual ~hole_manager();

      // common stuff!

      virtual bool mustRecycleManually() const {
        return false;
      }

      virtual bool firstSlotMustClearMSB() const {
        return true;
      }

      virtual bool lastSlotMustClearMSB() const {
        return true;
      }

#ifdef USE_SLOW_GET_CHUNK_ADDRESS
      virtual void* slowChunkAddress(node_address h) const {
        MEDDLY_DCASSERT(data);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1lu, h, data_alloc);
        return data + h;
      }
#endif

      virtual bool isValidHandle(node_address h) const {
        return (h >= 1) && (h<data_alloc);
      }

      virtual node_address getFirstAddress() const {
        return 1;
      }

      virtual bool isAddressInUse(node_address addr) const {
        if (0==addr) return false;
        if (addr > last_used_slot) return false;
        return !isHole(addr);
      }

      virtual node_address getNextAddress(node_address addr) const {
        if (0==addr) return 0;
        if (addr > last_used_slot) return 0;
        return addr + getHoleSize(addr);
      }

      void showInternal(output &s) const;
      void showInternalAddr(output &s, node_address addr, int slots) const;

    protected:
      /// Grab a hole from the end of the array
      node_address allocateFromArray(size_t numSlots);

      /** Try to merge the given hole with the end of the array.
          If this is not possible, do nothing.
            @param  addr      Address of hole.
            @param  numSlots  Size of hole.
            @return true      iff the hole was recycled.
      */
      inline bool recycleHoleInArray(node_address addr, size_t numSlots) {
        if (addr + numSlots - 1 == last_used_slot) {
#ifdef MEMORY_TRACE_DETAILS
          printf("\tMerging chunk with unused end bit\n");
#endif
          last_used_slot = addr-1;
          return true;
        } else {
          return false;
        }
      }

    protected:
      inline bool isHole(node_address h) const {
        MEDDLY_DCASSERT(data);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0lu, h, 1+last_used_slot);
        return data[h] & MSB;
        // Because we set data[0] to 0, this will work
        // correctly also for h=0 (which is not a hole).
      }
      inline INT getHoleSize(node_address h) const {
        MEDDLY_DCASSERT(isHole(h));
        return data[h] & (~MSB);
      }
      inline void setHoleSize(node_address h, INT hs) {
        MEDDLY_DCASSERT(data);
        MEDDLY_DCASSERT(hs>0);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1lu, h, 1+last_used_slot);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1lu, h+hs-1, 1+last_used_slot);
        data[h] = data[h+hs-1] = (hs | MSB);
      }
      inline void clearHole(node_address h, INT hs) const {
        MEDDLY_DCASSERT(isHole(h));
        data[h] = data[h+hs-1] = 0;
      }
      inline bool matchingHoleSizes(node_address h) const {
            const unsigned long hs = getHoleSize(h);
            return data[h] == data[h+ hs - 1];
      }

      inline INT readSlot(node_address h, const int slot) const {
        MEDDLY_DCASSERT(isHole(h));
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1lu, h+slot, 1+last_used_slot);
        return data[h+slot];
      }

      inline INT& refSlot(node_address h, const unsigned slot) {
        MEDDLY_DCASSERT(isHole(h));
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1lu, h+slot, 1+last_used_slot);
        return data[h+slot];
      }

      inline node_address max_handle() const {
        return (~MSB);
      }

      inline node_address getLastUsed() const {
        return last_used_slot;
      }

    private:
      /// @return true on success
      bool resize(size_t newalloc);

    private:
      INT* data;
      node_address data_alloc;
      node_address last_used_slot;

      INT MSB;

  }; // class hole_manager


};  // namespace MEDDLY


// ******************************************************************
// *                                                                *
// *                        holeman  methods                        *
// *                                                                *
// ******************************************************************

template <class INT>
MEDDLY::hole_manager<INT>::hole_manager(const char* n, memstats &stats)
  : memory_manager(n, stats)
{
  data = 0;
  data_alloc = 0;
  last_used_slot = 0;

  MSB = 1;
  MSB <<= (8*sizeof(INT) - 1);
#ifndef USE_SLOW_GET_CHUNK_ADDRESS
  setChunkBase(data);
  setChunkMultiplier(sizeof(INT));
#endif
}

// ******************************************************************

template <class INT>
MEDDLY::hole_manager<INT>::~hole_manager()
{
  decMemAlloc(data_alloc * sizeof(INT));
  free(data);
}

// ******************************************************************

template <class INT>
void MEDDLY::hole_manager<INT>::showInternal(output &s) const
{
  s << "  data pointer: ";
  s.put_hex((unsigned long)data);
  s << "\n";
  s << "  data_alloc: " << data_alloc << "\n";
  s << "  last_used_slot: " << last_used_slot << "\n\n";
}

// ******************************************************************

template <class INT>
void MEDDLY::hole_manager<INT>::showInternalAddr(output &s, node_address addr, int slots) const
{
  if (0==addr) return;
  if (addr > last_used_slot) {
    s << "free slots";
    return;
  }

  if (data[addr] & MSB) {
    s << "1:";
  } else {
    s << "0:";
    MEDDLY_DCASSERT(0);
  }
  s << getHoleSize(addr);
  for (int i=1; i<=slots; i++) {
    s << ", " << data[addr+i];
  }
  s << ", ..., ";
  addr += getHoleSize(addr)-1;
  if (data[addr] & MSB) {
    s << "1:";
  } else {
    s << "0:";
    MEDDLY_DCASSERT(0);
  }
  s << getHoleSize(addr);
}

// ******************************************************************

template <class INT>
MEDDLY::node_address MEDDLY::hole_manager<INT>::allocateFromArray(size_t numSlots)
{
  if (last_used_slot + numSlots >= data_alloc) {
    //
    // Expand.
    //

    bool ok = false;

    if (0==data_alloc) {
      if (numSlots < 512) ok = resize(1024);
      else                ok = resize(2*numSlots);
    } else {
      size_t want_size = last_used_slot + numSlots;
      want_size += want_size/2;
      ok = resize(want_size);
    }

    if (!ok) {
      //
      // Couldn't resize, fail cleanly
      //

      return 0;
    }
  }

  //
  // Grab node from the end
  //
  node_address h = last_used_slot + 1;
  if (h > max_handle()) {
    return 0;
  }
  last_used_slot += numSlots;
  return h;
}

// ******************************************************************

template <class INT>
bool MEDDLY::hole_manager<INT>::resize(size_t new_alloc)
{
#if 0
#ifdef HAVE_MALLOC_GOOD_SIZE
  size_t good_bytes = malloc_good_size(new_alloc * sizeof(INT));
  new_alloc = good_bytes / sizeof(INT);
#endif
#endif

#ifdef TRACE_REALLOCS
  if (new_alloc > data_alloc) printf("enlarging"); else printf("shrinking");
  printf(" data %lx, new size %ld\n", (unsigned long)data, new_alloc);
#endif

  INT* new_data = (INT*) realloc(data, new_alloc * sizeof(INT));

#ifdef TRACE_REALLOCS
  if (new_data != data) {
    printf("data moved to %lx\n", (unsigned long)new_data);
  }
#endif

  if (0==new_data && (new_alloc!=0)) {
    return false;
  }
  if ((0==data) && new_data) new_data[0] = 0;

  if (new_alloc > data_alloc) {
    incMemAlloc((new_alloc - data_alloc) * sizeof(INT));
  } else {
    decMemAlloc((data_alloc - new_alloc) * sizeof(INT));
  }

  data_alloc = new_alloc;
  data = new_data;

#ifndef USE_SLOW_GET_CHUNK_ADDRESS
  setChunkBase(data);
#endif

  return true;
}

// ******************************************************************

#endif  // #include guard

