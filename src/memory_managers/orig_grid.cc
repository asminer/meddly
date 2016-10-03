
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "orig_grid.h"

// new grid class here

namespace MEDDLY {

  // ******************************************************************
  // *                                                                *
  // *                 original_grid (template) class                 *
  // *                                                                *
  // ******************************************************************

  /**
    Grid structure for hole management.
    INT must be a signed integer type.

    Index hole:
    -------------------------------------------------------------
    [0] size with MSB set
    [1] up, >= 0
    [2] down, >= 0
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

  */
  template <class INT>
  class original_grid : public memory_manager {

    public:
      original_grid();
      virtual ~original_grid();

      virtual bool mustRecycleManually() const {
        return false;
      }

      virtual bool firstSlotMustClearMSB() const {
        return true;
      }

      virtual bool lastSlotMustClearMSB() const {
        return true;
      }

      virtual unsigned long requestChunk(size_t &numSlots);
      virtual void recycleChunk(unsigned long h, size_t numSlots);

      virtual void* getChunkAddress(unsigned long h) const {
        MEDDLY_DCASSERT(h < data_alloc);
        return data + h;
      }

      virtual void reportStats(output &s, const char* pad, bool details) const;
      virtual void dumpInternal(output &s) const;

    private:
      /// @return true on success
      bool resize(long newalloc);

    private:
      inline bool isHole(unsigned long h) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_DCASSERT(h <= last_used_slot);
        return data[h] & MSB;
      }
      inline INT getHoleSize(unsigned long h) const {
        MEDDLY_DCASSERT(data);
        MEDDLY_DCASSERT(h <= last_used_slot);
        return data[h] & (~MSB);
      }
      inline void setHoleSize(unsigned long h, INT hs) {
        MEDDLY_DCASSERT(data);
        MEDDLY_DCASSERT(h <= last_used_slot);
        MEDDLY_DCASSERT(hs>0);
        data[h] = data[h+hs-1] = (hs | MSB);
      }

      inline bool isIndexHole(unsigned long h) {
        MEDDLY_DCASSERT(data);
        MEDDLY_DCASSERT(h <= last_used_slot);
        return data[h+1] >= 0;
      }
      
      inline INT& Up(unsigned long h) {
        MEDDLY_DCASSERT(isHole(h));
        MEDDLY_DCASSERT(isIndexHole(h));
        return data[h+1];
      }

      inline INT& Down(unsigned long h) {
        MEDDLY_DCASSERT(isHole(h));
        MEDDLY_DCASSERT(isIndexHole(h));
        return data[h+2];
      }

      inline INT& Prev(unsigned long h) {
        MEDDLY_DCASSERT(isHole(h));
        MEDDLY_DCASSERT(!isIndexHole(h));
        return data[h+2];
      }

      inline INT& Next(unsigned long h) {
        MEDDLY_DCASSERT(isHole(h));
        return data[h+3];
      }


    private:
      INT* data; 
      long data_alloc;
      long last_used_slot;
      INT MSB;
  }; // class original_grid


  // ******************************************************************
  // *                                                                *
  // *                     original_grid  methods                     *
  // *                                                                *
  // ******************************************************************

  template <class INT>
  original_grid<INT>::original_grid()
  {
    data = 0;
    data_alloc = 0;
    last_used_slot = 0;

    MSB = 1;
    MSB <<= (8*sizeof(INT))-1;
  }

  template <class INT>
  original_grid<INT>::~original_grid()
  {
    // TBD - free data, update stats
  }

  template <class INT>
  unsigned long original_grid<INT>::requestChunk(size_t &numSlots)
  {
    // HUGE TBD : search the grid for a matching hole

    //
    // Still here?  We couldn't recycle a chunk.
    //

    if (last_used_slot + numSlots > data_alloc) {
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

        numSlots = 0;
        return 0;
      }
    }

    //
    // Grab node from the end
    //
    unsigned long h = last_used_slot + 1;
    last_used_slot += numSlots;
    return h;
  }

  template <class INT>
  void original_grid<INT>::recycleChunk(unsigned long h, size_t numSlots)
  {
    decMemUsed(sizeof(INT) * numSlots);

    // ANOTHER HUGE TBD HERE
  }

  template <class INT>
  void original_grid<INT>::reportStats(output &s, const char* pad, bool details) const
  {
    s << pad << "Report for original_grid memory manager:\n";
    // TBD - print stuff here
  }

  template <class INT>
  void original_grid<INT>::dumpInternal(output &s) const
  {
    s << "Internal storage for original_grid:\n";
    s << "  data pointer: ";
    s.put_hex((unsigned long)data);
    s << "\n";
    s << "  data_alloc: " << data_alloc << "\n";

    s << "Sanity check, MSB mask: ";
    s.put_hex(MSB);
    s << "\n";

    // TBD - print lots of stuff here
  }

  template <class INT>
  bool original_grid<INT>::resize(long new_alloc) {
    MEDDLY_DCASSERT(new_alloc >= 0);

#ifdef HAVE_MALLOC_GOOD_SIZE
/*
    size_t good_bytes = malloc_good_size(new_alloc * sizeof(INT));
    new_alloc = good_bytes / sizeof(INT);
*/
#endif

    INT* new_data = (INT*) realloc(data, new_alloc * sizeof(INT));

    if (0==new_data && (new_alloc!=0)) {
      return false;
    }
    if (0==data) new_data[0] = 0;

    if (new_alloc > data_alloc) {
      incMemAlloc((new_alloc - data_alloc) * sizeof(INT));
    } else {
      decMemAlloc((data_alloc - new_alloc) * sizeof(INT));
    }

    data_alloc = new_alloc;
    data = new_data;
  }


}; // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    orig_grid_style  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::orig_grid_style::orig_grid_style()
{
}

MEDDLY::orig_grid_style::~orig_grid_style()
{
}

MEDDLY::memory_manager*
MEDDLY::orig_grid_style::initManager(unsigned char granularity, unsigned char minsize) const
{
  if (sizeof(int) == granularity) {
    return new original_grid <int>;
  }

  if (sizeof(long) == granularity) {
    return new original_grid <long>;
  }

  if (sizeof(short) == granularity) {
    return new original_grid <short>;
  }

  // unsupported granularity

  return 0;
}

