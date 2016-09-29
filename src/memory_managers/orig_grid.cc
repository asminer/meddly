
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
        MEDDLY_CHECK_RANGE(0, h, data_alloc);
        return data + h;
      }

      virtual void reportStats(output &s, const char* pad, bool details) const;
      virtual void dumpInternal(output &s) const;

    private:
      INT* data; 
      long data_alloc;

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
  }

  template <class INT>
  original_grid<INT>::~original_grid()
  {
    // TBD - free data, update stats
  }

  template <class INT>
  unsigned long original_grid<INT>::requestChunk(size_t &numSlots)
  {
    // BIG FAT TBD HERE

    // FAIL
    numSlots = 0;
    return 0;
  }

  template <class INT>
  void original_grid<INT>::recycleChunk(unsigned long h, size_t numSlots)
  {
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

/*
    s << "Sanity check, MSB mask: ";
    INT x = 1;
    x <<= (8*sizeof(INT))-1;
    s.put_hex(x);
    s << "\n";
*/

    // TBD - print lots of stuff here
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

