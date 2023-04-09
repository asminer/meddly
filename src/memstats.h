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

#ifndef MEDDLY_MEMSTATS_H
#define MEDDLY_MEMSTATS_H

#include "defines.h"

#include <cstddef>  // for size_t

namespace MEDDLY {
    class memstats;
};


/**
    Class for memory statistics.
*/
class MEDDLY::memstats {
    public:
        memstats();

        void incMemUsed(size_t b);
        void decMemUsed(size_t b);
        void incMemAlloc(size_t b);
        void decMemAlloc(size_t b);

        void zeroMemUsed();
        void zeroMemAlloc();

        size_t getMemUsed() const;
        size_t getMemAlloc() const;
        size_t getPeakMemUsed() const;
        size_t getPeakMemAlloc() const;

        static size_t getGlobalMemUsed();
        static size_t getGlobalMemAlloc();
        static size_t getGlobalPeakMemUsed();
        static size_t getGlobalPeakMemAlloc();

    private:
        /// Current memory used
        size_t memory_used;
        /// Current memory allocated
        size_t memory_alloc;
        /// Peak memory used
        size_t peak_memory_used;
        /// Peak memory allocated
        size_t peak_memory_alloc;

        // global memory usage
        static size_t global_memory_used;
        static size_t global_memory_alloc;
        static size_t global_peak_used;
        static size_t global_peak_alloc;
};


// ******************************************************************
// *                                                                *
// *                    memstats inline methods                     *
// *                                                                *
// ******************************************************************

inline void MEDDLY::memstats::incMemUsed(size_t b)
{
  memory_used += b;
  if (memory_used > peak_memory_used) {
    peak_memory_used = memory_used;
  }
  global_memory_used += b;
  if (global_memory_used > global_peak_used) {
    global_peak_used = global_memory_used;
  }
}

inline void MEDDLY::memstats::decMemUsed(size_t b)
{
  MEDDLY_DCASSERT(memory_used >= b);
  memory_used -= b;
  MEDDLY_DCASSERT(global_memory_used >= b);
  global_memory_used -= b;
}

inline void MEDDLY::memstats::incMemAlloc(size_t b)
{
  memory_alloc += b;
  if (memory_alloc > peak_memory_alloc) {
    peak_memory_alloc = memory_alloc;
  }
  global_memory_alloc += b;
  if (global_memory_alloc > global_peak_alloc) {
    global_peak_alloc = global_memory_alloc;
  }
}

inline void MEDDLY::memstats::decMemAlloc(size_t b)
{
  MEDDLY_DCASSERT(memory_alloc >= b);
  memory_alloc -= b;
  MEDDLY_DCASSERT(global_memory_alloc >= b);
  global_memory_alloc -= b;
}

inline void MEDDLY::memstats::zeroMemUsed()
{
  MEDDLY_DCASSERT(global_memory_used >= memory_used);
  global_memory_used -= memory_used;
  memory_used = 0;
}

inline void MEDDLY::memstats::zeroMemAlloc()
{
  MEDDLY_DCASSERT(global_memory_alloc >= memory_alloc);
  global_memory_alloc -= memory_alloc;
  memory_alloc = 0;
}

inline size_t MEDDLY::memstats::getMemUsed() const
{
  return memory_used;
}

inline size_t MEDDLY::memstats::getMemAlloc() const
{
  return memory_alloc;
}

inline size_t MEDDLY::memstats::getPeakMemUsed() const
{
  return peak_memory_used;
}

inline size_t MEDDLY::memstats::getPeakMemAlloc() const
{
  return peak_memory_alloc;
}


inline size_t MEDDLY::memstats::getGlobalMemUsed()
{
  return global_memory_used;
}

inline size_t MEDDLY::memstats::getGlobalMemAlloc()
{
  return global_memory_alloc;
}

inline size_t MEDDLY::memstats::getGlobalPeakMemUsed()
{
  return global_peak_used;
}

inline size_t MEDDLY::memstats::getGlobalPeakMemAlloc()
{
  return global_peak_alloc;
}


#endif
