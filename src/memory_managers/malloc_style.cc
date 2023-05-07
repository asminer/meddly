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

#include "malloc_style.h"
#include "../io.h"

// #define MEMORY_TRACE

// new grid class here

// #define USE_SLOW_GET_CHUNK_ADDRESS

namespace MEDDLY {
  class malloc_manager;
};

// ******************************************************************
// *                                                                *
// *                      malloc_manager class                      *
// *                                                                *
// ******************************************************************

/**
    Simple malloc/free for hole management.
*/
class MEDDLY::malloc_manager : public memory_manager {
  public:
    malloc_manager(const char* n, memstats &stats, unsigned char gran);
    virtual ~malloc_manager();

    virtual bool mustRecycleManually() const {
      return true;
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
    unsigned long bytes_allocd_not_freed;
    unsigned char granularity;
}; // class malloc_manager


// ******************************************************************
// *                                                                *
// *                     malloc_manager methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::malloc_manager::malloc_manager(const char* n, memstats &stats,
  unsigned char gran) : memory_manager(n, stats)
{
  granularity = gran;
  bytes_allocd_not_freed = 0;
#ifndef USE_SLOW_GET_CHUNK_ADDRESS
  setChunkBase(0);
  setChunkMultiplier(1);
#endif
}

// ******************************************************************

MEDDLY::malloc_manager::~malloc_manager()
{
  // Nothing we can do
}

// ******************************************************************

MEDDLY::node_address MEDDLY::malloc_manager::requestChunk(size_t &numSlots)
{
  size_t bytes = numSlots * granularity;
  void* chunk = malloc(bytes);
  if (chunk) {
    bytes_allocd_not_freed += bytes;
    incMemUsed(bytes);
    incMemAlloc(bytes);
  }
#ifdef MEMORY_TRACE
  printf("Mallloc manager returning chunk %lx = %ld\n", chunk, node_address(chunk));
#endif
  return node_address(chunk);
}

// ******************************************************************

void MEDDLY::malloc_manager::recycleChunk(node_address h, size_t numSlots)
{
  size_t bytes = numSlots * granularity;
  void* chunk = (void*) h;
  free(chunk);
  bytes_allocd_not_freed -= bytes;
  decMemUsed(bytes);
  decMemAlloc(bytes);
}

// ******************************************************************

#ifdef USE_SLOW_GET_CHUNK_ADDRESS

void* MEDDLY::malloc_manager::slowChunkAddress(node_address h) const
{
  return (void*) h;
}

#endif

// ******************************************************************

bool MEDDLY::malloc_manager::isValidHandle(node_address h) const
{
  return h;
}

// ******************************************************************

void MEDDLY::malloc_manager::reportStats(output &s, const char* pad,
  bool human, bool details) const
{
}

// ******************************************************************

void MEDDLY::malloc_manager::dumpInternal(output &s) const
{
  s << "Malloc manager: uses malloc, nothing to display\n";
}

// ******************************************************************

MEDDLY::node_address MEDDLY::malloc_manager::getFirstAddress() const
{
  return 0;
}

// ******************************************************************

bool MEDDLY::malloc_manager::isAddressInUse(node_address addr) const
{
  return false;
  // we have no way of knowing
}

// ******************************************************************

MEDDLY::node_address MEDDLY::malloc_manager::getNextAddress(node_address addr) const
{
  return 0;
}

// ******************************************************************

void MEDDLY::malloc_manager::dumpInternalUnused(output &s, node_address addr) const
{
  s << "Malloc manager: uses malloc, nothing unused\n";
}


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      malloc_style methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::malloc_style::malloc_style(const char* n)
: memory_manager_style(n)
{
}

MEDDLY::malloc_style::~malloc_style()
{
}

MEDDLY::memory_manager*
MEDDLY::malloc_style::initManager(unsigned char granularity,
  unsigned char minsize, memstats &stats) const
{
  return new malloc_manager(getName(), stats, granularity);
}

