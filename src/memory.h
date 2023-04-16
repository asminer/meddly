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

#ifndef MEDDLY_MEMORY_H
#define MEDDLY_MEMORY_H

#include "memstats.h"

namespace MEDDLY {
    class memory_manager_style;
    class memory_manager;
    class output;
};

// ******************************************************************
// *                                                                *
// *                   memory_manager_style class                   *
// *                                                                *
// ******************************************************************

/** Abstract base class for memory manager factories.

    This allows us to use specialized implementations of
    memory managers (say, using templates) based on the granularity.
*/
class MEDDLY::memory_manager_style {
    const char* name;
  public:
    memory_manager_style(const char* n);
    virtual ~memory_manager_style();

    /**
        Build a new memory manager.

          @param  granularity   Unit of storage, in bytes.
                                Must be greater than 0.
                                All sizes specified to the memory manager,
                                for allocating and freeing chunks, are in
                                terms of the granularity.  For example, to
                                manage a collection of arrays of integers,
                                set the granularity to be sizeof(int) and
                                use sizes equal to the number of integers.
                                For behavior exactly the same as malloc,
                                use a granularity of 1.


          @param  minsize       The smallest size chunk that will ever be
                                requested.  Must be greater than 0.
                                This is specified here in case
                                that information can help the memory manager.

          @param  stats         Structure to use for updating memory stats.

          @return   A pointer to a new instance of a memory manager, or 0
                    if some error occurred, for example if the requested
                    granularity cannot be supported by this type of memory
                    manager.
    */
    virtual memory_manager* initManager(unsigned char granularity,
      unsigned char minsize, memstats& stats) const = 0;


    /**
        Human readable name.
        Used for debugging and reporting.
    */
    const char* getName() const;
};

// ******************************************************************
// *                                                                *
// *                      memory_manager class                      *
// *                                                                *
// ******************************************************************

/**
    Interface for memory managers.
*/
class MEDDLY::memory_manager {

  public:
    memory_manager(const char* sn, memstats& stats);
    virtual ~memory_manager();

    /**
        Is this memory manager unable to free everything on its own?

          @return   True, if destroying the memory manager DOES NOT
                    automatically recycle all non-freed requested chunks
                    (because the memory manager does not track them).

                    False, if destroying the memory manager DOES
                    automatically recycle all non-freed requested chunks.
    */
    virtual bool mustRecycleManually() const = 0;

    /**
        Does the memory manager require that the most significant bit
        of the first slot is cleared?

        (If true, the memory manager sets this bit for recycled chunks.)

          @return   True, if the first slot of a chunk must hold a value
                    such that the most significant bit is cleared.
                    If this cannot be guaranteed by what is stored there,
                    then you must allocate an extra slot and not use the
                    first one.

                    False, if there are no restrictions on what may be
                    stored in the first slot of a chunk.
    */
    virtual bool firstSlotMustClearMSB() const = 0;

    /**
        Does the memory manager require that the most significant bit
        of the last slot is cleared?

        (If true, the memory manager sets this bit for recycled chunks.)

          @return   True, if the last slot of a chunk must hold a value
                    such that the most significant bit is cleared.
                    If this cannot be guaranteed by what is stored there,
                    then you must allocate an extra slot and not use the
                    last one.

                    False, if there are no restrictions on what may be
                    stored in the last slot of a chunk.
    */
    virtual bool lastSlotMustClearMSB() const = 0;


    /**
        Request a chunk of memory.

          @param  numSlots    Number of slots.
                              INPUT: number of requested slots.
                              OUTPUT: number of slots in the given chunk,
                              might be larger than the number requested.
                              Will not be smaller than the number requested,
                              unless a failure occurred, in which case
                              it will be set to zero.

          @return   A non-zero handle for a new chunk of memory, containing
                    numSlots slots (and requiring numSlots * granularity
                    bytes), on success.
                    Zero, on failure.

    */
    virtual node_address requestChunk(size_t &numSlots) = 0;

    /**
        Recycle a chunk of memory.

          @param  h           Handle of the chunk, as returned by
                              method requestChunk().

          @param  numSlots    Total number of slots in the chunk.
    */
    virtual void recycleChunk(node_address h, size_t numSlots) = 0;

    /**
        Convert a handle to an actual pointer we can use.

        A default, fast, inlined implementation that will work for
        most memory managers is implemented here, based on
          address = base + m * h
        where m*h is the number of bytes to shift base by.
        This requires derived classes to maintain the pointer "base"
        and multiplier "m" by calling protected methods
            void setChunkBase(void* base)
            void setChunkMultiplier(unsigned int m)
        Note that if m is zero (its default value), this method
        will fall back to slowChunkAddress().

          @param  h     Non-null handle of the chunk, as returned by requestChunk().

          @return       If h is 0, or an invalid handle, then the result is
                        undefined.  Otherwise, we return a pointer to the
                        chunk given by handle h.
                        This pointer is guaranteed to be fixed until the next
                        call to requestChunk() or recycleChunk(); after that,
                        the pointer for a handle could change.
    */
    void* getChunkAddress(node_address h) const;


  protected:
    /**
        See getChunkAddress.
        This method (default behavior is assert(false))
        should be overridden in derived classes when
        the mapping from node address h to pointer
        does not follow the formula used by getChunkAddress().

        Convert a handle to an actual pointer we can use.

          @param  h     Handle of the chunk, as returned by requestChunk(), or 0.

          @return       If h is 0, then we return 0.  Otherwise, we return
                        a pointer to the chunk given by handle h.
                        This pointer is guaranteed to be fixed until the next
                        call to requestChunk() or recycleChunk(); after that,
                        the pointer for a handle could change.
    */
    virtual void* slowChunkAddress(node_address h) const;


  public:
    /**
        Check if a handle is valid (non-null or otherwise).
        Since this is not always possible,
        this method is conservative.

        @return   false if h is null or definitely invalid;
                  true otherwise.
    */
    virtual bool isValidHandle(node_address h) const = 0;

    /** Show various statistics.
          @param  s         Output stream to write to
          @param  pad       Padding string, written at the start of
                            each output line.
          @param  human     If false, just display raw byte counts.
                            If true, use units (e.g., Mbytes, Kbytes).
          @param  details   If false, just display basic statistics.
                            If true, display details.
    */
    virtual void reportStats(output &s, const char* pad,
      bool human, bool details) const = 0;


    /** Display manager-specific internals.
        For debugging.
          @param  s       Output stream to use
    */
    virtual void dumpInternal(output &s) const = 0;



    /** Get first address.
        Used to cycle through all addresses, if we can.
        If we cannot, then always return 0.
    */
    virtual node_address getFirstAddress() const = 0;

    /** Is a given address in use?
        Best effort answer only.
        If unsure, return false.
    */
    virtual bool isAddressInUse(node_address addr) const = 0;

    /** Get the next address.
        Used to cycle through all addresses, if we can.
        If we cannot, then always return 0.
          @param  addr    An address given by any of the methods that
                          "cycle through addresses", and one where
                          isAddressInUse() returned false.

          @return   Next address to check,  Presumably this address
                    is unused and is a hole of some kind.
    */
    virtual node_address getNextAddress(node_address addr) const = 0;

    /** Show information about an unused address.
        If we do not track unused addresses, then do nothing.
          @param  addr    An address given by any of the methods that
                          "cycle through addresses", and one where
                          isAddressInUse() returned false.
    */
    virtual void dumpInternalUnused(output &s, node_address addr) const = 0;

  protected:
    void incMemUsed(size_t b);
    void decMemUsed(size_t b);
    void incMemAlloc(size_t b);
    void decMemAlloc(size_t b);

    void zeroMemUsed();
    void zeroMemAlloc();

    /**
        Set base pointer used for fast getChunkAddress().
    */
    void setChunkBase(void* p);

    /*
        Set multiplier used for fast getChunkAddress().
        If zero, we call getSlowChunkAddress().
    */
    void setChunkMultiplier(unsigned int m);

  public:
    /**
        Return the name of the style that created us.
    */
    const char* getStyleName() const;

  private:
    /// Name of the style that invoked us
    const char* style_name;
    memstats &my_mem;

    /// Base pointer for getChunkAddress
    char* chunk_base;

    /// Handle multiplier for getChunkAddress; if zero must call virtual function
    unsigned int chunk_multiplier;
};

// ******************************************************************
// *                                                                *
// *              inlined memory_manager_style methods              *
// *                                                                *
// ******************************************************************

inline const char* MEDDLY::memory_manager_style::getName() const
{
  return name;
}

// ******************************************************************
// *                                                                *
// *                 inlined memory_manager methods                 *
// *                                                                *
// ******************************************************************

inline const char* MEDDLY::memory_manager::getStyleName() const
{
  return style_name;
}

inline void* MEDDLY::memory_manager::getChunkAddress(node_address h) const
{
  MEDDLY_DCASSERT(isValidHandle(h));

  return chunk_multiplier
    ?  chunk_base + chunk_multiplier * h
    :  slowChunkAddress(h);
}

inline void MEDDLY::memory_manager::incMemUsed(size_t b)
{
  my_mem.incMemUsed(b);
}

inline void MEDDLY::memory_manager::decMemUsed(size_t b)
{
  my_mem.decMemUsed(b);
}

inline void MEDDLY::memory_manager::incMemAlloc(size_t b)
{
  my_mem.incMemAlloc(b);
}

inline void MEDDLY::memory_manager::decMemAlloc(size_t b)
{
  my_mem.decMemAlloc(b);
}

inline void MEDDLY::memory_manager::zeroMemUsed()
{
  my_mem.zeroMemUsed();
}

inline void MEDDLY::memory_manager::zeroMemAlloc()
{
  my_mem.zeroMemAlloc();
}

inline void MEDDLY::memory_manager::setChunkBase(void* p)
{
  chunk_base = (char*) p;
}

inline void MEDDLY::memory_manager::setChunkMultiplier(unsigned int m)
{
  chunk_multiplier = m;
}



#endif
