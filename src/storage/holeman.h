
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


#ifndef HOLEMAN_H
#define HOLEMAN_H

#include "../defines.h"

namespace MEDDLY {
  class holeman;
};

/** Abstract base class for hole management.

    This class is responsible for allocating and
    recycling chunks of memory, for node storage.

    An active chunk of memory may be viewed as
    an array of node_handles where the first and last
    elements of the array are guaranteed to be non-negative:

        [non-negative]
        [     ??     ]
              ..
        [     ??     ]
        [non-negative]

    Also, note that in general, the hole manager does not know
    the size of the active chunks of memory, but may have a
    minimum size requirement.

    A recycled chunk of memory may also be viewed as
    an array of node_handles where the first and last
    elements of the array are guaranteed to be negative:

        [  negative  ]
        [     ??     ]
              ..
        [     ??     ]
        [  negative  ]

    Some of the other slots may be used to manage holes,
    and the hole manager can determine the size of each hole.
*/
class MEDDLY::holeman {
  public:
    holeman(int smallestHole, node_storage* p);
    virtual ~holeman();

  public:
    inline node_storage* getParent() const { 
      MEDDLY_DCASSERT(parent);
      return parent; 
    }

    inline forest* getForest() const {
      MEDDLY_DCASSERT(parent);
      return parent->getParent();
    }

    inline int numHoles() const {
      return num_holes;
    }

    inline node_address holeSlots() const {
      return hole_slots;
    }

    inline node_address fragmentSlots() const {
      return fragment_slots;
    }

    inline node_address lastSlot() const {
      return last_slot;
    };

    inline int smallestChunk() const {
      return smallest;
    }

    inline void releaseFragment(int slots) {
      fragment_slots -= slots;
    }

  public:
    /**
        Request a contiguous chunk of memory.

          @param  slots   Number of node_handles we want

          @return         Starting index of the chunk.
                          The first slot of the chunk
                          will be negative, specifically
                          -(sizeof chunk obtained).
                          The returned chunk may be larger
                          than what was requested.
    */
    virtual node_address requestChunk(int slots) = 0;


    /**
        Recycle a chunk of memory we obtained earlier.

          @param  addr    Starting index of the chunk
          @param  slots   Size of the chunk, as number of node_handles
    */
    virtual void recycleChunk(node_address addr, int slots) = 0;


    /**
        Get the address of the chunk right after this hole.
        Behavior is undefined if the address does not refer
        to a hole (normally - throw exception or assertion violation)

          @param  addr    Starting index of the hole.
    */
    virtual node_address chunkAfterHole(node_address addr) const;

    /**
        Dump information for debugging.

            @param  s   Output stream
    */
    virtual void dumpInternalInfo(FILE* s) const = 0;

    /**
        Dump the interesting contents of a hole.
        Useful for debugging.

            @param  s   Output stream
            @param  a   Starting index of the chunk that's a hole
    */
    virtual void dumpHole(FILE* s, node_address a) const = 0;

    /**
        Dump information for debugging.
        Called after we dump all the nodes and holes.

            @param  s   Output stream
    */
    virtual void dumpInternalTail(FILE* s) const = 0;

    /** 
        Print stats.
            @param  s       Output stream
            @param  pad     String printed at the start of each line
            @param  flats   Controls what is displayed

        Children can call this to save some code.
    */
    virtual void reportStats(FILE* s, const char* pad, unsigned flags) const;

    /**
        Clear hole data structure and maybe shrink.
        Called after garbage collection in our parent.

        Derived classes should call the parent's version.
    */
    virtual void clearHolesAndShrink(node_address new_last, bool shrink);
  protected:
    inline void incMemUsed(long delta) {
      MEDDLY_DCASSERT(parent);
      parent->incMemUsed(delta);
    }
    inline void decMemUsed(long delta) {
      MEDDLY_DCASSERT(parent);
      parent->decMemUsed(delta);
    }
    inline void incMemAlloc(long delta) {
      MEDDLY_DCASSERT(parent);
      parent->incMemAlloc(delta);
    }
    inline void decMemAlloc(long delta) {
      MEDDLY_DCASSERT(parent);
      parent->decMemAlloc(delta);
    }

    inline void newHole(int slots) {
      num_holes++;
      if (num_holes > max_holes) {
        max_holes = num_holes;
      }
      hole_slots += slots;
      if (hole_slots > max_hole_slots) {
        max_hole_slots = hole_slots;
      }
    }
    inline void useHole(int slots) {
      num_holes--;
      hole_slots -= slots;
    }
    
    inline void incFragments(int slots) {
      fragment_slots += slots;
      if (fragment_slots > max_fragment_slots) {
        max_fragment_slots = fragment_slots;
      }
    }
    inline void decFragments(int slots) {
      fragment_slots -= slots;
    }

    inline node_address allocFromEnd(int slots) {
#ifdef MEMORY_TRACE
      printf("No hole available\n");
#endif

      //
      // First -- try to compact if we need to expand
      //
      if (getForest()->getPolicies().compactBeforeExpand) {
        if (last_slot + slots >= size) getParent()->collectGarbage(false);
      }

      //
      // Do we need to expand?
      //
      if (last_slot + slots >= size) {
        // new size is 50% more than previous 
        node_address want_size = last_slot + slots;
        node_handle new_size = MAX(size, want_size) * 1.5;  // TBD: WTF?

        resize(new_size);
      }

      // 
      // Grab node from the end
      //
      node_handle h = last_slot + 1;
      last_slot += slots;
      data[h] = -slots;
      return h;
    }

    inline void releaseToEnd(node_address h, int slots) {
      last_slot -= slots;
      useHole(slots);
      if (size > min_size && (last_slot + 1) < size/2) {
        node_handle new_size = size/2;
        while (new_size > (last_slot + 1) * 2) new_size /= 2;
        if (new_size < min_size) new_size = min_size;
        resize(new_size);
      }
#ifdef MEMORY_TRACE
      printf("Made Last Hole %ld, last %ld\n", long(addr), long(last_slot));
#ifdef DEEP_MEMORY_TRACE
      getParent()->dumpInternal(stdout);
#endif
#endif
    }

  protected:
    /// data array; we manage this
    node_handle* data;

  private:
    void resize(node_address new_slots);

  private:
    /// Size of data array.
    node_address size;
    /// Last used data slot
    node_address last_slot;
    
    node_storage* parent;
    int smallest;

    /// Total number of holes
    int num_holes;
    /// Max number of holes
    int max_holes;
    
    /// Total slots wasted in holes
    node_address hole_slots;
    /// Maximum hole slots seen
    node_address max_hole_slots;

    /// Total slots wasted in fragments 
    node_address fragment_slots;
    /// Maximum fragment slots seen
    node_address max_fragment_slots;

  private:
    /// Don't shrink below this
    static const int min_size = 1024;
};

#endif

