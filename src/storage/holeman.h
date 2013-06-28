
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

    Let's see if this works :^)

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

    inline node_address holeSlots() const {
      return hole_slots;
    }

    inline node_address lastSlot() const {
      return last_slot;
    };

    inline int smallestChunk() const {
      return smallest;
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
    virtual node_address chunkAfterHole(node_address addr) const = 0;

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
        Report memory usage.

            @param  s     Output stream
            @param  pad   String printed at the start of each line
            @param  vL    verbosity level; higher means print more

    */
    virtual void reportMemoryUsage(FILE* s, const char* pad, int vL) const = 0;

    /**
        Clear hole data structure and maybe shrink.
        Called after garbage collection in our parent.
    */
    virtual void clearHolesAndShrink(node_address new_last, bool shrink) = 0;
  protected:
    inline void updateData(node_handle* d) {
      MEDDLY_DCASSERT(parent);
      parent->updateData(d);
    }
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

    /// Total slots wasted in holes
    node_address hole_slots;
    /// Last used data slot
    node_address last_slot;
  private:
    node_storage* parent;
    int smallest;
};

#endif

