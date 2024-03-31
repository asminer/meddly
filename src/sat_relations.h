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

#ifndef MEDDLY_SAT_RELATIONS_H
#define MEDDLY_SAT_RELATIONS_H

#include "opname.h"
#include "dd_edge.h"
#include "error.h"
#include "forest.h"

#include <unordered_map>
#include <map>
#include <set>

namespace MEDDLY {
    class pregen_relation;
};

//
// TBD: still need to combine these into one general-purpose relation
//

// ******************************************************************
// *                                                                *
// *                     pregen_relation  class                     *
// *                                                                *
// ******************************************************************

/** Class for a partitioned transition relation, already known
    The relation can be partitioned "by events" or "by levels".
    In the case of "by events", we can have more than one relation
    per level; otherwise, there is at most one relation per level.
*/
class MEDDLY::pregen_relation {
    public:
        /** Constructor, by events
              @param  mxd         MxD forest containing relations
              @param  num_events  Number of events; specifies the maximum
                                  number of calls to addToRelation().
        */
        pregen_relation(forest* mxd, unsigned num_events);

        /** Constructor, by levels
              @param  mxd         MxD forest containing relations
        */
        pregen_relation(forest* mxd);

        virtual ~pregen_relation();
        void addToRelation(const dd_edge &r);

        // Options for controlling the amount of processing performed by
        // \a finalize(splittingOption).
        enum splittingOption {
          // None.
          None,
          // Transitions from level K that do not effect level K,
          // are moved to a lower level.
          SplitOnly,
          // SplitOnly + duplicated transitions between adjacent levels
          // are removed from the higher level.
          SplitSubtract,
          // SplitOnly + all duplicate transitions are removed.
          SplitSubtractAll,
          // Same as SplitSubtractAll, but using an algorithm that
          // first combines all transitions before splitting it up per level.
          MonolithicSplit
        };

        /** To be called after all events have been added to
            the transition relation.
            This method modifies the decision diagrams stored at different
            levels, to reduce duplicated transitions.
              @param  split       This parameter only applies to "by levels",
                                  and it controls the amount of processing
                                  that is performed.
                                  Please refer to splittingOption for details.
        */
        void finalize(splittingOption split = SplitSubtract);

        inline bool isFinalized() const { return nullptr == next; }

        inline forest* getRelForest() const { return mxdF; }

        // the following methods assume the relation has been finalized.
        inline dd_edge* arrayForLevel(int k) const
        {
            MEDDLY_DCASSERT(isFinalized());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1u, (unsigned)k, K + 1);
            if (level_index) {
                // "by events"
                if (level_index[k - 1] > level_index[k]) {
                    return events + level_index[k];
                } else {
                    // empty list
                    return nullptr;
                }
            } else {
                // "by levels"
                return events+k;
            }
        }

        inline unsigned lengthForLevel(int k) const
        {
            MEDDLY_DCASSERT(isFinalized());
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1u, (unsigned)k, K + 1);
            if (level_index) {
                // "by events"
                return level_index[k - 1] - level_index[k];
            } else {
                // "by levels"
                return events[k].getNode() ? 1 : 0;
            }
        }


    private:
        // helper for finalize,
        // find intersection of diagonals of events[k],
        // subtracts the intersection of events[k] and adds it to events[k-1].
        void splitMxd(splittingOption split);
        // helper for finalize
        // adds all event[k]; sets all event[k] to 0;
        // sets events[level(sum)] = sum
        void unionLevels();

        forest* mxdF;
        unsigned K;
        // array of sub-relations
        dd_edge* events;
        // next pointers (plus one), unless we're finalized
        unsigned* next;
        // size of events array
        unsigned num_events;
        // one past last used element of events array
        unsigned last_event;

        // If null, then we are "by levels".  Otherwise, we are "by events",
        // and before we're finalized, level_index[k] "points" (index plus one)
        // to a linked-list of sub-relations that affect level k.
        // After we're finalized, the events array is sorted, so
        // level_index[k] is the (actual) index of the first event affecting level k.
        // Dimension is number of variables + 1.
        unsigned* level_index;
};



#endif // #include guard
