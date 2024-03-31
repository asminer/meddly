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

    class otf_subevent;
    class otf_event;
    class otf_relation;
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


// ******************************************************************
// *                                                                *
// *                       otf_subevent class                       *
// *                                                                *
// ******************************************************************

/**
    User must derive a subclass from this.
    Part of an enabling or updating function.
    It knows what variables it depends on, and how to build itself
    (provided by the user).
*/
class MEDDLY::otf_subevent {
    public:
        /// Constructor, specify variables that this function depends on,
        /// and if it is a firing or enabling event.
        otf_subevent(forest* f, int* v, int nv, bool firing);
        virtual ~otf_subevent();

        /// Get the forest to which this function belongs to.
        inline forest* getForest() const { return f; }

        /// Get number of variables this function depends on.
        inline int getNumVars() const { return num_vars; }

        /// Get array of variables this function depends on.
        inline const int* getVars() const { return vars; }

        /// Get the DD encoding of this function
        inline const dd_edge& getRoot() const { return root; }

        /// Get the "top" variable for this function
        inline int getTop() const { return top; }

        /// Is this a firing subevent?
        inline bool isFiring() const { return is_firing; }

        /// Is this an enabling subevent
        inline bool isEnabling() const { return !is_firing; }

        /**
          Rebuild the function to include the
          local state "index" for the variable "v".
          Updates root with the updated function.
          User MUST provide this method.
        */
        virtual void confirm(otf_relation &rel, int v, int index) = 0;

        /// If num_minterms > 0,
        ///   Add all minterms to the root
        ///   Delete all minterms.
        void buildRoot();

        /// Debugging info
        void showInfo(output& out) const;

        long mintermMemoryUsage() const;
        void clearMinterms();

    protected:
        bool addMinterm(const int* from, const int* to);
        inline bool usesExtensibleVariables() const {
            return uses_extensible_variables;
        }

        int* vars;
        int num_vars;
        dd_edge root;
        int top;
        forest* f;
        int** unpminterms;
        int** pminterms;
        int num_minterms;
        int size_minterms;
        bool is_firing;
        bool uses_extensible_variables;

};  // end of class otf_subevent


// ******************************************************************
// *                                                                *
// *                        otf_event  class                        *
// *                                                                *
// ******************************************************************

/**
    An "event".
    Produces part of the transition relation, from its sub-functions.

    TBD - do we need to split the enabling and updating sub-functions,
    or will one giant list work fine?
*/
class MEDDLY::otf_event {
        // TBD - put a list of events that have priority over this one

        // TBD - for priority - when is this event enabled?
    public:
        otf_event(otf_subevent** se, int nse);
        virtual ~otf_event();

        /// Get the forest to which the subevents belong to
        inline forest* getForest() { return f; }

        /// Get number of subevents
        inline int getNumOfSubevents() const { return num_subevents; }

        /// Get array of subevents
        inline otf_subevent** getSubevents() const { return subevents; }

        /// Get the "top" variable for this event
        inline int getTop() const { return top; }

        /// Get the number of variables that are effected by this event
        inline int getNumVars() const { return num_vars; }

        /// Get a (sorted) array of variables that are effected by this event
        inline const int* getVars() const { return vars; }

        inline const dd_edge& getRoot() const { return root; }

        inline bool isDisabled() const { return is_disabled; }

        inline bool needsRebuilding() const { return needs_rebuilding; }

        inline void markForRebuilding() { needs_rebuilding = true; }

        /**
            If this event has been marked for rebuilding:
              Build this event as a conjunction of its sub-events.

            @return               true, if the event needed rebuilding and
                                        the rebuilding modified the root.
                                  false, otherwise.
        */
        virtual bool rebuild();

        /// Enlarges the "from" variable to be the same size as the "to" variable
        void enlargeVariables();

        /// Debugging info
        void showInfo(output& out) const;

        long mintermMemoryUsage() const;

    protected:
        void buildEventMask();

    private:
        otf_subevent** subevents;
        int num_subevents;
        int top;
        int num_vars;
        int* vars;
        dd_edge root;
        bool needs_rebuilding;
        forest* f;

        bool is_disabled;
        int num_firing_vars;
        int* firing_vars;
        dd_edge event_mask;
        int* event_mask_from_minterm;
        int* event_mask_to_minterm;

};  // end of class event

// ******************************************************************
// *                                                                *
// *                       otf_relation class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::otf_relation {
    public:
        /** Constructor.
              @param  mxd         MxD forest containing relations
              @param  mdd         MDD forest containing result set
              @param  E           List of events
              @param  nE          Number of events
        */
        otf_relation(forest* mxd, forest* mdd, otf_event** E, int ne);

        virtual ~otf_relation();

        /// Returns the MXD forest that stores the events
        inline forest* getRelForest() const { return mxdF; }

        /// Returns the MDD forest that stores the result
        inline forest* getOutForest() const { return resF; }

        /// Returns true if the local state is already confirmed.
        inline bool isConfirmed(int level, int i) const
        {
            if (level < num_levels &&  i >= 0) {
                return (i < resF->getLevelSize(level) && confirmed[level][i]);
            }
            throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
        }

        /// Returns an array of local states for this level, such that
        /// result[i] == isConfirmed(level, i).
        inline const bool* getLocalStates(int level) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, level, num_levels);
            return confirmed[level];
        }

        /// Returns the number of confirmed states at this level
        inline int getNumConfirmed(int level) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, level, num_levels);
            return num_confirmed[level];
        }

        /// Confirms all local states enabled in the given MDD
        void confirm(const dd_edge& set);

        /** Confirm a variable's previously unconfirmed state.
            Any event that is dependent on this variable is marked
            as "stale" --- so that it is rebuilt before use.

            @param  level       variable's level
            @param  index       the state of the variable being confirmed.
            @return             false: if state was previously confirmed.
                                true: if state was previously unconfirmed.
         */
        bool confirm(int level, int index);

        /** Get the number of events at whose "top" is this level.

            @param  level       level for the events.
            @return             number of events whose "top" is this level.
         */
        inline int getNumOfEvents(int level) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, level, num_levels);
            return num_events_by_top_level[level];
        }

        /** Gets an event from the set of events whose "top" is this level.

            @param  level       level for the events.
            @param  i           index of the event.
            @return             if 0 <= i < getNumOfEvents(level),
                                the ith event at this level;
                                otherwise, 0.
         */
        inline const dd_edge& getEvent(int level, int i) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, getNumOfEvents(level));
            return events_by_top_level[level][i]->getRoot();
        }

        /** Rebuild an event.

            @param  i           index of the event.
            @return             true, if event was updated.
          */
        inline bool rebuildEvent(int level, int i) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, getNumOfEvents(level));
            return events_by_top_level[level][i]->rebuild();
        }

        /** Build a Monolithic Next State Function that is equivalent to
            the union of all events while restricting the size of each
            variable to that of the largest confirmed index.
        */
        void getBoundedMonolithicNSF(dd_edge &root) const;

        /** Bound all extensible variables
            using the maximum confirmed local state as the bound.
        */
        void bindExtensibleVariables();

        /** Get the number of arcs in the OTF relation
            restricted by the confirmed local states.

            Only works with non-extensible variables (call
            bindExtensibleVariables() prior to calling this).

            @param    count_duplicates  if false, counts arcs that are common
                                          to mutiple transitions as one.
            @return                     the number of arcs in the OTF relation.
        */
        double getArcCount(const dd_edge& mask, bool count_duplicates);

        /// For Debugging
        void showInfo(output &strm) const;

        long mintermMemoryUsage() const;

        void clearMinterms();

    protected:
        void findConfirmedStates(bool** confirmed, int* num_confirmed,
            node_handle mdd, int level, std::set<MEDDLY::node_handle>& visited);

        void enlargeConfirmedArrays(int level, int sz);
        // node_handle getBoundedMxd(node_handle mxd, const int* bounds_array, int sz,
            // std::unordered_map<node_handle, node_handle>& cache);

    private:
        forest* mxdF;
        forest* resF;
        int num_levels;

        // All events that begin at level i,
        // are listed in events_by_top_level[i].
        // An event will appear in only one list
        // (as there is only one top level per event).
        // num_events_by_top_level[i] gives the size of events_by_top_level[i]
        otf_event*** events_by_top_level;
        int *num_events_by_top_level;

        // All events that depend on a level i,
        // are listed in events_by_level[i]
        // Therefore, an event that depends on n levels,
        // will appear in n lists
        // num_events_by_level[i] gives the size of events_by_level[i]
        otf_event*** events_by_level;
        int *num_events_by_level;

        // All subevents that depend on a level i,
        // are listed in subevents_by_level[i]
        // Therefore, an subevent that depends on n levels,
        // will appear in n lists
        // num_subevents_by_level[i] gives the size of subevents_by_level[i]
        otf_subevent*** subevents_by_level;
        int *num_subevents_by_level;

        // List of confirmed local states at each level
        bool** confirmed;
        int* size_confirmed;
        int* num_confirmed;

};  // end of class otf_relation

#endif // #include guard
