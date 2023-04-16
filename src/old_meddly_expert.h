
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


/*! \file meddly_expert.h

    Low-level MDD library interface.

    This interface is for "expert" users who want to define new
    operations, or for library developers to define the built-in operations.
    Casual users probably only need the interface provided by "meddly.h".

    The first part of the interface describes the expert interface and the
    second part contains implementations of virtual functions in the interface.

    IMPORTANT: meddly.h must be included before including this file.
    TODO: Operations are not thread-safe.
*/

#ifndef MEDDLY_EXPERT_H
#define MEDDLY_EXPERT_H

#include <string.h>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <map>
#include <set>

#include "initializer.h"

// #define DEBUG_MARK_SWEEP
// #define DEBUG_BUILDLIST

// #define KEEP_LL_COMPUTES

namespace MEDDLY {

  // classes defined here

  class expert_variable;
  class expert_domain;

  // wrapper for temporary nodes
  class unpacked_node;  // replacement for node_reader, node_builder

  class relation_node;

  /*

    class op_initializer;

    Generalized to class initializer_list.
  */

  class initializer_list;

  /*
    class cleanup_procedure;

    Subsumed by class initializer_list.
  */

  // Memory managers, for node storage and compute tables
  class memory_manager_style;
  class memory_manager;

  // Node header storage
  class node_headers;

  // Actual node storage
  class node_storage_style;
  class node_storage;

  class expert_forest;

  class opname;
  class unary_opname;
  class binary_opname;
  class specialized_opname;
  class numerical_opname;
  class satpregen_opname;
  class satotf_opname;
  class satimpl_opname;
  class sathyb_opname;
  class constrained_opname;

  class ct_initializer;
  class compute_table_style;
  class compute_table;

  class operation;
  class unary_operation;
  class binary_operation;
  class specialized_operation;

  class global_rebuilder;

  // classes defined elsewhere
  class base_table;
  class unique_table;
  class impl_unique_table;

  class reordering_base;

  // ******************************************************************
  // *                                                                *
  // *                   Named numerical operations                   *
  // *                                                                *
  // ******************************************************************

  /** Computes y = y + xA.
      x and y are vectors, stored explicitly, and A is a matrix.
      x_ind and y_ind specify how minterms are mapped to indexes
      for vectors x and y, respectively.
  */
  extern const numerical_opname* EXPLVECT_MATR_MULT;
  // extern const numerical_opname* VECT_MATR_MULT; // renamed!

  /** Computes y = y + Ax.
      x and y are vectors, stored explicitly, and A is a matrix.
      x_ind and y_ind specify how minterms are mapped to indexes
      for vectors x and y, respectively.
  */
  extern const numerical_opname* MATR_EXPLVECT_MULT;
  // extern const numerical_opname* MATR_VECT_MULT; // renamed!

  // ******************************************************************
  // *                                                                *
  // *                  Named saturation operations                   *
  // *                                                                *
  // ******************************************************************

  /** Forward reachability using saturation.
      Transition relation is already known.
  */
  extern const satpregen_opname* SATURATION_FORWARD;

  /** Backward reachability using saturation.
      Transition relation is already known.
  */
  extern const satpregen_opname* SATURATION_BACKWARD;

  /** Forward reachability using saturation.
      Transition relation is not completely known,
      will be built along with reachability set.
  */
  extern const satotf_opname* SATURATION_OTF_FORWARD;

  /** Forward reachability using saturation.
      Transition relation is specified implicitly.
  */
  extern const satimpl_opname* SATURATION_IMPL_FORWARD;

  /** Forward reachability using saturation.
      Allows hybrid representation of transition relation.
  */
  extern const sathyb_opname* SATURATION_HYB_FORWARD;

  /** Minimum-witness operations.
  */
  extern const constrained_opname* CONSTRAINED_BACKWARD_BFS;
  extern const constrained_opname* CONSTRAINED_FORWARD_DFS;
  extern const constrained_opname* CONSTRAINED_BACKWARD_DFS;
  extern const constrained_opname* TRANSITIVE_CLOSURE_DFS;

  // ******************************************************************
  // *                                                                *
  // *                      Operation management                      *
  // *                                                                *
  // ******************************************************************

  /// Remove an existing operation from the operation cache.
  void removeOperationFromCache(operation* );

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest
        @param  res     Result forest
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    expert_forest* arg, expert_forest* res);

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest from this dd_edge
        @param  res     Result forest from this dd_edge
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    const dd_edge& arg, const dd_edge& res);

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest
        @param  res     Result type
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    expert_forest* arg, opnd_type result);

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest from this dd_edge
        @param  res     Result type
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    const dd_edge& arg, opnd_type result);


  /** Find, or build if necessary, a binary operation.
        @param  code    Operation we want
        @param  arg1    Argument 1 forest
        @param  arg2    Argument 2 forest
        @param  res     Result forest
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  binary_operation* getOperation(const binary_opname* code,
    expert_forest* arg1, expert_forest* arg2, expert_forest* res);

  /** Find, or build if necessary, a binary operation.
        @param  code    Operation we want
        @param  arg1    Argument 1 forest taken from this dd_edge
        @param  arg2    Argument 2 forest taken from this dd_edge
        @param  res     Result forest taken from this dd_edge
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  binary_operation* getOperation(const binary_opname* code,
    const dd_edge& arg1, const dd_edge& arg2, const dd_edge& res);

  /** Safely destroy the given unary operation.
      It should be unnecessary to call this directly.
  */
  void destroyOperation(unary_operation* &op);

  /** Safely destroy the given binary operation.
      It should be unnecessary to call this directly.
  */
  void destroyOperation(binary_operation* &op);

  /// Safely destroy the given numerical operation.
  void destroyOperation(specialized_operation* &op);

  /// Should not be called directly.
  void destroyOpInternal(operation* op);

  // ******************************************************************
  // *                                                                *
  // *                  library management functions                  *
  // *                                                                *
  // ******************************************************************

  /*
    /// Builds an initializer for MEDDLY's builtin operations.
    /// Use defaultInitializerList() instead
    op_initializer* makeBuiltinInitializer();

  */

}; // namespace MEDDLY


// ******************************************************************
// *                                                                *
// *                          opname class                          *
// *                                                                *
// ******************************************************************

/// Class for names of operations.
class MEDDLY::opname {
    const char* name;
    int index;
    static int next_index;

    friend void MEDDLY::initialize(initializer_list *);
    friend void MEDDLY::cleanup();
  public:
    opname(const char* n);
    virtual ~opname();

    int getIndex() const;
    const char* getName() const;
};

// ******************************************************************
// *                                                                *
// *                       unary_opname class                       *
// *                                                                *
// ******************************************************************

/// Unary operation names.
class MEDDLY::unary_opname : public opname {
public:
  unary_opname(const char* n);
  virtual ~unary_opname();

  virtual unary_operation*
  buildOperation(expert_forest* arg, expert_forest* res) const;

  virtual unary_operation*
  buildOperation(expert_forest* arg, opnd_type res) const;
};

// ******************************************************************
// *                                                                *
// *                      binary_opname  class                      *
// *                                                                *
// ******************************************************************

/// Binary operation names.
class MEDDLY::binary_opname : public opname {
  public:
    binary_opname(const char* n);
    virtual ~binary_opname();

    virtual binary_operation* buildOperation(expert_forest* arg1,
      expert_forest* arg2, expert_forest* res) const = 0;
};


// ******************************************************************
// *                                                                *
// *                    specialized_opname class                    *
// *                                                                *
// ******************************************************************

/// Specialized operation names.
class MEDDLY::specialized_opname : public opname {
  public:
    /**
      Abstract base class for arguments to buildOperation().
      Derived operation names must provide derived classes for arguments.
    */
    class arguments {
      public:
        arguments();
        virtual ~arguments();

        /**
            Specify if arguments should be destroyed or not.
            If yes (the default), the operation will destroy
            the arguments once they are no longer needed.
        */
        void setAutoDestroy(bool destroy);
        bool autoDestroy() const;

      private:
        bool destroyWhenDone;
    };

    specialized_opname(const char* n);
    virtual ~specialized_opname();

    /** Note - unlike the more general binary and unary ops,
        a specialized operation might be crafted to the specific
        arguments (passed as an abstract class).
        Examples are:
          - operations that will be called several times with
            several of the arguments unchanged, and we want
            to do some preprocessing on the unchanged arguments.

          - operations with very bizarre and/or user-defined
            parameters (like on-the-fly saturation).

        @param  a   Arguments.  Will be destroyed when we are finished,
                    if autoDestroy() is set for the arguments.
    */
    virtual specialized_operation* buildOperation(arguments* a) const = 0;
};

// ******************************************************************
// *                                                                *
// *                     numerical_opname class                     *
// *                                                                *
// ******************************************************************

/// Numerical operation names.
class MEDDLY::numerical_opname : public specialized_opname {
  public:
    class numerical_args : public specialized_opname::arguments {
      public:
        const dd_edge &x_ind;
        const dd_edge &A;
        const dd_edge &y_ind;

        numerical_args(const dd_edge &xi, const dd_edge &a, const dd_edge &yi);
        virtual ~numerical_args();
    };

    numerical_opname(const char* n);
    virtual ~numerical_opname();
    virtual specialized_operation* buildOperation(arguments* a) const = 0;

    /// For convenience, and backward compatability :^)
    specialized_operation* buildOperation(const dd_edge &x_ind,
      const dd_edge &A, const dd_edge &y_ind) const;
};


// ******************************************************************
// *                                                                *
// *                     satpregen_opname class                     *
// *                                                                *
// ******************************************************************

/** Saturation, with already generated transition relations, operation names.
    Implemented in operations/sat_pregen.cc
*/
class MEDDLY::satpregen_opname : public specialized_opname {
  public:
    satpregen_opname(const char* n);
    virtual ~satpregen_opname();

    /// Arguments should have type "pregen_relation".
    virtual specialized_operation* buildOperation(arguments* a) const = 0;

    /** Class for a partitioned transition relation, already known
        The relation can be partitioned "by events" or "by levels".
        In the case of "by events", we can have more than one relation
        per level; otherwise, there is at most one relation per level.
    */
    class pregen_relation : public specialized_opname::arguments {
      public:
        /** Constructor, by events
              @param  inmdd       MDD forest containing initial states
              @param  mxd         MxD forest containing relations
              @param  outmdd      MDD forest containing result
              @param  num_events  Number of events; specifies the maximum
                                  number of calls to addToRelation().
        */
        pregen_relation(forest* inmdd, forest* mxd, forest* outmdd,
          unsigned num_events);

        /** Constructor, by levels
              @param  inmdd       MDD forest containing initial states
              @param  mxd         MxD forest containing relations
              @param  outmdd      MDD forest containing result
        */
        pregen_relation(forest* inmdd, forest* mxd, forest* outmdd);

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

        bool isFinalized() const;

        forest* getInForest() const;
        forest* getRelForest() const;
        forest* getOutForest() const;

        // the following methods assume the relation has been finalized.
        dd_edge* arrayForLevel(int k) const;
        unsigned lengthForLevel(int k) const;

      private:
        // helper for constructors
        void setForests(forest* inf, forest* mxd, forest* outf);
        // helper for finalize,
        // find intersection of diagonals of events[k],
        // subtracts the intersection of events[k] and adds it to events[k-1].
        void splitMxd(splittingOption split);
        // helper for finalize
        // adds all event[k]; sets all event[k] to 0;
        // sets events[level(sum)] = sum
        void unionLevels();

        forest* insetF;
        expert_forest* mxdF;
        forest* outsetF;
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

};



// ******************************************************************
// *                                                                *
// *                      satotf_opname  class                      *
// *                                                                *
// ******************************************************************

#include "dd_edge.h"

/** Saturation, transition relations built on the fly, operation names.
    Implemented in operations/sat_otf.cc
*/
class MEDDLY::satotf_opname : public specialized_opname {
  public:
    satotf_opname(const char* n);
    virtual ~satotf_opname();

    /// Arguments should have type "otf_relation", below
    virtual specialized_operation* buildOperation(arguments* a) const = 0;

    class otf_relation;

    // ============================================================

    /**
        User must derive a subclass from this.
        Part of an enabling or updating function.
        It knows what variables it depends on, and how to build itself
        (provided by the user).
    */
    class subevent {
      public:
        /// Constructor, specify variables that this function depends on,
        /// and if it is a firing or enabling event.
        subevent(forest* f, int* v, int nv, bool firing);
        virtual ~subevent();

        /// Get the forest to which this function belongs to.
        expert_forest* getForest();

        /// Get number of variables this function depends on.
        int getNumVars() const;

        /// Get array of variables this function depends on.
        const int* getVars() const;

        /// Get the DD encoding of this function
        const dd_edge& getRoot() const;

        /// Get the "top" variable for this function
        int getTop() const;

        /// Is this a firing subevent?
        bool isFiring() const;

        /// Is this an enabling subevent
        bool isEnabling() const;

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
        bool usesExtensibleVariables() const;

        int* vars;
        int num_vars;
        dd_edge root;
        int top;
        expert_forest* f;
        int** unpminterms;
        int** pminterms;
        int num_minterms;
        int size_minterms;
        bool is_firing;
        bool uses_extensible_variables;

    };  // end of class subevent

    // ============================================================

    /**
        An "event".
        Produces part of the transition relation, from its sub-functions.

        TBD - do we need to split the enabling and updating sub-functions,
        or will one giant list work fine?
    */
    class event {
        // TBD - put a list of events that have priority over this one

        // TBD - for priority - when is this event enabled?
      public:
        event(subevent** se, int nse);
        virtual ~event();

        /// Get the forest to which the subevents belong to
        inline expert_forest* getForest() { return f; }

        /// Get number of subevents
        inline int getNumOfSubevents() const { return num_subevents; }

        /// Get array of subevents
        inline subevent** getSubevents() const { return subevents; }

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
        subevent** subevents;
        int num_subevents;
        int top;
        int num_vars;
        int* vars;
        dd_edge root;
        bool needs_rebuilding;
        expert_forest* f;

        bool is_disabled;
        int num_firing_vars;
        int* firing_vars;
        dd_edge event_mask;
        int* event_mask_from_minterm;
        int* event_mask_to_minterm;

    };  // end of class event

    // ============================================================

    /**
        Overall relation.
        This includes all events, and keeping track of which local
        variables are confirmed.

        TBD.
    */
    class otf_relation : public specialized_opname::arguments {
      public:
        /** Constructor.
              @param  inmdd       MDD forest containing initial states
              @param  mxd         MxD forest containing relations
              @param  outmdd      MDD forest containing result
              @param  E           List of events
              @param  nE          Number of events
        */
        otf_relation(forest* inmdd, forest* mxd, forest* outmdd,
          event** E, int ne);

        virtual ~otf_relation();

        /// Returns the MDD forest that stores the initial set of states
        expert_forest* getInForest() const;

        /// Returns the MXD forest that stores the events
        expert_forest* getRelForest() const;

        /// Returns the MDD forest that stores the resultant set of states
        expert_forest* getOutForest() const;

        /// Returns true if the local state is already confirmed.
        bool isConfirmed(int level, int index) const;

        /// Returns an array of local states for this level, such that
        /// result[i] == isConfirmed(level, i).
        const bool* getLocalStates(int level);

        /// Returns the number of confirmed states at this level
        int getNumConfirmed(int level) const;

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
        int getNumOfEvents(int level) const;

        /** Gets an event from the set of events whose "top" is this level.

            @param  level       level for the events.
            @param  i           index of the event.
            @return             if 0 <= i < getNumOfEvents(level),
                                the ith event at this level;
                                otherwise, 0.
         */
        const dd_edge& getEvent(int level, int i) const;

        /** Rebuild an event.

            @param  i           index of the event.
            @return             true, if event was updated.
          */
        bool rebuildEvent(int level, int i);

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
        void enlargeConfirmedArrays(int level, int sz);
        // node_handle getBoundedMxd(node_handle mxd, const int* bounds_array, int sz,
            // std::unordered_map<node_handle, node_handle>& cache);

      private:
        expert_forest* insetF;
        expert_forest* mxdF;
        expert_forest* outsetF;
        int num_levels;

        // All events that begin at level i,
        // are listed in events_by_top_level[i].
        // An event will appear in only one list
        // (as there is only one top level per event).
        // num_events_by_top_level[i] gives the size of events_by_top_level[i]
        event*** events_by_top_level;
        int *num_events_by_top_level;

        // All events that depend on a level i,
        // are listed in events_by_level[i]
        // Therefore, an event that depends on n levels,
        // will appear in n lists
        // num_events_by_level[i] gives the size of events_by_level[i]
        event*** events_by_level;
        int *num_events_by_level;

        // All subevents that depend on a level i,
        // are listed in subevents_by_level[i]
        // Therefore, an subevent that depends on n levels,
        // will appear in n lists
        // num_subevents_by_level[i] gives the size of subevents_by_level[i]
        subevent*** subevents_by_level;
        int *num_subevents_by_level;

        // List of confirmed local states at each level
        bool** confirmed;
        int* size_confirmed;
        int* num_confirmed;

    };  // end of class otf_relation

};  // end of class satotf_opname

// ******************************************************************
// *                                                                *
// *                      satimpl_opname class                      *
// *                                                                *
// ******************************************************************

/** Saturation, transition relations stored implcitly, operation names.
    Implemented in operations/sat_impl.cc
*/
class MEDDLY::satimpl_opname: public specialized_opname {
  public:

    satimpl_opname(const char* n);
    virtual ~satimpl_opname();

    /// Arguments should have type "implicit_relation", below
    virtual specialized_operation* buildOperation(arguments* a) const;

  public:

    /** An implicit relation, as a DAG of relation_nodes.

        The relation is partitioned by "events", where each event
        is the conjunction of local functions, and each local function
        is specified as a single relation_node.  The relation_nodes
        are chained together with at most one relation_node per state
        variable, and any skipped variables are taken to be unchanged
        by the event.

        If the bottom portion (suffix) of two events are identical,
        then they are merged.  This is done by "registering" nodes
        which assigns a unique ID to each node, not unlike an MDD forest.

        Note: node handles 0 and 1 are reserved.
        0 means null node.
        1 means special bottom-level "terminal" node
        (in case we need to distinguish 0 and 1).
    */
    class implicit_relation : public specialized_opname::arguments {
      public:
        /** Constructor.

            @param  inmdd       MDD forest containing initial states
            @param  outmdd      MDD forest containing result

            Not 100% sure we need these...
        */
        implicit_relation(forest* inmdd, forest* relmxd, forest* outmdd);
        virtual ~implicit_relation();

        /// Returns the Relation forest that stores the mix of relation nodes and mxd nodes
        expert_forest* getMixRelForest() const;


        /// Returns the MDD forest that stores the initial set of states
        expert_forest* getInForest() const;


        /// Returns the MDD forest that stores the resultant set of states
        expert_forest* getOutForest() const;

        /// Returns true iff a state in \a constraint is reachable
        /// from the states in \a initial_states
        /// Note: \a constraint can be an XDD
        bool isReachable(const dd_edge& initial_states, const dd_edge& constraint);

        /** Register a relation node.

            If we have seen an equivalent node before, then
            return its handle and destroy n.
            Otherwise, add n to the unique table, assign it a unique
            identifier, and return that identifier.

            @param  is_event_top    If true, this is also the top
                                    node of some event; register it
                                    in the list of events.

            @param  n               The relation node to register.

            @return Unique identifier to use to refer to n.
        */
        rel_node_handle registerNode(bool is_event_top, relation_node* n);

        /** Check if the relation node is unique
            @param n  The relation node.
            @return   If unique, 0
                      Else, existing node handle
        */
        rel_node_handle isUniqueNode(relation_node* n);


        /** Indicate that there will be no more registered nodes.
            Allows us to preprocess the events or cleanup or convert
            to a more useful representation for saturation.
        */

        //void finalizeNodes();

        /** Get the relation node associated with the given handle.

            Should be a fast, inlined implementation.
        */
        relation_node* nodeExists(rel_node_handle n);

        /** Get the relation node associated with the given handle.

            Should be a fast, inlined implementation.
        */
        bool isReserved(rel_node_handle n);

      private:
        expert_forest* insetF;
        expert_forest* outsetF;
        expert_forest* mixRelF;

        int num_levels;

      private:

        /// Last used ID of \a relation node.
        long last_in_node_array;

      private:
        // TBD - add a data structure for the "uniqueness table"
        // of relation_nodes, so if we register a node that
        // is already present in a node_array, we can detect it.

        std::unordered_map<rel_node_handle, relation_node*> impl_unique;

      private:
        // TBD - add a data structure for list of events with top level k,
        // for all possible k.
        // Possibly this data structure is built by method
        // finalizeNodes().

        rel_node_handle** event_list;
        long* event_list_alloc; // allocated space
        long* event_added; //how many events added so far

        long* confirm_states; //total no. of confirmed states of a level
        bool** confirmed; // stores whether a particular local state is confirmed
        long* confirmed_array_size; // stores size of confirmed array


      public:

        /// Get total number of events upto given level
        long getTotalEvent(int level);

        /// Resizes the Event List
        void resizeEventArray(int level);

        /// Returns the number of events that have this level as top
        long lengthForLevel(int level) const;

        /// Returns the array of events that have this level as top
        rel_node_handle* arrayForLevel(int level) const;

        /// Returns the number of confirmed states at a level
        long getConfirmedStates(int level) const;

        /// Confirms the local states at a level
        void setConfirmedStates(int level, int i);

        /// Confirms the local states in the given MDD
        void setConfirmedStates(const dd_edge &set);


        /// Checks if i is confirmed
        bool isConfirmedState(int level, int i);

        /// Expand confirm array
        void resizeConfirmedArray(int level, int index);

        /** Bound all extensible variables
            using the maximum confirmed local state as the bound.
        */
        void bindExtensibleVariables();

      public:
        /// Prints the implicit relation
        void show();

        /// Build mxd forest
        MEDDLY::node_handle buildMxdForest();

        /// Build each event_mxd
        dd_edge buildEventMxd(rel_node_handle event_top, forest *mxd);

        /// Get relation forest
        expert_forest* getRelForest() const;


      private:
        expert_forest* mxdF;

     public:

      /*
       Group the list of events at a given level by same next-of values
       @param level    level at which saturation is called
       @param i        number of tokens before event is fired
       @param R        array of events at level top

       Return the map.
       */
      std::unordered_map<long,std::vector<rel_node_handle>> getListOfNexts(int level, long i, relation_node **R);

      /*
       Returns whether there exist a possibility of doing union
       @param level    level at which saturation is called
       @param i        number of tokens before event is fired
       @param R        array of events at level top

       Return bool.
       */
      bool isUnionPossible(int level, long i, relation_node **R);

    };  // class implicit_relation



};
class MEDDLY::constrained_opname : public specialized_opname {
public:
	constrained_opname(const char* n);

	class constrained_args : public specialized_opname::arguments {
	public:
	  forest* consForest;
	  forest* inForest;
	  forest* relForest;
	  forest* outForest;

	  constrained_args(forest* consF, forest* inF, forest* relF, forest* outF)
	    : consForest(consF), inForest(inF), relForest(relF), outForest(outF)
	  {
	  }
	};
};

// ******************************************************************
// *                                                                *
// *                      sathyb_opname class                      *
// *                                                                *
// ******************************************************************

/** Saturation, transition relations stored implcitly, operation names.
    Implemented in operations/sat_impl.cc
*/
class MEDDLY::sathyb_opname: public specialized_opname {
  public:

    sathyb_opname(const char* n);
    virtual ~sathyb_opname();

    /// Arguments should have type "implicit_relation", below
    virtual specialized_operation* buildOperation(arguments* a) const;

  public:

    /** A hybrid relation, as a DAG of relation_nodes and MxD nodes.

        The relation is partitioned by "events", where each event
        is the "symbolic" conjunction of local functions, and each local function
        might be specified as a single relation_node or an MxD forest based on the cardinality of dependent variables of the local function.
        The relation_nodes are chained together.

        Ideally it behaves as a combination of MxD forest and implicit relation : covering the best of both worlds.
    */

    class hybrid_relation;


    class subevent { //: public semievent {
    public:
      /// Constructor, specify variables that this function depends on,
      /// and if it is a firing or enabling event.
      subevent(forest* f, int* v, int nv, bool firing);
      virtual ~subevent();

      //semievent* clone();

      /// Get the forest to which this function belongs to.
      expert_forest* getForest();

      /// Get number of variables this function depends on.
      int getNumVars() const; // Returns 1 for implicit node

      /// Get array of variables this function depends on.
      const int* getVars() const; // Returns only 1 value

      /// Get the DD encoding of this function
      const dd_edge& getRoot() const; // Non-existent for implicit node

      /// Get the node_handle of this function
      const node_handle getRootHandle() const; // getID() for implicit node

      /// Get the "top" variable for this function
      int getTop() const; // getLevel() for implicit node

      /// Is this a firing subevent?
      bool isFiring() const;

      /// Is this an enabling subevent
      bool isEnabling() const;

      /// Is this an implicit subevent
      bool isImplicit() const;

      /** Applicable if implicit node:
       A signature for this function.
       This helps class implicit_relation to detect duplicate
       functions (using a hash table where the signature
       is taken as the hash value).
       */
      //unsigned long getSignature() const;

      /** Applicable if implicit node :
       Get the ID of next piece, default = handleForValue(True).
       */
      node_handle getDown() const;

      /** Applicable if implicit node :
       Set the ID of next piece.
       */
      void setDown(node_handle down) ;

      /** Applicable if implicit node :
       Set the unique ID for this piece.
       */
      void setRootHandle(node_handle ID);

      //
      int getEnable() const;
      int getFire() const;

      /** Determine if this node is equal to another one.
       */
      //virtual bool equals(const relation_node* n) const;

      /**
       Rebuild the function to include the
       local state "index" for the variable "v".
       Updates root with the updated function.
       User MUST provide this method.
       */
      virtual void confirm(hybrid_relation &rel, int v, int index) = 0;

      /// If num_minterms > 0,
      /// Add all minterms to the root
      /// Delete all minterms.
      void buildRoot();

      /// Debugging info
      void showInfo(output& out) const;

      long mintermMemoryUsage() const;
      void clearMinterms();

    protected:
      bool addMinterm(const int* from, const int* to);
      bool usesExtensibleVariables() const;

      int* vars;
      int num_vars;
      dd_edge root; // NULL for implicit node
      node_handle root_handle;
      int top; // top = vars[0] for implicit node
      node_handle down;


      expert_forest* f;
      int** unpminterms; // unpminterms[0] for implicit node
      long enable;
      int** pminterms; // pminterms[0] for implicit node
      long fire;
      int num_minterms;
      int process_minterm_pos;
      int processed_minterm_pos;
      int size_minterms;
      bool is_firing;
      bool uses_extensible_variables;

    };  // end of class subevent

    // ============================================================

    /**
     An "event".
     Produces part of the transition relation, from its sub-functions.

     TBD - do we need to split the enabling and updating sub-functions,
     or will one giant list work fine?
     */
    class event {
      // TBD - put a list of events that have priority over this one

      // TBD - for priority - when is this event enabled?
    public:

      event(subevent** se, int nse, relation_node** r, int nr);
      virtual ~event();

      /// Get the forest to which the subevents belong to
      inline expert_forest* getForest() { return f; }

      /// Get number of subevents
      inline int getNumOfSubevents() const { return num_subevents; }

      /// Get array of subevents
      inline subevent** getSubevents() const { return subevents; }

      /// Get number of relation_nodes
      inline int getNumOfRelnodes() const { return num_relnodes; }

      /// Get array of relation_nodes
      inline relation_node** getRelNodes() const { return relnodes; }

      /// Get number of components
      inline int getNumOfComponents() const { return num_components; }

      ///Get the subevent handle or relation node handle whose top level is the given level
      inline std::set<node_handle> getComponentAt(int level)  { return level_component.find(level)->second; }

      inline node_handle getTopComponent()  { return *level_component[top].begin(); }

      ///Get the entire map of subevent-nodeHandles_by_topLevel
      inline std::map<int, std::set<node_handle>> getComponents()  { return level_component; }

      ///Get the entire map of subevent-nodeHandles_by_topLevel
      inline node_handle* getAllComponents()  {  return all_components; }

      ///Get the type of subevent-nodeHandles
      inline bool getSubeventType(node_handle nh)  { return component_se_type[nh]; }


      /// Get the "top" variable for this event
      inline int getTop() const { return top; }

      /// Get the number of variables that are effected by this event
      inline int getNumVars() const { return num_vars; }

      /// Get a (sorted) array of variables that are effected by this event
      inline const int* getVars() const { return vars; }

      inline const dd_edge& getRoot() const { return root; }

      inline const std::set<node_handle> getRootHandle() const { return root_handle; }

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

      /// Get down of a level of event
      int downLevel(int level) const;

      /// Enlarges the "from" variable to be the same size as the "to" variable
      void enlargeVariables();

      /// Debugging info
      void showInfo(output& out) const;

      long mintermMemoryUsage() const;

    protected:
      void buildEventMask();

    private:
      subevent** subevents;
      int num_subevents;
      relation_node** relnodes;
      int num_relnodes;
      int num_components;
      int top;
      int num_vars;
      int* vars;

      // set only if sub-events are conjuncted
      dd_edge partial_root;

      // set because multiple root_handles may exist if sub-events & relNodes are not conjuncted
      std::set<node_handle> root_handle;
      bool needs_rebuilding;
      expert_forest* f;
      dd_edge root;
      bool is_disabled;
      int num_firing_vars;
      int* firing_vars;
      dd_edge event_mask;
      int* event_mask_from_minterm;
      int* event_mask_to_minterm;
      bool first_time_build;


      int num_rel_vars;
      int* relNode_vars;

      std::map<int,std::set<node_handle>> level_component; // stores the set of subevent node_handles whose top is this level.
      std::map<node_handle,bool> component_se_type; //enabling:0 firing/rn:1
      node_handle* all_components; // set of subevent's top node_handles

    };  // end of class event



    /** An implicit relation, as a DAG of relation_nodes.

        The relation is partitioned by "events", where each event
        is the conjunction of local functions, and each local function
        is specified as a single relation_node.  The relation_nodes
        are chained together with at most one relation_node per state
        variable, and any skipped variables are taken to be unchanged
        by the event.

        If the bottom portion (suffix) of two events are identical,
        then they are merged.  This is done by "registering" nodes
        which assigns a unique ID to each node, not unlike an MDD forest.

        Note: node handles 0 and 1 are reserved.
        0 means null node.
        1 means special bottom-level "terminal" node
        (in case we need to distinguish 0 and 1).
    */
    class hybrid_relation : public specialized_opname::arguments {
      public:
        /** Constructor.

            @param  inmdd       MDD forest containing initial states
            @param  outmdd      MDD forest containing result

            Not 100% sure we need these...
        */
        hybrid_relation(forest* inmdd, forest* relmxd, forest* outmdd, event** E, int ne);
        virtual ~hybrid_relation();

        /// Returns the Relation forest that stores the mix of relation nodes and mxd nodes
        expert_forest* getHybridForest() const;


        /// Returns the MDD forest that stores the initial set of states
        expert_forest* getInForest() const;


        /// Returns the MDD forest that stores the resultant set of states
        expert_forest* getOutForest() const;


        /// If only relation_nodes are present
        /// return the vector of events pertaining to a given level
        /// that have the same net effect
        std::vector<node_handle> getRelNodeAtLevelWithEffect(int level, long effect);

      private:
        expert_forest* insetF;
        expert_forest* outsetF;
        expert_forest* hybRelF;

        int num_levels;

      private:

        /// Last used ID of \a relation node.
        long last_in_node_array;

      private:
        // TBD - add a data structure for the "uniqueness table"
        // of relation_nodes, so if we register a node that
        // is already present in a node_array, we can detect it.

        std::unordered_map<node_handle, relation_node*> impl_unique;

      private:
        // TBD - add a data structure for list of events with top level k,
        // for all possible k.
        // Possibly this data structure is built by method
        // finalizeNodes().

        node_handle** event_list;
        long* event_list_alloc; // allocated space
        long* event_added; //how many events added so far

        int* confirm_states; //total no. of confirmed states of a level
        bool** confirmed; // stores whether a particular local state is confirmed
        int* confirmed_array_size; // stores size of confirmed array


        // Obtained from OTF :

        // All events that begin at level i,
        // are listed in events_by_top_level[i].
        // An event will appear in only one list
        // (as there is only one top level per event).
        // num_events_by_top_level[i] gives the size of events_by_top_level[i]
        event*** events_by_top_level;
        int *num_events_by_top_level;

        // All events that depend on a level i,
        // are listed in events_by_level[i]
        // Therefore, an event that depends on n levels,
        // will appear in n lists
        // num_events_by_level[i] gives the size of events_by_level[i]
        event*** events_by_level;
        int *num_events_by_level;

        // All subevents that depend on a level i,
        // are listed in subevents_by_level[i]
        // Therefore, a subevent that depends on n levels,
        // will appear in n lists
        // num_subevents_by_level[i] gives the size of subevents_by_level[i]
        subevent*** subevents_by_level;
        int *num_subevents_by_level;

        // All relation_nodes that depend on a level i,
        // are listed in relnodes_by_level[i]
        // A relnode can only appear in one list
        // num_relnodes_by_level[i] gives the size of relnodes_by_level[i]
        relation_node*** relnodes_by_level;
        int *num_relnodes_by_level;

      public:

        /// Get total number of variables
        inline int getNumVariables() { return num_levels;}

        /** Rebuild an event.

         @param  i           index of the event.
         @return             true, if event was updated.
         */
        bool rebuildEvent(int level, int i);

        /// Get total number of events upto given level
        long getTotalEvent(int level);

        /// Resizes the Event List
        void resizeEventArray(int level);

        /// Returns the number of events that have this level as top
        long lengthForLevel(int level) const;

        /// Returns the array of events that have this level as top
        event** arrayForLevel(int level) const;

        /// Returns an array of local states for this level, such that
        /// result[i] == isConfirmed(level, i).
        const bool* getLocalStates(int level);


        /// Returns the number of confirmed states at a level
        int getConfirmedStates(int level) const;

        /// Confirms the local states at a level
        void setConfirmedStates(int level, int i);

        /// Confirms the local states in the given MDD
        void setConfirmedStates(const dd_edge &set);


        /// Checks if i is confirmed
        bool isConfirmedState(int level, int i);

        /// Expand confirm array
        void resizeConfirmedArray(int level, int index);

        /** Bound all extensible variables
            using the maximum confirmed local state as the bound.
        */
        void bindExtensibleVariables();

      public:
        /// Prints the hybrid relation
        void show();

     public:

      /*
       Group the list of events at a given level by same next-of values
       @param level    level at which saturation is called
       @param i        number of tokens before event is fired
       @param R        array of events at level top

       Return the map.
       */
      std::unordered_map<long,std::vector<node_handle>> getListOfNexts(int level, long i, relation_node **R);

      /*
       Returns whether there exist a possibility of doing union
       @param level    level at which saturation is called
       @param i        number of tokens before event is fired
       @param R        array of events at level top

       Return bool.
       */
      bool isUnionPossible(int level, long i, relation_node **R);

    };  // class hybrid_relation


};

// ******************************************************************
// *                                                                *
// *                        operation  class                        *
// *                                                                *
// ******************************************************************

/** Generic operation.
    Operations are tied to specific forests.
    Necessary for compute table entries.
*/
class MEDDLY::operation {
    const opname* theOpName;
    unsigned oplist_index;
    bool is_marked_for_deletion;

    // declared and initialized in meddly.cc
    static compute_table* Monolithic_CT;
    // declared and initialized in meddly.cc
    static operation** op_list;
    // declared and initialized in meddly.cc
    static unsigned* op_holes;
    // declared and initialized in meddly.cc
    static unsigned list_size;
    // declared and initialized in meddly.cc
    static unsigned list_alloc;
    // declared and initialized in meddly.cc
    static unsigned free_list;

    // should ONLY be called during library cleanup.
    static void destroyAllOps();

    /**
      Starting slot for entry_types, assigned
      by compute_table.
    */
    unsigned first_etid;

  protected:
    /// Compute table to use (for entry type 0), if any.
    compute_table* CT0;
    /// Array of compute tables, one per entry type.
    compute_table** CT;
    /** Array of entry types.
        Owned by the compute_table class; we have
        these pointers for convenience.
    */
    compute_table::entry_type** etype;
    /** Array of entry results.
        Use these during computation.
        We only ever need one result per entry type.
    */
    compute_table::entry_result* CTresult;

    /**
      Number of entry_types needed by this operation.
      This gives the dimension of arrays CT and etype.
    */
    unsigned num_etids;

    /// Struct for CT searches.
    // compute_table::entry_key* CTsrch;
    // for cache of operations.
    operation* next;

    virtual ~operation();
    void markForDeletion();
    void registerInForest(forest* f);
    void unregisterInForest(forest* f);
    void allocEntryForests(int nf);
    void addEntryForest(int index, expert_forest* f);
    void allocEntryObjects(int no);
    void addEntryObject(int index);

    virtual bool checkForestCompatibility() const = 0;

    void registerEntryType(unsigned slot, compute_table::entry_type* et);
    void buildCTs();

    friend class forest;
    friend void MEDDLY::destroyOpInternal(operation* op);
    friend void MEDDLY::cleanup();

    friend class ct_initializer;

  public:
    /**
        Constructor.
          @param  n         Operation "name"
          @param  et_slots  Number of different compute table entry types
                            used by this operation.
                            Derived class constructors must register
                            exactly this many entry types.
    */
    operation(const opname* n, unsigned et_slots);

    bool isMarkedForDeletion() const;
    void setNext(operation* n);
    operation* getNext();

    static bool usesMonolithicComputeTable();
    static void removeStalesFromMonolithic();
    static void removeAllFromMonolithic();

    /// Remove stale compute table entries for this operation.
    void removeStaleComputeTableEntries();

    /// Remove all compute table entries for this operation.
    void removeAllComputeTableEntries();

    // for compute tables.

    unsigned getIndex() const;
    static operation* getOpWithIndex(unsigned i);
    static unsigned getOpListSize();

    void setFirstETid(unsigned slot);
    unsigned getFirstETid() const;
    unsigned getNumETids() const;

    // for debugging:

    static void showMonolithicComputeTable(output &, int verbLevel);
    static void showAllComputeTables(output &, int verbLevel);
    static void countAllNodeEntries(const expert_forest* f, size_t* counts);

    void showComputeTable(output &, int verbLevel) const;
    void countCTEntries(const expert_forest* f, size_t* counts) const;

    // handy
    const char* getName() const;
    const opname* getOpName() const;
};

// ******************************************************************
// *                                                                *
// *                     unary_operation  class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a unary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::unary_operation : public operation {
  protected:
    expert_forest* argF;
    expert_forest* resF;
    opnd_type resultType;

    virtual ~unary_operation();

    virtual bool checkForestCompatibility() const;

  public:
    unary_operation(const unary_opname* code, unsigned et_slots,
      expert_forest* arg, expert_forest* res);

    unary_operation(const unary_opname* code, unsigned et_slots,
      expert_forest* arg, opnd_type res);

    bool matches(const expert_forest* arg, const expert_forest* res)
      const;

    bool matches(const expert_forest* arg, opnd_type res) const;

    // high-level front-ends

    /**
      Checks forest comatability and then calls computeDDEdge().
    */
    void compute(const dd_edge &arg, dd_edge &res);
    void computeTemp(const dd_edge &arg, dd_edge &res);

    virtual void compute(const dd_edge &arg, long &res);
    virtual void compute(const dd_edge &arg, double &res);
    virtual void compute(const dd_edge &arg, ct_object &c);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
};

// ******************************************************************
// *                                                                *
// *                     binary_operation class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a binary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::binary_operation : public operation {
  protected:
    bool can_commute;
    expert_forest* arg1F;
    expert_forest* arg2F;
    expert_forest* resF;
    opnd_type resultType;

    virtual ~binary_operation();
    void operationCommutes();

    // Check if the variables orders of relevant forests are compatible
    virtual bool checkForestCompatibility() const;

  public:
    binary_operation(const binary_opname* code, unsigned et_slots,
      expert_forest* arg1, expert_forest* arg2, expert_forest* res);

    bool matches(const expert_forest* arg1, const expert_forest* arg2,
      const expert_forest* res) const;

    // high-level front-end

    /**
      Checks forest comatability and then calls computeDDEdge().
    */
    void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);
    void computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

    virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res, bool userFlag)
      = 0;

    // low-level front ends.  TBD - REMOVE THESE BECAUSE THEY BREAK MARK & SWEEP

#ifdef KEEP_LL_COMPUTES

    /// Low-level compute on nodes a and b, return result.
    virtual node_handle compute(node_handle a, node_handle b);
    /// Low-level compute at level k on nodes a and b, return result.
    virtual node_handle compute(int k, node_handle a, node_handle b);

    /// Low-level compute on EV edges (av, ap) and (bv, bp), return result.
    virtual void compute(int av, node_handle ap, int bv, node_handle bp,
      int &cv, node_handle &cp);
    virtual void compute(long av, node_handle ap, long bv, node_handle bp,
      long &cv, node_handle &cp);
    virtual void compute(long av, node_handle ap, node_handle bp,
      long &cv, node_handle &cp);

    /// Low-level compute on EV edges (av, ap) and (bv, bp), return result.
    virtual void compute(float av, node_handle ap, float bv, node_handle bp,
      float &cv, node_handle &cp);

#endif

};

// ******************************************************************
// *                                                                *
// *                  specialized_operation  class                  *
// *                                                                *
// ******************************************************************

/** Mechanism to apply specialized operations.
*/
class MEDDLY::specialized_operation : public operation {
  public:
    specialized_operation(const specialized_opname* code, unsigned et_slots);
  protected:
    virtual ~specialized_operation();
  public:

    /** For unary (like) operations.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &arg, dd_edge &res);

    /** For unary (like) operations with boolean results.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &arg, bool &res);

    /** For binary (like) operations.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

    /** For numerical operations.
        compute y += some function of x, depending on the operation.
        Default behavior is to throw an exception.
    */
    virtual void compute(double* y, const double* x);

    /** For tenary (like) operations.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &ar1, const dd_edge &ar2, const dd_edge &ar3, dd_edge &res);
};

// ******************************************************************
// *                                                                *
// *                    global_rebuilder  class                     *
// *                                                                *
// ******************************************************************

/** Rebuild the dd_edge from the source forest in the target forest.
    The source and target forests may have different variable orders.
    While rebuilding, extra nodes may be created in the source forest
    because of the restrict operation.
*/

class MEDDLY::global_rebuilder {
private:
  struct RestrictKey {
    node_handle p;
    int var;
    int idx;

    bool operator==(const RestrictKey &other) const {
      return (p == other.p && var == other.var && idx == other.idx);
    }
  };

  struct RestrictKeyHasher {
    size_t operator()(const RestrictKey &key) const;
  };

  struct TransformKey {
    int sig;
//    int var;

    bool operator==(const TransformKey &other) const {
      return sig == other.sig;
//      return (sig == other.sig && var == other.var);
    }
  };

  struct TransformEntry {
    // Partial assignment
    std::vector<int> pa;
    node_handle p;

    bool operator==(const TransformEntry &other) const {
      return p == other.p;
    }
  };

  struct TransformKeyHasher {
    size_t operator()(const TransformKey &key) const;
  };

  class SignatureGenerator {
  protected:
    global_rebuilder &_gr;

  public:
    SignatureGenerator(global_rebuilder& gr);
    virtual ~SignatureGenerator();
    virtual void precompute() = 0;
    virtual int signature(node_handle p) = 0;
  };

  class TopDownSignatureGenerator: public SignatureGenerator {
  public:
    TopDownSignatureGenerator(global_rebuilder& gr);
    void precompute() override;
    int signature(node_handle p) override;
  };

  class BottomUpSignatureGenerator: public SignatureGenerator {
  private:
    std::unordered_map<node_handle, int> _cache_sig;
    std::unordered_map<node_handle, int> _cache_rec_sig;

    int rec_signature(node_handle p);

  public:
    BottomUpSignatureGenerator(global_rebuilder& gr);
    void precompute() override;
    int signature(node_handle p) override;
  };

  std::unordered_map<RestrictKey, node_handle, RestrictKeyHasher> _computed_restrict;
  std::unordered_multimap<TransformKey, TransformEntry, TransformKeyHasher> _computed_transform;
  SignatureGenerator* _sg;

  expert_forest* _source;
  expert_forest* _target;
  node_handle _root;
  int _hit;
  int _total;

  node_handle transform(node_handle p, int target_level, std::vector<int>& pa);
  node_handle restrict(node_handle p, std::vector<int>& pa);

  bool restrict_exist(node_handle p, const std::vector<int>& pa,
      unsigned start, node_handle& result);
  int signature(node_handle p) const;

  // Return the top variable in the sub-order of the target variable order
  // starting from 0 to target_level
  // such that the given decision diagram depends on it.
  int check_dependency(node_handle p, int target_level) const;

public:
  friend class SignatureGenerator;

  global_rebuilder(expert_forest* source, expert_forest* target);
  ~global_rebuilder();

  dd_edge rebuild(const dd_edge& e);
  void clearCache();
  double hitRate() const;
};

#endif

