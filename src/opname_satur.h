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

#ifndef MEDDLY_OPNAME_SATUR_H
#define MEDDLY_OPNAME_SATUR_H

#include "opname.h"
#include "dd_edge.h"
#include "error.h"
#include "forest.h"

#include <unordered_map>
#include <map>
#include <set>

namespace MEDDLY {
    class satimpl_opname;
    class sathyb_opname;
    class constrained_opname;
};


// ******************************************************************
// *                                                                *
// *                      satimpl_opname class                      *
// *                                                                *
// ******************************************************************

#if 0

/** Saturation, transition relations stored implcitly, operation names.
    Implemented in operations/sat_impl.cc
*/
class MEDDLY::satimpl_opname: public specialized_opname {
  public:

    satimpl_opname(const char* n);
    virtual ~satimpl_opname();

    /// Arguments should have type "implicit_relation", below
    virtual specialized_operation* buildOperation(arguments* a);

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
        forest* getMixRelForest() const;


        /// Returns the MDD forest that stores the initial set of states
        forest* getInForest() const;


        /// Returns the MDD forest that stores the resultant set of states
        forest* getOutForest() const;

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
        forest* insetF;
        forest* outsetF;
        forest* mixRelF;

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
        forest* getRelForest() const;


      private:
        forest* mxdF;

     public:

      /*
       Group the list of events at a given level by same next-of values
       @param level    level at which saturation is called
       @param i        number of tokens before event is fired
       @param R        array of events at level top

       Return the map.
       */
      std::unordered_map<long,std::vector<rel_node_handle> > getListOfNexts(int level, long i, relation_node **R);

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




// ******************************************************************
// *                                                                *
// *                 inlined  satimpl_opname methods                *
// *                                                                *
// ******************************************************************


inline MEDDLY::relation_node*
MEDDLY::satimpl_opname::implicit_relation::nodeExists(rel_node_handle n)
{
  std::unordered_map<rel_node_handle, relation_node*>::iterator finder = impl_unique.find(n);
  if(finder!=impl_unique.end())
    return finder->second;
  else
    return NULL;
}

inline bool
MEDDLY::satimpl_opname::implicit_relation::isReserved(rel_node_handle n)
{
  return (n==1);
}

//************************************************************************

inline MEDDLY::forest*
MEDDLY::satimpl_opname::implicit_relation::getInForest() const
{
  return insetF;
}

inline MEDDLY::forest*
MEDDLY::satimpl_opname::implicit_relation::getOutForest() const
{
  return outsetF;
}

inline MEDDLY::forest*
MEDDLY::satimpl_opname::implicit_relation::getMixRelForest() const
{
  return mixRelF;
}

// ***********************************************************************

inline long
MEDDLY::satimpl_opname::implicit_relation::getTotalEvent(int level)
{
  int total_event = 0;
  for(int i=1;i<=level;i++)
    total_event +=  lengthForLevel(i);

  return total_event;
}

inline long
MEDDLY::satimpl_opname::implicit_relation::lengthForLevel(int level) const
{
  return event_added[level];
}

inline MEDDLY::rel_node_handle*
MEDDLY::satimpl_opname::implicit_relation::arrayForLevel(int level) const
{
  return event_list[level];
}

// ****************************************************************************

inline MEDDLY::forest*
MEDDLY::satimpl_opname::implicit_relation::getRelForest() const
{
  return mxdF;
}


inline long
MEDDLY::satimpl_opname::implicit_relation::getConfirmedStates(int level) const
{
  return confirm_states[level];
}

inline void
MEDDLY::satimpl_opname::implicit_relation::setConfirmedStates(int level,int i)
{
  resizeConfirmedArray(level,i);
  MEDDLY_DCASSERT(confirmed_array_size[level]>i);
  if(!isConfirmedState(level,i))
    {
      confirmed[level][i]=true;
      confirm_states[level]++;
    }
}

inline bool
MEDDLY::satimpl_opname::implicit_relation::isConfirmedState(int level,int i)
{

  return (i < insetF->getLevelSize(level) && confirmed[level][i]);
}

#endif

// ******************************************************************
// *                                                                *
// *                    constrained_opname class                    *
// *                                                                *
// ******************************************************************

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
    virtual specialized_operation* buildOperation(arguments* a);

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
      forest* getForest();

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


      forest* f;
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
      inline forest* getForest() { return f; }

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
      inline std::map<int, std::set<node_handle> > getComponents()  { return level_component; }

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
      forest* f;
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

      std::map<int,std::set<node_handle> > level_component; // stores the set of subevent node_handles whose top is this level.
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
        forest* getHybridForest() const;


        /// Returns the MDD forest that stores the initial set of states
        forest* getInForest() const;


        /// Returns the MDD forest that stores the resultant set of states
        forest* getOutForest() const;


        /// If only relation_nodes are present
        /// return the vector of events pertaining to a given level
        /// that have the same net effect
        std::vector<node_handle> getRelNodeAtLevelWithEffect(int level, long effect);

      private:
        forest* insetF;
        forest* outsetF;
        forest* hybRelF;

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
      std::unordered_map<long,std::vector<node_handle> > getListOfNexts(int level, long i, relation_node **R);

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
// *                 inlined  sathyb_opname methods                 *
// *                                                                *
// ******************************************************************

inline MEDDLY::forest*
MEDDLY::sathyb_opname::hybrid_relation::getInForest() const
{
  return insetF;
}

inline MEDDLY::forest*
MEDDLY::sathyb_opname::hybrid_relation::getOutForest() const
{
  return outsetF;
}

inline MEDDLY::forest*
MEDDLY::sathyb_opname::hybrid_relation::getHybridForest() const
{
  return hybRelF;
}

// ***********************************************************************

inline long
MEDDLY::sathyb_opname::hybrid_relation::getTotalEvent(int level)
{
  int total_event = 0;
  for(int i=1;i<=level;i++)
    total_event +=  lengthForLevel(i);

  return total_event;
}

inline long
MEDDLY::sathyb_opname::hybrid_relation::lengthForLevel(int level) const
{
  return num_events_by_top_level[level];
}

inline MEDDLY::sathyb_opname::event**
MEDDLY::sathyb_opname::hybrid_relation::arrayForLevel(int level) const
{
  return events_by_top_level[level];
}

// ****************************************************************************

inline int
MEDDLY::sathyb_opname::hybrid_relation::getConfirmedStates(int level) const
{
  return confirm_states[level];
}


inline bool
MEDDLY::sathyb_opname::hybrid_relation::isConfirmedState(int level,int i)
{
  return (i < insetF->getLevelSize(level) && (i < confirmed_array_size[level]) && confirmed[level][i]);
}

// ******************************************************************


inline MEDDLY::forest*
MEDDLY::sathyb_opname::subevent::getForest() {
  return f;
}

inline int
MEDDLY::sathyb_opname::subevent::getNumVars() const {
  return num_vars;
}

inline const int*
MEDDLY::sathyb_opname::subevent::getVars() const {
  return vars;
}

inline int
MEDDLY::sathyb_opname::subevent::getTop() const {
  return top;
}

inline bool
MEDDLY::sathyb_opname::subevent::isFiring() const {
  return is_firing;
}

inline bool
MEDDLY::sathyb_opname::subevent::isEnabling() const {
  return !is_firing;
}

inline bool
MEDDLY::sathyb_opname::subevent::isImplicit() const {
  return num_vars == 1;
}

inline const MEDDLY::dd_edge&
MEDDLY::sathyb_opname::subevent::getRoot() const {
  return root;
}

inline const MEDDLY::node_handle
MEDDLY::sathyb_opname::subevent::getRootHandle() const {
  return root_handle;
}

inline void
MEDDLY::sathyb_opname::subevent::setRootHandle( node_handle ID ) {
  root_handle = ID;
}

inline MEDDLY::node_handle
MEDDLY::sathyb_opname::subevent::getDown() const {
  return down;
}

inline void
MEDDLY::sathyb_opname::subevent::setDown( node_handle d_ID ) {
   down = d_ID;
}

inline int
MEDDLY::sathyb_opname::subevent::getEnable() const {
  return enable;
}

inline int
MEDDLY::sathyb_opname::subevent::getFire() const {
  return fire;
}

inline bool
MEDDLY::sathyb_opname::subevent::usesExtensibleVariables() const {
  return uses_extensible_variables;
}


#endif // #include guard
