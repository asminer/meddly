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

#ifndef MEDDLY_POLICIES_H
#define MEDDLY_POLICIES_H

#include "defines.h"

namespace MEDDLY {
    struct policies;
    class memory_manager_style;
    class node_storage_style;

    // ******************************************************************
    // *                    Memory management schemes                   *
    // ******************************************************************

    extern const memory_manager_style* ORIGINAL_GRID;
    extern const memory_manager_style* ARRAY_PLUS_GRID;
    extern const memory_manager_style* MALLOC_MANAGER;
    extern const memory_manager_style* HEAP_MANAGER;
    extern const memory_manager_style* FREELISTS;   // used for compute tables

    // ******************************************************************
    // *                     Node storage mechanisms                    *
    // ******************************************************************

    /**
        New, "simple" style with memory manager removed.
        This node storage mechanism relies on the
        memory_manager_style for memory management.
    */
    extern const node_storage_style* SIMPLE_STORAGE;
    // extern const node_storage_style* PATTERN_STORAGE;
    // extern const node_storage_style* BEST_STORAGE;

    //
    // From here are "old" mechanisms for node storage,
    // with built-in memory managers.
    //

    /** "Classic" node storage mechanism.
        The node storage mechanism from early versions of this library.
    */
    extern const node_storage_style* CLASSIC_STORAGE;

    /** Similar to classic.
        Differences are in class design, so we can measure overhead
        (if any) of class design differences.
    */
    extern const node_storage_style* SIMPLE_GRID;

    /// Like classic, but use an array of lists for hole management.
    extern const node_storage_style* SIMPLE_ARRAY;

    /// Like classic, but use heaps for hole management.
    extern const node_storage_style* SIMPLE_HEAP;

    /// Classic node storage but no hole management.
    extern const node_storage_style* SIMPLE_NONE;

    /** Nodes are stored in a compressed form.
        Holes are managed using the original grid structure.
    */
    extern const node_storage_style* COMPACT_GRID;

    // ******************************************************************
    // *                         reduction rules                        *
    // ******************************************************************

    /** Supported node reduction rules.
        Currently, the following reduction rules are allowed:
            - Quasi-reduced.
                All zero nodes are eliminated.
                Pointers from a node are always directly to the level
                below, or directly to terminal node zero.

            - Fully-reduced.
                All redundant nodes are eliminated.

            - Identity-reduced.
                For relations only.
                At unprimed levels, all redundant nodes are eliminated.
                At primed levels, an i-singleton node is one that
                is all zero, execept for index i. An i-singleton node
                cannot have an incoming pointer from any level except
                the unprimed level directly above it, and cannot have
                an incoming pointer that is from the i-th child.

        In all cases, duplicate nodes are eliminated.
    */
    enum class reduction_rule {
        /// Nodes are fully reduced.
        FULLY_REDUCED,
        /// Nodes are quasi-reduced.
        QUASI_REDUCED,
        /// Nodes are identity-reduced.
        IDENTITY_REDUCED,
        /// Nodes are user-defined reduced
        USER_DEFINED
    };

    inline const char* nameOf(reduction_rule rr)
    {
        switch (rr) {
            case reduction_rule::FULLY_REDUCED:     return "fully-reduced";
            case reduction_rule::QUASI_REDUCED:     return "quasi-reduced";
            case reduction_rule::IDENTITY_REDUCED:  return "identity-reduced";
            case reduction_rule::USER_DEFINED:      return "userdef-reduced";
        }
    }

    // ******************************************************************
    // *                       node storage  rules                      *
    // ******************************************************************

    /// Flags for node storage.
    typedef unsigned char node_storage_flags;

    const node_storage_flags FULL_ONLY      = 0x01;
    const node_storage_flags SPARSE_ONLY    = 0x02;
    const node_storage_flags FULL_OR_SPARSE = 0x03;

    // ******************************************************************
    // *                         set vs relation                        *
    // ******************************************************************

    /// Set vs relation identifiers, for improved readability.
    typedef bool set_or_rel;

    const set_or_rel    SET             = false;
    const set_or_rel    VECTOR          = false;
    const set_or_rel    UNPRIMED_ONLY   = false;

    const set_or_rel    RELATION        = true;
    const set_or_rel    MATRIX          = true;
    const set_or_rel    WITH_PRIMED     = true;

    // ******************************************************************
    // *                         edge  labelings                        *
    // ******************************************************************

    /// Edge annotation mechanism.
    enum class edge_labeling {
        /// Edges unlabeled, all values stored in distinct terminals.
        MULTI_TERMINAL,
        /// Edges labeled, values summed along path.
        EVPLUS,
        /// Special case of EVPLUS for indexed sets.
        INDEX_SET,
        /// Edges labeled, values multiplied along path.
        EVTIMES
        // TBD: there may be others in the future :^)
    };

    inline const char* nameOf(edge_labeling el) {
        switch (el) {
            case edge_labeling::MULTI_TERMINAL: return "MT";

            case edge_labeling::EVPLUS:         return "EV+";
            case edge_labeling::INDEX_SET:      return "EV+";

            case edge_labeling::EVTIMES:        return "EV*";
        }
    }

    // ******************************************************************
    // *                         reporting flags                        *
    // ******************************************************************

    /// Flags for forest reports
    typedef unsigned reporting_flags;

    const reporting_flags HUMAN_READABLE_MEMORY   = 0x0001;
    const reporting_flags BASIC_STATS             = 0x0002;
    const reporting_flags EXTRA_STATS             = 0x0004;
    const reporting_flags FOREST_STATS            = 0x0008;
    const reporting_flags STORAGE_STATS           = 0x0010;
    const reporting_flags STORAGE_DETAILED        = 0x0020;
    const reporting_flags UNIQUE_TABLE_STATS      = 0x0040;
    const reporting_flags UNIQUE_TABLE_DETAILED   = 0x0080;
    const reporting_flags HOLE_MANAGER_STATS      = 0x0100;
    const reporting_flags HOLE_MANAGER_DETAILED   = 0x0200;


    // ******************************************************************
    // *                          display flags                         *
    // ******************************************************************

    typedef unsigned display_flags;

    const display_flags SHOW_TERMINALS    = 0x01;
    const display_flags SHOW_INDEX        = 0x02;
    const display_flags SHOW_DETAILS      = 0x04;
    const display_flags SHOW_UNREACHABLE  = 0x08;
    const display_flags SHOW_DELETED      = 0x10;

};

// ******************************************************************
// *                         policies struct                        *
// ******************************************************************

/// Collection of forest policies.
struct MEDDLY::policies {

    /// Supported node deletion policies.
    enum class node_deletion {
        /** Never delete nodes.
            Useful for debugging the garbage collector.
        */
        NEVER,
        /** Nodes are marked for deletion.
            We assume that a deleted node might be used again,
            so we keep it until it no longer appears in any compute table.
        */
        OPTIMISTIC,
        /// Nodes are removed as soon as they become disconnected.
        PESSIMISTIC
    };

    /// Supported variable swap strategies.
    /// Work for relations only.
    enum class variable_swap_type {
  	    /// Swap adjacent variables in one go.
  	    VAR,
  	    /// Swap adjacent variables through swapping levels 4 time.
  	    /// Do not work for fully-identity reduced relations.
  	    LEVEL
    };

    /// Schedule heuristic for reordering via swaps.
    enum class reordering_type {
        /// Always choose the lowest swappable inversion
        LOWEST_INVERSION,
        /// Always choose the highest swappable inversion
        HIGHEST_INVERSION,
        /// Sink down the "heaviest" variables
        SINK_DOWN,
        /// Bubble up the "lightest" variables
        BRING_UP,
        /// Always choose the swappabel inversion where the variable on the top
        /// has the fewest associated nodes
        LOWEST_COST,
        /// Always choose the swappable inversion that will result in
        /// the lowest memory consumption
        LOWEST_MEMORY,
        /// Choose the swappable inversion randomly
        RANDOM,
        /// Always choose the swappable inversion with the
        /// lowest average reference count
        LARC
    };

    /// Defaults: how may we store nodes for all levels in the forest.
    node_storage_flags storage_flags;
    /// Default reduction rule for all levels in the forest.
    reduction_rule reduction;
    /// Default deletion policy for all levels in the forest.
    node_deletion deletion;
    // Default variable reorder strategy.
    reordering_type reorder;
    // Default variable swap strategy.
    variable_swap_type swap;

    /// Backend memory management mechanism for nodes.
    const memory_manager_style* nodemm;

    /// Backend storage mechanism for nodes.
    const node_storage_style* nodestor;

    /// Memory compactor: never run if fewer than this many unused slots.
    // int compact_min;
    /// Memory compactor: always run if more than this many unused slots.
    // int compact_max;
    /// Memory compactor: fraction of unused slots to trigger.
    // int compact_frac;

    /// Number of zombie nodes to trigger garbage collection
    // int zombieTrigger;
    /// Number of orphan nodes to trigger garbage collection
    // int orphanTrigger;
    /// Should we run the memory compactor after garbage collection
    // bool compactAfterGC;
    /// Should we run the memory compactor before trying to expand
    // bool compactBeforeExpand;

    /// Use reference counts to know when nodes can be recycled.
    /// Otherwise, use mark and sweep to recycle disconnected nodes.
    bool useReferenceCounts;

    /// Empty constructor, for setting up defaults later
    policies();

    /// Constructor; sets reasonable defaults
    policies(set_or_rel rel);

    /// Set to hard-wired defaults
    void useDefaults(set_or_rel rel);

    //
    // Node storage helper methods
    //

    inline void setFullStorage() {
        storage_flags = FULL_ONLY;
    }

    inline void setSparseStorage() {
        storage_flags = SPARSE_ONLY;
    }

    inline void setFullOrSparse() {
        storage_flags = FULL_OR_SPARSE;
    }

    inline bool allowsSparse() const {
        return SPARSE_ONLY & storage_flags;
    }

    inline bool allowsFull() const {
        return FULL_ONLY & storage_flags;
    }

    //
    // Reduction rule helper methods
    //

    inline void setFullyReduced() {
        reduction = reduction_rule::FULLY_REDUCED;
    }
    inline bool isFullyReduced() const {
        return reduction_rule::FULLY_REDUCED == reduction;
    }

    inline void setQuasiReduced() {
        reduction = reduction_rule::QUASI_REDUCED;
    }
    inline bool isQuasiReduced() const {
        return reduction_rule::QUASI_REDUCED == reduction;
    }

    inline void setIdentityReduced() {
        reduction = reduction_rule::IDENTITY_REDUCED;
    }
    inline bool isIdentityReduced() const {
        return reduction_rule::IDENTITY_REDUCED == reduction;
    }

    inline void setUserDefinedReduced() {
        reduction = reduction_rule::USER_DEFINED;
    }
    inline bool isUserDefinedReduced() const {
        return reduction_rule::USER_DEFINED == reduction;
    }

    //
    // Node deletion policy helper methods
    //

    inline void setNeverDelete() {
        deletion = node_deletion::NEVER;
    }
    inline bool isNeverDelete() const {
        return node_deletion::NEVER == deletion;
    }

    inline void setOptimistic() {
        deletion = node_deletion::OPTIMISTIC;
    }
    inline bool isOptimistic() const {
        return node_deletion::OPTIMISTIC == deletion;
    }

    inline void setPessimistic() {
        deletion = node_deletion::PESSIMISTIC;
    }
    inline bool isPessimistic() const {
        return node_deletion::PESSIMISTIC == deletion;
    }

    //
    // Reordering policy helper methods
    //

    inline void setLowestInversion() {
        reorder = reordering_type::LOWEST_INVERSION;
    }

    inline void setHighestInversion() {
        reorder = reordering_type::HIGHEST_INVERSION;
    }

    inline void setSinkDown() {
        reorder = reordering_type::SINK_DOWN;
    }

    inline void setBringUp() {
        reorder = reordering_type::BRING_UP;
    }

    inline void setLowestCost() {
        reorder = reordering_type::LOWEST_COST;
    }

    inline void setLowestMemory() {
        reorder = reordering_type::LOWEST_MEMORY;
    }

    inline void setRandom() {
        reorder = reordering_type::RANDOM;
    }

    inline void setLARC() {
        reorder = reordering_type::LARC;
    }

    inline void setVarSwap() {
        swap = variable_swap_type::VAR;
    }

    inline bool isVarSwap() const {
        return variable_swap_type::VAR == swap;
    }

    inline void setLevelSwap() {
        swap = variable_swap_type::LEVEL;
    }

    inline bool isLevelSwap() const {
        return variable_swap_type::VAR == swap;
    }

};

#endif
