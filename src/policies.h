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
    extern const node_storage_style* PATTERN_STORAGE;
    extern const node_storage_style* BEST_STORAGE;

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
          - Fully reduced, meaning that duplicate and redundant nodes are
            eliminated.
          - Quasi reduced, meaning that duplicate nodes are eliminated.
          - Identity reduced, for relations only, meaning that duplicate
            nodes are eliminated, as are "identity" pairs of primed, unprimed
            variables.
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


    // ******************************************************************
    // *                       node storage  rules                      *
    // ******************************************************************

    /// Flags for node storage.
    typedef unsigned char node_storage_flags;

    const node_storage_flags FULL_ONLY      = 0x01;
    const node_storage_flags SPARSE_ONLY    = 0x02;
    const node_storage_flags FULL_OR_SPARSE = 0x03;


    // ******************************************************************
    // *                           range types                          *
    // ******************************************************************

    /** Types of values that we can currently store in forests.
        I.e., if every node in a forest is a function,
        these are the possible ranges for a function.
    */
    enum class range_type {
        /// boolean-valued functions.
        BOOLEAN,
        /// integer-valued functions.
        INTEGER,
        /// real-valued functions.
        REAL
    };

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
    policies(bool rel);

    /// Set to hard-wired defaults
    void useDefaults(bool rel);

    inline void setFullStorage() {
        storage_flags = FULL_ONLY;
    }

    inline void setSparseStorage() {
        storage_flags = SPARSE_ONLY;
    }

    inline void setFullOrSparse() {
        storage_flags = FULL_OR_SPARSE;
    }

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

    inline void setLevelSwap() {
        swap = variable_swap_type::LEVEL;
    }

};

#endif
