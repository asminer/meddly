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

#ifndef MEDDLY_FOREST_H
#define MEDDLY_FOREST_H

#include "defines.h"
#include "memstats.h"
#include "varorder.h"
#include "operators.h"
#include "node_headers.h"
#include "node_storage.h"
#include "unpacked_node.h"
#include "policies.h"
#include "domain.h"
#include "enumerator.h"

#include <memory>
#include <map>

namespace MEDDLY {
    class domain;

    class operation;
    class forest;
    class expert_forest;
    class node_marker;
    class node_storage;
    class unpacked_node;
    class unique_table;
    class impl_unique_table;
    class relation_node;
    class global_rebuilder;

    /** Front-end function to destroy a forest.
    */
    void destroyForest(forest* &f);
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                          forest class                          *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Forest class.
    Abstract base class.

    A data structure for managing collections of functions
    (or sets, or vectors, or relations, or matrices,
    depending on your conceptual view) represented in a single
    decision diagram forest over a common domain.

    Each forest is assigned a unique identifier "forever"
    (actually, until the library is re-initialized).

    TBD: node header / node storage split

    TBD: discussion of garbage collection.

*/
class MEDDLY::forest {

    // ------------------------------------------------------------
    public: // static "constructors"
    // ------------------------------------------------------------

        /** Create a forest.
            Conceptually, a forest is a structure used to represent a
            collection of functions over a common domain. For simplicity
            (although it is a slight abuse of notation) a forest may represent
            "vectors" or "sets" over a domain, or "matrices" or "relations".
            In case of matrices / relations, the forest uses primed and
            unprimed versions of every variable in the domain.

            @param  sr      Is this a relation / matrix, versus a set / vector.

            @param  t       Range type of the functions, namely,
                            booleans, integers, or reals.

            @param  ev      Edge labeling mechanism, i.e., should this be a
                            Multi-terminal decision diagram forest,
                            edge-valued with plus/times decision diagram forest.

            @param  p       Policies to use within the forest.

            @param  level_reduction_rule
                            Rules for reduction on different levels.

            @param  tv      Transparent value.

            @return nullptr if an error occurs, a new forest otherwise.
        */
        static forest* create(domain* d, set_or_rel sr, range_type t,
            edge_labeling ev, const policies &p,
            int* level_reduction_rule=nullptr, int tv=0);

        /// Create a forest using the library default policies.
        inline static forest* create(domain* d, set_or_rel sr, range_type t,
            edge_labeling ev)
        {
            return create(d, sr, t, ev,
                sr ? getDefaultPoliciesMXDs() : getDefaultPoliciesMDDs(),
                nullptr, 0
            );
        }

    // ------------------------------------------------------------
    public: // Retrieve a node from the forest
    // ------------------------------------------------------------

        /// Copy into an unpacked node.
        ///     @param  un      Unpacked node to copy into
        ///     @param  p       Handle of node to retrieve
        ///     @param  st2     Storage: full only, sparse only, or
        ///                     whatever is easier.
        void unpackNode(unpacked_node* un, node_handle p,
                node_storage_flags st2) const;

        /// Copy into a new unpacked node.
        inline unpacked_node* newUnpacked(node_handle p,
                node_storage_flags f) const
        {
            unpacked_node* un = unpacked_node::New();
            unpackNode(un, p, f);
            return un;
        }


    // ------------------------------------------------------------
    public: // Getters related to the domain / levels
    // ------------------------------------------------------------

        /// Returns a non-modifiable pointer to this forest's domain.
        inline const domain* getDomain() const {
            return d;
        }

        /// Returns a pointer to this forest's domain.
        inline domain* getDomain() {
            return d;
        }

#ifdef ALLOW_DEPRECATED_0_17_3
        /// Returns a pointer to this forest's domain.
        /// Deprecated; use getDomain() instead.
        inline domain* useDomain() {
            return getDomain();
        }
#endif

        inline int getNumVariables() const {
            return d->getNumVariables();
        }

        /**
            Negative values are used for primed levels or variables.
        */
        inline int getVarByLevel(int level) const {
            return level > 0
                ? var_order->getVarByLevel(level)
                : -var_order->getVarByLevel(-level);
        }
        inline int getLevelByVar(int var) const {
            return var > 0
                ? var_order->getLevelByVar(var)
                : -var_order->getLevelByVar(-var);
        }

        /// returns 0 or -K, depending if it's a relation
        inline int getMinLevelIndex() const {
            return isForRelations() ? -getNumVariables() : 0;
        }

        /// Check if the given level is valid
        inline bool isValidLevel(int k) const {
            return (k >= getMinLevelIndex()) && (k <= getNumVariables());
        }

        /// Can we have extensible nodes at level k?
        inline bool isExtensibleLevel(int k) const {
            MEDDLY_DCASSERT(isValidLevel(k));
            return d->getVar(k < 0? -k: k)->isExtensible();
        }

        /// The maximum size (number of indices) a node at this level can have
        inline int getLevelSize(int k) const {
            MEDDLY_DCASSERT(isValidLevel(k));
            int var=getVarByLevel(k);
            if (var < 0) {
                return getDomain()->getVariableBound(-var, true);
            } else {
                return getDomain()->getVariableBound(var, false);
            }
        }

        /// The maximum size (number of indices) a variable can have.
        inline int getVariableSize(int var) const {
            return getDomain()->getVariableBound(var, false);
        }


    // ------------------------------------------------------------
    private: // Domain info
    // ------------------------------------------------------------
        domain* d;
        friend class domain;

    protected:
        std::shared_ptr<const variable_order> var_order;



    // ------------------------------------------------------------
    public: // interpreting edges in the forest
    // ------------------------------------------------------------

        /**
            Get the transparent node.
            This is the default node value for "skipped" edges in sparse nodes.
        */
        inline node_handle getTransparentNode() const {
            return transparent_node;
        }

        /**
            Get the transparent edge value.
            This is the default node value for "skipped" edges in sparse nodes.
        */
        inline const edge_value& getTransparentEdge() const {
            return transparent_edge;
        }

        /**
            Is the given edge transparent?
            If so it may be "skipped" in a sparse node.
                @param  ep    Node part of the edge to check
                @param  ev    Value part of the edge to check
        */
        inline bool isTransparentEdge(node_handle ep, const edge_value &ev) const {
            if (ep != transparent_node) return false;
            return ev == transparent_edge;
        }

        /**
            Get the transparent edge value.
            This is the default edge value for "skipped" edges in sparse nodes.
            Copy the transparent edge value into the parameters.
                @param  ep      Put node part of the edge here.
                @param  ev      Put value part of the edge here.
        */
        inline void getTransparentEdge(node_handle &ep, edge_value &ptr) const {
            ep = transparent_node;
            ptr = transparent_edge;
        }

        /**
            Make a transparent edge
        */
        inline void getTransparentEdge(dd_edge &e) const {
            e.set(transparent_node, transparent_edge);
        }

        /// Get the edge type.
        inline edge_type getEdgeType() const {
            return the_edge_type;
        }

        /// Are edge values included when computing the hash.
        inline bool areEdgeValuesHashed() const {
            return hash_edge_values;
        }

        /// Get the terminal node type.
        inline terminal_type getTerminalType() const {
            return the_terminal_type;
        }

        /**
            Convenience function.
            Based on the forest type, convert the desired value
            into a terminal node handle.
                @param  v   Value to encode
                @return     Handle for terminal node
        */
        template <typename T>
        inline node_handle handleForValue(T v) const {
            terminal t(v, the_terminal_type);
            return t.getHandle();
        }

        /**
            Convenience function.
            Based on the forest type, convert the terminal node handle
            into its encoded value.
                @param  n   Node handle
                @param  v   Output: encoded value
        */
        template <typename T>
        inline void getValueFromHandle(node_handle n, T& v) const {
            MEDDLY_DCASSERT(n <= 0);
            terminal t(the_terminal_type, n);
            t.getValue(v);
        }

        inline bool getBooleanFromHandle(MEDDLY::node_handle n) const {
            bool v;
            getValueFromHandle(n, v);
            return v;
        }

        inline int getIntegerFromHandle(MEDDLY::node_handle n) const {
            int v;
            getValueFromHandle(n, v);
            return v;
        }

        inline float getRealFromHandle(MEDDLY::node_handle n) const {
            float v;
            getValueFromHandle(n, v);
            return v;
        }

    // ------------------------------------------------------------
    public: // node header getters/setters
    // ------------------------------------------------------------
        /// Return the level number for a given node
        inline int getNodeLevel(node_handle p) const {
            if (isTerminalNode(p)) return 0;
            return nodeHeaders.getNodeLevel(p);
        }

        /// Is a given node at a primed level?
        inline bool isPrimedNode(node_handle p) const {
            return getNodeLevel(p) < 0;
        }

        /// Is a given node at an unprimed level?
        inline bool isUnprimedNode(node_handle p) const {
            return getNodeLevel(p) > 0;
        }

        /// Is a given node still active (not recycled)?
        inline bool isActiveNode(node_handle p) const {
            return nodeHeaders.isActive(p);
        }

        /// Is a given node deleted?
        inline bool isDeletedNode(node_handle p) const {
            return nodeHeaders.isDeleted(p);
        }

        /// Is a given node a terminal node?
        inline static bool isTerminalNode(node_handle p) {
            return (p < 1);
        }

#ifdef ALLOW_DEPRECATED_0_17_3
        inline bool isValidNonterminalIndex(MEDDLY::node_handle p) const {
            return (p > 0) && (p <= nodeHeaders.lastUsedHandle());
        }

        inline bool isValidNodeIndex(MEDDLY::node_handle p) const {
            return p <= nodeHeaders.lastUsedHandle();
        }
#endif

        inline node_handle getLastNode() const {
            return nodeHeaders.lastUsedHandle();
        }

        /// Returns the in-count for a node.
        inline unsigned long getNodeInCount(node_handle p) const {
            MEDDLY_DCASSERT(deflt.useReferenceCounts);
            return nodeHeaders.getIncomingCount(p);
        }


        /** Increase the link count to this node.
            Call this when another node is made to point to this node.
                @return p, for convenience.
        */
        inline node_handle linkNode(node_handle p) {
            if (deflt.useReferenceCounts) {
                return nodeHeaders.linkNode(p);
            } else {
                return p;
            }
        }

        /** Increase the link count to this node.
            Call this when another node is made to point to this node.
                @return the node handle, for convenience.
        */
        inline node_handle linkNode(const dd_edge &p) {
            MEDDLY_DCASSERT(p.isAttachedTo(this));
            if (deflt.useReferenceCounts) {
                return nodeHeaders.linkNode(p.getNode());
            } else {
                return p.getNode();
            }
        }

        /** Decrease the link count to this node.
            If link count reduces to 0, this node may get marked for deletion.
            Call this when another node releases its connection to this node.
        */
        inline void unlinkNode(node_handle p) {
            if (deflt.useReferenceCounts) {
                nodeHeaders.unlinkNode(p);
            }
        }


        /** Increase the cache count for this node.
            Call this whenever this node is added to a cache.
            Do nothing if we are not using reference counts for caches.
                @param  p     Node we care about.
        */
        inline void cacheNode(node_handle p) {
            if (deflt.useReferenceCounts) {
                nodeHeaders.cacheNode(p);
            }
        }

        /** Decrease the cache count for this node.
            Call this whenever this node is added to a cache.
            Do nothing if we are not using reference counts for caches.
                @param  p     Node we care about.
        */
        inline void uncacheNode(node_handle p) {
            if (deflt.useReferenceCounts) {
                nodeHeaders.uncacheNode(p);
            }
        }

        /** Mark the node as belonging to some cache entry.
            Do nothing if we are using reference counts for caches.
                @param p    Node we care about.
        */
        inline void setCacheBit(node_handle p) {
            if (!deflt.useReferenceCounts) {
                nodeHeaders.setInCacheBit(p);
            }
        }

        /** Clear all cache bits.
            Do nothing if we are using reference counts for caches.
        */
        inline void clearAllCacheBits() {
            if (!deflt.useReferenceCounts) {
#ifdef DEBUG_MARK_SWEEP
                std::cerr << "Clearing cache bits for forest "
                          << FID() << "\n";
#endif
                nodeHeaders.clearAllInCacheBits();
            }
        }


        /** Sweep all cache bits.
            Called by CT to indicate we can begin sweep phase
            on cache bits.
        */
        inline void sweepAllCacheBits() {
            if (!deflt.useReferenceCounts) {
#ifdef DEBUG_MARK_SWEEP
                std::cerr << "Sweeping cache bits for forest "
                          << FID() << "\n";
#endif
                nodeHeaders.sweepAllInCacheBits();
            }
        }

        /// Is the implicit bit set for the given node?
        inline bool isImplicit(node_handle p) const {
            return nodeHeaders.getNodeImplicitFlag(p);
        }

        /// Get a node's address
        inline node_address getNodeAddress(node_handle p) const {
            return nodeHeaders.getNodeAddress(p);
        }

    // ------------------------------------------------------------
    protected: // protected node header getters/setters
    // ------------------------------------------------------------

        /// Set the level of a node
        inline void setNodeLevel(node_handle p, int level) {
            nodeHeaders.setNodeLevel(p, level);
        }

        inline void setNodeAddress(node_handle p, node_address a) {
            nodeHeaders.setNodeAddress(p, a);
        }

    // ------------------------------------------------------------
    protected:  // Eventually: make this private?
    // ------------------------------------------------------------
        /// Node header information
        node_headers nodeHeaders;

        /// Used for mark & sweep
        node_marker* reachable;


    // ------------------------------------------------------------
    public: // node storage getters/setters
    // ------------------------------------------------------------

        inline const node_storage* getNodeManager() const {
            return nodeMan;
        }

        /// Get next node in the unique table chain
        inline node_handle getNext(node_handle p) const {
            return nodeMan->getNextOf(getNodeAddress(p));
        }

        /// Set the next node in the unique table chain
        inline void setNext(node_handle p, node_handle n) {
            nodeMan->setNextOf(getNodeAddress(p), n);
        }

        /// Compute a hash for a node.
        inline unsigned hashNode(node_handle p) const {
            return nodeMan->hashNode(getNodeLevel(p), getNodeAddress(p));
        }


        /** Check and find the index of a single downward pointer.

            @param  node    Node we care about
            @param  down    Output:
                            The singleton downward pointer, or undefined.

            @return     If the node has only one non-zero downward pointer,
                        then return the index for that pointer.
                        Otherwise, return a negative value.
        */
        inline int getSingletonIndex(node_handle p, node_handle &down) const {
            return nodeMan->getSingletonIndex(getNodeAddress(p), down);
        }


        /** Check and get a single downward pointer.

            @param  node    Node we care about
            @param  index   Index we're trying to match

            @return     If the only non-zero downward pointer for
                        this node happens at \a index, then return the pointer.
                        Otherwise, return 0.
        */
        inline node_handle getSingletonDown(node_handle node, int index) const {
            MEDDLY::node_handle down;
            if (getSingletonIndex(node, down) == index) return down;
            return 0;
        }


        /** For a given node, get a specified downward pointer.

            This is designed to be used for one or two indexes only.
            For reading all or several downward pointers, an
            unpacked_node should be used instead.

            @param  p       Node to look at
            @param  index   Index of the pointer we want.

            @return         The downward pointer at that index.
        */
        inline node_handle getDownPtr(node_handle p, int index) const {
            return nodeMan->getDownPtr(getNodeAddress(p), index);
        }


        /** For a given node, get a specified downward pointer.

            This is designed to be used for one or two indexes only.
            For reading all or several downward pointers, an
            unpacked_node should be used instead.

            @param  p       Node to look at
            @param  index   Index of the pointer we want.

            @param  ev      Output: edge value at that index.
            @param  dn      Output: downward pointer at that index.
        */
        inline void getDownPtr(node_handle p, int index, int& ev,
                node_handle& dn) const
        {
            nodeMan->getDownPtr(getNodeAddress(p), index, ev, dn);
        }


        /** For a given node, get a specified downward pointer.

            This is designed to be used for one or two indexes only.
            For reading all or several downward pointers, an
            unpacked_node should be used instead.

            @param  p       Node to look at
            @param  index   Index of the pointer we want.

            @param  ev      Output: edge value at that index.
            @param  dn      Output: downward pointer at that index.
        */
        inline void getDownPtr(node_handle p, int index, long& ev,
                node_handle& dn) const
        {
            nodeMan->getDownPtr(getNodeAddress(p), index, ev, dn);
        }


        /** For a given node, get a specified downward pointer.

            This is designed to be used for one or two indexes only.
            For reading all or several downward pointers, a
            unpacked_node should be used instead.

            @param  p       Node to look at
            @param  index   Index of the pointer we want.

            @param  ev      Output: edge value at that index.
            @param  dn      Output: downward pointer at that index.
        */
        inline void getDownPtr(node_handle p, int index, float& ev,
                node_handle& dn) const
        {
            nodeMan->getDownPtr(getNodeAddress(p), index, ev, dn);
        }


        /** Check if an unpacked node duplicates one in the forest.
                @param  p       Handle to a node in the forest.
                @param  nr      Unpacked node to check.

                @return   true, iff the nodes are duplicates.
        */
        inline bool areDuplicates(node_handle p, const unpacked_node &nr) const
        {
            MEDDLY_DCASSERT(p > 0);
            if (nodeHeaders.getNodeLevel(p) != nr.getLevel()) {
                return false;
            }
            return nodeMan->areDuplicates(nodeHeaders.getNodeAddress(p), nr);
        }

    // ------------------------------------------------------------
    protected: // node storage members
    // ------------------------------------------------------------

        /// Class that stores nodes.
        node_storage* nodeMan;

    // ------------------------------------------------------------
    public: // Getting/setting the global default policies
    // ------------------------------------------------------------

        /**
            Get the default policies for all MDD forests.
        */
        inline static const policies& getDefaultPoliciesMDDs() {
            return mddDefaults;
        }

        /**
            Get the default policies for all MxD forests.
        */
        inline static const policies& getDefaultPoliciesMXDs() {
            return mxdDefaults;
        }

        /**
            Set the default policies for all MDD forests.
        */
        inline static void setDefaultPoliciesMDDs(const policies &p) {
            mddDefaults = p;
        }

        /**
            Set the default policies for all MxD forests.
        */
        inline static void setDefaultPoliciesMXDs(const policies &p) {
            mxdDefaults = p;
        }

    // ------------------------------------------------------------
    public: // Getting/setting policies for this forest
    // ------------------------------------------------------------

        /// Returns the current policies used by this forest.
        inline const policies& getPolicies() const {
            return deflt;
        }

        /// Returns the current policies used by this forest.
        inline policies& getPolicies() {
            return deflt;
        }

#ifdef ALLOW_DEPRECATED_0_17_3
        inline bool isVarSwap() const {
    	    return deflt.isVarSwap();
        }

        inline bool isLevelSwap() const {
    	    return deflt.isLevelSwap();
        }

        /// Returns the storage mechanism used by this forest.
        inline node_storage_flags getNodeStorage() const {
            return deflt.storage_flags;
        }

        /// Returns the node deletion policy used by this forest.
        inline policies::node_deletion getNodeDeletion() const {
            return deflt.deletion;
        }

        /// Are we using pessimistic deletion
        inline bool isPessimistic() const {
            return deflt.isPessimistic();
        }

        /// Can we store nodes sparsely
        inline bool areSparseNodesEnabled() const {
            return deflt.allowsSparse();
        }

        /// Can we store nodes fully
        inline bool areFullNodesEnabled() const {
            return deflt.allowsFull();
        }
#endif

    // ------------------------------------------------------------
    protected: // methods to set edge info, for derived classes
    // ------------------------------------------------------------

        // Call one of these first...

        inline void setVoidEdges() {
            switch (rangeType) {
                case range_type::BOOLEAN:
                        the_terminal_type = terminal_type::BOOLEAN;
                        break;
                case range_type::INTEGER:
                        the_terminal_type = terminal_type::INTEGER;
                        break;
                case range_type::REAL:
                        the_terminal_type = terminal_type::REAL;
                        break;
                default:
                        MEDDLY_DCASSERT(false);
            }
            the_edge_type = edge_type::VOID;
            hash_edge_values = false;
        }
        inline void setIntEdges(bool hashed = true) {
            the_terminal_type = terminal_type::OMEGA;
            the_edge_type = edge_type::INT;
            hash_edge_values = hashed;
        }
        inline void setLongEdges(bool hashed = true) {
            the_terminal_type = terminal_type::OMEGA;
            the_edge_type = edge_type::LONG;
            hash_edge_values = hashed;
        }
        inline void setFloatEdges(bool hashed = false) {
            the_terminal_type = terminal_type::OMEGA;
            the_edge_type = edge_type::FLOAT;
            hash_edge_values = hashed;
        }
        inline void setDoubleEdges(bool hashed = false) {
            the_terminal_type = terminal_type::OMEGA;
            the_edge_type = edge_type::DOUBLE;
            hash_edge_values = hashed;
        }

        // ... then one of these

        inline void setTransparentEdge(node_handle p) {
            MEDDLY_DCASSERT(edge_type::VOID == the_edge_type);
            transparent_node = p;
            transparent_edge.set();
        }
        inline void setTransparentEdge(node_handle p, int v) {
            MEDDLY_DCASSERT(edge_type::INT == the_edge_type);
            transparent_node = p;
            transparent_edge.set(v);
        }
        inline void setTransparentEdge(node_handle p, long v) {
            MEDDLY_DCASSERT(edge_type::LONG == the_edge_type);
            transparent_node = p;
            transparent_edge.set(v);
        }
        inline void setTransparentEdge(node_handle p, float v) {
            MEDDLY_DCASSERT(edge_type::FLOAT == the_edge_type);
            transparent_node = p;
            transparent_edge.set(v);
        }
        inline void setTransparentEdge(node_handle p, double v) {
            MEDDLY_DCASSERT(edge_type::DOUBLE == the_edge_type);
            transparent_node = p;
            transparent_edge.set(v);
        }


    // ------------------------------------------------------------
    private: // transparent edge info
    // ------------------------------------------------------------

        /// Transparent node value
        node_handle transparent_node;

        /// Transparent edge value
        edge_value  transparent_edge;

        /// Edge type.
        edge_type the_edge_type;

        /// Terminal node type.
        terminal_type the_terminal_type;

        /// Are edge values hashed?
        bool hash_edge_values;



    public:
        /// Returns the forest identifier, a unique positive integer per forest.
        /// FID 0 can safely be used to indicate "no forest".
        inline unsigned FID() const { return fid; }

        /// Returns the largest forest identifier ever seen.
        static inline unsigned MaxFID() { return gfid; }

        static inline forest* getForestWithID(unsigned id) {
            if (id >= max_forests) return nullptr;
            return all_forests[id];
        }

    /*
     *  Methods and members for keeping track of root edges.
     *  They should normally not be called directly.
     *
     */

    public:
        /// Register a dd_edge with this forest.
        /// Called automatically in dd_edge.
        void registerEdge(dd_edge& e);

        /// Unregister a dd_edge with this forest.
        /// Called automatically in dd_edge.
        void unregisterEdge(dd_edge& e);

        /// For debugging: count number of registered dd_edges.
        unsigned countRegisteredEdges() const;

        /// Mark all registered dd_edges.
        void markAllRoots();

    private:
        void unregisterDDEdges();

    private:
        /// Registry of dd_edges
        dd_edge* roots;


    /*
     *  Methods and members for the global forest registry.
     *
     */

    private:
        static void initStatics();
        static void freeStatics();
        static void registerForest(forest* f);
        static void unregisterForest(forest* f);

    private:
        // All forests
        static forest** all_forests;
        // Size of forests array
        static unsigned max_forests;

        // our ID
        unsigned fid;

        // global ID
        static unsigned gfid;

        friend class initializer_list;

// ===================================================================
//
// To be cleaned up still, below here.
//
// ===================================================================

  public:
    int *level_reduction_rule;

    /// Status indicators for nodes.
    enum node_status {
      /// Node is active: it can be used without issue.
      ACTIVE,
      /// Node is recoverable: it has been marked for garbage collection
      /// but it can still be used without issue.
      RECOVERABLE,
      /// Node is not-recoverable: it has been marked for garbage collection
      /// and it cannot be used.
      DEAD
    };



    /// Collection of various stats for performance measurement
    struct statset {
      /// Number of times we scanned for reachable nodes
      long reachable_scans;
      /// Number of times a dead node was resurrected.
      long reclaimed_nodes;
      /// Number of times the forest storage array was compacted
      long num_compactions;
      /// Number of times the garbage collector ran.
      long garbage_collections;

#ifdef TRACK_UNREACHABLE_NODES
      /// Current number of unreachable (disconnected) nodes
      long unreachable_nodes;
#endif

      /// Current number of connected nodes
      long active_nodes;
      /// Current number of temporary nodes
      // long temp_nodes;

      /// Peak number of active nodes
      long peak_active;

      statset();

      // nice helpers

      void incActive(long b);
      void decActive(long b);
    };

    /**
        Abstract base class for logging of forest stats.
        The idea is to allow an external tool to read this information,
        and display stats in real time.  Potentially much more detailed
        that the stats collected above.
    */
    class logger {
        bool nfix;
        bool node_counts;
        bool time_stamps;

        long startsec;
        long startusec;
      public:
        logger();
        virtual ~logger();

      /*
          Settings.
          Cannot be changed once we attach to a forest.
      */
      public:
        bool recordingNodeCounts() const;
        void recordNodeCounts();
        void ignoreNodeCounts();

        bool recordingTimeStamps() const;
        void recordTimeStamps();
        void ignoreTimeStamps();

      /*
          Hooks, used in various places.
          Must be overridden in derived classes.
      */
      public:
        /**
            Insert a comment string.
            May be ignored depending on the file format.

              @param  comment   String to insert.
        */
        virtual void addComment(const char* comment) = 0;

        /**
            Start a new phase of computation.
            May be ignored depending on the file format.

              @param  f         Forest this applies to.
              @param  comment   Info to display about this phase.
        */
        virtual void newPhase(const forest* f, const char* comment) = 0;

        /**
            Called once, when the logger is attached to a forest.
            Must call method fixLogger().

              @param  f       Forest info to log
              @param  name    Forest name; better for displaying if we
                              have multiple forests.
        */
        virtual void logForestInfo(const forest* f, const char* name) = 0;

        /**
            Change node counts by specified amounts.
        */
        virtual void addToActiveNodeCount(const forest* f, int level, long delta) = 0;

      /*
          Helpers for derived classes
      */
      protected:
        /* Use this for generating timestamps. */
        void currentTime(long &sec, long &usec);

        /* Call this inside logForestInfo() */
        void fixLogger();
    };


  protected:
    /** Constructor -- this class cannot be instantiated.
      @param  d       domain to which this forest belongs to.
      @param  rel     does this forest represent a relation.
      @param  t       the range of the functions represented in this forest.
      @param  ev      edge annotation.
      @param  p       Polcies for reduction, storage, deletion.
      @param  level_reduction_rule       Rules for reduction on different levels.
    */
    forest(domain* d, bool rel, range_type t, edge_labeling ev,
      const policies &p, int* level_reduction_rule);

    /// Destructor.
    virtual ~forest();

  // ------------------------------------------------------------
  // static inlines.
  public:
    /** Go "down a level" in a relation.
        Safest to use this, in case in later versions
        the level numbering changes, or becomes forest dependent.
            @param  k   Current level
            @return Level immediately below level k
    */
    static inline int downLevel(int k) {
        return (k>0) ? (-k) : (-k-1);
    }


    /** Go "up a level" in a relation.
        Safest to use this, in case in later versions
        the level numbering changes, or becomes forest dependent.
            @param  k   Current level
            @return Level immediately above level k
    */
    static inline int upLevel(int k) {
        return (k<0) ? (-k) : (-k-1);
    }


  // ------------------------------------------------------------
  // inlines.
  public:

    /// Does this forest represent relations or matrices?
    bool isForRelations() const;

    /// Returns the range type.
    range_type getRangeType() const;

    /// Query the range type.
    bool isRangeType(range_type r) const;

    /// Returns the edge labeling mechanism.
    edge_labeling getEdgeLabeling() const;

    /// Is the edge labeling "MULTI_TERMINAL".
    bool isMultiTerminal() const;

    /// Is the edge labeling "EV_PLUS".
    bool isEVPlus() const;

    /// Is the edge labeling "INDEX_SET".
    bool isIndexSet() const;

    /// Is the edge labeling "EV_TIMES".
    bool isEVTimes() const;

    /// Check if we match a specific type of forest
    bool matches(bool isR, range_type rt, edge_labeling el) const;

    /// Returns the reduction rule used by this forest.
    reduction_rule getReductionRule() const;

    /// Returns true if the forest is fully reduced.
    bool isFullyReduced() const;

    /// Returns true if the forest is quasi reduced.
    bool isQuasiReduced() const;

    /// Returns true if the forest is identity reduced.
    bool isIdentityReduced() const;

     /// Returns true if the forest is user_defined reduced.
    bool isUserDefinedReduced() const;

    /// Returns true if the level is fully reduced.
    bool isFullyReduced(int k) const;

    /// Returns true if the level is quasi reduced.
    bool isQuasiReduced(int k) const;

    /// Returns true if the level is identity reduced.
    bool isIdentityReduced(int k) const;

    int* getLevelReductionRule() const;




    /// Get forest performance stats.
    const statset& getStats() const;

    /// Get forest memory stats.
    const memstats& getMemoryStats() const;

    /** Get the current number of nodes in the forest, at all levels.
        @return     The current number of nodes, not counting deleted or
                    marked for deletion nodes.
    */
    long getCurrentNumNodes() const;

    /** Get the current total memory used by the forest.
        This should be equal to summing getMemoryUsedForVariable()
        over all variables.
        @return     Current memory used by the forest.
    */
    size_t getCurrentMemoryUsed() const;

    /** Get the current total memory allocated by the forest.
        This should be equal to summing getMemoryAllocatedForVariable()
        over all variables.
        @return     Current total memory allocated by the forest.
    */
    size_t getCurrentMemoryAllocated() const;

    /** Get the peak number of nodes in the forest, at all levels.
        This will be at least as large as calling getNumNodes() after
        every operation and maintaining the maximum.
        @return     The peak number of nodes that existed at one time,
                    in the forest.
    */
    long getPeakNumNodes() const;

    /** Set the peak number of nodes to the number current number of nodes.
    */
    inline void resetPeakNumNodes() {
      stats.peak_active = stats.active_nodes;
    }

    /** Get the peak memory used by the forest.
        @return     Peak total memory used by the forest.
    */
    size_t getPeakMemoryUsed() const;

	/** Set the peak memory to the current memory.
	*/
  /*
	inline void resetPeakMemoryUsed() {
	  stats.peak_memory_used = stats.memory_used;
	}
  */

    /** Get the peak memory allocated by the forest.
        @return     Peak memory allocated by the forest.
    */
    size_t getPeakMemoryAllocated() const;

    /// Are we about to be deleted?
    bool isMarkedForDeletion() const;

  // ------------------------------------------------------------
  // non-virtual.
  public:
    /// Remove any stale compute table entries associated with this forest.
    void removeStaleComputeTableEntries();

    /// Remove all compute table entries associated with this forest.
    void removeAllComputeTableEntries();

  // ------------------------------------------------------------
  // virtual, with default implementation.
  public:

   virtual node_handle unionOneMinterm(node_handle a,  int* from,  int* to, int level);

    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = terms[i] for 0 <= i < size(vh).

        For example, in a forest with range_type BOOLEAN, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
        X   0   X  terms[0]
        X   1   X  terms[1]
        X   2   X  terms[2]
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  vp
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  terms Array of boolean terminal values.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = terms[i]
                      for 0 <= i < size(vh).

        @throws       TYPE_MISMATCH, if the forest's range is not BOOLEAN.
    */
    virtual void createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a);

    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = terms[i] for 0 <= i < size(vh).

        For example, in a forest with range_type INTEGER, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
        X   0   X  terms[0]
        X   1   X  terms[1]
        X   2   X  terms[2]
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  vp
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  terms Array of boolean terminal values.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = terms[i]
                      for 0 <= i < size(vh).

        @throws       TYPE_MISMATCH, if the forest's range is not INTEGER.
    */
    virtual void createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a);

    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = terms[i] for 0 <= i < size(vh).

        For example, in a forest with range_type REAL, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
        X   0   X  terms[0]
        X   1   X  terms[1]
        X   2   X  terms[2]
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  vp
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  terms Array of boolean terminal values.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = terms[i]
                      for 0 <= i < size(vh).

        @throws       TYPE_MISMATCH, if the forest's range is not REAL.
    */
    virtual void createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a);

    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = i for 0 <= i < size(vh).

        For example, in a forest with range_type INTEGER, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
        X   0   X  0
        X   1   X  1
        X   2   X  2
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  pr
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = i for 0 <= i < size(vh).
    */
    void createEdgeForVar(int vh, bool pr, dd_edge& a);

    /** Create an edge as the union of several explicit vectors.
        @param  vlist Array of vectors. Each vector has dimension equal
                      to one plus the largest variable handle in the domain.
                      A vector \a x indicates a set of variable assignments,
                      where x[vh] may have a value of DONT_CARE;
                      otherwise x[vh] gives the variable assignment for vh.
        @param  N     Number of vectors (dimension of \a vlist).
        @param  e     returns a handle to a node in the forest, such that
                      f(v_1, ..., v_n) = 1, iff there is a vector
                      x in \a vlist corresponding to the variable assignments
                      v_1, ..., v_n; f(v_1, ..., v_n) = 0, otherwise.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not BOOLEAN,
                        or the forest is for relations.
    */
    virtual void createEdge(const int* const* vlist, int N, dd_edge &e);

    /** Create an edge as the union of several vectors and return values.
        @param  vlist Array of vectors. Each vector has dimension equal
                      to one plus the largest variable handle in the domain.
                      A vector \a x indicates a set of variable assignments,
                      where x[vh] may have a value of DONT_CARE;
                      otherwise x[vh] gives the variable assignment for vh.
        @param  terms Array of return values, same dimension as \a vlist.
        @param  N     Number of vectors (dimension of \a vlist).
        @param  e     returns a handle to a node in the forest that encodes
                      the function f: f(v_1, ..., v_n) = terms[j]
                      iff j is the smallest integer such that vector vlist[j]
                      corresponds to the variable assignments v_1, ..., v_n;
                      f(v_1, ..., v_n) = 0, otherwise.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not INTEGER,
                        or the forest is for relations.
    */
    virtual void createEdge(const int* const* vlist, const long* terms, int N, dd_edge &e);

    /** Create an edge as the union of several vectors and return values.
        @param  vlist Array of vectors. Each vector has dimension equal
                      to one plus the largest variable handle in the domain.
                      A vector \a x indicates a set of variable assignments,
                      where x[vh] may have a value of DONT_CARE;
                      otherwise x[vh] gives the variable assignment for vh.
        @param  terms Array of return values, same dimension as \a vlist.
        @param  N     Number of vectors (dimension of \a vlist).
        @param  e     returns a handle to a node in the forest that encodes
                      the function f: f(v_1, ..., v_n) = terms[j]
                      iff j is the smallest integer such that vector vlist[j]
                      corresponds to the variable assignments v_1, ..., v_n;
                      f(v_1, ..., v_n) = 0, otherwise.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not REAL,
                        or the forest is for relations.
    */
    virtual void createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e);


    /** Create an edge as the union of several explicit matrices.
        @param  vlist   Array of vectors. Each vector has dimension equal to
                        one plus the largest variable handle in the domain.
                        Vector vlist indicates the set of unprimed variable
                        assignments, and vector vplist indicates the set of
                        primed variable assignments.  Both vlist and vplist
                        may have elements equal to DONT_CARE, and vplist
                        may have elements equal to DONT_CHANGE.
        @param  vplist  Array of vectors, same dimension as \a vlist.
                        A vector \a x = vplist[i] indicates a set of primed
                        variable assignments, where x[vh] less than 0 means
                        "don't change for variable vh" iff we are not
                        changing the unprimed variable, and
                        "don't care for primed variable vh" otherwise.
                        Otherwise x[vh] gives the variable assignment for
                        primed variable vh.
        @param  N       Number of vectors (dimension of \a vlist).
        @param  e       returns a handle to a node in the forest, such that
                        f(v_1, v'_1, ..., v_n, v'_n) = 1, iff there is a
                        vector vlist[i] corresponding to the variable
                        assignments v_1, ..., v_n, and vplist[i] corresponds
                        to variable assignments v'_1, ..., v'_n.
                        f(v_1, v'_1, ..., v_n, v'_n) = 0, otherwise.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not BOOLEAN,
                        or the forest is not for relations.
    */
    virtual void createEdge(const int* const* vlist, const int* const* vplist, int N, dd_edge &e);


    /** Create an edge as the union of several explicit matrices.
        @param  vlist   Array of vectors. Each vector has dimension equal to
                        one plus the largest variable handle in the domain.
                        Vector vlist indicates the set of unprimed variable
                        assignments, and vector vplist indicates the set of
                        primed variable assignments.  Both vlist and vplist
                        may have elements equal to DONT_CARE, and vplist
                        may have elements equal to DONT_CHANGE.
        @param  vplist  Array of vectors, same dimension as \a vlist.
                        A vector \a x = vplist[i] indicates a set of primed
                        variable assignments, where x[vh] less than 0 means
                        "don't change for variable vh" iff we are not
                        changing the unprimed variable, and
                        "don't care for primed variable vh" otherwise.
                        Otherwise x[vh] gives the variable assignment for
                        primed variable vh.
        @param  terms   Array of return values, same dimension as \a vlist.
        @param  N       Number of vectors (dimension of \a vlist).
        @param  e       returns a handle to a node in the forest, such that
                        f(v_1, v'_1, ..., v_n, v'_n) = term[i], iff there is
                        a vector vlist[i] corresponding to the variable
                        assignments v_1, ..., v_n, and vplist[i] corresponds
                        to variable assignments v'_1, ..., v'_n.
                        f(v_1, v'_1, ..., v_n, v'_n) = 0, otherwise.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not INTEGER,
                        or the forest is not for relations.
    */
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const long* terms, int N, dd_edge &e);


    /** Create an edge as the union of several explicit matrices.
        @param  vlist   Array of vectors. Each vector has dimension equal to
                        one plus the largest variable handle in the domain.
                        Vector vlist indicates the set of unprimed variable
                        assignments, and vector vplist indicates the set of
                        primed variable assignments.  Both vlist and vplist
                        may have elements equal to DONT_CARE, and vplist
                        may have elements equal to DONT_CHANGE.
        @param  vplist  Array of vectors, same dimension as \a vlist.
                        A vector \a x = vplist[i] indicates a set of primed
                        variable assignments, where x[vh] less than 0 means
                        "don't change for variable vh" iff we are not
                        changing the unprimed variable, and
                        "don't care for primed variable vh" otherwise.
                        Otherwise x[vh] gives the variable assignment for
                        primed variable vh.
        @param  terms   Array of return values, same dimension as \a vlist.
        @param  N       Number of vectors (dimension of \a vlist).
        @param  e       returns a handle to a node in the forest, such that
                        f(v_1, v'_1, ..., v_n, v'_n) = term[i], iff there is
                        a vector vlist[i] corresponding to the variable
                        assignments v_1, ..., v_n, and vplist[i] corresponds
                        to variable assignments v'_1, ..., v'_n.
                        f(v_1, v'_1, ..., v_n, v'_n) = 0, otherwise.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not REAL,
                        or the forest is not for relations.
    */
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge &e);


    /** Create an edge for a boolean constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not BOOLEAN.
    */
    virtual void createEdge(bool val, dd_edge &e);

    /** Create an edge for an integer constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not BOOLEAN.
    */
    virtual void createEdge(long val, dd_edge &e);

    /** Create an edge for an integer constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not BOOLEAN.
    */
    virtual void createEdge(float val, dd_edge &e);

    /** Evaluate the function encoded by an edge.
        @param  f     Edge (function) to evaluate.
        @param  vlist List of variable assignments, of dimension one higher
                      than the largest variable handle.
        @param  term  Output parameter, will be set to
                      f(v1 = vlist[v1], ..., vn = vlist[vn]).

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not BOOLEAN,
                        or the forest is for relations.
    */
    virtual void evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;

    /** Evaluate the function encoded by an edge.
        @param  f     Edge (function) to evaluate.
        @param  vlist List of variable assignments, of dimension one higher
                      than the largest variable handle.
        @param  term  Output parameter, will be set to
                      f(v1 = vlist[v1], ..., vn = vlist[vn]).

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not INTEGER,
                        or the forest is for relations.
    */
    virtual void evaluate(const dd_edge &f, const int* vlist, long &term)
      const;

    /** Evaluate the function encoded by an edge.
        @param  f     Edge (function) to evaluate.
        @param  vlist List of variable assignments, of dimension one higher
                      than the largest variable handle.
        @param  term  Output parameter, will be set to
                      f(v1 = vlist[v1], ..., vn = vlist[vn]).

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not REAL,
                        or the forest is for relations.
    */
    virtual void evaluate(const dd_edge &f, const int* vlist, float &term)
      const;

    /** Evaluate the function encoded by an edge.
        @param  f       Edge (function) to evaluate.
        @param  vlist   List of variable assignments for unprimed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  vplist  List of variable assignments for primed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  term    Output parameter, will be set to
                        f(v1 = vlist[v1], v'1 = vplist[v1], ...,
                        vn = vlist[vn], v'n = vplist[vn]).

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not BOOLEAN,
                        or the forest is not for relations.
    */
    virtual void evaluate(const dd_edge& f, const int* vlist,
      const int* vplist, bool &term) const;

    /** Evaluate the function encoded by an edge.
        @param  f       Edge (function) to evaluate.
        @param  vlist   List of variable assignments for unprimed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  vplist  List of variable assignments for primed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  term    Output parameter, will be set to
                        f(v1 = vlist[v1], v'1 = vplist[v1], ...,
                        vn = vlist[vn], v'n = vplist[vn]).

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not INTEGER,
                        or the forest is not for relations.
    */
    virtual void evaluate(const dd_edge& f, const int* vlist,
      const int* vplist, long &term) const;

    /** Evaluate the function encoded by an edge.
        @param  f       Edge (function) to evaluate.
        @param  vlist   List of variable assignments for unprimed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  vplist  List of variable assignments for primed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  term    Output parameter, will be set to
                        f(v1 = vlist[v1], v'1 = vplist[v1], ...,
                        vn = vlist[vn], v'n = vplist[vn]).

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not REAL,
                        or the forest is not for relations.
    */
    virtual void evaluate(const dd_edge& f, const int* vlist,
      const int* vplist, float &term) const;

    /** Returns element \a e at index \a i from an Index Set EV+MDD.

        size(e) = number of variables in the forest + 1 (for terminals).
        TODO: complete this description

        on return: e[0] will be 1 if the element could be found, 0 otherwise.

        @throws       INVALID_OPERATION, if this is not an Index Set EV+MDD.
    */
    virtual void getElement(const dd_edge& a, int index, int* e);
    virtual void getElement(const dd_edge& a, long index, int* e);

  // ------------------------------------------------------------
  // abstract virtual.
  public:
    /** Write edges to a file in a format that can be read back later.
        This implies that all nodes "below" those edges are also
        written to the file.
          @param  s   Stream to write to
          @param  E   Array of edges
          @param  n   Dimension of the edge array

          @throws     COULDNT_WRITE, if writing failed
    */
    virtual void writeEdges(output &s, const dd_edge* E, unsigned n) const = 0;

    /** Read edges from a file.
        Allows reconstruction of edges that we
        saved using \a writeEdges().
        The forest does not need to be empty;
        the edges are added to the forest as necessary.
          @param  s   Stream to read from
          @param  E   Array of edges
          @param  n   Dimension of the edge array

          @throws     INVALID_FILE, if the file does not match what we expect,
                      including different number of edges specified.
    */
    virtual void readEdges(input &s, dd_edge* E, unsigned n) = 0;

//        virtual void readEdgeValue(input &s, dd_edge &E) const = 0;
        // virtual void writeEdgeValue(output &s, const dd_edge &E) const = 0;
        // virtual void showEdgeValue(output &s, const dd_edge &E) const = 0;

        // void readEdge(input &s, dd_edge &E, const node_handle* map);
        // void writeEdge(output &s, const dd_edge &E, const node_handle* map) const;

        /**
            Show an edge, compactly.
            Called for example for each child when displaying an entire node.
                @param  s       Stream to write to.
                @param  ev      Edge value
                @param  d       Down pointer
        */
        virtual void showEdge(output &s, const edge_value &ev, node_handle d) const = 0;

        /** Write an edge value in machine-readable format.
            @param  s       Stream to write to.
            @param  edge    Edge value
        */
        // virtual void writeEdgeValue(output &s, const edge_value &edge) const;

        /** Read an edge value in machine-readable format.
            @param  s       Stream to read from.
            @param  edge    Edge value
        */
        // virtual void readEdgeValue(input &s, edge_value &edge);

private:
    /**
        TBD: virtual for now;
        eventually the implementation will move here.

        Disconnects all downward pointers from p,
        and removes p from the unique table.
    */
    virtual void deleteNode(node_handle p) = 0;

    friend class node_headers;

public:

    /** Force garbage collection.
        All disconnected nodes in this forest are discarded along with any
        compute table entries that may include them.
    */
    // virtual void garbageCollect() = 0;

    /** Compact the memory for all variables in this forest.
        This is not the same as garbage collection.
    */
    // virtual void compactMemory() = 0;

    /** Create an edge representing the subset of a Matrix Diagram.

        size(vlist) == number of variables in the forest + 1 (for terminals)
        size(vlist[i]) == size (i.e. bound) of variable i
        size(vlist) == size(vplist)
        size(vlist[i]) == size(vplist[i])
        If vlist[i][j] is true, that index is included in the mask
        If vlist[i][j] is false, that index is excluded from the mask.

        TODO: write a better description (an example might also help).
    */
    /*
    virtual void createSubMatrix(const bool* const* vlist,
      const bool* const* vplist, const dd_edge a, dd_edge& b) = 0;
      */




    /** Display all active (i.e., connected) nodes in the forest.
        This is primarily for aid in debugging.
        @param  strm      Stream to write to.
        @param  verbosity How much information to display.
                          0 : just statistics.
                          1 : all forest nodes + statistics.
                          2 : internal forest + statistics.
    */
    virtual void showInfo(output &strm, int verbosity=0) = 0;


    /** Start logging stats.
          @param  L     Logger to use; if 0, we don't log anything.
                        Will overwrite the old logger.
                        The logger WILL NOT be deleted by the forest
                        (in case we want to log for multiple forests).

          @param  name  Name to use for the forest (for display only).
    */
    void setLogger(logger* L, const char* name);


  // ------------------------------------------------------------
  // For derived classes.
  protected:
    // for debugging:
    void showComputeTable(output &s, int verbLevel) const;

  protected:
    policies deflt;
    statset stats;
    memstats mstats;
    logger *theLogger;

  protected:
    /// uniqueness table, still used by derived classes.
    unique_table* unique;

    /// uniqueness table for relation nodes.
    impl_unique_table* implUT;

    /// Should a terminal node be considered a stale entry in the compute table.
    /// per-forest policy, derived classes may change as appropriate.
    MEDDLY::forest::node_status terminalNodesStatus;



  // ------------------------------------------------------------
  // Ugly details from here down.
  private:  // Defaults
    friend class forest_initializer;
    static policies mddDefaults;
    static policies mxdDefaults;

  private:  // For operation registration
    friend class operation;

    unsigned* opCount;
    unsigned szOpCount;

    /// Register an operation with this forest.
    /// Called only within operation.
    void registerOperation(const operation* op);

    /// Unregister an operation.
    /// Called only within operation.
    void unregisterOperation(const operation* op);

  private:
    bool isRelation;
    bool is_marked_for_deletion;
    range_type rangeType;
    edge_labeling edgeLabel;

    /// Mark for deletion.
    void markForDeletion();

    friend void MEDDLY::destroyForest(MEDDLY::forest* &f);


    // We should be able to remove this after updating unpacked_node
    //
    friend class unpacked_node;
};


// ******************************************************************
// *                                                                *
// *                     inlined forest methods                     *
// *                                                                *
// ******************************************************************



// forest::statset::
inline void MEDDLY::forest::statset::incActive(long b) {
  active_nodes += b;
  if (active_nodes > peak_active)
    peak_active = active_nodes;
  MEDDLY_DCASSERT(active_nodes >= 0);
}
inline void MEDDLY::forest::statset::decActive(long b) {
  active_nodes -= b;
  MEDDLY_DCASSERT(active_nodes >= 0);
}
// end of forest::statset

// forest::logger::
inline bool MEDDLY::forest::logger::recordingNodeCounts() const {
  return node_counts;
}

inline void MEDDLY::forest::logger::recordNodeCounts() {
  if (nfix) node_counts = true;
}

inline void MEDDLY::forest::logger::ignoreNodeCounts() {
  if (nfix) node_counts = false;
}

inline bool MEDDLY::forest::logger::recordingTimeStamps() const {
  return time_stamps;
}

inline void MEDDLY::forest::logger::recordTimeStamps() {
  if (nfix) time_stamps = true;
}

inline void MEDDLY::forest::logger::ignoreTimeStamps() {
  if (nfix) time_stamps = false;
}

inline void MEDDLY::forest::logger::fixLogger() {
  nfix = false;
}

// end of forest::logger

// forest::

inline bool MEDDLY::forest::isForRelations() const {
  return isRelation;
}

inline MEDDLY::range_type MEDDLY::forest::getRangeType() const {
  return rangeType;
}

inline bool MEDDLY::forest::isRangeType(MEDDLY::range_type r) const {
  return r == rangeType;
}

inline MEDDLY::edge_labeling MEDDLY::forest::getEdgeLabeling() const {
  return edgeLabel;
}

inline bool MEDDLY::forest::isMultiTerminal() const {
  return edge_labeling::MULTI_TERMINAL == edgeLabel;
}

inline bool MEDDLY::forest::isEVPlus() const {
  return edge_labeling::EVPLUS == edgeLabel;
}

inline bool MEDDLY::forest::isIndexSet() const {
  return edge_labeling::INDEX_SET == edgeLabel;
}

inline bool MEDDLY::forest::isEVTimes() const {
  return edge_labeling::EVTIMES == edgeLabel;
}

inline bool MEDDLY::forest::matches(bool isR, MEDDLY::range_type rt,
  MEDDLY::edge_labeling el) const {
  return (isRelation == isR) && (rangeType == rt) && (edgeLabel == el);
}

inline MEDDLY::reduction_rule MEDDLY::forest::getReductionRule() const {
  return deflt.reduction;
}


inline bool MEDDLY::forest::isFullyReduced() const {
    return deflt.isFullyReduced();
}

inline bool MEDDLY::forest::isQuasiReduced() const {
    return deflt.isQuasiReduced();
}

inline bool MEDDLY::forest::isIdentityReduced() const {
    return deflt.isIdentityReduced();
}

inline bool MEDDLY::forest::isUserDefinedReduced() const {
    return deflt.isUserDefinedReduced();
}

inline int* MEDDLY::forest::getLevelReductionRule() const{
    return level_reduction_rule;

}

inline bool MEDDLY::forest::isFullyReduced(int k) const {
    if(k<0)
        return getLevelReductionRule()[2*(-k)]==-1;
    else
        return getLevelReductionRule()[2*(k) - 1]==-1;

}

inline bool MEDDLY::forest::isQuasiReduced(int k) const {
    if(k<0)
        return getLevelReductionRule()[2*(-k)]==-2;
    else
        return getLevelReductionRule()[2*(k) - 1]==-2;
}

inline bool MEDDLY::forest::isIdentityReduced(int k) const {
    if(k<0)
        return getLevelReductionRule()[2*(-k)]==-3;
    else
        return false;
}



inline const MEDDLY::forest::statset& MEDDLY::forest::getStats() const {
  return stats;
}

inline const MEDDLY::memstats& MEDDLY::forest::getMemoryStats() const {
  return mstats;
}

inline long MEDDLY::forest::getCurrentNumNodes() const {
  return stats.active_nodes;
}

inline size_t MEDDLY::forest::getCurrentMemoryUsed() const {
  return mstats.getMemUsed();
}

inline size_t MEDDLY::forest::getCurrentMemoryAllocated() const {
  return mstats.getMemAlloc();
}

inline long MEDDLY::forest::getPeakNumNodes() const {
  return stats.peak_active;
}

inline size_t MEDDLY::forest::getPeakMemoryUsed() const {
  return mstats.getPeakMemUsed();
}

inline size_t MEDDLY::forest::getPeakMemoryAllocated() const {
  return mstats.getPeakMemAlloc();
}

inline bool MEDDLY::forest::isMarkedForDeletion() const {
  return is_marked_for_deletion;
}

inline void MEDDLY::forest::createEdgeForVar(int vh, bool pr, dd_edge& a)
{
    switch (rangeType) {
        case range_type::BOOLEAN:
                createEdgeForVar(vh, pr, (bool*)  0, a);
                break;

        case range_type::INTEGER:
                createEdgeForVar(vh, pr, (long*)  0, a);
                break;

        case range_type::REAL:
                createEdgeForVar(vh, pr, (float*) 0, a);
                break;

        default:
                throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
};

inline void MEDDLY::forest::setLogger(logger* L, const char* name) {
  theLogger = L;
  if (theLogger) theLogger->logForestInfo(this, name);
}

// end of class forest


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      expert_forest  class                      *
// *                                                                *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_forest: public MEDDLY::forest
{
  public:
    friend class reordering_base;

    /** Constructor.
      @param  d       domain to which this forest belongs to.
      @param  rel     does this forest represent a relation.
      @param  t       the range of the functions represented in this forest.
      @param  ev      edge annotation.
      @param  p       Polcies for reduction, storage, deletion.
      @param  level_reduction_rule       Rules for reduction on different levels.
    */
    expert_forest(domain *d, bool rel, range_type t,
                  edge_labeling ev, const policies &p, int* level_reduction_rule);

  // ------------------------------------------------------------
  // inlined helpers.


    memstats& changeMemStats();

    /// Extra bytes per node, not hashed.
    unsigned char unhashedHeaderBytes() const;
    /// Extra bytes per node, hashed.
    unsigned char hashedHeaderBytes() const;



  // --------------------------------------------------
  // Node status
  // --------------------------------------------------
  public:
    /*
    bool isActiveNode(node_handle p) const;
    // bool isZombieNode(node_handle p) const;
    bool isDeletedNode(node_handle p) const;
    static bool isTerminalNode(node_handle p);
    /// Sanity check: is this a valid nonterminal node index.
    bool isValidNonterminalIndex(node_handle p) const;
    /// Sanity check: is this a valid node index.
    bool isValidNodeIndex(node_handle p) const;
    node_handle getLastNode() const;
    */

    // --------------------------------------------------
    // Extensible Node Information:
    // --------------------------------------------------
    /// Is this an extensible node
    bool isExtensible(node_handle p) const;
    /// If \a p is an extensible node, this returns an index
    /// that can be used by getDownPtr(...) to access the
    /// extensible portion of this node.
    int getExtensibleIndex(node_handle p) const;

  public:

    /// Get the cardinality of an Index Set.
    int getIndexSetCardinality(node_handle node) const;



    /// A node can be discarded once it goes stale. Whether a node is
    /// considered stale depends on the forest's deletion policy.
    /// Optimistic deletion: A node is said to be stale only when both the
    ///   in-count and cache-count are zero.
    /// Pessimistic deletion: A node is said to be stale when the in-count
    ///  is zero regardless of the cache-count.
    /// If we don't use reference counts and instead mark and sweep,
    ///  then a node cannot be recovered once it is "unreachable"
    ///  because its children might have been recycled
    MEDDLY::forest::node_status getNodeStatus(node_handle node) const;

  // ------------------------------------------------------------
  // non-virtual, handy methods for debugging or logging.

    /**
        Display all nodes in the forest.
          @param  s       File stream to write to
          @param  flags   Switches to control output;
                          see constants "SHOW_DETAILED", etc.
    */
    void dump(output &s, unsigned int flags) const;
    void dumpInternal(output &s) const;
    void dumpUniqueTable(output &s) const;
    void validateIncounts(bool exact);
    void validateCacheCounts() const;
    void countNodesByLevel(long* active) const;



  // ------------------------------------------------------------
  // non-virtual, handy methods.

#ifdef ALLOW_DEPRECATED_0_17_3

    /** Build a list of nodes in the subgraph below the given node.
        DEPRECATED; use a node_marker object for this instead.

        This for example is used to determine which nodes must
        be printed to display a subgraph.
        Terminal nodes are NOT included.

          @param  roots   Array of root nodes in the forest.
                          Each root node will be included in the list,
                          except for terminal nodes.

          @param  N       Dimension of \a roots array.

          @param  sort    If true, the list will be in increasing order.
                          Otherwise, the list will be in some convenient order
                          (currently, it is the order that nodes
                          are discovered).

          @return   A malloc'd array of non-terminal nodes, terminated by 0.
                    Or, a null pointer, if the list is empty.
    */
    node_handle*
    markNodesInSubgraph(const node_handle* roots, int N, bool sort) const;

    /** Count and return the number of non-terminal nodes
        in the subgraph below the given node.
        DEPRECATED; use a node_marker object for this.
    */
    unsigned long getNodeCount(node_handle p) const;

    /** Count and return the number of non-terminal nodes
        in the subgraph below the given nodes.
        DEPRECATED; use a node_marker object for this.
    */
    unsigned long getNodeCount(const node_handle* roots, int N) const;

    /** Count and return the number of edges
        in the subgraph below the given node.
        DEPRECATED; use a node_marker object for this.
    */
    unsigned long getEdgeCount(node_handle p, bool countZeroes) const;

    /// Show all the nodes in the subgraph below the given nodes.
    void showNodeGraph(output &s, const node_handle* node, int n) const;

    /// DEPRECATED; use dot_maker instead (see io_dot.h).
    /// Write all the nodes in the subgraph below the given nodes
    /// in a graphical format specified by the extension.
    void writeNodeGraphPicture(const char* filename, const char *extension,
        const node_handle* nodes, const char* const* labels, int n);

#endif

    /** Display the contents of a single node.
          @param  s       File stream to write to.
          @param  node    Node to display.
          @param  flags   Switches to control output;
                          see constants "SHOW_DETAILED", etc.

          @return true, iff we displayed anything
    */
    bool showNode(output &s, node_handle node, unsigned int flags = 0) const;


    /** Show various stats for this forest.
          @param  s       Output stream to write to
          @param  pad     Padding string, written at the start of
                          each output line.
          @param  flags   Which stats to display, as "flags";
                          use bitwise or to combine values.
                          For example, BASIC_STATS | FOREST_STATS.
    */
    void reportStats(output &s, const char* pad, unsigned flags) const;





  /**   Return a forest node equal to the one given.
        The node is constructed as necessary.
        This version should be used only for
        multi terminal forests.
        The unpacked node is recycled.
          @param  in    Incoming pointer index;
                        used for identity reductions.
          @param  un    Temporary node; will be recycled.

          @return       A node handle equivalent
                        to \a un, taking into account
                        the forest reduction rules
                        and if a duplicate node exists.

  */
  node_handle createReducedNode(int in, unpacked_node *un);

  /** Return a forest node for implicit node equal to the one given.
   The implicit node is already constructed inside satimpl_opname.
   This version should be used only for
   multi terminal forests.
   @param  un    Implicit relation node

   @return       A node handle equivalent
   to \a un.

   */
  node_handle createRelationNode(MEDDLY::relation_node *un);

  unsigned getImplicitTableCount() const;
  inline MEDDLY::relation_node* buildImplicitNode(node_handle rnh);
  inline MEDDLY::node_handle getImplicitTerminalNode() const;


  /** Return a forest node equal to the one given.
      The node is constructed as necessary.
      This version should be used only for
      edge valuded forests.
      The node builder nb is recycled.
        @param  in    Incoming pointer index;
                      used for identity reductions.
        @param  nb    Constructed node.
        @param  ev    Output: edge value
        @param  node  Output: node handle.
                      On exit, the edge value and the node
                      handle together are equivalent to nb;
                      taking into account the forest reduction rules
                      and if a duplicate node exists.
  */
  template <class T>
  void createReducedNode(int in, unpacked_node* nb, T& ev, node_handle& node);


    /** Swap the content of nodes.
        Do not update their parents and inCount.
    */
    void swapNodes(node_handle p, node_handle q);

    /*
     * Modify a node in place.
     * Does not check if the modified node is duplicate or redundant.
     * The level of the node may change.
     * Keep the reference number and the cache count of the node.
     */
    node_handle modifyReducedNodeInPlace(unpacked_node* un, node_handle p);

  // ------------------------------------------------------------
  // virtual in the base class, but implemented here.
  // See meddly.h for descriptions of these methods.

    virtual void writeEdges(output &s, const dd_edge* E, unsigned n) const;
    virtual void readEdges(input &s, dd_edge* E, unsigned n);
    // virtual void garbageCollect();
    // virtual void compactMemory();
    virtual void showInfo(output &strm, int verbosity);

  // ------------------------------------------------------------
  // abstract virtual, must be overridden.
  //


    virtual void normalize(unpacked_node &nb, long& ev) const;


    /** Is this a redundant node that can be eliminated?
        Must be implemented in derived forests
        to deal with the default edge value.
          @param  nb    Node we're trying to build.

          @return   True, if nr is a redundant node
                          AND it should be eliminated.
    */
    virtual bool isRedundant(const unpacked_node &nb) const = 0;

    /** Is the specified edge an identity edge, that can be eliminated?
        Must be implemented in derived forests
        to deal with the default edge value.

          @param  nr    Node we're trying to build.
                        We know there is a single non-zero downward pointer.

          @param  i     Candidate edge (or edge index for sparse nodes).

          @return True, if nr[i] is an identity edge, and the
                        identity node should be eliminated.
    */
    virtual bool isIdentityEdge(const unpacked_node &nb, int i) const = 0;


    /**
        Build an iterator.
        Used by class enumerator.
    */
    virtual enumerator::iterator* makeFullIter() const = 0;

    /**
        Build an iterator with a fixed row.
        Default behavior - throw an "INVALID_FOREST" error.
    */
    virtual enumerator::iterator* makeFixedRowIter() const;

    /**
        Build an iterator with a fixed column.
        Default behavior - throw an "INVALID_FOREST" error.
    */
    virtual enumerator::iterator* makeFixedColumnIter() const;

    /*
     * Reorganize the variables in a certain order.
     */
    void reorderVariables(const int* level2var);

    void getVariableOrder(int* level2var) const;

    std::shared_ptr<const variable_order> variableOrder() const;

    /*
     * Swap the variables at level and level+1.
     * This method should only be called by domain.
     */
    virtual void swapAdjacentVariables(int level) = 0;

    /*
     * Move the variable at level high down to level low.
     * The variables from level low to level high-1 will be moved one level up.
     */
    virtual void moveDownVariable(int high, int low) = 0;

    /*
     * Move the variable at level low up to level high.
     * The variables from level low+1 to level high will be moved one level down.
     */
    virtual void moveUpVariable(int low, int high) = 0;

    virtual void dynamicReorderVariables(int top, int bottom) {
    	throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    }

  public:
    /** Show the header information.
          @param  s       Stream to write to.
          @param  nr      Unpacked node containing header data.
    */
    virtual void showHeaderInfo(output &s, const unpacked_node &nr) const;

    /** Write the header information in machine-readable format.
          @param  s       Stream to write to.
          @param  nr      Unpacked node containing header data.
    */
    virtual void writeHeaderInfo(output &s, const unpacked_node &nr) const;

    /** Read the header information in machine-readable format.
          @param  s       Stream to read from.
          @param  nb      Node we're building.
    */
    virtual void readHeaderInfo(input &s, unpacked_node &nb) const;




  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                       protected  methods                       |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:

    /// Destructor.
    virtual ~expert_forest();

    /** Initialize data.
        Should be called in the child class constructors.
        Allows us to use class properties to initialize the data.
    */
    void initializeForest();


  // ------------------------------------------------------------
  // inlined setters for derived classes to use.

    void setUnhashedSize(unsigned char ubytes);
    void setHashedSize(unsigned char hbytes);

  // ------------------------------------------------------------
  // virtual, with default implementation.
  // Should be overridden in appropriate derived classes.

    /// Character sequence used when writing forests to files.
    virtual const char* codeChars() const;


    /** Normalize a node.
        Used only for "edge valued" DDs with range type: integer.
        Different forest types will have different normalization rules,
        so the default behavior given here (throw an error) will need
        to be overridden by all edge-valued forests.

          @param  nb    Array of downward pointers and edge values;
                        may be modified.
          @param  ev    The incoming edge value, may be modified
                        as appropriate to normalize the node.
    */
    virtual void normalize(unpacked_node &nb, int& ev) const;

    /** Normalize a node.
        Used only for "edge valued" DDs with range type: real.
        Different forest types will have different normalization rules,
        so the default behavior given here (throw an error) will need
        to be overridden by all edge-valued forests.

          @param  nb    Array of downward pointers and edge values;
                        may be modified.
          @param  ev    The incoming edge value, may be modified
                        as appropriate to normalize the node.
    */
    virtual void normalize(unpacked_node &nb, float& ev) const;

    /** Show forest-specific stats.
          @param  s     Output stream to use
          @param  pad   String to display at the beginning of each line.
    */
    virtual void reportForestStats(output &s, const char* pad) const;



  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                        private  methods                        |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:

  // ------------------------------------------------------------
  // inlined helpers for this class

    // bool isTimeToGc() const;

    /** Change the location of a node.
        Used by node_storage during compaction.
        Should not be called by anything else.
          @param  node        Node we're moving
          @param  old_addr    Current address of node, for sanity check
          @param  new_addr    Where we're moving the node to
    */
    void moveNodeOffset(node_handle node, node_address old_addr, node_address new_addr);
    friend class MEDDLY::node_storage;
    friend class MEDDLY::node_headers;  // calls deleteNode().

    friend class MEDDLY::global_rebuilder;

  // ------------------------------------------------------------
  // helpers for this class

    /**
        Disconnects all downward pointers from p,
        and removes p from the unique table.
    */
    virtual void deleteNode(node_handle p);

    /** Apply reduction rule to the temporary node and finalize it.
        Once a node is reduced, its contents cannot be modified.
          @param  in    Incoming index, used only for identity reduction;
                        Or -1.
          @param  un    Unpacked node.
          @return       Handle to a node that encodes the same thing.
    */
    node_handle createReducedHelper(int in, unpacked_node &nb);

    /** Create implicit node in the forest. Just add a handle and it points to original location
     @param  un    Relation node.
     @return       Handle to a node that encodes the same thing.
     */
    node_handle createImplicitNode(MEDDLY::relation_node &nb);

    unsigned getImplTableCount() const;
    MEDDLY::relation_node* buildImplNode(node_handle rnh);
    MEDDLY::node_handle getImplTerminalNode() const;


    /** Apply reduction rule to the temporary extensible node and finalize it.
        Once a node is reduced, its contents cannot be modified.
          @param  in    Incoming index, used only for identity reduction;
                        Or -1.
          @param  un    Unpacked extensible node. Must be sorted by indices if sparse.
          @return       Handle to a node that encodes the same thing.
    */
    node_handle createReducedExtensibleNodeHelper(int in, unpacked_node &nb);

    // Sanity check; used in development code.
    void validateDownPointers(const unpacked_node &nb) const;

    static void recycle(unpacked_node *n);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                              Data                              |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  private:
    // Garbage collection in progress
    bool performing_gc;

    // memory for validating incounts
    node_handle* in_validate;
    int  in_val_size;
    // depth of delete/zombie stack; validate when 0
    int delete_depth;

    /// Number of bytes of unhashed header
    unsigned char unhashed_bytes;
    /// Number of bytes of hashed header
    unsigned char hashed_bytes;



    class nodecounter;
    class nodemarker;
};
// end of expert_forest class.


// ******************************************************************
// *                                                                *
// *                 inlined  expert_forest methods                 *
// *                                                                *
// ******************************************************************


inline MEDDLY::memstats&
MEDDLY::expert_forest::changeMemStats()
{
  return mstats;
}

inline unsigned char
MEDDLY::expert_forest::unhashedHeaderBytes() const
{
  return unhashed_bytes;
}

inline unsigned char
MEDDLY::expert_forest::hashedHeaderBytes() const
{
  return hashed_bytes;
}


// --------------------------------------------------
// Node status
// --------------------------------------------------

// --------------------------------------------------
// Extensible Node Information:
// --------------------------------------------------
inline bool
MEDDLY::expert_forest::isExtensible(node_handle p) const
{
  return nodeMan->isExtensible(getNodeAddress(p));
}

//
// Unorganized from here
//

inline int
MEDDLY::expert_forest::getIndexSetCardinality(MEDDLY::node_handle node) const
{
  MEDDLY_DCASSERT(isIndexSet());
  if (isTerminalNode(node)) return (node != 0) ? 1 : 0;
  // yes iff the unhashed extra header is non-zero.
  const int* uhh = (const int*) nodeMan->getUnhashedHeaderOf(getNodeAddress(node));
  MEDDLY_DCASSERT(*uhh > 0);
  return *uhh;
}

inline MEDDLY::forest::node_status
MEDDLY::expert_forest::getNodeStatus(MEDDLY::node_handle node) const
{
  if (isMarkedForDeletion()) {
    return MEDDLY::forest::DEAD;
  }
  if (isTerminalNode(node)) {
    return terminalNodesStatus;
  }
  if (isDeletedNode(node)) {
    return MEDDLY::forest::DEAD;
  }
  // Active node.

  // If we're using reference counts,
  // and the incoming count is zero,
  // then we must be using optimistic
  // and the node is stale but recoverable.

  // If we're NOT using reference counts,
  // since we're not a deleted node,
  // assume we are still active.

  if (deflt.useReferenceCounts) {
    if (getNodeInCount(node) == 0) {
      return MEDDLY::forest::RECOVERABLE;
    }
  }

  return MEDDLY::forest::ACTIVE;
}


inline unsigned
MEDDLY::expert_forest::getImplicitTableCount() const
{
  unsigned q = getImplTableCount();
  return q;
}

inline MEDDLY::relation_node*
MEDDLY::expert_forest::buildImplicitNode(node_handle rnh)
{
  return buildImplNode(rnh);
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::getImplicitTerminalNode() const
{
  return getImplTerminalNode();
}


inline MEDDLY::node_handle
MEDDLY::expert_forest::createRelationNode(MEDDLY::relation_node *un)
{
  MEDDLY_DCASSERT(un);
  MEDDLY::node_handle q = createImplicitNode(*un);
#ifdef TRACK_DELETIONS
  printf("Created node %d\n", q);
#endif
  return q;
}

inline MEDDLY::node_handle
MEDDLY::expert_forest::createReducedNode(int in, MEDDLY::unpacked_node *un)
{
  MEDDLY_DCASSERT(un);
  MEDDLY::node_handle q = createReducedHelper(in, *un);
#ifdef TRACK_DELETIONS
  printf("Created node %d\n", q);
#endif
  recycle(un);
  return q;
}

template<class T>
inline void
MEDDLY::expert_forest::createReducedNode(int in, MEDDLY::unpacked_node *un, T& ev,
      MEDDLY::node_handle& node)
{
  MEDDLY_DCASSERT(un);
  // MEDDLY_DCASSERT(un->isBuildNode());
  normalize(*un, ev);
  MEDDLY_DCASSERT(ev >= 0);
  node = createReducedHelper(in, *un);
#ifdef TRACK_DELETIONS
  printf("Created node %d\n", node);
#endif
  recycle(un);
}


inline void
MEDDLY::expert_forest::setUnhashedSize(unsigned char ubytes)
{
  MEDDLY_DCASSERT(0 == unhashed_bytes);
  unhashed_bytes = ubytes;
}

inline void
MEDDLY::expert_forest::setHashedSize(unsigned char hbytes)
{
  MEDDLY_DCASSERT(0 == hashed_bytes);
  hashed_bytes = hbytes;
}

/*
inline bool
MEDDLY::expert_forest::isTimeToGc() const
{
  return isPessimistic() ? (stats.zombie_nodes > deflt.zombieTrigger)
      : (stats.orphan_nodes > deflt.orphanTrigger);
}
*/

inline void
MEDDLY::expert_forest::moveNodeOffset(MEDDLY::node_handle node, node_address old_addr,
    node_address new_addr)
{
  nodeHeaders.moveNodeAddress(node, old_addr, new_addr);
}

inline void
MEDDLY::expert_forest::getVariableOrder(int* level2var) const
{
  // Assume sufficient space has been allocated for order
  level2var[0] = 0;
  for (int i = 1; i < getNumVariables() + 1; i++) {
    level2var[i] = var_order->getVarByLevel(i);
  }
}

inline std::shared_ptr<const MEDDLY::variable_order>
MEDDLY::expert_forest::variableOrder() const
{
  return var_order;
}


#endif
