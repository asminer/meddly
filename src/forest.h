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
#include "node_headers.h"
#include "node_storage.h"
#include "unpacked_node.h"
#include "policies.h"
#include "rangeval.h"
#include "domain.h"
#include "enumerator.h"

#include "statset.h"

#include <vector>

namespace MEDDLY {
    class forest;

    class MDD_levels;
    class MXD_levels;

    class domain;

    class node_marker;
    class node_storage;
    class unpacked_node;
    class unique_table;
    class impl_unique_table;
    class relation_node;
    class global_rebuilder;

    class logger;
};

// ******************************************************************
// *                                                                *
// *                        MDD_levels class                        *
// *                                                                *
// ******************************************************************

/**
    Small class for MDD level behavior,
    used for template parameters in operations.
 */
class MEDDLY::MDD_levels {
    public:

        /// Go "down a level"
        static inline int downLevel(int k) {
            return k-1;
        }

        ///  Go "up a level"
        static inline int upLevel(int k) {
            return k+1;
        }

        /// Determine the top of two levels
        static inline int topLevel(int k1, int k2) {
            return MAX(k1, k2);
        }
};


// ******************************************************************
// *                                                                *
// *                        MXD_levels class                        *
// *                                                                *
// ******************************************************************

/**
    Small class for MXD level behavior,
    used for template parameters in operations.
 */
class MEDDLY::MXD_levels {
    public:

        /// Go "down a level"
        static inline int downLevel(int k) {
            return (k>0) ? (-k) : (-k-1);
        }

        ///  Go "up a level"
        static inline int upLevel(int k) {
            return (k<0) ? (-k) : (-k-1);
        }

        /// Determine the top of two levels
        static inline int topLevel(int k1, int k2) {
            if (ABS(k1) == ABS(k2)) {
                return MAX(k1, k2);
            }
            return (ABS(k1) > ABS(k2)) ? k1 : k2;
        }

        /// Determine the top unprimed level of two levels.
        /// This is ABS(topLevel(k1, k2)) but computed more efficiently.
        static inline int topUnprimed(int k1, int k2) {
            return MAX(ABS(k1), ABS(k2));
        }
};


// ******************************************************************
// *                                                                *
// *                       EdgeOp_none  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    /**
        Small class for MT behavior on edges.
    */
    class EdgeOp_none {
        public:
            /// Do we have edge values
            static inline bool hasEdgeValues()
            {
                return false;
            }
            /// Edge value type; for CT entries
            static inline char edgeValueTypeLetter()
            {
                return 0;
            }
            /// Clear an edge value.
            static inline void clear(edge_value &v)
            {
                v.set();
            }
            /// Return a + b
            static inline edge_value applyOp(const edge_value &a,
                    const edge_value &b)
            {
                return edge_value();
            }
            /// Accumulate a += b
            static inline void accumulateOp(edge_value &a,
                    const edge_value &b)
            {
            }
            /// Accumulate a += -b
            static inline void accumulateInverse(edge_value &a,
                    const edge_value &b)
            {
            }
            static inline terminal buildTerm(
                            const edge_value &val, node_handle p)
            {
                MEDDLY_DCASSERT(false);
                // Keep compiler happy:
                return terminal(p, terminal_type::OMEGA);
            }
            static inline bool isZeroFunction(const edge_value &val,
                    node_handle p)
            {
                return (0==p);
            }

    };
};


// ******************************************************************
// *                                                                *
// *                       EdgeOp_plus  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    /**
        Small template class for EV+ behavior on edges.
    */
    template <class TYPE>
    class EdgeOp_plus {
        public:
            /// Do we have edge values
            static inline bool hasEdgeValues()
            {
                return true;
            }
            /// Edge value type; for CT entries
            static char edgeValueTypeLetter();
            /// Clear an edge value.
            static inline void clear(edge_value &v)
            {
                v.set(TYPE(0));
            }
            /// Accumulate two edge values. Returns a + b
            static inline edge_value applyOp(const edge_value &a,
                    const edge_value &b)
            {
                TYPE av, bv;
                a.get(av);
                b.get(bv);
                return edge_value(av+bv);
            }
            /// Accumulate one edge value into another: a += b
            static inline void accumulateOp(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.add(bv);
            }
            /// Accumulate an inverse: a += -b
            static inline void accumulateInverse(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.subtract(bv);
            }
            /// Build terminal value for edge <val, p>
            static inline terminal buildTerm(
                            const edge_value &val, node_handle p)
            {
                if (OMEGA_NORMAL == p) {
                    TYPE v;
                    val.get(v);
                    return terminal(v);
                }
                return terminal(p, terminal_type::OMEGA);
            }
            static inline bool isZeroFunction(const edge_value &val,
                    node_handle p)
            {
                if (p != OMEGA_NORMAL) return false;
                TYPE bv;
                val.get(bv);
                return (0==bv);
            }
    };
    template <>
    inline char EdgeOp_plus<int>::edgeValueTypeLetter() {
        return 'I';
    }
    template <>
    inline char EdgeOp_plus<long>::edgeValueTypeLetter() {
        return 'L';
    }
    template <>
    inline char EdgeOp_plus<float>::edgeValueTypeLetter() {
        return 'F';
    }
    template <>
    inline char EdgeOp_plus<double>::edgeValueTypeLetter() {
        return 'D';
    }
};

// ******************************************************************
// *                                                                *
// *                       EdgeOp_times class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    /**
      Small template class for EV+ behavior on edges.
    */
    template <class TYPE>
    class EdgeOp_times {
        public:
            /// Do we have edge values
            static inline bool hasEdgeValues()
            {
                return true;
            }
            /// Edge value type; for CT entries
            static char edgeValueTypeLetter();
            /// Clear an edge value.
            static inline void clear(edge_value &v)
            {
                v.set(TYPE(1));
            }
            /// Accumulate two edge values. Returns a * b
            static inline edge_value applyOp(const edge_value &a,
                    const edge_value &b)
            {
                TYPE av, bv;
                a.get(av);
                b.get(bv);
                return edge_value(av*bv);
            }
            /// Accumulate one edge value into another: a *= b
            static inline void accumulateOp(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.multiply(bv);
            }
            /// Accumulate an inverse: a *= 1/b
            static inline void accumulateInverse(edge_value &a,
                    const edge_value &b)
            {
                TYPE bv;
                b.get(bv);
                a.divide(bv);
            }
            /// Build terminal value for edge <val, p>
            static inline terminal buildTerm(
                            const edge_value &val, node_handle p)
            {
                if (OMEGA_NORMAL == p) {
                    TYPE v;
                    val.get(v);
                    return terminal(v);
                }
                return terminal(p, terminal_type::OMEGA);
            }
            static inline bool isZeroFunction(const edge_value &val,
                    node_handle p)
            {
                return (OMEGA_ZERO == p);
            }
    };
    template <>
    inline char EdgeOp_times<float>::edgeValueTypeLetter() {
        return 'F';
    }
    template <>
    inline char EdgeOp_times<double>::edgeValueTypeLetter() {
        return 'D';
    }
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

    TBD: fix and discuss the by-level reduction rules
*/
class MEDDLY::forest {

// ===================================================================
//
// Static methods to construct, destroy forests
//
// ===================================================================

    // ------------------------------------------------------------
    public:
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

        /// Front-end function to destroy a forest.
        static void destroy(forest* &f);


// ===================================================================
//
// Methods to manage nodes in the forest
//
// ===================================================================


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
            unpacked_node* un = unpacked_node::New(this);
            unpackNode(un, p, f);
            return un;
        }

    // ------------------------------------------------------------
    public: // Add /remove a node to the forest
    // ------------------------------------------------------------

        /** Return a forest node equal to the one given.
            The node is constructed as necessary.
            The unpacked node is recycled.
                @param  un  Unpacked node; will be recycled.

                @param  ev  Output: edge value from normalizing
                            the node. Will be void for multi-terminal.

                @param  nh  Output: A node handle equivalent to
                            \a un (when using edge value ev),
                            taking into account the forest reduction rules
                            and if a duplicate node exists.

                @param  in  Incoming pointer index;
                            used for identity reductions.
                            A value of -1 may be used to not
                            attempt identity reductions.
        */
        void createReducedNode(unpacked_node *un, edge_value &ev,
                node_handle &node, int in=-1);

#ifdef ALLOW_DEPRECATED_0_17_6
        /** Return a forest node equal to the one given.
            The node is constructed as necessary.
            This version should be used only for
            multi terminal forests.
            The unpacked node is recycled.
                @param  un  Unpacked node; will be recycled.

                @param  in  Incoming pointer index;
                            used for identity reductions.
                            A value of -1 may be used to not
                            attempt identity reductions.

                @return     A node handle equivalent
                            to \a un, taking into account
                            the forest reduction rules
                            and if a duplicate node exists.

        */
        inline node_handle createReducedNode(int in, unpacked_node *un) {
            MEDDLY_DCASSERT(un);
            edge_value ev;
            node_handle node;
            createReducedNode(un, ev, node, in);
            MEDDLY_DCASSERT(ev.isVoid());
            return node;
        }


        /** Return a forest node equal to the one given.
            The node is constructed as necessary.
            This version should be used only for
            edge valuded forests.
            The node builder nb is recycled.
                @param  in      Incoming pointer index;
                                used for identity reductions.

                @param  un      Constructed node.

                @param  ev      Output: edge value

                @param  node    Output: node handle.
                                On exit, the edge value and the node
                                handle together are equivalent to nb;
                                taking into account the forest reduction rules
                                and if a duplicate node exists.

            TBD: rewrite this to use edge_value instead
        */
        template <class T>
        inline void createReducedNode(int in, unpacked_node* un,
                T& ev, node_handle& node)
        {
            MEDDLY_DCASSERT(un);
            edge_value _ev;
            createReducedNode(un, _ev, node, in);
            _ev.get(ev);
        }
#endif

        /** Return a forest node for implicit node equal to the one given.
            The implicit node is already constructed inside satimpl_opname.
            This version should be used only for multi terminal forests.
                @param  un    Implicit relation node

                @return       A node handle equivalent to \a un.

        */
        inline node_handle createRelationNode(MEDDLY::relation_node *un) {
            MEDDLY_DCASSERT(un);
            node_handle q = createImplicitNode(*un);
#ifdef TRACK_DELETIONS
            std::cerr << "Created relation node " << q << "\n";
#endif
            return q;
        }

        relation_node* buildImplicitNode(node_handle rnh);

        /** Remove node p.
            Disconnects all downward pointers from p,
            and removes p from the unique table.
            Unless you are implementing a garbage collector,
            you shouldn't call this yourself.
        */
        void deleteNode(node_handle p);


        /**
            Useful helper function.
            Build a chain of redundant nodes from node p up to
            a node at level L.
            More precisely, build a chain of nodes indicating
            that the function does not depend on levels from
            the one above K to level L.
            If the forest is fully-reduced, do nothing instead.
                @param  p   Bottom node
                @param  K   Start redundant nodes above this level.
                            Usually this is the same level as p.
                @param  L   Stop redundant nodes here
                @return     A node that can be referenced at level L
                            (i.e., definitely at level L if we're quasi
                            reduced) with redundant nodes added down
                            to level K.

            Illustration:

            Level:      input:              output:

            L               |                  [.]  (new redundant node)
                            |                   |
                            |                   |
            ...             |                  [.]  (new redundant node)
                            |                   |
                            |                   |
            K+1             |                  [.]  (new redundant node)
                            |                   |
                            |                   |
            K               |                   |
                            |                   |
                            v                   v
            p.lvl           p                   p

         */
        inline node_handle makeRedundantsTo(node_handle p, int K, int L)
        {
            if (0==L) return p;
            if (K==L) return p;
            MEDDLY_DCASSERT( ABS(L) >= ABS(K) );
            if (isFullyReduced()) return p;
            if (0==p) return p;

            return _makeRedundantsTo(p, K, L);
        }

        /**
            Useful helper function.
            Build a chain of redundant/identity nodes from node p
            up to a node at level L.
            If the forest is identity-reduced, do nothing instead.

                @param  p   Bottom node
                @param  K   Start adding nodes above this level
                @param  L   Stop adding nodes here
                @param  in  Index where root edge originates.
                            Only needed if L is a primed level.

                @return     A node that can be referenced at level L
                            (i.e., definitely at level L if we're quasi
                            reduced) with redundant/identity nodes added down
                            to level K.

            Illustration:

            Level:      input:              output:

            L               |                  [.]  (new redundant node)
                            |                   |
            -L              |                  [.]  (new identity nodes)
                            |                   |
            ...             |                   |
                            |                   |
            K+1             |                  [.]  (new redundant node)
                            |                   |
            -(K+1)          |                  [.]  (new identity nodes)
                            |                   |
            K               |                   |
                            |                   |
            ...             |                   |
                            v                   v
                            p                   p
        */
        inline node_handle makeIdentitiesTo(node_handle p, int K, int L, int in)
        {
            if (0==L) return p;
            if (K==L) return p;
            MEDDLY_DCASSERT(MXD_levels::topLevel(L, K) == L);
            if (0==p) return p;
            if (isIdentityReduced()) {
                //
                // A bit extra in case p is a singleton node
                //
                if (K<0) p = _makeRedundantsTo(p, K, -K);
                return p;
            }
            if (isFullyReduced() && (K<0) && (-K == L)) {
                // Just adding one redundant node;
                // that's a no-op for fully reduced
                return p;
            }
            return _makeIdentitiesTo(p, K, L, in);
        }

        /**
            Redirect illegal i-singleton edges if we're identity reduced.

                @param  i   Which child edge are we.

                @param  p   Target node.

                @return p   If we are allowed to point to p.
                            This is always the case for forests that
                            are not identity reduced and
                            for nodes p that are not at a primed level.
        */
        inline node_handle redirectSingleton(unsigned i, node_handle p)
        {
            if (!isIdentityReduced()) return p;
            if (getNodeLevel(p) >= 0) return p;

            // check if p is a singleton node
            unsigned pi;
            node_handle pd;
            if (isSingletonNode(p, pi, pd)) {
                // p is a singleton node.
                //  child pi points to pd.
                if (i != pi) {
                    return p;
                }
                // redirect
                // and fix link counts
#ifdef REFCOUNTS_ON
                if (deflt.useReferenceCounts) {
                    nodeHeaders.linkNode(pd);
                    nodeHeaders.unlinkNode(p);
                }
#endif
                return pd;
            } else {
                // p is not a singleton node.
                return p;
            }
        }


    // ------------------------------------------------------------
    protected:  // helpers for reducing nodes
    // ------------------------------------------------------------

        /** Apply reduction rule to the temporary extensible node and finalize it.
            Once a node is reduced, its contents cannot be modified.
                @param  in      Incoming index, used only for identity reduction;
                                Or -1.
                @param  un      Unpacked extensible node. Must be sorted by
                                indices if sparse.
                @return         Handle to a node that encodes the same thing.
        */
#ifdef ALLOW_EXTENSIBLE
        node_handle createReducedExtensibleNodeHelper(int in, unpacked_node &nb);
#endif


        /** Create implicit node in the forest.
            Just add a handle and it points to original location
                @param  un    Relation node.
                @return       Handle to a node that encodes the same thing.
        */
        node_handle createImplicitNode(MEDDLY::relation_node &nb);


        /**
            Heavy implementation for makeRedundantsTo().
                @param  p   Bottom node
                @param  K   Level of node p
                @param  L   Level of node we want
         */
        node_handle _makeRedundantsTo(node_handle p, int K, int L);

        /**
            Heavy implementation for makeIdentitiesTo().
                @param  p   Bottom node
                @param  K   Level of node p
                @param  L   Level of node we want
                @param  in  Index of edge above level L
         */
        node_handle _makeIdentitiesTo(node_handle p, int K, int L, int in);

    // ------------------------------------------------------------
    protected: // Moving nodes around; called in derived classes for reordering
    // ------------------------------------------------------------

        /**
            Swap the content of nodes.
            Do not update their parents and inCount.
        */
        void swapNodes(node_handle p, node_handle q);

        /**
            Modify a node in place.
            Does not check if the modified node is duplicate or redundant.
            The level of the node may change.
            Keep the reference number and the cache count of the node.
        */
        node_handle modifyReducedNodeInPlace(unpacked_node* un, node_handle p);


// ===================================================================
//
// Methods for getting/setting node information
//
// ===================================================================

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

        /// Is a CT entry containing this node, dead?
        /// I.e., cannot be used at all?
        inline bool isDeadEntry(node_handle p) const {
            if (isMarkedForDeletion()) return true;
            if (isTerminalNode(p)) return false;
            return isDeletedNode(p);
        }

        /// Is a CT entry containing this node, stale?
        /// Still usable if matched, otherwise should be
        /// removed from the CT?
        inline bool isStaleEntry(node_handle p) const {
            if (isMarkedForDeletion()) return true;
            if (isTerminalNode(p)) return false;
            if (isDeletedNode(p)) return true;

            if (deflt.useReferenceCounts) {
                if (getNodeInCount(p) == 0) {
                    return true;
                }
            }
            return false;
        }

        /// Is the forest about to be deleted?
        inline bool isMarkedForDeletion() const {
            return is_marked_for_deletion;
        }

        inline node_handle getLastNode() const {
            return nodeHeaders.lastUsedHandle();
        }

        /// Returns the in-count for a node.
        inline unsigned long getNodeInCount(node_handle p) const {
#ifdef REFCOUNTS_ON
            MEDDLY_DCASSERT(deflt.useReferenceCounts);
            return nodeHeaders.getIncomingCount(p);
#else
            return 0;
#endif
        }


        /** Increase the link count to this node.
            Call this when another node is made to point to this node.
                @return p, for convenience.
        */
        inline node_handle linkNode(node_handle p) {
#ifdef REFCOUNTS_ON
            if (deflt.useReferenceCounts) {
                return nodeHeaders.linkNode(p);
            } else {
                return p;
            }
#else
            return p;
#endif
        }

        /** Increase the link count to this node.
            Call this when another node is made to point to this node.
                @return the node handle, for convenience.
        */
        inline node_handle linkNode(const dd_edge &p) {
#ifdef REFCOUNTS_ON
            MEDDLY_DCASSERT(p.isAttachedTo(this));
            if (deflt.useReferenceCounts) {
                return nodeHeaders.linkNode(p.getNode());
            } else {
                return p.getNode();
            }
#else
            return p.getNode();
#endif
        }

        /** Decrease the link count to this node.
            If link count reduces to 0, this node may get marked for deletion.
            Call this when another node releases its connection to this node.
        */
        inline void unlinkNode(node_handle p) {
#ifdef REFCOUNTS_ON
            if (deflt.useReferenceCounts) {
                nodeHeaders.unlinkNode(p);
            }
#endif
        }


        /** Increase the cache count for this node.
            Call this whenever this node is added to a cache.
            Do nothing if we are not using reference counts for caches.
                @param  p     Node we care about.
        */
        inline void cacheNode(node_handle p) {
#ifdef REFCOUNTS_ON
            if (deflt.useReferenceCounts) {
                nodeHeaders.cacheNode(p);
            }
#endif
        }

        /** Decrease the cache count for this node.
            Call this whenever this node is added to a cache.
            Do nothing if we are not using reference counts for caches.
                @param  p     Node we care about.
        */
        inline void uncacheNode(node_handle p) {
#ifdef REFCOUNTS_ON
            if (deflt.useReferenceCounts) {
                nodeHeaders.uncacheNode(p);
            }
#endif
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

        /** Change the address of a node.
            Used by node_storage during compaction.
            Should not be called by anything else.
                @param  node        Node we're moving
                @param  old_addr    Current address of node, for sanity check
                @param  new_addr    Where we're moving the node to
        */
        inline void moveNodeAddress(node_handle node,
                node_address old_addr, node_address new_addr)
        {
            nodeHeaders.moveNodeAddress(node, old_addr, new_addr);
        }


        /// Unlink down pointers in an  unpacked node.
        inline void unlinkAllDown(const unpacked_node &un, unsigned i=0) {
#ifdef REFCOUNTS_ON
            if (deflt.useReferenceCounts) {
                for ( ; i<un.getSize(); i++) {
                    nodeHeaders.unlinkNode(un.down(i));
                }
            }
#endif
        }

        /// Link down pointers in an  unpacked node.
        inline void linkAllDown(unpacked_node &un, unsigned i=0) {
#ifdef REFCOUNTS_ON
            if (deflt.useReferenceCounts) {
                for ( ; i<un.getSize(); i++) {
                    nodeHeaders.linkNode(un.down(i));
                }
            }
#endif
        }



    // ------------------------------------------------------------
    private:  // private node header information
    // ------------------------------------------------------------
        /// Node header information
        node_headers nodeHeaders;

        /// Used for mark & sweep
        node_marker* reachable;

        /// Is the forest marked for deletion
        bool is_marked_for_deletion;

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


        /** Determine if this is a singleton node.
            Used for identity reductions.
            @param  p       Handle of node to look at

            @param  index   On output:
                            if we're a singleton node, the index of the
                            only non-zero pointer;
                            otherwise undefined.

            @param  down    On output:
                            if we're a singleton node, the only non-zero
                            downward pointer;
                            otherwise undefined.

            @return     If the node is a singleton node
                        (only one non-zero downward pointer),
                        return true; otherwise return false.
        */
        inline bool isSingletonNode(node_handle p, unsigned &index,
                node_handle &down) const
        {
            if (p<=0) return false;
            return nodeMan->isSingletonNode(getNodeAddress(p), index, down);
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
            For reading all or several downward pointers, a
            unpacked_node should be used instead.

            @param  p       Node to look at
            @param  index   Index of the pointer we want.

            @param  ev      Output: edge value at that index.
            @param  dn      Output: downward pointer at that index.
        */
        inline void getDownPtr(node_handle p, int index, edge_value& ev,
                node_handle& dn) const
        {
            ev = transparent_edge;
            nodeMan->getDownPtr(getNodeAddress(p), index, ev, dn);
        }

        // TBD: remove these older 'getDownPtr' methods

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
            edge_value v = transparent_edge;
            nodeMan->getDownPtr(getNodeAddress(p), index, v, dn);
            ev = v.getInt();
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
            edge_value v = transparent_edge;
            nodeMan->getDownPtr(getNodeAddress(p), index, v, dn);
            ev = v.getLong();
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
            edge_value v = transparent_edge;
            nodeMan->getDownPtr(getNodeAddress(p), index, v, dn);
            ev = v.getFloat();
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

        /// Is this an extensible node
        inline bool isExtensible(node_handle p) const {
            return nodeMan->isExtensible(getNodeAddress(p));
        }

        /// Get the cardinality of an Index Set.
        inline int getIndexSetCardinality(node_handle node) const {
            if (!isIndexSet()) {
                throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
            }
            if (isTerminalNode(node)) return (node != 0) ? 1 : 0;
            // yes iff the unhashed extra header is non-zero.
            const int* uhh = (const int*) nodeMan->getUnhashedHeaderOf(
                    getNodeAddress(node)
            );
            MEDDLY_DCASSERT(*uhh > 0);
            return *uhh;
        }

    // ------------------------------------------------------------
    private: // node storage members
    // ------------------------------------------------------------

        /// Class that stores nodes.
        node_storage* nodeMan;

// ===================================================================
//
// Methods for terminals, edge values in the forest
//
// ===================================================================

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
            Get the no-op edge value.
        */
        inline const edge_value& getNoOpEdge() const {
            return noop_edge;
        }

        /**
            Is the given edge transparent?
            If so it may be "skipped" in a sparse node.
                @param  ev    Value part of the edge to check
                @param  ep    Node part of the edge to check
        */
        inline bool isTransparentEdge(const edge_value &ev, node_handle ep)
                    const
        {
            if (ep != transparent_node) return false;
            return ev == transparent_edge;
        }

        /**
            Get the transparent edge value.
            This is the default edge value for "skipped" edges in sparse nodes.
            Copy the transparent edge value into the parameters.
                @param  ev      Put value part of the edge here.
                @param  ep      Put node part of the edge here.
        */
        inline void getTransparentEdge(edge_value &ptr, node_handle &ep)
                    const
        {
            ptr = transparent_edge;
            ep = transparent_node;
        }

        /**
            Make a transparent edge
        */
        inline void getTransparentEdge(dd_edge &e) const {
            e.set(transparent_edge, transparent_node);
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
            MEDDLY_DCASSERT(isMultiTerminal());
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
            MEDDLY_DCASSERT(isMultiTerminal());
            MEDDLY_DCASSERT(n <= 0);
            terminal t(the_terminal_type, n);
            t.getValue(v);
        }

        inline bool getBooleanFromHandle(MEDDLY::node_handle n) const {
            MEDDLY_DCASSERT(isMultiTerminal());
            bool v;
            getValueFromHandle(n, v);
            return v;
        }

        inline int getIntegerFromHandle(MEDDLY::node_handle n) const {
            MEDDLY_DCASSERT(isMultiTerminal());
            int v;
            getValueFromHandle(n, v);
            return v;
        }

        inline float getRealFromHandle(MEDDLY::node_handle n) const {
            MEDDLY_DCASSERT(isMultiTerminal());
            float v;
            getValueFromHandle(n, v);
            return v;
        }


        /**
            Build an edge at level 0 (i.e., always to a terminal node)
            for a constant value.
         */
        void getEdgeForValue(rangeval T, edge_value &v, node_handle &p)
            const;

        /**
            Get the value of a level 0 edge (to a terminal node).
        */
        void getValueForEdge(const edge_value &v, node_handle p, rangeval &T)
            const;

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
            MEDDLY_DCASSERT(isTerminalNode(p));
            transparent_node = p;
            transparent_edge.set();
            noop_edge.set();
        }
        inline void setTransparentEdge(node_handle p, int v) {
            MEDDLY_DCASSERT(edge_type::INT == the_edge_type);
            MEDDLY_DCASSERT(isTerminalNode(p));
            transparent_node = p;
            transparent_edge.set(v);
            setNoOpEdgeval();
        }
        inline void setTransparentEdge(node_handle p, long v) {
            MEDDLY_DCASSERT(edge_type::LONG == the_edge_type);
            MEDDLY_DCASSERT(isTerminalNode(p));
            transparent_node = p;
            transparent_edge.set(v);
            setNoOpEdgeval();
        }
        inline void setTransparentEdge(node_handle p, float v) {
            MEDDLY_DCASSERT(edge_type::FLOAT == the_edge_type);
            MEDDLY_DCASSERT(isTerminalNode(p));
            transparent_node = p;
            transparent_edge.set(v);
            setNoOpEdgeval();
        }
        inline void setTransparentEdge(node_handle p, double v) {
            MEDDLY_DCASSERT(edge_type::DOUBLE == the_edge_type);
            MEDDLY_DCASSERT(isTerminalNode(p));
            transparent_node = p;
            transparent_edge.set(v);
            setNoOpEdgeval();
        }

    // ------------------------------------------------------------
    private:
    // ------------------------------------------------------------

        inline void setNoOpEdgeval() {
            MEDDLY_DCASSERT(!isMultiTerminal());
            if (isEVTimes()) {
                noop_edge.setTempl(transparent_edge.getType(), 1);
            } else {
                noop_edge.setTempl(transparent_edge.getType(), 0);
            }
        }

    // ------------------------------------------------------------
    private: // transparent edge info
    // ------------------------------------------------------------

        /// Transparent node value
        node_handle transparent_node;

        /// Transparent edge value
        edge_value  transparent_edge;

        /// Edge value for no-op (e.g., 0 for +, 1 for *)
        edge_value  noop_edge;

        /// Edge type.
        edge_type the_edge_type;

        /// Terminal node type.
        terminal_type the_terminal_type;

        /// Are edge values hashed?
        bool hash_edge_values;


// ===================================================================
//
// Methods for extra "header" information in a node (typically none)
//
// ===================================================================

    // ------------------------------------------------------------
    public: // methods to get header size
    // ------------------------------------------------------------

        /// Extra bytes per node, not hashed.
        inline unsigned unhashedHeaderBytes() const {
            return unhashed_bytes;
        }

        /// Extra bytes per node, hashed.
        inline unsigned hashedHeaderBytes() const {
            return hashed_bytes;
        }


    // ------------------------------------------------------------
    protected: // methods to set header size, for derived classes
    // ------------------------------------------------------------

        /// Set the unhashed extra bytes per node.
        /// Call as needed in derived classes;
        /// otherwise we assume 0 bytes.
        inline void setUnhashedSize(unsigned ubytes) {
            MEDDLY_DCASSERT(0 == unhashed_bytes);
            unhashed_bytes = ubytes;
        }

        /// Set the hashed extra bytes per node.
        /// Call as needed in derived classes;
        /// otherwise we assume 0 bytes.
        inline void setHashedSize(unsigned hbytes) {
            MEDDLY_DCASSERT(0 == hashed_bytes);
            hashed_bytes = hbytes;
        }


    // ------------------------------------------------------------
    private: // members for header size
    // ------------------------------------------------------------

        /// Number of bytes of unhashed header
        unsigned unhashed_bytes;

        /// Number of bytes of hashed header
        unsigned hashed_bytes;


    // ------------------------------------------------------------
    protected: // must be called in derived class constructors
    // ------------------------------------------------------------

        /** Initialize node storage.
            Should be called in the child class constructors,
            after calling setHashedSize() and setUnhashedSize().
            Allows us to use class properties to initialize the data.
        */
        void initializeStorage();


// ===================================================================
//
// Methods for the forest's variables, and ordering
//
// ===================================================================

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

        inline unsigned getNumVariables() const {
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
            return isForRelations() ? -int(getNumVariables()) : 0;
        }

        /// returns K
        inline int getMaxLevelIndex() const {
            return int(getNumVariables());
        }

        /// Check if the given level is valid
        inline bool isValidLevel(int k) const {
            return (k >= getMinLevelIndex()) && (k <= getMaxLevelIndex());
        }


        /** Go "down a level" in a relation.
            Safest to use this, in case in later versions
            the level numbering changes, or becomes forest dependent.
                @param  k   Current level
                @return Level immediately below level k
        */
        // static inline int downLevel(int k) {
            // return (k>0) ? (-k) : (-k-1);
        // }

        /** Go "up a level" in a relation.
            Safest to use this, in case in later versions
            the level numbering changes, or becomes forest dependent.
                @param  k   Current level
                @return Level immediately above level k
        */
        // static inline int upLevel(int k) {
            // return (k<0) ? (-k) : (-k-1);
        // }


        /// Can we have extensible nodes at level k?
        inline bool isExtensibleLevel(int k) const {
            MEDDLY_DCASSERT(isValidLevel(k));
            return d->getVar(unsigned(k < 0? -k: k))->isExtensible();
        }

        /// The maximum size (number of indices) a node at this level can have
        inline int getLevelSize(int k) const {
            MEDDLY_DCASSERT(isValidLevel(k));
            int var=getVarByLevel(k);
            if (var < 0) {
                return getDomain()->getVariableBound(unsigned(-var), true);
            } else {
                return getDomain()->getVariableBound(unsigned(var), false);
            }
        }

        /// The maximum size (number of indices) a variable can have.
        inline int getVariableSize(int var) const {
            return getDomain()->getVariableBound(unsigned(var), false);
        }


    // ------------------------------------------------------------
    private: // Domain info
    // ------------------------------------------------------------
        domain* d;
        friend class domain;

    // ------------------------------------------------------------
    public: // public methods for variable re/ordering
    // ------------------------------------------------------------

        /*
         * Reorganize the variables in a certain order.
         */
        void reorderVariables(const int* level2var);

        inline void getVariableOrder(int* level2var) const {
            // Assume sufficient space has been allocated for order
            level2var[0] = 0;
            for (unsigned i = 1; i < getNumVariables() + 1; i++) {
                level2var[i] = var_order->getVarByLevel((int) i);
            }
        }

        inline std::shared_ptr<const variable_order> variableOrder() const {
            return var_order;
        }


        /*
         * Swap the variables at level and level+1.
         * This method should only be called by domain.
         */
        virtual void swapAdjacentVariables(int level) = 0;

        /*
         * Move the variable at level high down to level low.
         * The variables from level low to level high-1 will be moved
         * one level up.
         */
        virtual void moveDownVariable(int high, int low) = 0;

        /*
         * Move the variable at level low up to level high.
         * The variables from level low+1 to level high will be moved
         * one level down.
         */
        virtual void moveUpVariable(int low, int high) = 0;

        virtual void dynamicReorderVariables(int top, int bottom) {
    	    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }

        // needed by reordering_base.h
        inline const unique_table* getUT() {
            return unique;
        }

    protected:
        std::shared_ptr<const variable_order> var_order;

        /// uniqueness table, still used by derived classes.
        unique_table* unique;

        /// uniqueness table for relation nodes.
        impl_unique_table* implUT;


// ===================================================================
//
// Methods for forest settings
//
// ===================================================================


    // ------------------------------------------------------------
    public: // Getters related to the mdd type and other policies
    // ------------------------------------------------------------

        /// Does this forest represent relations or matrices?
        inline bool isForRelations() const {
            return isRelation;
        }

        /// Returns the range type.
        inline range_type getRangeType() const {
            return rangeType;
        }

        /// Query the range type.
        inline bool isRangeType(range_type r) const {
            return r == rangeType;
        }

        /// Returns the edge labeling mechanism.
        inline edge_labeling getEdgeLabeling() const {
            return edgeLabel;
        }

        /// Is the edge labeling "MULTI_TERMINAL".
        inline bool isMultiTerminal() const {
            return edge_labeling::MULTI_TERMINAL == edgeLabel;
        }

        /// Is the edge labeling "EV_PLUS".
        inline bool isEVPlus() const {
            return edge_labeling::EVPLUS == edgeLabel;
        }

        /// Is the edge labeling "INDEX_SET".
        inline bool isIndexSet() const {
            return edge_labeling::INDEX_SET == edgeLabel;
        }

        /// Is the edge labeling "EV_TIMES".
        inline bool isEVTimes() const {
            return edge_labeling::EVTIMES == edgeLabel;
        }

        /// Check if we match a specific type of forest
        inline bool matches(bool isR, range_type rt, edge_labeling el) const {
            return (isRelation == isR) && (rangeType == rt) && (edgeLabel == el);
        }


        /// Returns the current policies used by this forest.
        inline const policies& getPolicies() const {
            return deflt;
        }

        /// Returns the current policies used by this forest.
        inline policies& getPolicies() {
            return deflt;
        }

        /// Returns the reduction rule used by this forest.
        inline reduction_rule getReductionRule() const {
            return deflt.reduction;
        }

        /// Returns true if the forest is fully reduced.
        inline bool isFullyReduced() const {
            return deflt.isFullyReduced();
        }

        /// Returns true if the forest is quasi reduced.
        inline bool isQuasiReduced() const {
            return deflt.isQuasiReduced();
        }

        /// Returns true if the forest is identity reduced.
        inline bool isIdentityReduced() const {
            return deflt.isIdentityReduced();
        }

        /// Returns true if the forest is user_defined reduced.
        inline bool isUserDefinedReduced() const {
            return deflt.isUserDefinedReduced();
        }

        /// Returns true if the level is fully reduced.
        inline bool isFullyReduced(int k) const {
            if(k<0) return level_reduction_rule[2*(-k)   ] == -1;
            else    return level_reduction_rule[2*(k) - 1] == -1;
        }

        /// Returns true if the level is quasi reduced.
        inline bool isQuasiReduced(int k) const {
            if(k<0) return level_reduction_rule[2*(-k)   ] == -2;
            else    return level_reduction_rule[2*(k) - 1] == -2;
        }

        /// Returns true if the level is identity reduced.
        inline bool isIdentityReduced(int k) const {
            if(k<0) return level_reduction_rule[2*(-k)   ] == -3;
            else    return level_reduction_rule[2*(k) - 1] == -3;
        }

        inline const int* getLevelReductionRule() const {
            return level_reduction_rule;
        }


    // ------------------------------------------------------------
    private: // private members for mdd type / reductions
    // ------------------------------------------------------------

        policies deflt;
        int *level_reduction_rule;

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


// ===================================================================
//
// Methods for I/O
//
// ===================================================================


    // ------------------------------------------------------------
    public: // Public I/O methods
    // ------------------------------------------------------------

        /**
            Show an edge, compactly.
            Called for example for each child when displaying an entire node.
                @param  s       Stream to write to.
                @param  ev      Edge value
                @param  d       Down pointer
        */
        virtual void showEdge(output &s, const edge_value &ev, node_handle d)
            const = 0;

        /** Show the header information in human-readable format.
            Default behavior does nothing; override in derived
            forests if there is any extra header information.
                @param  s       Stream to write to.
                @param  nr      Unpacked node containing header data.
        */
        virtual void showHeaderInfo(output &s, const unpacked_node &nr) const;

        /** Write the header information in machine-readable format.
            Default behavior does nothing; override in derived
            forests if there is any extra header information.
                @param  s       Stream to write to.
                @param  nr      Unpacked node containing header data.
        */
        virtual void writeHeaderInfo(output &s, const unpacked_node &nr) const;

        /** Read the header information in machine-readable format.
            Default behavior does nothing; override in derived
            forests if there is any extra header information.
                @param  s       Stream to read from.
                @param  nb      Node we're building.
        */
        virtual void readHeaderInfo(input &s, unpacked_node &nb) const;


// ===================================================================
//
// Managing root edges (using a registry)
//
// ===================================================================


    // ------------------------------------------------------------
    public: // Public methods for root edge registry
    // ------------------------------------------------------------

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

    // ------------------------------------------------------------
    private: // Private methods for root edge registry
    // ------------------------------------------------------------
        void unregisterDDEdges();

    // ------------------------------------------------------------
    private: // Private members for root edge registry
    // ------------------------------------------------------------
        /// Registry of dd_edges (doubly-linked list)
        dd_edge* roots;


// ===================================================================
//
// Managing operations on this forest
//
// ===================================================================

  public:
    /// Remove all compute table entries associated with this forest.
    void removeAllComputeTableEntries();


// ===================================================================
//
// Managing all forests in a registry
//
// ===================================================================

    // ------------------------------------------------------------
    public: // public methods for the forest registry
    // ------------------------------------------------------------

        /// Returns the forest identifier, a unique positive integer per forest.
        /// FID 0 can safely be used to indicate "no forest".
        inline unsigned FID() const { return fid; }

        /// Returns the largest forest identifier ever seen.
        static inline unsigned MaxFID() {
            return all_forests.size()-1;
        }

        /// Find the forest with the given FID.
        /// If the forest has been deleted, this will be null.
        static inline forest* getForestWithID(unsigned id) {
            if (id >= all_forests.size()) return nullptr;
            return all_forests[id];
        }


    // ------------------------------------------------------------
    private: // private methods for the forest registry
    // ------------------------------------------------------------
        static void initStatics();
        static void freeStatics();
        static void registerForest(forest* f);
        static void unregisterForest(forest* f);

        friend class forest_initializer;

    // ------------------------------------------------------------
    private: // private members for the forest registry
    // ------------------------------------------------------------
        /// "Registry" of all forests, ever.
        static std::vector <forest*> all_forests;

        /// our ID
        unsigned fid;


        friend class initializer_list;

// ===================================================================
//
// Misc. debugging or logging
//
// ===================================================================

    // ------------------------------------------------------------
    public: // public methods for debugging or logging
    // ------------------------------------------------------------

        /** Display the contents of a single node.
                @param  s       File stream to write to.
                @param  node    Node to display.
                @param  flags   Switches to control output;
                                see constants "SHOW_DETAILS", etc.
                @return true, iff we displayed anything
        */
        bool showNode(output &s, node_handle node, display_flags flags = 0) const;


        /** Show various stats for this forest.
                @param  s       Output stream to write to
                @param  pad     Padding string, written at the start of
                                each output line.
                @param  flags   Which stats to display, as "flags";
                                use bitwise or to combine values.
                                For example, BASIC_STATS | FOREST_STATS.
        */
        void reportStats(output &s, const char* pad, unsigned flags) const;


        /**
            Display all nodes in the forest.
                @param  s       File stream to write to
                @param  flags   Switches to control output;
                                see constants "SHOW_DETAILS", etc.
        */
        void dump(output &s, display_flags flags) const;
        void dumpInternal(output &s) const;
        void dumpUniqueTable(output &s) const;
        void validateIncounts(bool exact, const char* fn, unsigned ln,
                const char* opname=nullptr) const;
        void validateCacheCounts() const;
        void countNodesByLevel(long* active) const;

        // Sanity check; used in development code.
        void validateDownPointers(const unpacked_node &nb) const;


        /** Display all active (i.e., connected) nodes in the forest.
            This is primarily for aid in debugging.
            @param  strm      Stream to write to.
            @param  verbosity How much information to display.
                                0 : just statistics.
                                1 : all forest nodes + statistics.
                                2 : internal forest + statistics.
        */
        virtual void showInfo(output &strm, int verbosity=0);


        /** Start logging stats.
                @param  L       Logger to use; if 0, we don't log anything.
                                Will overwrite the old logger.
                                The logger WILL NOT be deleted by the forest
                                (in case we want to log for multiple forests).

                @param  name    Name to use for the forest (for display only).
        */
        void setLogger(logger* L, const char* name);

    // ------------------------------------------------------------
    protected: // protected methods for debugging
    // ------------------------------------------------------------

        /** Show forest-specific stats.
            Default does nothing; override in derived classes
            if there are other stats to display.
                @param  s     Output stream to use
                @param  pad   String to display at the beginning of each line.
        */
        virtual void reportForestStats(output &s, const char* pad) const;

    // ------------------------------------------------------------
    private: // private members for debugging
    // ------------------------------------------------------------

        logger *theLogger;
        int delete_depth;

// ===================================================================
//
// Various forest stats
//
// ===================================================================

    // ------------------------------------------------------------
    public: // public methods for performance stats
    // ------------------------------------------------------------

        /// Get forest performance stats.
        inline const statset& getStats() const {
            return stats;
        }

        /** Get the current number of nodes in the forest, at all levels.
            @return     The current number of nodes, not counting deleted or
                        marked for deletion nodes.
        */
        inline long getCurrentNumNodes() const {
            return stats.active_nodes;
        }

        /** Get the peak number of nodes in the forest, at all levels.
            This will be at least as large as calling getNumNodes() after
            every operation and maintaining the maximum.
            @return     The peak number of nodes that existed at one time,
                        in the forest.
        */
        inline long getPeakNumNodes() const {
            return stats.peak_active;
        }

        /** Set the peak number of nodes to the number current number of nodes.
        */
        inline void resetPeakNumNodes() {
            stats.peak_active = stats.active_nodes;
        }

    // ------------------------------------------------------------
    public: // public methods for memory stats
    // ------------------------------------------------------------

        /// Get forest memory stats.
        inline const memstats& getMemoryStats() const {
            return mstats;
        }

        /** Get the current total memory used by the forest.
            This should be equal to summing getMemoryUsedForVariable()
            over all variables.
            @return     Current memory used by the forest.
        */
        inline size_t getCurrentMemoryUsed() const {
            return mstats.getMemUsed();
        }

        /** Get the current total memory allocated by the forest.
            This should be equal to summing getMemoryAllocatedForVariable()
            over all variables.
            @return     Current total memory allocated by the forest.
        */
        inline size_t getCurrentMemoryAllocated() const {
            return mstats.getMemAlloc();
        }

        /** Get the peak memory used by the forest.
            @return     Peak total memory used by the forest.
        */
        inline size_t getPeakMemoryUsed() const {
            return mstats.getPeakMemUsed();
        }

        /** Get the peak memory allocated by the forest.
            @return     Peak memory allocated by the forest.
        */
        inline size_t getPeakMemoryAllocated() const {
            return mstats.getPeakMemAlloc();
        }

    protected:
        statset stats;
        memstats mstats;


// ===================================================================
//
// To be cleaned up still, below here.
//
// ===================================================================


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


// ===================================================================
//
// Building methods
//
// ===================================================================

  public:

    /** Create an edge for a constant function.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.
    */
    inline void createConstant(rangeval val, dd_edge &e)
    {
        if (!e.isAttachedTo(this)) {
            throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
        }
        edge_value v;
        node_handle p;
        getEdgeForValue(val, v, p);
        p = makeRedundantsTo(p, 0, getNumVariables());
        e.set(v, p);
    }


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
    inline void createEdgeForVar(int vh, bool pr, dd_edge& a) {
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




    //
    // misc.
    //

  // ------------------------------------------------------------
  // Ugly details from here down.
  private:  // Defaults
    static policies mddDefaults;
    static policies mxdDefaults;

  private:
    bool isRelation;
    range_type rangeType;
    edge_labeling edgeLabel;

    /// Mark for deletion.
    void markForDeletion();


// ===================================================================
//
// Deprecated as of version 0.17.7
//
// ===================================================================

#ifdef ALLOW_DEPRECATED_0_17_7
    public:

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
                        the range type of the forest is not INTEGER.
    */
    inline  void createEdge(int val, dd_edge &e) {
        createEdge(long(val), e);
    }

    /** Create an edge for an integer constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not INTEGER.
    */
    virtual void createEdge(long val, dd_edge &e);

    /** Create an edge for a real constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not REAL.
    */
    virtual void createEdge(float val, dd_edge &e);

    /** Create an edge for a real constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.

        @throws       TYPE_MISMATCH, if
                        the range type of the forest is not REAL.
    */
    virtual void createEdge(double val, dd_edge &e);



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

#endif

};


#endif
