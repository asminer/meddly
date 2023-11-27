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

#ifndef MEDDLY_NODE_STORAGE
#define MEDDLY_NODE_STORAGE

#include "policies.h"

namespace MEDDLY {
    class node_storage_style;
    class node_storage;
    class node_marker;

    class expert_forest;
    class memory_manager_style;
    class output;
    class unpacked_node;
};

// ******************************************************************
// *                                                                *
// *                    node_storage_style class                    *
// *                                                                *
// ******************************************************************

/** Abstract base class for node storage factories.

    The base class is implemented in node_storage.cc;
    various backends are implemented in directory storage/.

*/
class MEDDLY::node_storage_style {
        const char* name;
    public:
        node_storage_style(const char* n);
        virtual ~node_storage_style();


        /** Build a new node storage mechanism, bound to the given forest.

            @param  f   Forest to bind to
            @param  mm  Memory manager style to use

            @return     A pointer to a node storage class,
                        initialized for forest f.
        */
        virtual node_storage* createForForest(expert_forest* f,
            const memory_manager_style* mmst) const = 0;

        inline const char* getName() const { return name; }
};


// ******************************************************************
// *                                                                *
// *                       node_storage class                       *
// *                                                                *
// ******************************************************************

/** Abstract base class for node storage.

    The base class is implemented in node_storage.cc;
    various backends are implemented in directory storage/,
    or you may implement your own scheme :^)

    Nodes are represented by an address, which for a valid
    node, will be greater than 0.  The first valid address
    must be 1.

    Whatever scheme is used to store nodes internally,
    it must be possible to set pointers \a count and \a next
    such that count[addr] and next[addr] give the incoming
    count and next pointer for the given node addr.
    Derived classes are responsible for setting up
    and maintaining these pointers.
*/
class MEDDLY::node_storage {
    public:
        node_storage(const char* sn, expert_forest* f);
        virtual ~node_storage();

        /** Go through and collect any garbage.

            @param  shrink  If true, we will shrink data structures
                            as we free space; otherwise, we won't.
        */
        virtual void collectGarbage(bool shrink) = 0;

        /** Show various stats.
            @param  s       Output stream to write to
            @param  pad     Padding string, written at the start of
                            each output line.
            @param  flags   Controls what is displayed.
        */
        virtual void reportStats(output &s, const char* pad, unsigned flags)
            const = 0;

        /** Dump the internal storage details.
            Primarily used for debugging.

            @param  s       Output stream to use
            @param  flags   What to show.
                            0x01  Show active memory
                            0x02  Show memory "holes"

            TBD - remove this
        */
        void dumpInternal(output &s, unsigned flags) const;

        /** Allocate space for, and store, a node.
            I.e., create a new node that is a copy of the given one.
            The node might be "compressed" in various ways to reduce
            storage requirements.  (Indeed, that is the whole point
            of the node_storage class.)
                @param  p     Node handle number, in case it is used
                @param  nb    Node data is copied from here.
                @param  opt   Ways we can store the node.
                @return       The "address" of the new node.
        */
        virtual node_address makeNode(node_handle p, const unpacked_node &nb,
                                  node_storage_flags opt) = 0;

        /** Destroy a node.
            Unlink the downward pointers, and recycle the memory
            used by the node.
                @param  addr    Address of the node.
        */
        virtual void unlinkDownAndRecycle(node_address addr) = 0;


        /** Schedule down pointers for exploration.
                @param  mark    Node marking mechanism.
                @param  addr    Address of the node.
        */
        virtual void addDownToQueue(node_marker &mark, node_address addr) const = 0;

    // various ways to read a node

        /** Check for duplicates.
            @param  addr    Node address in this structure
            @param  nr      Node to compare against

            @return true    iff the nodes are duplicates
        */
        virtual bool areDuplicates(node_address addr, const unpacked_node &nr)
            const = 0;

        /**
            Copy the node at the specified address, into an unpacked node.
            Useful for reading an entire node.
            @param  un      Result will be stored here.  Will be resized
                            if needed.
            @param  addr    Node address in this structure.
            @param  opt     Ways we can fill the unpacked node.
        */
        virtual void fillUnpacked(unpacked_node &un, node_address addr,
            node_storage_flags opt) const = 0;

        /** Compute the hash value for a node.
            Should give the same answer as filling a unpacked_node
            and computing the hash on the unpacked_node.

            @param  levl  Level of the node of interest
            @param  addr  Address of the node of interest
        */
        virtual unsigned hashNode(int level, node_address addr) const = 0;

        /** Determine if this is an extensible node.
            @param  addr    Node Address
            @return         True if the node stores an extensible edge,
                            False otherwise.
        */
        virtual bool isExtensible(node_address addr) const = 0;

        /** Determine if this is a singleton node.
            Used for identity reductions.
            @param  addr    Address of the node we care about
            @param  down    Output:
                            The singleton downward pointer, or undefined.

            @return     If the node has only one non-zero downward pointer,
                        then return the index for that pointer.
                        Otherwise, return a negative value.
        */
        virtual int getSingletonIndex(node_address addr, node_handle &down)
            const = 0;


        /** Get the specified downward pointer for a node.
            Fast if we just want one.
            @param  addr    Address of the node we care about
            @param  index   Index of downward pointer
            @return         Desired pointer
            @throw          INVALID_VARIABLE, if index is negative.
        */
        virtual node_handle getDownPtr(node_address addr, int index) const = 0;

        /** Get the specified outgoing edge for a node.
            Fast if we just want one.

            @param  addr    Address of the node we care about
            @param  ind     Index of the pointer we want.
            @param  ev      Output: edge value at that index.
            @param  dn      Output: downward pointer at that index.
        */
        virtual void getDownPtr(node_address addr, int ind, int& ev,
            node_handle& dn) const = 0;
        virtual void getDownPtr(node_address addr, int ind, long& ev,
            node_handle& dn) const = 0;

        /** Get the specified outgoing edge for a node.
            Fast if we just want one.

            @param  addr    Address of the node we care about
            @param  ind     Index of the pointer we want.
            @param  ev      Output: edge value at that index.
            @param  dn      Output: downward pointer at that index.
        */
        virtual void getDownPtr(node_address addr, int ind, float& ev,
            node_handle& dn) const = 0;


        /** Read the unhashed header portion of a node.

            @param  addr    Address of the node we care about
        */
        virtual const void* getUnhashedHeaderOf(node_address addr) const = 0;

        /** Read the hashed header portion of a node.

            @param  addr    Address of the node we care about
        */
        virtual const void* getHashedHeaderOf(node_address addr) const = 0;


        /**
            Get next pointer for node at this address.
            Used by unique table in case of chaining.
            @param  addr    Address of node
            @return         Handle of next node, or 0.
        */
        virtual node_handle getNextOf(node_address addr) const = 0;


        /**
            Set next pointer for node at this address.
            Used by unique table in case of chaining.
            @param  addr    Address of node
            @param  n       Non-negative.  Handle of
                            next node in a unique table chain,
                            or 0 if none.
        */
        virtual void setNextOf(node_address addr, node_handle n) = 0;


        /**
            Return the name of the style that created us.
        */
        inline const char* getStyleName() const {
            return style_name;
        }

    protected:
        /// TBD - remove this
        /// Dump information not related to individual nodes.
        virtual void dumpInternalInfo(output &s) const = 0;

        /** Get the first interesting address.
            If this cannot be determined, return 0.
            TBD - remove this
        */
        virtual node_address firstNodeAddress() const = 0;

        /** Dump the node/hole information at the given address.
            @param  s       Output stream to use
            @param  addr    Address
            @param  flags   What chunks should be displayed

            @return     Next interesting address, if we can determine this;
                        otherwise 0.

            TBD - make this public
        */
        virtual node_address dumpInternalNode(output &s, node_address addr,
            unsigned flags) const = 0;

        /// TBD - remove this
        /// Dump final info (after node info)
        virtual void dumpInternalTail(output &s) const = 0;

        // Hooks from other classes, so we don't need to make
        // all the derived classes "friends".


        //
        // Methods for derived classes to deal with
        // members owned by the base class
        //

        inline const expert_forest* getParent() const {
            MEDDLY_DCASSERT(parent);
            return parent;
        }
        inline expert_forest* getParent() {
            MEDDLY_DCASSERT(parent);
            return parent;
        }

    protected:
        /// Parent forest.
        expert_forest* parent;

    private:
        /// Name of the style that invoked us
        const char* style_name;
};

#endif
