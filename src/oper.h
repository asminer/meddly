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

#ifndef MEDDLY_OPER_H
#define MEDDLY_OPER_H

#include "defines.h"
#include "opname.h"
#include "io.h"

namespace MEDDLY {
    class opname;
    class operation;
    class compute_table;
    class forest;

    class initializer_list;
    class ct_initializer;
    class ct_entry_type;
    class ct_entry_result;

    void cleanup();

    // ******************************************************************
    // *                      Operation management                      *
    // ******************************************************************

    /// Safely destroy the given operation.
    /// TBD: see if this can go in the operation destructor?
    void destroyOperation(operation* &op);
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
        friend class initializer_list;
        friend void destroyOperation(operation* &op);
    public:
        /** Constructor.
                @param  n           Operation "name"
                @param  et_slots    Number of different compute table entry
                                    types used by this operation.
                                    Derived class constructors must register
                                    exactly this many entry types.
        */
        operation(opname* n, unsigned et_slots);

    protected:
        virtual ~operation();

    public:

        /// Get the name of this operation; for display
        inline const char* getName() const { return theOpName->getName(); }

        /// Get the opname parent class that built us
        inline const opname* getOpName() const { return theOpName; }



        /** Are we marked for deletion?
            If so, all compute table entries for this operation
            will be removed.
        */
        inline bool isMarkedForDeletion() const {
            return is_marked_for_deletion;
        }

        void markForDeletion();

        //
        // List of operations; used as operation cache in opname classes
        //

        inline void setNext(operation* n) { next = n; }
        inline operation* getNext() const { return next; }

        //
        // Used primarily by compute tables.
        //

        /// Unique index for the operation
        inline unsigned getIndex() const { return oplist_index; }

        /// Find the operation with the given index.
        inline static operation* getOpWithIndex(unsigned i) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i, list_size);
            return op_list[i];
        }

        /// How many operations are active?
        inline static unsigned getOpListSize() { return list_size; }

        /// Set the first entry type ID;
        /// the remaining will come right after.
        inline void setFirstETid(unsigned slot) { first_etid = slot; }

        /// Get the first entry type ID
        inline unsigned getFirstETid() const { return first_etid; }

        /// Get the number of entry type IDs
        inline unsigned getNumETids() const { return num_etids; }

        //
        // CT operations
        //

        /// Remove stale compute table entries for this operation.
        void removeStaleComputeTableEntries();

        /// Remove all compute table entries for this operation.
        void removeAllComputeTableEntries();


        inline static bool usesMonolithicComputeTable() {
            return Monolithic_CT;
        }
        static void removeStalesFromMonolithic();
        static void removeAllFromMonolithic();

        static void countAllNodeEntries(const forest* f, size_t* counts);
        void countCTEntries(const forest* f, size_t* counts) const;

        //
        // Display, for debugging
        //

        static void showMonolithicComputeTable(output &, int verbLevel);
        static void showAllComputeTables(output &, int verbLevel);
        void showComputeTable(output &, int verbLevel) const;


        // TBD MOVE THIS HERE
        static void purgeAllMarked();

    protected:
        void registerInForest(forest* f);
        void unregisterInForest(forest* f);

        virtual bool checkForestCompatibility() const = 0;

        void registerEntryType(unsigned slot, ct_entry_type* et);
        void buildCTs();

        inline opname* getParent() { return theOpName; }

    private:
        // should be called during library init.
        static void initializeStatics();

        // should ONLY be called during library cleanup.
        static void destroyAllOps();


    //
    // Members
    //

    protected:
        /// Compute table to use (for entry type 0), if any.
        compute_table* CT0;
        /// Array of compute tables, one per entry type.
        compute_table** CT;
        /** Array of entry types.
            Owned by the compute_table class; we have
            these pointers for convenience.
        */
        ct_entry_type** etype;
        /** Array of entry results.
            Use these during computation.
            We only ever need one result per entry type.
        */
        ct_entry_result* CTresult;

        /**
            Number of entry_types needed by this operation.
            This gives the dimension of arrays CT and etype.
        */
        unsigned num_etids;

        /// Struct for CT searches.
        // compute_table::entry_key* CTsrch;
        // for cache of operations.
        operation* next;

    private:

    private:

        opname* theOpName;
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

        /**
            Starting slot for entry_types, assigned
            by compute_table.
        */
        unsigned first_etid;

    friend class forest;
    friend class initializer_list;

    friend class ct_initializer;

};


#endif // #include guard
