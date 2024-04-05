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
#include "io.h"

#include <vector>

#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif

namespace MEDDLY {
    class operation;
    class compute_table;
    class forest;

    class initializer_list;
    class ct_initializer;
    class ct_entry_type;
    class ct_entry_result;
    class ct_object;

    void cleanup();


    /// Argument and result types for apply operations.
    enum class opnd_type {
        FOREST      = 0,
        BOOLEAN     = 1,
        INTEGER     = 2,
        REAL        = 3,
        HUGEINT     = 4,
        FLOATVECT   = 5,
        DOUBLEVECT  = 6
    };

    // ******************************************************************
    // *                    Wrapper for GMP integers                    *
    // ******************************************************************

#ifdef HAVE_LIBGMP
    ct_object& get_mpz_wrapper();
    void unwrap(const ct_object &, mpz_t &value);
#endif

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
    public:
        /** Constructor.
                @param  n           Operation name, for debugging
                @param  et_slots    Number of different compute table entry
                                    types used by this operation.
                                    Derived class constructors must register
                                    exactly this many entry types.
        */
        operation(const char* n, unsigned et_slots);

        /// Safely destroy the given operation.
        static void destroy(operation* op);

    protected:
        virtual ~operation();

    public:
        /// Get the name of this operation; for display
        inline const char* getName() const { return name; }

        /** Are we marked for deletion?
            If so, all compute table entries for this operation
            will be removed.
        */
        inline bool isMarkedForDeletion() const {
            return is_marked_for_deletion;
        }

        void markForDeletion();

        //
        // Used primarily by compute tables.
        //

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


    // ------------------------------------------------------------
    public:  // public  methods for the operation registry
    // ------------------------------------------------------------
        /// Unique ID > 0 for the operation.
        inline unsigned getID() const { return oplist_index; }

        /// Find the operation with the given ID.
        inline static operation* getOpWithID(unsigned i) {
#ifdef DEVELOPMENT_CODE
            return op_list.at(i);
#else
            return op_list[i];
#endif
        }

        /// How many operations are active?
        inline static unsigned getOpListSize() {
            return op_list.size();
        }


    // ------------------------------------------------------------
    private: // private methods for the operation registry
    // ------------------------------------------------------------
        // should be called during library init.
        static void initializeStatics();

        // should ONLY be called during library cleanup.
        static void destroyAllOps();

        static void registerOperation(operation &o);
        static void unregisterOperation(operation &o);

    // ------------------------------------------------------------
    private: // private members for the operation registry
    // ------------------------------------------------------------
        unsigned oplist_index;
        static std::vector <operation*> op_list;
        static std::vector <unsigned>   free_list;

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

    private:
        const char* name;
        bool is_marked_for_deletion;

        // declared and initialized in meddly.cc
        static compute_table* Monolithic_CT;

        /**
            Starting slot for entry_types, assigned
            by compute_table.
        */
        unsigned first_etid;

    friend class forest;
    friend class initializer_list;

    friend class ct_initializer;

};

#ifdef ALLOW_DEPRECATED_0_17_5
namespace MEDDLY {
    inline void destroyOperation(operation* op)
    {
        operation::destroy(op);
    }
};
#endif

#endif // #include guard
