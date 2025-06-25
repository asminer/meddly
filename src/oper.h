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

namespace MEDDLY {
    class operation;
    class compute_table;
    class forest;

    class initializer_list;
    class ct_initializer;
    class ct_entry_type;
    class ct_entry_result;

    /// Argument and result types for apply operations.
    enum class opnd_type {
        FOREST      = 0,
        BOOLEAN     = 1,
        INTEGER     = 2,
        REAL        = 3,
        HUGEINT     = 4
    };
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
        /** Constructor (OLD)
                @param  n           Operation name, for debugging
                @param  et_slots    Number of different compute table entry
                                    types used by this operation.
                                    Derived class constructors must register
                                    exactly this many entry types.
        */
        operation(const char* n, unsigned et_slots=0);

        /// New constructor.
        operation();


        /** Destructor.
            Safe to call directly now :)
        */
        virtual ~operation();

#ifdef ALLOW_DEPRECATED_0_17_6
        /// Safely destroy the given operation.
        inline static void destroy(operation* op) {
            delete op;
        }
#endif

        /// Safely destroy all operations associated with the given forest.
        /// Called in forest destructor.
        static void destroyAllWithForest(const forest* f);

        /// Get the name of this operation; for display / debugging
        inline const char* getName() const { return name; }

    protected:
        /// Set the name after constructing.
        inline void setName(const char* n) {
            if (!n) return;
            MEDDLY_DCASSERT(!name);
            name = n;
        }

        void registerInForest(forest* f);
        void unregisterInForest(forest* f);

#ifdef ALLOW_DEPRECATED_0_17_6
        void registerEntryType(unsigned slot, ct_entry_type* et);
        void buildCTs();
#endif


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
        /// Called during library initialization.
        static void initializeStatics();

        /// Called during library cleanup.
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
#ifdef ALLOW_DEPRECATED_0_17_6
        /// Compute table to use (for entry type 0), if any.
        compute_table* CT0;
        /// Array of compute tables, one per entry type.
        compute_table** CT;

        //
        // TBD: remove all of these
        //

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
#endif

    private:
        const char* name;
        /// List of forest IDs associated with this operation.
        std::vector <unsigned> FList;
};

#endif // #include guard
