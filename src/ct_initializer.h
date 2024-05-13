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

#ifndef MEDDLY_CT_INITIALIZER_H
#define MEDDLY_CT_INITIALIZER_H

#include "initializer.h"
#include <cstddef>

namespace MEDDLY {
    struct ct_settings;
    class ct_initializer;
    class compute_table;
    class ct_entry_type;

    class memory_manager_style;
    class compute_table_style;


    enum class staleRemovalOption {
        /// Whenever we see a stale entry, remove it.
        Aggressive,
        /// Only remove stales when we need to expand the table
        Moderate,
        /// Only remove stales during Garbage Collection.
        Lazy
    };

    enum class compressionOption {
        /// No compression at all
        None,
        /// Compression based on item type
        TypeBased
        // TBD - others
    };

};

// ******************************************************************
// *                                                                *
// *                       ct_settings  struct                      *
// *                                                                *
// ******************************************************************

struct MEDDLY::ct_settings {
    public:
        /// Memory manager to use for compute table entries
        const memory_manager_style *MMS;
        /// Maximum compute table size
        unsigned long maxSize;
        /// Stale removal policy
        staleRemovalOption staleRemoval;
        /// Compression policy
        compressionOption compression;
        /// Allow for huge (more than 4 billion entries) tables?
        bool allowHugeTables;
    public:
        ct_settings() {
            MMS = nullptr;
            allowHugeTables = false;
        }
};


// ******************************************************************
// *                                                                *
// *                      ct_initializer  class                     *
// *                                                                *
// ******************************************************************

/** Interface for initializing Meddly's compute table(s).
    Note - this is a singleton class but this is not enforced.

    This is exposed here because it allows us to avoid a
    "chicken and egg" problem:  to initialize the library, we want to
    set the compute table style, but we cannot guarantee that those
    pointers are set up before we initialize the library.
    So, settings for compute tables should be made as follows.

    (1) call defaultInitializerList(), to build an instance of this class,
        and save the result.  That will set up the default settings.

    (2) change settings using static members

    (3) initialize Meddly using the saved initializer list.

*/
class MEDDLY::ct_initializer : public initializer_list {
    public:
        enum builtinCTstyle {
            /// One huge hash table that uses chaining.
            MonolithicChainedHash,

            /// One huge hash table that does not use chaining.
            MonolithicUnchainedHash,

            /// A hash table (with chaining) for each operation.
            OperationChainedHash,

            /// A hash table (no chaining) for each operation.
            OperationUnchainedHash,
        };

    public:
        ct_initializer(initializer_list* previous);
        virtual ~ct_initializer();

    protected:
        virtual void setup();
        virtual void cleanup();
        static void setMemoryManager(const memory_manager_style*);

    // use these to change defaults, before library initialization
    public:
        static void setStaleRemoval(staleRemovalOption sro);
        static void setMaxSize(unsigned long ms);
        static void setBuiltinStyle(builtinCTstyle cts);
        static void setUserStyle(const compute_table_style*);
        static void setCompression(compressionOption co);
        static void setHugeTables(bool on);

        // for convenience
        static compute_table* createForOp(const ct_entry_type* et);

    private:
        static ct_settings the_settings;
        static const compute_table_style* ct_factory;
        static compute_table_style* builtin_ct_factory;
};


#endif // #include guard
