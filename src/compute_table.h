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

#ifndef MEDDLY_COMPUTE_TABLE_H
#define MEDDLY_COMPUTE_TABLE_H

#include "ct_entry_type.h"
#include "ct_entry_key.h"

#include "forest.h"
#include "oper.h"

namespace MEDDLY {
    class operation;
    class ct_entry_type;
    class ct_entry_key;
    class ct_entry_result;

    class ct_object;
    class ct_initializer;
    struct ct_settings;
    class compute_table_style;
    class compute_table;
};

// ******************************************************************
// *                                                                *
// *                    compute_table_style class                   *
// *                                                                *
// ******************************************************************

/** Interface for building compute tables.
*/
class MEDDLY::compute_table_style {
    public:
        compute_table_style();
        virtual ~compute_table_style();

        /** Build a new, monolithic table.
            Monolithic means that the table stores entries for several
            (ideally, all) operations.

            Default throws an error.
        */
        virtual compute_table* create(const ct_settings &s)
            const;


        /**
            Build a new table for a single operation.
            Default throws an error.
        */
        virtual compute_table* create(const ct_settings &s,
            operation* op, unsigned slot) const;


        /**
            Does this style build monolithic CTs?
        */
        virtual bool usesMonolithic() const = 0;
};

// ******************************************************************
// *                                                                *
// *                      compute_table  class                      *
// *                                                                *
// ******************************************************************

/** Interface for compute tables.
    Anyone implementing an operation (see below) will
    probably want to use this.
    Implementation is in compute_table.cc.
*/
class MEDDLY::compute_table {
    public:
        struct stats {
            unsigned long numEntries;
            unsigned long hits;
            unsigned long pings;
            static const unsigned searchHistogramSize = 256;
            unsigned long searchHistogram[searchHistogramSize];
            unsigned long numLargeSearches;
            unsigned maxSearchLength;
            unsigned long resizeScans;
        };

    public:

        //
        // convenience methods, for grabbing edge values
        //

        inline static void readEV(const MEDDLY::node_handle* p, int &ev) {
            ev = p[0];
        }
        inline static void readEV(const MEDDLY::node_handle* p, long &ev) {
            ev = reinterpret_cast<const long*>(p)[0];
        }
        inline static void readEV(const MEDDLY::node_handle* p, float &ev) {
            ev = reinterpret_cast<const float*>(p)[0];
        }

    public:

        /** Constructor.
                @param  s       Settings for compute table.
                @param  op      For MONOLITHIC tables, this should be 0;
                                otherwise, a pointer to the operation specific
                                to this table.
                @param  slot    Ignored for MONOLITHIC tables.  For
                                operation-specific tables, the (op, slot) pair
                                completely identifies the kinds of entries
                                in the table.
        */
        compute_table(const ct_settings &s, operation* op,
                unsigned slot);

        /** Destructor.
            Does NOT properly discard all table entries;
            use \a removeAll() for this.
        */
        virtual ~compute_table();

        /// Is this a per-operation compute table?
        inline bool isOperationTable() const { return global_et; }

        /// Get performance stats for the table.
        inline const stats& getStats() { return perf; }

// ********************************************************************
//
//   centralized key recycling; obsolete with ct_vectors.
//
// ********************************************************************

    public:
        /**
            Start using an ct_entry_key for the given operation.
        */
        inline static ct_entry_key* useEntryKey(const ct_entry_type* et,
                unsigned repeats)
        {
            if (!et) return nullptr;
            MEDDLY_DCASSERT( (0==repeats) || et->isRepeating() );

            ct_entry_key* k;
            if (free_keys) {
                k = free_keys;
                free_keys = free_keys->next;
            } else {
                k = new ct_entry_key();
            }
            k->setup(et, repeats);
            return k;
        }

        /**
            Done using an ct_entry_key.
        */
        inline static void recycle(ct_entry_key* k) {
            if (k) {
                k->next = free_keys;
                free_keys = k;
            }
        }

    private:
        static ct_entry_key* free_keys;


// ********************************************************************
//
//   operation registry, and entry type registry.
//
// ********************************************************************

    public:
        /// Find entry_type for operation and slot number.
        inline static const ct_entry_type* getEntryType(operation* op,
                unsigned slot)
        {
            MEDDLY_DCASSERT(op);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, slot, op->getNumETids());
            unsigned etid = op->getFirstETid() + slot;
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, etid, entryInfoSize);
            return entryInfo[etid];
        }

        /// Find entry type for given entryID
        inline static const ct_entry_type* getEntryType(unsigned etid)
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, etid, entryInfoSize);
            return entryInfo[etid];
        }

    protected:
        friend class ct_initializer;
        friend class operation;

        /// Initialize the entry registry
        static void initStatics();
        /// Destroy the entry registry
        static void doneStatics();

        /** Register an operation.
            Sets aside a number of entry_type slots for the operation.
        */
        static void registerOp(operation* op, unsigned num_ids);

        /// Register an entry_type.
        static void registerEntryType(unsigned etid, ct_entry_type* et);

        /** Unregister an operation.
            Frees the entry_type slots for the operation.
        */
        static void unregisterOp(operation* op, unsigned num_ids);

    private:
        static ct_entry_type** entryInfo;
        static unsigned entryInfoAlloc;
        static unsigned entryInfoSize;


// ********************************************************************
//
//   Compute table operations
//
// ********************************************************************
    public:
        //
        // Overridden in different compute table implementations
        //

        /** Find an entry in the compute table based on the key provided.
                @param  key   Key to search for.
                @param  res   Where to store the result, if any.
        */
        virtual void find(ct_entry_key* key, ct_entry_result &res) = 0;

        /**
            Add an entry (key plus result) to the compute table.
                @param  key   Key portion of the entry.  Will be recycled.
                @param  res   Result portion of the entry.
        */
        virtual void addEntry(ct_entry_key* key, const ct_entry_result &res)
            = 0;

        /**
            Update an existing entry in the compute table.
                @param  key   Key portion of the entry.  Will be recycled.
                @param  res   Updated result portion of the entry.
        */
        virtual void updateEntry(ct_entry_key* key, const ct_entry_result &res)
            = 0;

        /** Remove all stale entries.
            Scans the table for entries that are no longer valid (i.e. they are
            stale, according to operation::isEntryStale) and removes them. This
            can be a time-consuming process (proportional to the number of
            cached entries).
        */
        virtual void removeStales() = 0;

        /// Removes all entries.
        virtual void removeAll() = 0;

        /// For debugging.
        virtual void show(output &s, int verbLevel = 0) = 0;

        /** Also for debugging.
            Examine all entries, and for each pointer to forest f node p,
            increment counts[p].
        */
        virtual void countNodeEntries(const forest* f, size_t* counts)
            const = 0;


    protected:

        inline static void setHash(ct_entry_key *k, unsigned h) {
            MEDDLY_DCASSERT(k);
            k->setHash(h);
        }

        /// Clear CT Bits in forests that could have entries in this table.
        void clearForestCTBits(bool* skipF, unsigned n) const;

        /// Start sweep phase for forests that could have entries in this table.
        void sweepForestCTBits(bool* whichF, unsigned n) const;

    protected:
        /// The maximum size of the hash table.
        unsigned maxSize;
        /// Do we try to eliminate stales during a "find" operation
        bool checkStalesOnFind;
        /// Do we try to eliminate stales during a "resize" operation
        bool checkStalesOnResize;
        /// Global entry type, if we're an operation cache; otherwise null.
        const ct_entry_type* global_et;
        /// Performance statistics
        stats perf;

};

#endif
