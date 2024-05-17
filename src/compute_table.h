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

#include "defines.h"
#include "io.h"

namespace MEDDLY {
    class operation;
    class forest;

#ifdef ALLOW_DEPRECATED_0_17_6
    class ct_entry_key;
    class ct_entry_result;
#endif

    class ct_entry_type;
    class ct_vector;
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
        bool uses_mono;
    public:
        compute_table_style(bool um);
        virtual ~compute_table_style();

        /** Build a new, monolithic table.
            Monolithic means that the table stores entries for several
            (ideally, all) operations.

            Default throws an error.
        */
        virtual compute_table* create(const ct_settings &s) const;

        /**
            Build a new table for a single operation.
            Default throws an error.
        */
        virtual compute_table* create(const ct_settings &s,
                unsigned etid) const;


        /**
            Does this style build monolithic CTs?
        */
        inline bool usesMonolithic() const {
            return uses_mono;
        }
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
                @param  etid    Which entry type we are tied to.
                                For MONOLITHIC tables, this should be 0;
                                otherwise, it is the unique ID for the
                                ct_entry_type.
        */
        compute_table(const ct_settings &s, unsigned etid);

        /** Destructor.
            Does NOT properly discard all table entries;
            use \a removeAll() for this.
        */
        virtual ~compute_table();

        /// Is this a per-operation compute table?
        inline bool isOperationTable() const { return global_etid; }

        /// Get performance stats for the table.
        inline const stats& getStats() { return perf; }

// ********************************************************************
//
//   centralized key recycling; obsolete with ct_vectors.
//
// ********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

    public:
        /**
            Start using an ct_entry_key for the given operation.
        */
        static ct_entry_key* useEntryKey(const ct_entry_type* et,
                unsigned repeats);

        static void recycle(ct_entry_key* k);

    private:
        static ct_entry_key* free_keys;

#endif


// ********************************************************************
//
//   Compute table operations
//
// ********************************************************************

    public:
        //
        // Overridden in different compute table implementations
        //

        // **********************************************************
        //  OLD IFACE
        // **********************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

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
        virtual void updateEntry(ct_entry_key* key, const ct_entry_result &res);


        /**
            Indicate that we are done with a key.
            Might be needed if we do a find but then do NOT use addEntry.
         */
        virtual void doneKey(ct_entry_key* key) { };
#endif

        // **********************************************************
        //  NEW IFACE
        // **********************************************************

        /** Find an entry in the compute table based on the key provided.
                @param  ET    Entry type we're looking for.
                @param  key   Key to search for.
                @param  res   Where to store the result, if found.
                @return true, if a result was found; false otherwise.
        */
        virtual bool find(const ct_entry_type& ET, ct_vector &key,
                ct_vector &res) = 0;

        /**
            Add an entry (key plus result) to the compute table.
                @param  ET    Entry type we're adding.
                @param  key   Key portion of the entry.
                @param  res   Result portion of the entry.
        */
        virtual void addEntry(const ct_entry_type& ET, ct_vector &key,
                const ct_vector &res) = 0;


        // **********************************************************

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
        virtual void countNodeEntries(const forest* f,
                std::vector <unsigned long> &counts) const = 0;



        /// Clear CT Bits in forests that could have entries in this table.
        void clearForestCTBits(std::vector <bool> &skipF) const;

        /// Start sweep phase for forests that could have entries in this table.
        void sweepForestCTBits(std::vector <bool> &whichF) const;

    protected:
        /// The maximum size of the hash table.
        unsigned long maxSize;
        /// Do we try to eliminate stales during a "find" operation
        bool checkStalesOnFind;
        /// Do we try to eliminate stales during a "resize" operation
        bool checkStalesOnResize;
        /// Global entry type id, if we're an operation cache; otherwise 0.
        const unsigned global_etid;
        /// Performance statistics
        stats perf;

    // ===================================================================
    // Monolithic compute table members/methods
    // ===================================================================

    private:
        /// Monolithic compute table
        static compute_table* Monolithic_CT;

    public:
        /// Initialize the monolithic compute table.
        /// Normally done automatically in ct_initializer::setup().
        static void initStatics(compute_table* mct);

        /// Destroy the monolithic compute table.
        /// Normally done automatically in ct_initializer::cleanup().
        static void doneStatics();

        static inline compute_table* Monolithic() {
            return Monolithic_CT;
        }

        inline static bool removeStalesFromMonolithic() {
            if (Monolithic_CT) {
                Monolithic_CT->removeStales();
                return true;
            }
            return false;
        }
        inline static bool removeAllFromMonolithic() {
            if (Monolithic_CT) {
                Monolithic_CT->removeAll();
                return true;
            }
            return false;
        }
        inline static bool showMonolithicComputeTable(output &s, int verbLevel)
        {
            if (Monolithic_CT) {
                Monolithic_CT->show(s, verbLevel);
                return true;
            }
            return false;
        }
        inline static bool countMonolithicNodeEntries(const forest* f,
                std::vector<unsigned long> &counts)
        {
            if (Monolithic_CT) {
                Monolithic_CT->countNodeEntries(f, counts);
                return true;
            }
            return false;
        }
};

#endif
