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


#include "old_meddly_expert.h"   // for operation
#include "forest.h"

namespace MEDDLY {
    class operation;
    class ct_entry_type;
    class ct_entry_key;
    class ct_entry_result;

    class ct_object;
    class ct_initializer;
    class compute_table_style;
    class compute_table;
};

// ******************************************************************
// *                                                                *
// *                      ct_initializer  class                     *
// *                                                                *
// ******************************************************************

/** Interface for initializing Meddly's compute table(s).
    Implemented in compute_table.cc.
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
    enum staleRemovalOption {
      /// Whenever we see a stale entry, remove it.
      Aggressive,
      /// Only remove stales when we need to expand the table
      Moderate,
      /// Only remove stales during Garbage Collection.
      Lazy
    };

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

    enum compressionOption {
      /// No compression at all
      None,
      /// Compression based on item type
      TypeBased
      // TBD - others
    };

    struct settings {
      public:
        /// Memory manager to use for compute table entries
        const memory_manager_style *MMS;
        /// Maximum compute table size
        size_t maxSize;
        /// Stale removal policy
        staleRemovalOption staleRemoval;
        /// Compression policy
        compressionOption compression;
      public:
        settings() {
          MMS = 0;
        }
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
    static void setMaxSize(unsigned ms);
    static void setBuiltinStyle(builtinCTstyle cts);
    static void setUserStyle(const compute_table_style*);
    static void setCompression(compressionOption co);

    // for convenience
    static compute_table* createForOp(operation* op, unsigned slot);

  private:
    static settings the_settings;
    static const compute_table_style* ct_factory;
    static compute_table_style* builtin_ct_factory;
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
    virtual compute_table* create(const ct_initializer::settings &s)
      const;


    /**
        Build a new table for a single operation.
        Default throws an error.
    */
    virtual compute_table* create(const ct_initializer::settings &s,
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

      //
      // ******************************************************************
      //

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

      //
      // ******************************************************************
      //


      //
      // ******************************************************************
      //

      // convenience methods, for grabbing edge values
      static void readEV(const node_handle* p, int &ev);
      static void readEV(const node_handle* p, long &ev);
      static void readEV(const node_handle* p, float &ev);

      /** Constructor.
            @param  s   Settings for compute table.
            @param  op  For MONOLITHIC tables, this should be 0;
                        otherwise, a pointer to the operation specific to
                        this table.
            @param  slot  Ignored for MONOLITHIC tables.  For operation-specific
                          tables, the (op, slot) pair completely identifies
                          the kinds of entries in the table.
      */
      compute_table(const ct_initializer::settings &s, operation* op, unsigned slot);

      /** Destructor.
          Does NOT properly discard all table entries;
          use \a removeAll() for this.
      */
      virtual ~compute_table();

      /**
          Start using an ct_entry_key for the given operation.
      */
      static ct_entry_key* useEntryKey(const ct_entry_type* et, unsigned repeats);

      /**
          Done using an ct_entry_key.
      */
      static void recycle(ct_entry_key* k);


      /// Is this a per-operation compute table?
      bool isOperationTable() const;

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
      virtual void addEntry(ct_entry_key* key, const ct_entry_result &res) = 0;

      /**
          Update an existing entry in the compute table.
            @param  key   Key portion of the entry.  Will be recycled.
            @param  res   Updated result portion of the entry.
      */
      virtual void updateEntry(ct_entry_key* key, const ct_entry_result &res) = 0;

      /** Remove all stale entries.
          Scans the table for entries that are no longer valid (i.e. they are
          stale, according to operation::isEntryStale) and removes them. This
          can be a time-consuming process (proportional to the number of cached
          entries).
      */
      virtual void removeStales() = 0;

      /** Removes all entries.
      */
      virtual void removeAll() = 0;

      /// Get performance stats for the table.
      const stats& getStats();

      /// For debugging.
      virtual void show(output &s, int verbLevel = 0) = 0;

      /** Also for debugging.
          Examine all entries, and for each pointer to forest f node p,
          increment counts[p].
      */
      virtual void countNodeEntries(const expert_forest* f, size_t* counts) const = 0;


      static void initialize();
      static void destroy();

    protected:
      /// Clear CT Bits in forests that could have entries in this table.
      void clearForestCTBits(bool* skipF, unsigned n) const;

      /// Start sweep phase for forests that could have entries in this table.
      void sweepForestCTBits(bool* whichF, unsigned n) const;

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

    public:
      /// Find entry_type for operation and slot number.
      static const ct_entry_type* getEntryType(operation* op, unsigned slot);

      /// Find entry type for given entryID
      static const ct_entry_type* getEntryType(unsigned etID);

    protected:
      void setHash(ct_entry_key *k, unsigned h);

    protected:
      /// The maximum size of the hash table.
      unsigned maxSize;
      /// Do we try to eliminate stales during a "find" operation
      bool checkStalesOnFind;
      /// Do we try to eliminate stales during a "resize" operation
      bool checkStalesOnResize;
      /// Global entry type, if we're an operation cache; otherwise 0.
      const ct_entry_type* global_et;
      /// Performance statistics
      stats perf;

    private:
      static ct_entry_type** entryInfo;
      static unsigned entryInfoAlloc;
      static unsigned entryInfoSize;

    private:
      static ct_entry_key* free_keys;

    friend class operation;
};

// ******************************************************************
// *                                                                *
// *                 inlined  compute_table methods                 *
// *                                                                *
// ******************************************************************

// ******************************************************************

// ******************************************************************

inline bool
MEDDLY::compute_table::isOperationTable() const
{
  return global_et;
}

// convenience methods, for grabbing edge values
inline void
MEDDLY::compute_table::readEV(const MEDDLY::node_handle* p, int &ev)
{
  ev = p[0];
}
inline void
MEDDLY::compute_table::readEV(const MEDDLY::node_handle* p, long &ev)
{
  long* l = (long*) p;
  ev = l[0];
}
inline void
MEDDLY::compute_table::readEV(const MEDDLY::node_handle* p, float &ev)
{
  float* f = (float*) p;
  ev = f[0];
}

inline MEDDLY::ct_entry_key*
MEDDLY::compute_table::useEntryKey(const ct_entry_type* et, unsigned repeats)
{
  if (0==et) return 0;
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

inline void
MEDDLY::compute_table::recycle(ct_entry_key* k)
{
  if (k) {
    k->next = free_keys;
    free_keys = k;
  }
}

inline const MEDDLY::compute_table::stats&
MEDDLY::compute_table::getStats()
{
  return perf;
}

inline const MEDDLY::ct_entry_type*
MEDDLY::compute_table::getEntryType(operation* op, unsigned slot)
{
  MEDDLY_DCASSERT(op);
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, slot, op->getNumETids());
  unsigned etid = op->getFirstETid() + slot;
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, etid, entryInfoSize);
  return entryInfo[etid];
}

inline const MEDDLY::ct_entry_type*
MEDDLY::compute_table::getEntryType(unsigned etid)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, etid, entryInfoSize);
  return entryInfo[etid];
}

inline void
MEDDLY::compute_table::setHash(ct_entry_key *k, unsigned h)
{
  MEDDLY_DCASSERT(k);
  k->setHash(h);
}

#endif
