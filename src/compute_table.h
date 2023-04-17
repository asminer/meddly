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
#include "old_meddly_expert.h"   // for operation
#include "forest.h"

namespace MEDDLY {
    class operation;
    class ct_entry_type;

    class ct_object;
    class ct_initializer;
    class compute_table_style;
    class compute_table;
};

// ******************************************************************
// *                                                                *
// *                         ct_object class                        *
// *                                                                *
// ******************************************************************

/** Generic objects in compute tables.
    Used for things other than dd_edges and simple types.
    Defined in compute_table.cc
*/
class MEDDLY::ct_object {
  public:
    ct_object();
    virtual ~ct_object();
    virtual opnd_type getType() = 0;
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

      /*
       * Moved to ct_entry_type.h, with name
       * enum class ct_typeID
       *
      enum typeID {
        ERROR = 0,
        NODE = 1,
        INTEGER = 2,
        LONG = 3,
        FLOAT = 4,
        DOUBLE = 5,
        GENERIC = 6
      };
      */

      //
      // ******************************************************************
      //

      union entry_item {
        int I;
        unsigned int U;
        long L;
        unsigned long UL;
        node_handle N;
        float F;
        double D;
        ct_object* G;
      };

      //
      // ******************************************************************
      //

      /**
        The key portion of an entry.
        Internally, in the compute table, we may store
        entries differently.  This class is used to build
        keys for searching and to construct CT entries.
      */
      class entry_key {
        public:
          entry_key();
          ~entry_key();

        protected:
          /// Start using for this operation
          void setup(const ct_entry_type* et, unsigned repeats);

        public:
          const ct_entry_type* getET() const;

          // interface, for operations.  All inlined in meddly_expert.hh
          void writeN(node_handle nh);
          void writeI(int i);
          void writeL(long i);
          void writeF(float f);
          // For templates
          inline void write_ev(long i)  { writeL(i); }
          inline void write_ev(float f) { writeF(f); }

        public:
          // interface, for compute_table.  All inlined in meddly_expert.hh
          const entry_item* rawData() const;
          unsigned dataLength() const;
          unsigned numRepeats() const;

          const void* readTempData() const;
          unsigned numTempBytes() const;
          void* allocTempData(unsigned bytes);
          /// Increase cache counters for nodes in this portion of the entry.
          void cacheNodes() const;
          unsigned getHash() const;

        protected:
          // protected interface, for compute_table.  All inlined in meddly_expert.hh
          void setHash(unsigned h);

        private:
          ct_typeID theSlotType() const;

        private:
          const ct_entry_type* etype;
          entry_item* data;
          void* temp_data;
          unsigned temp_bytes;
          unsigned temp_alloc;
          unsigned num_repeats;
          unsigned hash_value;
          unsigned data_alloc;

          unsigned currslot;
          unsigned total_slots;
#ifdef DEVELOPMENT_CODE
          bool has_hash;
#endif
        protected:
          /// Used for linked-list of recycled search keys in compute_table
          entry_key* next;

        friend class compute_table;
      };

      //
      // ******************************************************************
      //

      /**
        The result portion of an entry.
        Internally, in the compute table, we may store
        entries differently.  This class is used to return
        results from searches and to construct CT entries.
      */
      class entry_result {
        public:
          entry_result();
          ~entry_result();

        public:
          // For delayed construction
          void initialize(const ct_entry_type* et);

          // interface, for operations (reading).
          node_handle readN();
          int readI();
          float readF();
          long readL();
          double readD();
          ct_object* readG();
          // for templates
          void read_ev(long &l)   { l = readL(); }
          void read_ev(float &f)  { f = readF(); }

          // interface, for operations (building).
          void reset();
          void writeN(node_handle nh);
          void writeI(int i);
          void writeF(float f);
          void writeL(long L);
          void writeD(double D);
          void writeG(ct_object* G);

          // interface, for compute tables.
          void setValid();
          void setValid(const entry_item* d);
          void setInvalid();
          operator bool() const;
          /// Increase cache counters for nodes in this portion of the entry.
          void cacheNodes() const;

          const entry_item* rawData() const;
          unsigned dataLength() const;


        private:
          const ct_entry_type* etype;
          entry_item* build;
          const entry_item* data;
          bool is_valid;
          unsigned currslot;
      };

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
          Start using an entry_key for the given operation.
      */
      static entry_key* useEntryKey(const ct_entry_type* et, unsigned repeats);

      /**
          Done using an entry_key.
      */
      static void recycle(entry_key* k);


      /// Is this a per-operation compute table?
      bool isOperationTable() const;

      /** Find an entry in the compute table based on the key provided.
          @param  key   Key to search for.
          @param  res   Where to store the result, if any.
      */
      virtual void find(entry_key* key, entry_result &res) = 0;

      /**
          Add an entry (key plus result) to the compute table.
            @param  key   Key portion of the entry.  Will be recycled.
            @param  res   Result portion of the entry.
      */
      virtual void addEntry(entry_key* key, const entry_result &res) = 0;

      /**
          Update an existing entry in the compute table.
            @param  key   Key portion of the entry.  Will be recycled.
            @param  res   Updated result portion of the entry.
      */
      virtual void updateEntry(entry_key* key, const entry_result &res) = 0;

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
      void setHash(entry_key *k, unsigned h);

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
      static entry_key* free_keys;

    friend class operation;
};

// ******************************************************************
// *                                                                *
// *                   inlined ct_object  methods                   *
// *                                                                *
// ******************************************************************


// ******************************************************************
// *                                                                *
// *                 inlined  compute_table methods                 *
// *                                                                *
// ******************************************************************

inline void
MEDDLY::compute_table::entry_key::setup(const ct_entry_type* et, unsigned repeats)
{
  MEDDLY_DCASSERT(et);
  etype = et;
  num_repeats = repeats;
  MEDDLY_DCASSERT( 0==repeats || et->isRepeating() );
  total_slots = et->getKeySize(repeats);
  if (total_slots > data_alloc) {
    data_alloc = (1+(data_alloc / 8)) * 8;   // allocate in chunks of size 8
    data = (entry_item*) realloc(data, data_alloc*sizeof(entry_item));
    if (0==data) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(data, 0, total_slots * sizeof(entry_item));
  currslot = 0;
#ifdef DEVELOPMENT_CODE
  has_hash = false;
#endif
}

inline const MEDDLY::ct_entry_type*
MEDDLY::compute_table::entry_key::getET() const
{
  return etype;
}

inline void MEDDLY::compute_table::entry_key::writeN(node_handle nh)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::NODE == theSlotType());
  data[currslot++].N = nh;
}

inline void MEDDLY::compute_table::entry_key::writeI(int i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::INTEGER == theSlotType());
  data[currslot++].I = i;
}

inline void MEDDLY::compute_table::entry_key::writeL(long i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::LONG == theSlotType());
  data[currslot++].L = i;
}

inline void MEDDLY::compute_table::entry_key::writeF(float f)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::FLOAT == theSlotType());
  data[currslot++].F = f;
}

inline const MEDDLY::compute_table::entry_item*
MEDDLY::compute_table::entry_key::rawData() const
{
  return data;
}

inline unsigned MEDDLY::compute_table::entry_key::dataLength() const
{
  return total_slots;
}

inline unsigned MEDDLY::compute_table::entry_key::numRepeats() const
{
  return num_repeats;
}

inline const void*
MEDDLY::compute_table::entry_key::readTempData() const
{
  return temp_data;
}

inline unsigned
MEDDLY::compute_table::entry_key::numTempBytes() const
{
  return temp_bytes;
}

inline void*
MEDDLY::compute_table::entry_key::allocTempData(unsigned bytes)
{
  temp_bytes = bytes;
  if (bytes > temp_alloc) {
    temp_alloc = (1+(temp_bytes/64)) * 64;    // allocate in chunks of 64 bytes
    temp_data = realloc(temp_data, temp_alloc);
    if (0==temp_data) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  return temp_data;
}

inline void
MEDDLY::compute_table::entry_key::cacheNodes() const
{
  for (unsigned i=0; i<total_slots; i++) {
    expert_forest* f = etype->getKeyForest(i);
    if (f) {
      f->cacheNode(data[i].N);
    }
  }
}

inline unsigned MEDDLY::compute_table::entry_key::getHash() const
{
  MEDDLY_DCASSERT(has_hash);
  return hash_value;
}

inline void MEDDLY::compute_table::entry_key::setHash(unsigned h)
{
  hash_value = h;
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}

inline MEDDLY::ct_typeID MEDDLY::compute_table::entry_key::theSlotType() const
{
  //
  // Adjust currslot for OP entry, and number of repeats entry
  //
  // return etype->getKeyType(currslot - (etype->isRepeating() ? 2 : 1) );
  return etype->getKeyType(currslot);
}

// ******************************************************************

inline MEDDLY::node_handle MEDDLY::compute_table::entry_result::readN()
{
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(ct_typeID::NODE == etype->getResultType(currslot));
  return data[currslot++].N;
}

inline int MEDDLY::compute_table::entry_result::readI()
{
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(ct_typeID::INTEGER == etype->getResultType(currslot));
  return data[currslot++].I;
}

inline float MEDDLY::compute_table::entry_result::readF()
{
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(ct_typeID::FLOAT == etype->getResultType(currslot));
  return data[currslot++].F;
}

inline long MEDDLY::compute_table::entry_result::readL()
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::LONG == etype->getResultType(currslot));
  return data[currslot++].L;
}

inline double MEDDLY::compute_table::entry_result::readD()
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::DOUBLE == etype->getResultType(currslot));
  return data[currslot++].D;
}

inline MEDDLY::ct_object* MEDDLY::compute_table::entry_result::readG()
{
  MEDDLY_DCASSERT(data);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::GENERIC == etype->getResultType(currslot));
  return data[currslot++].G;
}


inline void MEDDLY::compute_table::entry_result::reset()
{
  currslot = 0;
}

inline void MEDDLY::compute_table::entry_result::writeN(node_handle nh)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::NODE == etype->getResultType(currslot));
  build[currslot++].N = nh;
}

inline void MEDDLY::compute_table::entry_result::writeI(int i)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::INTEGER == etype->getResultType(currslot));
  build[currslot++].I = i;
}

inline void MEDDLY::compute_table::entry_result::writeF(float f)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::FLOAT == etype->getResultType(currslot));
  build[currslot++].F = f;
}

inline void MEDDLY::compute_table::entry_result::writeL(long L)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::LONG == etype->getResultType(currslot));
  build[currslot++].L = L;
}

inline void MEDDLY::compute_table::entry_result::writeD(double D)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::DOUBLE == etype->getResultType(currslot));
  build[currslot++].D = D;
}

inline void MEDDLY::compute_table::entry_result::writeG(ct_object* G)
{
  MEDDLY_DCASSERT(build);
  MEDDLY_DCASSERT(currslot < dataLength());
  MEDDLY_DCASSERT(ct_typeID::GENERIC == etype->getResultType(currslot));
  build[currslot++].G = G;
}

inline void
MEDDLY::compute_table::entry_result::setValid()
{
  is_valid = true;
  data = build;
}

inline void
MEDDLY::compute_table::entry_result::setValid(const entry_item* d)
{
  is_valid = true;
  data = d;
}

inline void
MEDDLY::compute_table::entry_result::setInvalid()
{
  is_valid = false;
}

inline
MEDDLY::compute_table::entry_result::operator bool() const
{
  return is_valid;
}

inline void
MEDDLY::compute_table::entry_result::cacheNodes() const
{
  for (unsigned i=0; i<etype->getResultSize(); i++) {
    expert_forest* f = etype->getResultForest(i);
    if (f) {
      f->cacheNode(build[i].N);
    }
  }
}

inline const MEDDLY::compute_table::entry_item*
MEDDLY::compute_table::entry_result
::rawData() const
{
  return build;
}

inline unsigned MEDDLY::compute_table::entry_result
::dataLength() const
{
  return etype->getResultSize();
}


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

inline MEDDLY::compute_table::entry_key*
MEDDLY::compute_table::useEntryKey(const ct_entry_type* et, unsigned repeats)
{
  if (0==et) return 0;
  MEDDLY_DCASSERT( (0==repeats) || et->isRepeating() );

  entry_key* k;
  if (free_keys) {
    k = free_keys;
    free_keys = free_keys->next;
  } else {
    k = new entry_key();
  }
  k->setup(et, repeats);
  return k;
}

inline void
MEDDLY::compute_table::recycle(entry_key* k)
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
MEDDLY::compute_table::setHash(entry_key *k, unsigned h)
{
  MEDDLY_DCASSERT(k);
  k->setHash(h);
}

#endif
