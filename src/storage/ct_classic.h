
// $Id$

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

#ifndef CT_CLASSIC_H
#define CT_CLASSIC_H

#include "../defines.h"

#include <map>  // for operation_map
#include <limits.h>

namespace MEDDLY {
  /// base class for all compute tables here;
  /// handles allocation of entries :^)
  class base_table;

  /// base class for hash tables (gives hash function)
  class base_hash;

  class base_chained;
  class base_unchained;

  class monolithic_chained;
  class operation_chained;

  class monolithic_unchained;
  class operation_unchained;

  class base_map;
}


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                          base_table class                          *
// *                                                                    *
// *                                                                    *
// **********************************************************************

class MEDDLY::base_table : public compute_table {
  protected:
    class old_search_key : public compute_table::search_key {
        friend class MEDDLY::base_table;
        unsigned hash_value;
        node_handle* data;
        bool killData;
        node_handle* key_data;
        int currslot;
#ifdef DEVELOPMENT_CODE
        /// used only for range checking during "development".
        int keyLength;  
        bool has_hash;
#endif
      public:
        old_search_key(operation* op);
        virtual ~old_search_key();

        virtual void reset() {
          currslot = 0;
#ifdef DEVELOPMENT_CODE
          has_hash = false;
#endif
        }

        virtual void writeNH(node_handle nh) {
          MEDDLY_CHECK_RANGE(0, currslot, keyLength);
          key_data[currslot++] = nh;
        }
        virtual void write(int i) {
          MEDDLY_CHECK_RANGE(0, currslot, keyLength);
          key_data[currslot++] = i;
        }
        virtual void write(float f) {
          MEDDLY_CHECK_RANGE(0, currslot, keyLength);
          float* x = (float*) (key_data+currslot);
          x[0] = f;
          currslot++;
        }

        inline node_handle* rawData() const { return data; }
        // inline int dataLength() const { return hashLength; }
        inline int dataLength() const { return currslot + (key_data - data);}

        inline void setHash(unsigned h) {
          hash_value = h;
#ifdef DEVELOPMENT_CODE
          has_hash = true;
#endif
        }
        inline unsigned getHash() const {
          MEDDLY_DCASSERT(has_hash);
          return hash_value;
        }
    };

    class old_search_result : public compute_table::search_result {
        const node_handle* data;
        unsigned currslot;
#ifdef DEVELOPMENT_CODE
        unsigned ansLength;
#endif
      public:
        old_search_result() {
          data = 0;
        }
        virtual ~old_search_result() {
        }

        virtual node_handle readNH() {
          MEDDLY_CHECK_RANGE(0, currslot, ansLength);
          return data[currslot++];
        }
        virtual void read(int &i) {
          MEDDLY_CHECK_RANGE(0, currslot, ansLength);
          i = data[currslot++];
        }
        virtual void read(float &f) {
          MEDDLY_CHECK_RANGE(0, currslot, ansLength);
          f = ((float*)(data + currslot))[0];
          currslot++;
        }
        virtual void read(long &L) {
          MEDDLY_CHECK_RANGE(0, 
            currslot+sizeof(long)/sizeof(node_handle), ansLength+1);
          memcpy(&L, data+currslot, sizeof(long));
          currslot += sizeof(long) / sizeof(node_handle);
        }
        virtual void read(double &D) {
          MEDDLY_CHECK_RANGE(0, 
            currslot+sizeof(double)/sizeof(node_handle), ansLength+1);
          memcpy(&D, data+currslot, sizeof(double));
          currslot += sizeof(double) / sizeof(node_handle);
        }
        virtual void read(void* &P) {
          MEDDLY_CHECK_RANGE(0, 
            currslot+sizeof(void*)/sizeof(node_handle), ansLength+1);
          memcpy(&P, data+currslot, sizeof(void*));
          currslot += sizeof(void*) / sizeof(node_handle);
        }

        inline void setResult(const node_handle* d, int sz) {
          setValid();
          data = d;
          currslot = 0;
#ifdef DEVELOPMENT_CODE
          ansLength = sz;
#endif
        }

    };



    class old_temp_entry : public compute_table::entry_builder {
          friend class MEDDLY::base_table;
          int handle;
          unsigned hash_value;
          node_handle* entry;
          node_handle* res_entry;
          unsigned resSlot;
          // The remaining entries are used only in development code
#ifdef DEVELOPMENT_CODE
          unsigned resLength;
#endif
        public:
          old_temp_entry() { }
          virtual ~old_temp_entry() { }

          virtual void writeResultNH(node_handle nh)
          {
            MEDDLY_CHECK_RANGE(0, resSlot, resLength);
            res_entry[resSlot++] = nh;
          }
          virtual void writeResult(int i)
          {
            MEDDLY_CHECK_RANGE(0, resSlot, resLength);
            res_entry[resSlot++] = i;
          }
          virtual void writeResult(float f)
          {
            MEDDLY_CHECK_RANGE(0, resSlot, resLength);
            float* x = (float*) res_entry + resSlot;
            x[0] = f;
            resSlot++;
          }
        protected:
          inline void writeResult(void* data, size_t slots)
          {
            MEDDLY_DCASSERT(slots>0);
            MEDDLY_DCASSERT(resSlot>=0);
            MEDDLY_DCASSERT(resSlot+slots<=resLength);
            memcpy(res_entry+resSlot, data, slots * sizeof(node_handle));
            resSlot += slots;
          }
        public:
          virtual void writeResult(long L)
          {
            writeResult(&L, sizeof(long) / sizeof(node_handle));
          }
          virtual void writeResult(double D)
          {
            writeResult(&D, sizeof(double) / sizeof(node_handle));
          }
          virtual void writeResult(void* P)
          {
            writeResult(&P, sizeof(void*) / sizeof(node_handle));
          }



          // The following are used by the compute table.
          inline const node_handle* readEntry(int off) const { return entry+off; }
          inline int readHandle() const { return handle; }
          inline unsigned getHash() const { return hash_value; }
          inline node_handle& data(int i) {
            return entry[i];
          }
    };

  public:
    base_table(const settings::computeTableSettings &s);
    virtual ~base_table();

  protected:
    old_temp_entry currEntry;
    

    inline void sawSearch(int c) {
      if (c>=stats::searchHistogramSize) {
        perf.numLargeSearches++;
      } else {
        perf.searchHistogram[c]++;
      }
      if (c>perf.maxSearchLength) perf.maxSearchLength = c;
    }
    int newEntry(int size);
    inline void recycleEntry(int h, int size) {
#ifdef DEBUG_CTALLOC
      fprintf(stderr, "Recycling entry %d size %d\n", h, size);
#endif
      entries[h] = freeList[size];
      freeList[size] = h;
      perf.numEntries--;
    }
    void dumpInternal(FILE* s, int verbLevel) const;
    void report(FILE* s, int indent, int &level) const;

    static inline bool equal_sw(const int* a, const int* b, int N) {
      switch (N) {  // note: cases 8 - 2 fall through
        case  8:    if (a[7] != b[7]) return false;
        case  7:    if (a[6] != b[6]) return false;
        case  6:    if (a[5] != b[5]) return false;
        case  5:    if (a[4] != b[4]) return false;
        case  4:    if (a[3] != b[3]) return false;
        case  3:    if (a[2] != b[2]) return false;
        case  2:    if (a[1] != b[1]) return false;
        case  1:    return a[0] == b[0];
        case  0:    return true;
        default:    return (0==memcmp(a, b, N*sizeof(int)));
      };
    }

    // M is 1 if we need a slot for the operation index, 0 otherwise.
    static inline search_key* init(operation* op, int M) {
      old_search_key* key = new old_search_key(op);
      MEDDLY_DCASSERT(0==key->data);
      key->data = new int[op->getKeyLength()+M];
      key->killData = true;
      key->key_data = key->data+M;
      if (M) key->data[0] = op->getIndex();
#ifdef DEVELOPMENT_CODE
      key->keyLength = op->getKeyLength();
#endif
      return key;
    }

    // M is 1 if we need a slot for the operation index, 0 otherwise.
    static inline search_key* init(int* data, operation* op, int M) {
      old_search_key* key = new old_search_key(op);
      MEDDLY_DCASSERT(0==key->data);
      key->data = data;
      key->killData = false;
      key->key_data = key->data+M;
      if (M) key->data[0] = op->getIndex();
#ifdef DEVELOPMENT_CODE
      key->keyLength = op->getKeyLength();
#endif
      return key;
    }
    
    // C is the number of slots we need for "chaining".
    // M is 1 if we need a slot for the operation index, 0 otherwise.
    inline void startIndexedEntry(old_search_key* Key, int C, int M) {
      MEDDLY_DCASSERT(Key);
      operation* op = Key->getOp();
      MEDDLY_DCASSERT(op);
      currEntry.hash_value = Key->getHash();
      currEntry.resSlot = 0;
      currEntry.handle = newEntry(op->getCacheEntryLength()+C+M);
      currEntry.entry = entries + currEntry.handle;
      memcpy(currEntry.entry+C, Key->rawData(), Key->dataLength()*sizeof(node_handle));
      currEntry.res_entry = currEntry.entry + Key->dataLength() + C;
#ifdef DEVELOPMENT_CODE
      currEntry.resLength = op->getAnsLength();
#endif
      op->doneCTkey(Key);
    }

    // C is the number of slots we need for "chaining".
    // M is 1 if we need a slot for the operation index, 0 otherwise.
    inline int* startPtrEntry(old_search_key* Key, int C, int M) {
      MEDDLY_DCASSERT(Key);
      operation* op = Key->getOp();
      MEDDLY_DCASSERT(op);
      currEntry.hash_value = Key->getHash();
      currEntry.resSlot = 0;
      int* data = new int[op->getKeyLength()+M+C];
      currEntry.entry = data;
      memcpy(currEntry.entry+C, Key->rawData(), Key->dataLength()*sizeof(node_handle));
      currEntry.res_entry = currEntry.entry + Key->dataLength() + C;
#ifdef DEVELOPMENT_CODE
      currEntry.resLength = op->getAnsLength();
#endif
      op->doneCTkey(Key);
      return data;
    }
    
  protected:
    int*  entries;
    int entriesSize;
    int entriesAlloc;

    unsigned long currMemory;
    unsigned long peakMemory;
  private:
    static const int maxEntrySize = 15;
    static const int maxEntryBytes = sizeof(int) * maxEntrySize;
    int* freeList;
};



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                          base_hash  class                          *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/** Abstract base class for hashing compute tables.
  * All use the same hash function:
  *
  *   Bob Jenkin's Hash
  *   Free to use for educational or commerical purposes
  *   http://burtleburtle.net/bob/hash/doobs.html
  */
class MEDDLY::base_hash : public base_table {
  public:
    base_hash(const settings::computeTableSettings &s, int initTsz, 
      int initTex);
    virtual ~base_hash();

  protected:
    inline unsigned hash(const int* k, int length) const {
      return raw_hash(k, length) % tableSize;
    }
    inline unsigned hash(search_key* k) const {
      old_search_key* key = smart_cast<old_search_key*>(k);
      MEDDLY_DCASSERT(key);
      key->setHash(raw_hash(key->rawData(), key->dataLength()));
      return key->getHash() % tableSize;
    }
    
    void dumpInternal(FILE* s, int verbLevel) const;
    void report(FILE* s, int indent, int &level) const;
  protected:
    int*  table;
    unsigned int tableSize;
    unsigned int tableExpand;
    unsigned int tableShrink;
  private:
    static inline unsigned rot(unsigned x, int k) {
        return (((x)<<(k)) | ((x)>>(32-(k))));
    }
    static inline void mix(unsigned &a, unsigned &b, unsigned &c) {
        a -= c;  a ^= rot(c, 4);  c += b; 
        b -= a;  b ^= rot(a, 6);  a += c;
        c -= b;  c ^= rot(b, 8);  b += a;
        a -= c;  a ^= rot(c,16);  c += b;
        b -= a;  b ^= rot(a,19);  a += c;
        c -= b;  c ^= rot(b, 4);  b += a;
    }
    static inline void final(unsigned &a, unsigned &b, unsigned &c) {
        c ^= b; c -= rot(b,14);
        a ^= c; a -= rot(c,11);
        b ^= a; b -= rot(a,25);
        c ^= b; c -= rot(b,16);
        a ^= c; a -= rot(c,4); 
        b ^= a; b -= rot(a,14);
        c ^= b; c -= rot(b,24);
    }
    static unsigned raw_hash(const int* k, int length);
};



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                         base_chained class                         *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/// Abstract base class for hash tables with chaining.
class MEDDLY::base_chained : public base_hash {
  public:
    base_chained(const settings::computeTableSettings &s);
    virtual ~base_chained();

    virtual void addEntry();
    virtual void removeStales();
  protected:
    void dumpInternal(FILE* s, int verbLevel) const;
    virtual int convertToList(bool removeStales) = 0;
    virtual void listToTable(int h) = 0;
    virtual void showEntry(FILE* s, int h) const = 0;
};



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                      monolithic_chained class                      *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/*
    Anatomy of an entry:

      entries[h]      : next pointer in chain
      entries[h+1]    : operation index
      entries[h+2]    : first "payload" item
      ...
      entries[h+L+1]  : last "payload" item, L is cache entry length.
*/
class MEDDLY::monolithic_chained : public base_chained {
  public:
    monolithic_chained(const settings::computeTableSettings &s);
    virtual ~monolithic_chained();

    // required functions

    virtual bool isOperationTable() const   { return false; }
    virtual search_key* initializeSearchKey(operation* op);
    virtual search_result& find(search_key *key);
    virtual entry_builder& startNewEntry(search_key *key);
    virtual void removeAll();

    virtual void show(FILE *s, int verbLevel);

  protected:
    virtual int convertToList(bool removeStales);
    virtual void listToTable(int h);
    virtual void showEntry(FILE* s, int h) const;

    inline bool checkStale(unsigned h, int prev, int &curr) {
        operation* currop = operation::getOpWithIndex(entries[curr+1]);
        MEDDLY_DCASSERT(currop);
        if (currop->isEntryStale(entries+curr+2)) {
          currop->discardEntry(entries+curr+2);
          int next = entries[curr];
          if (prev) entries[prev] = next;
          else      table[h] = next;
          int length = currop->getCacheEntryLength();
          recycleEntry(curr, length+2);
          curr = next;
          return true;
        }
        return false;
    }
};



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                      operation_chained  class                      *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/*
    Anatomy of an entry:

      entries[h]      : next pointer in chain
      entries[h+1]    : first "payload" item
      ...
      entries[h+L]    : last "payload" item, L is cache entry length.
*/
class MEDDLY::operation_chained : public base_chained {
  public:
    operation_chained(const settings::computeTableSettings &s, operation* op);
    virtual ~operation_chained();

    // required functions

    virtual bool isOperationTable() const   { return true; }
    virtual search_key* initializeSearchKey(operation* op);
    virtual entry_builder& startNewEntry(search_key *key);
    virtual void removeAll();

    virtual void show(FILE *s, int verbLevel);

  protected:
    virtual int convertToList(bool removeStales);
    virtual void listToTable(int h);
    virtual void showEntry(FILE* s, int h) const { 
      global_op->showEntry(s, entries + h + 1);
    }

    inline bool checkStale(unsigned h, int prev, int &curr) {
        if (global_op->isEntryStale(entries+curr+1)) {
          global_op->discardEntry(entries+curr+1);
          int next = entries[curr];
          if (prev) entries[prev] = next;
          else      table[h] = next;
          int length = global_op->getCacheEntryLength();
          recycleEntry(curr, length+1);
          curr = next;
          return true;
        }
        return false;
    }
  protected:
    operation* global_op;
};


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                    operation_chained_fast class                    *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  template <int N>
  class operation_chained_fast : public operation_chained {
    public:
      operation_chained_fast(const settings::computeTableSettings &s, 
        operation* op) : operation_chained(s, op) { }
      virtual ~operation_chained_fast() { }
      virtual search_result& find(search_key *key);
  };
};


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                        base_unchained class                        *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/// Abstract base class for hash tables without chaining.
class MEDDLY::base_unchained : public base_hash {
  public:
    base_unchained(const settings::computeTableSettings &s);
    virtual ~base_unchained();

    virtual void show(FILE *s, int verbLevel);
  protected:
    virtual void showTitle(FILE* s) const = 0;
    virtual void showEntry(FILE* s, int h) const = 0;

    inline void incMod(unsigned &h) {
      h++;
      if (h>=tableSize) h=0;
    }
    template <int M>
    inline void remove(int curr) {
        operation* currop = M 
          ?   operation::getOpWithIndex(entries[curr])
          :   global_op;
        MEDDLY_DCASSERT(currop);
#ifdef DEBUG_CT
        printf("Removing CT entry ");
        currop->showEntry(stdout, entries+curr+M);
        printf("\n");
#endif  
        currop->discardEntry(entries+curr+M);
        int length = currop->getCacheEntryLength();
        recycleEntry(curr, length+M);
    }
    template <int M>
    inline void setTable(unsigned h, int curr) {
        unsigned hfree = h;
        for (int i=maxCollisionSearch; i>=0; i--, incMod(hfree)) {
          // find a free slot
          if (0==table[hfree]) {
            table[hfree] = curr;
            return;
          }
        }
        // full; remove entry at our slot.
        collisions++;    
        remove<M>(table[h]);
        table[h] = curr;
    }
    template <int M>
    inline bool checkStale(unsigned h, int curr) {
        operation* currop = M 
          ?   operation::getOpWithIndex(entries[curr])
          :   global_op;
        MEDDLY_DCASSERT(currop);
        if (currop->isEntryStale(entries+curr+M)) {
#ifdef DEBUG_CT
          printf("Removing CT stale entry ");
          currop->showEntry(stdout, entries+curr+M);
          printf("\n");
#endif  
          currop->discardEntry(entries+curr+M);
          table[h] = 0;
          int length = currop->getCacheEntryLength();
          recycleEntry(curr, length+M);
          return true;
        }
        return false;
    }
    template <int M>
    inline void scanForStales() {
      for (unsigned i=0; i<tableSize; i++) {
        if (0==table[i]) continue;
        checkStale<M>(i, table[i]);
      }
    }
    template <int M>
    inline void rehashTable(int* oldT, unsigned oldS) {
        for (unsigned i=0; i<oldS; i++) {
          int curr = oldT[i];
          if (0==curr) continue;
          operation* currop = M 
            ?   operation::getOpWithIndex(entries[curr])
            :   global_op;
          MEDDLY_DCASSERT(currop);
          int hashlength = M+currop->getKeyLength();
          unsigned h = hash(entries + curr, hashlength);
          setTable<M>(h, curr);
        }
    }
    template <int M>
    inline void addEntryT() {
      unsigned h = currEntry.getHash() % tableSize;

#ifdef DEBUG_CT
      printf("Adding CT entry ");
      showEntry(stdout, currEntry.readHandle());
      // fprintf(stderr, " to slot %u", h);
      printf("\n");
#endif

      setTable<M>(h, currEntry.readHandle());
  
      if (perf.numEntries < tableExpand) return;

#ifdef DEBUG_SLOW
      fprintf(stdout, "Running GC in compute table (size %d, entries %ld)\n", 
        tableSize, perf.numEntries
      );
#endif

      if (checkStalesOnResize) {
        scanForStales<M>();
        if (perf.numEntries < tableExpand / 4) {
#ifdef DEBUG_SLOW
          fprintf(stdout, "Done CT GC, no resizing (now entries %ld)\n", 
            perf.numEntries
          );
#endif
          return;
        }
      }

      unsigned newsize = tableSize*2;
      if (newsize > maxSize) newsize = maxSize;
      if (tableSize == newsize) return;

      int* oldT = table;
      unsigned oldSize = tableSize;
      tableSize = newsize;
      table = (int*) malloc(newsize * sizeof(int));
      if (0==table) {
        table = oldT;
        tableSize = oldSize;
        throw error(error::INSUFFICIENT_MEMORY);
      }
      for (unsigned i=0; i<newsize; i++) table[i] = 0;

      currMemory += newsize * sizeof(int);

      rehashTable<M>(oldT, oldSize);
      free(oldT);

      currMemory -= oldSize * sizeof(int);
      if (currMemory > peakMemory) peakMemory = currMemory;

      if (tableSize == maxSize) {
        tableExpand = INT_MAX;
      } else {
        tableExpand = tableSize / 2;
      }
      tableShrink = tableSize / 8;

#ifdef DEBUG_SLOW
      fprintf(stdout, "CT enlarged to size %d\n", tableSize);
#endif
    }
    template <int M>
    inline void removeStalesT() {
#ifdef DEBUG_SLOW
      fprintf(stdout, "Removing stales in CT (size %d, entries %ld)\n", 
        tableSize, perf.numEntries
      );
#endif

      scanForStales<M>();

      if (perf.numEntries < tableShrink) {
        // shrink table
        unsigned newsize = tableSize / 2;
        if (newsize < 1024) newsize = 1024;
        if (newsize < tableSize) {
          int* oldT = table;
          unsigned oldSize = tableSize;
          tableSize = newsize;
          table = (int*) malloc(newsize * sizeof(int));
          if (0==table) {
            table = oldT;
            tableSize = oldSize;
            throw error(error::INSUFFICIENT_MEMORY);
          }
          for (unsigned i=0; i<newsize; i++) table[i] = 0;
          currMemory += newsize * sizeof(int);
    
          rehashTable<M>(oldT, oldSize);
          free(oldT);
      
          currMemory -= oldSize * sizeof(int);
          if (currMemory > peakMemory) peakMemory = currMemory;
    
          tableExpand = tableSize / 2;
          if (1024 == tableSize) {
            tableShrink = 0;
          } else {
            tableShrink = tableSize / 8;
          }
        } // if different size
      }
    
#ifdef DEBUG_SLOW
      fprintf(stdout, "Done removing CT stales (size %d, entries %ld)\n", 
        tableSize, perf.numEntries
      );
#endif
    }


    
  protected:
    static const int maxCollisionSearch = 2;
    long collisions;
    operation* global_op;
};



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                     monolithic_unchained class                     *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/*
    Anatomy of an entry:

      entries[h]      : operation index
      entries[h+1]    : first "payload" item
      ...
      entries[h+L]    : last "payload" item, L is cache entry length.
*/
class MEDDLY::monolithic_unchained : public base_unchained {
  public:
    monolithic_unchained(const settings::computeTableSettings &s);
    virtual ~monolithic_unchained();

    // required functions

    virtual bool isOperationTable() const   { return false; }
    virtual search_key* initializeSearchKey(operation* op);
    virtual search_result& find(search_key *key);
    virtual entry_builder& startNewEntry(search_key *key);
    virtual void addEntry();
    virtual void removeStales();
    virtual void removeAll();
  protected:
    virtual void showTitle(FILE *s) const;
    virtual void showEntry(FILE *s, int curr) const;
};


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                     operation_unchained  class                     *
// *                                                                    *
// *                                                                    *
// **********************************************************************

/*
    Anatomy of an entry:

      entries[h]      : first "payload" item
      ...
      entries[h+L-1]  : last "payload" item, L is cache entry length.
*/
class MEDDLY::operation_unchained : public base_unchained {
  public:
    operation_unchained(const settings::computeTableSettings &s, operation*);
    virtual ~operation_unchained();

    // required functions

    virtual bool isOperationTable() const   { return true; }
    virtual search_key* initializeSearchKey(operation* op);
    virtual entry_builder& startNewEntry(search_key *key);
    virtual void addEntry();
    virtual void removeStales();
    virtual void removeAll();
  protected:
    virtual void showTitle(FILE *s) const;
    virtual void showEntry(FILE *s, int curr) const;
};



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                   operation_unchained_fast class                   *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  template <int N>
  class operation_unchained_fast : public operation_unchained {
    public:
      operation_unchained_fast(const settings::computeTableSettings &s, 
        operation* op) : operation_unchained(s, op) { }
      virtual ~operation_unchained_fast() { }
      virtual search_result& find(search_key *key);
  };
};


// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                           base_map class                           *
// *                                                                    *
// *                                                                    *
// **********************************************************************

class MEDDLY::base_map : public base_table {
  protected:
    template <int K>
    class less {
      public:
        inline bool operator() (const int* p, const int* q) const {
          switch (K) {
              case 8:   if (p[7] < q[7]) return true;
                        if (p[7] > q[7]) return false;
              case 7:   if (p[6] < q[6]) return true;
                        if (p[6] > q[6]) return false;
              case 6:   if (p[5] < q[5]) return true;
                        if (p[5] > q[5]) return false;
              case 5:   if (p[4] < q[4]) return true;
                        if (p[4] > q[4]) return false;
              case 4:   if (p[3] < q[3]) return true;
                        if (p[3] > q[3]) return false;
              case 3:   if (p[2] < q[2]) return true;
                        if (p[2] > q[2]) return false;
              case 2:   if (p[1] < q[1]) return true;
                        if (p[1] > q[1]) return false;
              case 1:   return p[0] < q[0];
              case 0:   return false;
              default:  assert(0);
          }
        }
    };
  public:
    base_map(const settings::computeTableSettings &s, operation *op);
    virtual ~base_map();
    virtual bool isOperationTable() const   { return true; }
    virtual search_key* initializeSearchKey(operation* op);
    virtual entry_builder& startNewEntry(search_key *key);
  protected:
    inline void showEntry(FILE* s, int* h) const {
      global_op->showEntry(s, h);
    }
    inline bool isStale(const int* entry) {
      return global_op->isEntryStale(entry);
    }
    inline void removeEntry(int* h) {
      global_op->discardEntry(h);
      delete[] h;
    }
  protected:
    operation* global_op;
    int* search;
    int* current;
};



// **********************************************************************
// *                                                                    *
// *                                                                    *
// *                        operation_map  class                        *
// *                                                                    *
// *                                                                    *
// **********************************************************************

namespace MEDDLY {
  template <int K>
  class operation_map : public base_map {
    public:
      operation_map(const settings::computeTableSettings &s, operation* op);
      virtual ~operation_map();

      virtual search_result& find(search_key *key);

      virtual void addEntry();
      virtual void removeStales();
      virtual void removeAll();

      virtual void show(FILE *s, int verbLevel = 0);
    protected:
      std::map<int*, int*, less<K> > ct;
  };
};

#endif
