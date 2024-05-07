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


#include "ct_styles.h"

#include "../memory.h"
#include "../hash_stream.h"
#include "../node_storage.h"
#include "../compute_table.h"
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../ct_vector.h"
#include "../ct_initializer.h"
#include "../operators.h"

// #define DEBUG_FIND
// #define DEBUG_REHASH
// #define DEBUG_CT

#define USE_NEW_TEMPLATE

#ifndef USE_NEW_TEMPLATE
#include "ct_typebased.h"
#include "ct_none.h"
#endif

// #define HASH_EVERYTHING

// #define OLD_DEADCHECK

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***                                                                     ***
// ***                    New templatized compute table                    ***
// ***                                                                     ***
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

#define maxCollisionSearch 2

namespace MEDDLY {
    /**
        Template-based compute table implementation.
        The template parameters:
            @param  MONOLITHIC      If true, entries for several different
                                    operations will be stored in one table.
                                    That means each entry needs to store
                                    the operation index.

            @param  CHAINED         If true, we use hashing with chaining,
                                    and each entry needs a next pointer.

            @param  INTSLOTS        If true, we use unsigned integer slots
                                    instead of unsigned long slots for
                                    entries, and some items will require
                                    two slots instead of one.

        An entry is stored as follows.

            +---------------------------------------+
            |                                       |
            |     next pointer  (unsigned long)     |
            |                                       |
            | Only needed when CHAINED is true.     |
            |                                       |
            +---------------------------------------+
            |                                       |
            |      operation index  (unsigned)      |
            |                                       |
            | Only needed when MONOLITHIC is true.  |
            |                                       |
            +---------------------------------------+
            |                                       |
            |          key size (unsigned)          |
            |                                       |
            | Only needed when the key has a        |
            | repeatable portion; otherwise the key |
            | size is fixed for the operation.      |
            |                                       |
            +---------------------------------------+
            |                                       |
            |               key items               |
            |                                       |
            | INTSLOTS true: each item requires one |
            |   or two slots depending on its type. |
            | INTSLOTS false: each item requires    |
            |   one (larger) slot.                  |
            |                                       |
            +---------------------------------------+
            |                                       |
            |             result  items             |
            |                                       |
            | INTSLOTS true: each item requires one |
            |   or two slots depending on its type. |
            | INTSLOTS false: each item requires    |
            |   one (larger) slot.                  |
            |                                       |
            +---------------------------------------+


    */
    template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
    class ct_tmpl : public compute_table {
        public:
            ct_tmpl(const ct_settings &s, operation* op, unsigned slot);
            virtual ~ct_tmpl();

            // Required functions

#ifdef ALLOW_DEPRECATED_0_17_6
            virtual void find(ct_entry_key* key, ct_entry_result &res);
            virtual void addEntry(ct_entry_key* key,
                const ct_entry_result& res);
            virtual void doneKey(ct_entry_key* key);
#endif

            virtual bool find(const ct_entry_type &ET, ct_vector &key,
                    ct_vector &res)
            {
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            }
            virtual void addEntry(const ct_entry_type &ET, ct_vector &key,
                    const ct_vector &res)
            {
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
            }

            // calls removeStaleEntries, then maybe shrinks table
            virtual void removeStales();
            virtual void removeAll();
            virtual void show(output &s, int verbLevel = 0);
            virtual void countNodeEntries(const forest* f, size_t* counts)
                const;

        protected:
            // Helper functions

            void removeStaleEntries();

#ifdef ALLOW_DEPRECATED_0_17_6
            void* key2entry(const ct_entry_key& key, void* e);
            static void result2entry(const ct_entry_result& res, void* e);
#endif
            /*
            inline static ct_entry_item* vector2entry(const ct_vector &v,
                    ct_entry_item* e)
            {
                MEDDLY_DCASSERT(!INTSLOTS);
                for (unsigned i=0; i<v.size(); i++)
                {
                    v[i].get(*e);
                    e++;
                }
                return e;
            }
            inline static ct_entry_item* vector2entry(const ct_vector &v,
                    ct_entry_item* e, hash_stream &H)
            {
                MEDDLY_DCASSERT(!INTSLOTS);
                for (unsigned i=0; i<v.size(); i++)
                {
                    v[i].get(*e, H);
                    e++;
                }
                return e;
            }
            inline static unsigned* vector2entry(const ct_vector &v,
                    unsigned* e)
            {
                MEDDLY_DCASSERT(INTSLOTS);
                for (unsigned i=0; i<v.size(); i++)
                {
                    e += v[i].getRaw(e);
                }
                return e;
            }
            inline static unsigned* vector2entry(const ct_vector &v,
                    unsigned* e, hash_stream &H)
            {
                MEDDLY_DCASSERT(INTSLOTS);
                for (unsigned i=0; i<v.size(); i++)
                {
                    e += v[i].getRaw(e, H);
                }
                return e;
            }
            */

            //
            //  Short vector equality check
            //
            inline static bool equal_sw(const unsigned* a, const unsigned* b,
                    unsigned N)
            {
                switch (N) {  // note: cases 12 - 2 fall through
                    case 12:    if (a[11] != b[11]) return false;
                    case 11:    if (a[10] != b[10]) return false;
                    case 10:    if (a[9] != b[9]) return false;
                    case  9:    if (a[8] != b[8]) return false;
                    case  8:    if (a[7] != b[7]) return false;
                    case  7:    if (a[6] != b[6]) return false;
                    case  6:    if (a[5] != b[5]) return false;
                    case  5:    if (a[4] != b[4]) return false;
                    case  4:    if (a[3] != b[3]) return false;
                    case  3:    if (a[2] != b[2]) return false;
                    case  2:    if (a[1] != b[1]) return false;
                    case  1:    return a[0] == b[0];
                    case  0:    return true;
                    default:    return (0==memcmp(a, b, N*sizeof(unsigned)));
                };
            }

            //
            //  Short vector equality check
            //
            inline static bool equal_sw(const ct_entry_item* a,
                    const ct_entry_item* b, unsigned N)
            {
                switch (N) {  // note: cases 12 - 2 fall through
                    case 12:    if (a[11].UL != b[11].UL) return false;
                    case 11:    if (a[10].UL != b[10].UL) return false;
                    case 10:    if (a[9].UL != b[9].UL) return false;
                    case  9:    if (a[8].UL != b[8].UL) return false;
                    case  8:    if (a[7].UL != b[7].UL) return false;
                    case  7:    if (a[6].UL != b[6].UL) return false;
                    case  6:    if (a[5].UL != b[5].UL) return false;
                    case  5:    if (a[4].UL != b[4].UL) return false;
                    case  4:    if (a[3].UL != b[3].UL) return false;
                    case  3:    if (a[2].UL != b[2].UL) return false;
                    case  2:    if (a[1].UL != b[1].UL) return false;
                    case  1:    return a[0].UL == b[0].UL;
                    case  0:    return true;
                    default:    return (0==memcmp(a, b, N*sizeof(ct_entry_item)));
                };
            }

            //
            // Chain management
            //

            inline static unsigned long getNext(const void* raw)
            {
                MEDDLY_DCASSERT(CHAINED);
                if (INTSLOTS) {
                    const unsigned* u = (const unsigned*) raw;
                    unsigned long n = u[0];
                    n <<= 32;
                    n |= u[1];
                    return n;
                } else {
                    const ct_entry_item* ct = (const ct_entry_item*) raw;
                    return ct[0].UL;
                }
            }
            inline static void setNext(void* raw, unsigned long n)
            {
                MEDDLY_DCASSERT(CHAINED);
                if (INTSLOTS) {
                    unsigned* u = (unsigned*) raw;
                    u[0] = (n >> 32);
                    u[1] = n & 0xffffffff;
                } else {
                    ct_entry_item* ct = (ct_entry_item*) raw;
                    ct[0].UL = n;
                }
            }

            /// Update stats: we just searched through c items
            inline void sawSearch(unsigned c) {
                if (c>=stats::searchHistogramSize) {
                    perf.numLargeSearches++;
                } else {
                    perf.searchHistogram[c]++;
                }
                if (c>perf.maxSearchLength) perf.maxSearchLength = c;
            }

            //
            // Should we discard a CT hit?
            //

#ifdef ALLOW_DEPRECATED_0_17_6

            /**
                Copy the result portion of an entry, and check if it is dead.
                    @param  ET      Entry type information
                    @param  sres    Pointer to result portion of entry
                    @param  dres    Where to copy the result
             */
            static bool isDead(const ct_entry_type &ET, const void* sres,
                    ct_entry_result &dres);

#endif

            //
            // Should we discard a CT miss?
            //
            /**
                See if an entry is stale.
                    @param  e       Pointer to the entire entry
                    @param  mark    If true, mark nodes in the entry
             */
            bool isStale(const void* e, bool mark) const;

            /**
                Check if the key portion of an entry equals Key.
                Check if the entry should be discarded (either dead or stale).

                @param  ET      Entry type information
                @param  entry   Complete entry in CT to check.
                @param  key     Key to compare against.
                @param  keylen  Length of key portion, in 'slots'.
                @param  dres    Where to copy the result
                @param  discard On output: should this entry be discarded

                @return true if equal, false if not equal.
            */
            /*
            inline bool equalAndUsable(const ct_entry_type &ET,
                    void* entry, void* key, unsigned keylen,
                    ct_entry_result &dres, bool &discard)
            {
                // TBD THIS?
            }
            */


            /**
                Try to set table[h] to curr.
                If the slot is full, check ahead the next couple for empty.
                If none are empty, then recycle current table[h]
                and overwrite it.
            */
            inline void setTable(unsigned long h, unsigned long curr)
            {
                MEDDLY_DCASSERT(!CHAINED);
                unsigned long hfree = h;
                //
                // Look for a free slot
                //
                unsigned i=0;
                for (;;) {
                    if (!table[hfree]) {
                        table[hfree] = curr;
                        return;
                    }
                    ++i;
                    if (i > maxCollisionSearch) break;
                    hfree = (1+hfree) % table.size();
                }

                //
                // Nothing free; remove entry at our slot.
                //
                collisions++;
                MEDDLY_DCASSERT(table[h]);
                deleteEntry(table[h]);
                table[h] = curr;
            }



            //
            // batch deletion of table entries in array toDelete.
            // These have already been removed from the table,
            // but now we update reference counts and recycle
            // the memory, and other bookkeeping tasks.
            //
            inline void batchDelete()
            {
                while (toDelete.size()) {
                    // std::cout << "batchdel " << toDelete.back();

                    deleteEntry(toDelete.back());

                    // std::cout << " now " << perf.numEntries << " entries\n";

                    toDelete.pop_back();
                }
            }

            /**
               Delete an entry.
                @param  hnd         'pointer' (in the memory manager)
                                    to the entire entry. Will be set to 0.
            */
            void deleteEntry(unsigned long &hnd);

            /**
                Resize the table.
                    @param  newsz   New size to use
            */
            void resizeTable(unsigned long newsz);

            /**
                Compute the hash of an entry.
                    @param  entry   Raw pointer to entire entry

                    @return     Hashed value; still needs to be modded
                                by hash table size.

                TBD - switch to 64-bit hash
             */
            unsigned hashEntry(const void* entry);


            /**
                Display an entry.
                Used for debugging, and by method(s) to display the entire CT.
                    @param  s         Stream to write to.
                    @param  h         Handle of the entry.
                    @param  keyonly     true: we don't have a result yet
                                        false: it's a proper entry with both
            */
            void showEntry(output &s, unsigned long h, bool keyonly) const;

            /// Display a chain
            inline void showChain(output &s, unsigned long L) const
            {
                s << L;
                if (CHAINED) {
                    while (L) {
                        L = getNext(MMAN->getChunkAddress(L));
                        s << "->" << L;
                    }
                }
                s << "\n";
            }

            /// Convert unsigned pointer to node
            static inline node_handle u2n(const unsigned* ptr)
            {
                return *((const node_handle*)ptr);
            }

        private:
            //
            // Abstraction layer for hashing array
            //
            inline void HSreset() {
#ifndef HASH_EVERYTHING
                hsi = 0;
#endif
            }
            inline void HSreserve(unsigned sz) {
#ifndef HASH_EVERYTHING
                if (sz > HS.size()) {
                    HS.resize(sz);
                }
                hsi = 0;
#endif
            }
            inline unsigned HSpush(unsigned x) {
#ifdef HASH_EVERYTHING
                return x;
#else
#ifdef DEVELOPMENT_CODE
                return HS.at(hsi++) = x;
#else
                return HS[hsi++] = x;
#endif
#endif
            }
#ifndef HASH_EVERYTHING
            inline unsigned hash32() {
                return hash_stream::raw_hash(HS.data(), hsi);
            }
#endif

        private:
#ifndef HASH_EVERYTHING
            /// Hashing array
            std::vector <unsigned> HS;

            /// Current index in hashing array
            unsigned hsi;
#endif

            /// List of entries to delete
            std::vector <unsigned long> toDelete;

            /// Hash table
            std::vector <unsigned long> table;

            /// When to next expand the table
            unsigned long tableExpand;

            /// When to next shrink the table
            unsigned long tableShrink;

            /// Manager for the entries
            memory_manager* MMAN;

            /// Our memory stats
            memstats mstats;

            /// Stats: count collisions
            unsigned long collisions;
    }; // class ct_tmpl

};

// **********************************************************************
// *                                                                    *
// *                          ct_tmpl  methods                          *
// *                                                                    *
// **********************************************************************

template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
MEDDLY::ct_tmpl<MONOLITHIC, CHAINED, INTSLOTS>::ct_tmpl(
    const ct_settings &s, operation* op, unsigned slot)
    : compute_table(s, op, slot)
{
    if (MONOLITHIC) {
        MEDDLY_DCASSERT(0==op);
        MEDDLY_DCASSERT(0==slot);
    } else {
        MEDDLY_DCASSERT(op);
    }

    /*
        Initialize memory management for entries.
    */
    MEDDLY_DCASSERT(s.MMS);
    if (INTSLOTS) {
        MMAN = s.MMS->initManager(sizeof(unsigned), 2, mstats);
    } else {
        MMAN = s.MMS->initManager(sizeof(ct_entry_item), 2, mstats);
    }

    /*
        Initialize hash table
    */
    tableExpand = CHAINED ? 4*1024 : 512;
    tableShrink = 0;
    table.resize(1024, 0);

    mstats.incMemUsed(table.size() * sizeof(unsigned long));
    mstats.incMemAlloc(table.size() * sizeof(unsigned long));

    collisions = 0;
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
MEDDLY::ct_tmpl<MONOLITHIC, CHAINED, INTSLOTS>::~ct_tmpl()
{
    /*
        Clean up memory manager for entries
    */
    delete MMAN;

    /*
        Update stats: important for global usage
    */
    mstats.zeroMemUsed();
    mstats.zeroMemAlloc();

    //
    // table will be destroyed automatically
}

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<MONOLITHIC,CHAINED,INTSLOTS>::find(ct_entry_key* key,
        ct_entry_result &res)
{
    MEDDLY_DCASSERT(key);
    const ct_entry_type* et = key->getET();
    MEDDLY_DCASSERT(et);
    if (!MONOLITHIC) {
        if (et != global_et)
            throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }

    //
    // Go ahead and start building an entry (the key portion anyway)
    // as this makes searching for matches easier.
    //

    union genentry {
        unsigned*       uptr;
        ct_entry_item*  ctptr;
        void*           vptr;
    };

    //
    // (1) determine entry size, and allocate
    //

    const unsigned chain_slots =
        CHAINED ? (INTSLOTS ? sizeof(unsigned long)/sizeof(unsigned) : 1) : 0;

    const unsigned op_slots = MONOLITHIC ? 1 : 0;

    const unsigned key_slots =
        INTSLOTS  ? et->getKeyIntslots(key->numRepeats())
                  : et->getKeySize(key->numRepeats());

    const unsigned res_slots =
        INTSLOTS  ? et->getResultIntslots()
                  : et->getResultSize();

    const unsigned entry_slots =
          op_slots
        + (et->isRepeating() ? 1 : 0)
        + key_slots
    ;

    size_t num_slots =
          chain_slots
        + entry_slots
        + res_slots
    ;

    key->entry_slots = num_slots;
    key->my_entry = MMAN->requestChunk(num_slots);
    if (num_slots < key->entry_slots) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    MEDDLY_DCASSERT(key->my_entry);

    genentry keyentry, keyafterchain;
    keyentry.vptr = MMAN->getChunkAddress(key->my_entry);

    //
    // (2) build the entry, except for the result.
    //     we hash it as we go.
    //

    HSreserve(num_slots);
    if (INTSLOTS) {
        unsigned* e = keyentry.uptr;
        if (CHAINED) {
            e += sizeof(unsigned long)/sizeof(unsigned);
            // skip over NEXT pointer slot(s)
        }
        keyafterchain.uptr = e;
        if (MONOLITHIC) {
            *e = HSpush(et->getID());
            ++e;
        }
        if (et->isRepeating()) {
            *e = HSpush(key->numRepeats());
            ++e;
        }
        key->result_shift = (unsigned*) key2entry(*key, e) - keyentry.uptr;
    } else {
        ct_entry_item* e = keyentry.ctptr;
        if (CHAINED) {
            ++e;    // skip over NEXT pointer
        }
        keyafterchain.ctptr = e;
        if (MONOLITHIC) {
            e->raw[0] = HSpush(et->getID());
            e->raw[1] = 0;
            ++e;
        }
        if (et->isRepeating()) {
            e->raw[0] = HSpush(key->numRepeats());
            e->raw[1] = 0;
            ++e;
        }
        key->result_shift = (ct_entry_item*) key2entry(*key, e) - keyentry.ctptr;
    }

    //
    // (3) Save hash value
    //

#ifdef HASH_EVERYTHING
    if (INTSLOTS) {
        setHash(key, hash_stream::raw_hash(keyafterchain.uptr, entry_slots));
    } else {
        setHash(key, hash_stream::raw_hash(keyafterchain.uptr, 2*entry_slots));
    }
#else
    setHash(key, hash32());
#endif
    MEDDLY_DCASSERT(key->getHash() == hashEntry(keyentry.vptr));
    const unsigned long hslot = key->getHash() % table.size();
    unsigned long hcurr = hslot;
    // TBD: use a 64-bit hash

#ifdef DEBUG_FIND
    ostream_output out(std::cout);
    out << "Searching for entry ";
    showEntry(out, key->my_entry, true);
    out << " at slot " << hcurr << "...\n";
#endif


    //
    // (4) Search the table with the following policy.
    //
    //  equal:
    //      dead:   discard and stop searching, return not found
    //      !dead:  return as found
    //
    //  !equal:
    //      stale:  discard and continue search
    //      !stale: keep and continue search
    //
    //  We use one loop that does both linked list traversal
    //  (for chained) or check next hash slots (for not chained).
    //
    perf.pings++;
    unsigned chainlen;
#ifdef DEVELOPMENT_CODE
    unsigned long currhnd = table.at(hcurr);
#else
    unsigned long currhnd = table[hcurr];
#endif
    unsigned long prevhnd = 0;
    for (chainlen = 0; ; )
    {
        if (!currhnd) {
            if (CHAINED) break;
            // ^ we got to the end of the list

            //
            // Advance
            //
            ++chainlen;
            if (chainlen > maxCollisionSearch) break;
            // ^ we're not chained, only check the next few slots
            // and break once we exceed that.
            hcurr = (1+hcurr) % table.size();
#ifdef DEVELOPMENT_CODE
            currhnd = table.at(hcurr);
#else
            currhnd = table[hcurr];
#endif
            continue;
        } else {
            if (CHAINED) chainlen++;
        }

        //
        // Get the next entry
        //
        genentry currentry, currafterchain, currres;
        currentry.vptr = MMAN->getChunkAddress(currhnd);
        if (INTSLOTS) {
            currafterchain.uptr = currentry.uptr + chain_slots;
            currres.uptr = currafterchain.uptr + entry_slots;
        } else {
            currafterchain.ctptr = currentry.ctptr + chain_slots;
            currres.ctptr = currafterchain.ctptr + entry_slots;
        }

#ifdef DEBUG_FIND
        out << "        trying node " << currhnd << "\n";
#endif

        //
        // Check if it is equal to the entry we built,
        // and if we need to remove the entry.
        //
        bool equal, remove;

        if (INTSLOTS) {
            if (equal_sw(currafterchain.uptr, keyafterchain.uptr, entry_slots))
            {
                equal = true;
                remove = isDead(*et, currres.vptr, res);
            } else {
                equal = false;
                if (checkStalesOnFind) {
                    remove = isStale(currentry.vptr, false);
                } else {
                    remove = false;
                }
            }
        } else {
            if (equal_sw(currafterchain.ctptr, keyafterchain.ctptr, entry_slots))
            {
                equal = true;
                remove = isDead(*et, currres.vptr, res);
            } else {
                equal = false;
                if (checkStalesOnFind) {
                    remove = isStale(currentry.vptr, false);
                } else {
                    remove = false;
                }
            }
        }

#ifdef DEBUG_FIND
        out << "            equal: " << (equal ? "yes" : "no")
            << ", remove: " << (remove ? "yes\n" : "no\n");
#endif

        if (remove) {
            //
            // Remove the entry
            //
            if (CHAINED) {
                unsigned long next = getNext(currentry.vptr);
                if (prevhnd) {
                    void* prev = MMAN->getChunkAddress(prevhnd);
                    setNext(prev, next);
                } else {
#ifdef DEVELOPMENT_CODE
                    table.at(hcurr) = next;
#else
                    table[hcurr] = next;
#endif
                }
            } else {
#ifdef DEVELOPMENT_CODE
                table.at(hcurr) = 0;
#else
                table[hcurr] = 0;
#endif
            }
            toDelete.push_back(currhnd);

            if (equal) {
                //
                // There's no hope of finding another equal entry,
                // so stop searching.
                //
                break;
            }
        } // if remove

        if (equal) {
            //
            // Move to front if we're chained
            //
            if (CHAINED) {
                if (prevhnd) {
                    void* prev = MMAN->getChunkAddress(prevhnd);
                    setNext(prev, getNext(currentry.vptr));
                    setNext(currentry.vptr, table[hcurr]);
                    table[hcurr] = currhnd;
#ifdef DEBUG_FIND
                    out << "        moved to front\n";
#endif
                }
                // if prev is null, then
                // we are already at the front.
#ifdef DEBUG_FIND
                out << "        new chain: ";
                showChain(out, table[hcurr]);
#endif
            }

            //
            // Set up the result portion
            //
            res.reset();    // rewind
            res.setValid();
            sawSearch(chainlen);
            perf.hits++;
            batchDelete();
            MMAN->recycleChunk(key->my_entry, key->entry_slots);
            key->entry_slots = 0;
            return;
        }

        //
        // Advance to the next entry we need to check
        //
        if (CHAINED) {
            //
            // Advance previous, unless we removed the current entry
            //
            if (!remove) {
                prevhnd = currhnd;
            }
            currhnd = getNext(currentry.vptr);
        } else {
            ++chainlen;
            if (chainlen > maxCollisionSearch) break;
            hcurr = (1+hcurr) % table.size();
#ifdef DEVELOPMENT_CODE
            currhnd = table.at(hcurr);
#else
            currhnd = table[hcurr];
#endif
        }


    } // for chainlen

    //
    // not found
    //

    sawSearch(chainlen);
    res.setInvalid();
    batchDelete();
#ifdef DEBUG_FIND
    out << "        not found\n";
#endif
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<MONOLITHIC,CHAINED,INTSLOTS>::addEntry(ct_entry_key* key,
        const ct_entry_result& res)
{
    MEDDLY_DCASSERT(key);
    const ct_entry_type* et = key->getET();
    MEDDLY_DCASSERT(et);
    if (!MONOLITHIC) {
        if (et != global_et)
            throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }

    unsigned h = key->getHash() % table.size();

    //
    // Increment cache counters for nodes in the key
    //
    key->cacheNodes();

    //
    // Copy the result portion, and increment cache counters in the result
    //
    MEDDLY_DCASSERT(key->my_entry);
    void* rawentry = MMAN->getChunkAddress(key->my_entry);
    if (INTSLOTS) {
        unsigned* entry = (unsigned*) rawentry;
        result2entry(res, entry + key->result_shift);
    } else {
        ct_entry_item* entry = (ct_entry_item*) rawentry;
        result2entry(res, entry + key->result_shift);
    }

    //
    // Add to the CT
    //
    if (CHAINED) {
        // Add this to front of chain
        setNext(rawentry, table[h]);
        table[h] = key->my_entry;
    } else {
        setTable(h, key->my_entry);
    }

    //
    // Recycle key
    //
    key->my_entry = 0;
    recycle(key);

    //
    // Is it time to GC / resize the table?
    //
    perf.numEntries++;
    if (perf.numEntries < tableExpand) return;
    perf.resizeScans++;

    //
    // Clear CT bits for all forests we affect.
    // Should be a no-op for forests that use reference counts.
    //
    const unsigned NF = forest::MaxFID()+1;
    bool* skipF = new bool[NF];
    for (unsigned i=0; i<NF; i++) {
        skipF[i] = false;
    }
    clearForestCTBits(skipF, NF);
    removeStaleEntries();
    sweepForestCTBits(skipF, NF);
    delete[] skipF;

    //
    // Do we still need to enlarge the table?
    //
    if (CHAINED) {
        if (perf.numEntries < table.size()) return;
    } else {
        if (perf.numEntries < tableExpand / 4) return;
    }
    unsigned long newsize = table.size() * 2;
    if (newsize > maxSize) newsize = maxSize;
    if (newsize == table.size()) return;    // already max size

    //
    // Enlarge table
    //
    resizeTable(newsize);

    //
    // Update expand/shrink sizes
    //
    if (CHAINED) {
        if (table.size() == maxSize) {
            tableExpand = std::numeric_limits<int>::max();
        } else {
            tableExpand = 4*table.size();
        }
        tableShrink = table.size() / 2;
    } else {
        if (table.size() == maxSize) {
            tableExpand = std::numeric_limits<int>::max();
        } else {
            tableExpand = table.size() / 2;
        }
        tableShrink = table.size() / 8;
    }
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
void MEDDLY::ct_tmpl<M,C,I>::doneKey(ct_entry_key* key)
{
    if (key) {
        if (key->my_entry) {
            MMAN->recycleChunk(key->my_entry, key->entry_slots);
            key->entry_slots = 0;
        }
    }
}

#endif

// **********************************************************************

template <bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<M, CHAINED, I>::removeStales()
{
    removeStaleEntries();

    //
    // Is it time to shrink the table?
    //
    if (perf.numEntries >= tableShrink) return;
    unsigned long newsize = table.size() / 2;
    if (newsize < 1024) newsize = 1024;
    if (newsize == table.size()) return;

    //
    // Shrink table
    //
    resizeTable(newsize);

    //
    // Update expand/shrink sizes
    //
    if (CHAINED) {
        tableExpand = 4*table.size();
        if (1024 == table.size()) {
            tableShrink = 0;
        } else {
            tableShrink = table.size() / 2;
        }
    } else {
        tableExpand = table.size() / 2;
        if (1024 == table.size()) {
            tableShrink = 0;
        } else {
            tableShrink = table.size() / 8;
        }
    }
}

// **********************************************************************

template <bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<M,CHAINED,I>::removeAll()
{
    for (unsigned long i=0; i<table.size(); i++) {
        while (table[i]) {
            unsigned long curr = table[i];
            if (CHAINED) {
                table[i] = getNext(MMAN->getChunkAddress(curr));
            } else {
                table[i] = 0;
            }
            deleteEntry(curr);
        } // while
    } // for
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<MONOLITHIC,CHAINED,INTSLOTS>
    ::show(output &s, int verbLevel)
{
    if (verbLevel < 1) return;

    if (MONOLITHIC) {
        s << "Monolithic compute table\n";
    } else {
        s << "Compute table for " << global_et->getName() << " (index "
          << long(global_et->getID()) << ")\n";
    }

    s.put("", 6);
    s << "Current CT memory   :\t" << mstats.getMemUsed() << " bytes\n";
    s.put("", 6);
    s << "Peak    CT memory   :\t" << mstats.getPeakMemUsed() << " bytes\n";
    s.put("", 6);
    s << "Current CT alloc'd  :\t" << mstats.getMemAlloc() << " bytes\n";
    s.put("", 6);
    s << "Peak    CT alloc'd  :\t" << mstats.getPeakMemAlloc() << " bytes\n";
    if (!CHAINED) {
        s.put("", 6);
        s << "Collisions          :\t" << long(collisions) << "\n";
    }
    s.put("", 6);
    s << "Hash table size     :\t" << long(table.size()) << "\n";
    s.put("", 6);
    s << "Number of entries   :\t" << long(perf.numEntries) << "\n";

    if (--verbLevel < 1) return;

    s.put("", 6);
    s << "Pings               :\t" << long(perf.pings) << "\n";
    s.put("", 6);
    s << "Hits                :\t" << long(perf.hits) << "\n";
    s.put("", 6);
    s << "Resize (GC) scans   :\t" << long(perf.resizeScans) << "\n";

    if (--verbLevel < 1) return;

    s.put("", 6);
    s << "Search length histogram:\n";
    for (unsigned i=0; i<stats::searchHistogramSize; i++) {
        if (perf.searchHistogram[i]) {
            s.put("", 10);
            s.put(long(i), 3);
            s << ": " << long(perf.searchHistogram[i]) << "\n";
        }
    }
    if (perf.numLargeSearches) {
        s.put("", 6);
        s << "Searches longer than " << long(stats::searchHistogramSize-1)
          << ": " << long(perf.numLargeSearches) << "\n";
    }
    s.put("", 6);
    s << "Max search length: " << long(perf.maxSearchLength) << "\n";


    if (--verbLevel < 1) return;

    s << "Hash table:\n";

    for (unsigned long i=0; i<table.size(); i++) {
        if (0==table[i]) continue;
        s << "table[";
        s.put(long(i), 9);
        s << "]: ";
        showChain(s, table[i]);
    }

    if (--verbLevel < 1) return;

    s << "\nHash table nodes:\n";

    for (unsigned long i=0; i<table.size(); i++) {
        unsigned long curr = table[i];
        while (curr) {
            s << "\tNode ";
            s.put(long(curr), 9);
            s << ":  ";
            showEntry(s, curr, false);
            s.put('\n');
            if (CHAINED) {
                curr = getNext(MMAN->getChunkAddress(curr));
            } else {
                curr = 0;
            }
        } // while curr
    } // for i
    s.put('\n');

    if (--verbLevel < 1) return;

    if (MMAN) MMAN->dumpInternal(s);
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<MONOLITHIC,CHAINED,INTSLOTS>
    ::countNodeEntries(const forest* f, size_t* counts) const
{
    for (unsigned long h=0; h<table.size(); h++) {
        unsigned long curr = table[h];
        while (curr) {
            const void* entry = MMAN->getChunkAddress(curr);
            const unsigned* ue = (const unsigned*) entry;
            const ct_entry_item* cte = (const ct_entry_item*) entry;

            //
            // Skip over chained portion
            //
            if (CHAINED) {
                if (INTSLOTS) {
                    ue += sizeof(unsigned long) / sizeof(unsigned);
                } else {
                    ++cte;
                }
            }

            //
            // Get operation and advance
            //
            const ct_entry_type* et = MONOLITHIC
                ?   getEntryType( INTSLOTS ? *ue : cte->U )
                :   global_et;
            if (MONOLITHIC) {
                if (INTSLOTS) {
                    ++ue;
                } else {
                    ++cte;
                }
            }

            //
            // Get key size and advance
            //
            unsigned reps;
            if (et->isRepeating()) {
                if (INTSLOTS) {
                    reps = *ue;
                    ++ue;
                } else {
                    reps = cte->U;
                    ++cte;
                }
            } else {
                reps = 0;
            }

            //
            // Go through key portion
            //
            ct_item x;
            const unsigned klen = et->getKeySize(reps);
            for (unsigned i=0; i<klen; i++) {
                const ct_itemtype &item = et->getKeyType(i);
                if (INTSLOTS) {
                    ue = x.set(item.getType(), ue);
                } else {
                    x.set(item.getType(), *cte);
                    ++cte;
                }

                if (item.hasForest(f)) {
                    if (x.getN() > 0) {
                        ++counts[ x.getN() ];
                    }
                }
            } // for i

            //
            // Go through result portion
            //
            for (unsigned i=0; i<et->getResultSize(); i++) {
                const ct_itemtype &item = et->getResultType(i);
                if (INTSLOTS) {
                    ue = x.set(item.getType(), ue);
                } else {
                    x.set(item.getType(), *cte);
                    ++cte;
                }

                if (item.hasForest(f)) {
                    if (x.getN() > 0) {
                        ++counts[ x.getN() ];
                    }
                }
            } // for i

            //
            // Advance in chain
            //
            if (CHAINED) {
                curr = getNext(entry);
            } else {
                curr = 0;
            }
        } // while curr
    } // for h

}

// **********************************************************************
//
//      Helper methods
//
// **********************************************************************

template <bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<M, CHAINED, I>::removeStaleEntries()
{
    if (CHAINED) {
        // Chained
        ostream_output out(std::cout);
        for (unsigned long i=0; i<table.size(); i++) {
            if (!table[i]) continue;

            void* preventry = nullptr;
            unsigned long curr;
            for (curr = table[i]; curr; ) {
                void* entry = MMAN->getChunkAddress(curr);
                unsigned long next = getNext(entry);

                if (isStale(entry, true)) {
                    // remove from list
                    if (preventry) {
                        setNext(preventry, next);
                    } else {
                        table[i] = next;
                    }
                    // add to list of entries to delete
                    toDelete.push_back(curr);
                } else {
                    // Not stale.
                    // Keep this entry, which means
                    // we can advance the previous entry.
                    preventry = entry;
                }

                // advance
                curr = next;
            } // for curr

            batchDelete();
        } // for i

        // End of chained
    } else {
        // Not chained
        for (unsigned long i=0; i<table.size(); i++) {
            if (!table[i]) continue;

            if (!isStale(MMAN->getChunkAddress(table[i]), true))
            {
                continue;
            }

            deleteEntry(table[i]);
            table[i] = 0;
        } // for i
        // End of not chained
    }
}

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
void* MEDDLY::ct_tmpl<M,C,I>::key2entry(const ct_entry_key& key, void* e)
{
    const ct_entry_type* et = key.getET();
    MEDDLY_DCASSERT(et);
    const unsigned keylen = et->getKeySize(key.numRepeats());
    const ct_entry_item* data = key.rawData();

    unsigned* ue = (unsigned*) e;
    ct_entry_item* cte = (ct_entry_item*) e;

    for (unsigned i=0; i<keylen; i++) {
        const ct_itemtype &it = et->getKeyType(i);

#ifndef HASH_EVERYTHING
        if (it.shouldBeHashed()) {
#endif
            if (it.requiresTwoSlots()) {
                if (I) {
                    ue[0] = HSpush(data[i].raw[0]);
                    ue[1] = HSpush(data[i].raw[1]);
                    ue+=2;
                } else {
                    cte->UL = data[i].UL;
                    HSpush(cte->raw[0]);
                    HSpush(cte->raw[1]);
                    cte++;
                }
                continue;
            } else {
                if (I) {
                    *ue = HSpush(data[i].U);
                    ue++;
                } else {
                    cte->raw[0] = HSpush(data[i].U);
                    cte->raw[1] = 0;
                    // zero out the entire slot, makes comparisons easier
                    cte++;
                }
                continue;
            }
#ifndef HASH_EVERYTHING
        } else {
            if (it.requiresTwoSlots()) {
                if (I) {
                    ue[0] = data[i].raw[0];
                    ue[1] = data[i].raw[1];
                    ue+=2;
                } else {
                    cte->UL = data[i].UL;
                    cte++;
                }
                continue;
            } else {
                if (I) {
                    *ue = data[i].U;
                    ue++;
                } else {
                    cte->raw[0] = data[i].U;
                    cte->raw[1] = 0;
                    // zero out the entire slot, makes comparisons easier
                    cte++;
                }
                continue;
            }
        }
#endif

    } // for i

    if (I)  return ue;
    else    return cte;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
void MEDDLY::ct_tmpl<M,C,I>::result2entry(const ct_entry_result& res, void* e)
{
    const ct_entry_type* et = res.getET();
    MEDDLY_DCASSERT(et);

    //
    // Copy the result
    //
    unsigned* ue = (unsigned*) e;
    ct_entry_item* cte = (ct_entry_item*) e;

    const ct_entry_item* data = res.rawData();
    const unsigned datalen = res.dataLength();
    for (unsigned i=0; i<datalen; i++) {
        const ct_itemtype &it = et->getResultType(i);
        switch (it.getType())
        {
            case ct_typeID::NODE:
                    it.cacheNode(data[i].N);
                    // FALL THROUGH

            case ct_typeID::INTEGER:
            case ct_typeID::FLOAT:
                    MEDDLY_DCASSERT(sizeof(data[i].N) == sizeof(data[i].U));
                    MEDDLY_DCASSERT(sizeof(data[i].I) == sizeof(data[i].U));
                    MEDDLY_DCASSERT(sizeof(data[i].F) == sizeof(data[i].U));
                    if (I) {
                        *ue = data[i].U;
                        ue++;
                    } else {
                        cte->U = data[i].U;
                        cte++;
                    }
                    continue;

            case ct_typeID::DOUBLE:
            case ct_typeID::GENERIC:
            case ct_typeID::LONG:
                    MEDDLY_DCASSERT(sizeof(data[i].G) == sizeof(data[i].UL));
                    MEDDLY_DCASSERT(sizeof(data[i].D) == sizeof(data[i].UL));
                    if (I) {
                        ue[0] = data[i].raw[0];
                        ue[1] = data[i].raw[1];
                        ue+=2;
                    } else {
                        cte->UL = data[i].UL;
                        cte++;
                    }
                    continue;

            default:
                  MEDDLY_DCASSERT(0);
        } // switch t
    } // for i
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<M,C,I>::isDead(const ct_entry_type &ET,
        const void* sres, ct_entry_result &dres)
{
    MEDDLY_DCASSERT(sres);
    ct_entry_item* data = dres.rawData();
    MEDDLY_DCASSERT(dres.dataLength() == ET.getResultSize());

    const unsigned* ures = (const unsigned*) sres;
    const ct_entry_item* ctres = (const ct_entry_item*) sres;

    //
    // Copy the entry result from sres into dres,
    // and make sure it's not a dead entry while scanning it.
    //
    for (unsigned i=0; i<ET.getResultSize(); i++) {
        const ct_itemtype &item = ET.getResultType(i);
        if (I) {
            if (item.requiresTwoSlots()) {
                data[i].raw[0] = ures[0];
                data[i].raw[1] = ures[1];
                ures += 2;
            } else {
                data[i].raw[0] = ures[0];
                ures++;
            }
        } else {
            data[i] = ctres[i];
        }
        if (item.hasNodeType()) {
            if (item.isDeadEntry(data[i].N)) {
                return true;
            }
        }
    }
    return false;
}

#endif

// **********************************************************************

template <bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<M,C,I>::isStale(const void* entry, bool mark) const
{
    const unsigned* uptr = (const unsigned*) entry;
    const ct_entry_item* ctptr = (const ct_entry_item*) entry;

    //
    // Skip over chained portion
    //
    if (C) {
        if (I) {
            uptr += sizeof(unsigned long) / sizeof(unsigned);
        } else {
            ++ctptr;
        }
    }

    //
    // Get entry type and check if those are marked
    // If monolithic, advance entry pointer
    //
    const ct_entry_type* et = M
        ?   getEntryType( I ? *uptr : ctptr->U )
        :   global_et;
    if (M) {
        if (I) {
            ++uptr;
        } else {
            ++ctptr;
        }
    }
    MEDDLY_DCASSERT(et);
    if (et->isMarkedForDeletion()) return true;

    //
    // Get key size and advance if needed
    //
    unsigned reps;
    if (et->isRepeating()) {
        if (I) {
            reps = *uptr;
            ++uptr;
        } else {
            reps = ctptr->U;
            ++ctptr;
        }
    } else {
        reps = 0;
    }
    const unsigned klen = et->getKeySize(reps);

    //
    // Check the key portion of the entry
    //
    for (unsigned i=0; i<klen; i++) {
        const ct_itemtype &item = et->getKeyType(i);
        if (item.hasNodeType()) {
            //
            // Check the node
            //
            if (I) {
                if (item.isStaleEntry(u2n(uptr), mark)) {
                    return true;
                }
                uptr++;
            } else {
                if (item.isStaleEntry(ctptr->N, mark)) {
                    return true;
                }
                ctptr++;
            }
        } else {
            //
            // Some other type; skip
            //
            if (I) {
                uptr += item.intslots();
            } else {
                ctptr++;
            }
        }
    } // for i

    //
    // Check result portion of the entry
    //
    for (unsigned i=0; i<et->getResultSize(); i++) {
        const ct_itemtype &item = et->getResultType(i);
        if (item.hasNodeType()) {
            //
            // Check the node
            //
            if (I) {
                if (item.isStaleEntry(u2n(uptr), mark)) {
                    return true;
                }
                uptr++;
            } else {
                if (item.isStaleEntry(ctptr->N, mark)) {
                    return true;
                }
                ctptr++;
            }
        } else {
            //
            // Some other type; skip
            //
            if (I) {
                uptr += item.intslots();
            } else {
                ctptr++;
            }
        }
    } // for i

    return false;
}

// **********************************************************************

template <bool M, bool C, bool I>
void MEDDLY::ct_tmpl<M,C,I>::deleteEntry(unsigned long &h)
{
    MEDDLY_DCASSERT(h);
    const void* hptr = MMAN->getChunkAddress(h);
    const unsigned* uptr = (unsigned*) hptr;
    const ct_entry_item* ctptr = (ct_entry_item*) hptr;

    //
    // Skip over chained portion
    //
    if (C) {
        if (I) {
            uptr += sizeof(unsigned long) / sizeof(unsigned);
        } else {
            ++ctptr;
        }
    }

    //
    // Get operation and advance
    //
    const ct_entry_type* et = M
        ?   getEntryType( I ? *uptr : ctptr->U )
        :   global_et;
    if (M) {
        if (I) {
            ++uptr;
        } else {
            ++ctptr;
        }
    }

    //
    // Get key size and advance
    //
    unsigned reps;
    if (et->isRepeating()) {
        if (I) {
            reps = *uptr;
            ++uptr;
        } else {
            reps = ctptr->U;
            ++ctptr;
        }
    } else {
        reps = 0;
    }

    //
    // Go through key portion
    //
    const unsigned klen = et->getKeySize(reps);
    for (unsigned i=0; i<klen; i++) {
        const ct_itemtype &item = et->getKeyType(i);
        if (item.hasNodeType()) {
            if (I) {
                item.uncacheNode(u2n(uptr));
                uptr += item.intslots();
            } else {
                item.uncacheNode( ctptr->N );
                ++ctptr;
            }
        } else {
            if (I) {
                uptr += item.intslots();
            } else {
                ++ctptr;
            }
        }
    } // for i

    //
    // Go through result portion
    //
    for (unsigned i=0; i<et->getResultSize(); i++) {
        const ct_itemtype &item = et->getResultType(i);
        if (item.hasNodeType()) {
            if (I) {
                item.uncacheNode(u2n(uptr));
                uptr += item.intslots();
            } else {
                item.uncacheNode( ctptr->N );
                ++ctptr;
            }
        } else {
            if (I) {
                uptr += item.intslots();
            } else {
                ++ctptr;
            }
        }
    } // for i

    //
    // Recycle the actual entry
    //
    const unsigned slots = I    ? (uptr - (unsigned*) hptr)
                                : (ctptr - (ct_entry_item*) hptr);

    MMAN->recycleChunk(h, slots);
    perf.numEntries--;

    h = 0;
}

// **********************************************************************

template <bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<M,CHAINED,I>::resizeTable(unsigned long newsz)
{
#ifdef DEBUG_REHASH
    std::cout << "Resizing table, old size: " << table.size() << ", new size: " << newsz << "\n";
#endif
    //
    // Allocate the new table
    //
    std::vector <unsigned long> oldtab(newsz, 0);
    mstats.incMemUsed(oldtab.size() * sizeof(unsigned long));
    mstats.incMemAlloc(oldtab.size() * sizeof(unsigned long));
    table.swap(oldtab);

    //
    // Copy entries; requires re-hashing
    //
    if (CHAINED) {
        for (unsigned long i=0; i<oldtab.size(); i++) {
            while (oldtab[i]) {
                // Remove from old
                const unsigned long curr = oldtab[i];
                void* ptr = MMAN->getChunkAddress(curr);
                oldtab[i] = getNext(ptr);
                // Add to new
                const unsigned long h = hashEntry(ptr) % table.size();
                setNext(ptr, table[h]);
                table[h] = curr;
            }
        } // for i
    } else {
        for (unsigned long i=0; i<oldtab.size(); i++) {
            if (!oldtab[i]) continue;
            const void* ptr = MMAN->getChunkAddress(oldtab[i]);
            const unsigned long h = hashEntry(ptr) % table.size();
            setTable(h, oldtab[i]);
        } // for i
    }

    //
    // account for deletion of old table
    //
    mstats.decMemUsed(oldtab.size() * sizeof(unsigned long));
    mstats.decMemAlloc(oldtab.size() * sizeof(unsigned long));
}

// **********************************************************************

template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
unsigned MEDDLY::ct_tmpl<MONOLITHIC, CHAINED, INTSLOTS>
            ::hashEntry(const void* entry)
{
    HSreset();

    const unsigned* ue = (const unsigned*) entry;
    const ct_entry_item* cte = (const ct_entry_item*) entry;

    //
    // Skip next pointer if chained
    //
    if (CHAINED) {
        if (INTSLOTS) {
            ue += sizeof(unsigned long) / sizeof(unsigned);
        } else {
            ++cte;
        }
    }
    const ct_entry_type* et;

#ifdef HASH_EVERYTHING

    //
    // Determine entry type
    //
    unsigned eptr = 0;
    if (MONOLITHIC) {
        if (INTSLOTS) {
            et = getEntryType(ue[eptr++]);
        } else {
            et = getEntryType(cte[eptr++].U);
        }
    } else {
        et = global_et;
    }
    //
    // Determine entry size
    //
    unsigned repeats;
    if (et->isRepeating()) {
        repeats = INTSLOTS ? ue[eptr++] : cte[eptr++].U;
    } else {
        repeats = 0;
    }
    const unsigned key_slots = INTSLOTS
                                ? et->getKeyIntslots(repeats)
                                : et->getKeySize(repeats);

    //
    // Hash the whole key
    //
    if (INTSLOTS) {
        return hash_stream::raw_hash(ue, key_slots + eptr);
    } else {
        return hash_stream::raw_hash((unsigned*) cte, 2*(key_slots + eptr));
    }


    if (INTSLOTS) {
        //
        // Determine entry size, and hash away
        //
        unsigned ueptr = 0;

        return hash_stream::raw_hash(ue, key_slots + ueptr);
    }



#else

    //
    // Determine entry type
    //
    if (MONOLITHIC) {
        if (INTSLOTS) {
            et = getEntryType(HSpush(*ue));
            ++ue;
        } else {
            et = getEntryType(HSpush(cte->U));
            ++cte;
        }
    } else {
        et = global_et;
    }
    MEDDLY_DCASSERT(et);

    //
    // Determine key size
    //
    unsigned repeats;
    if (et->isRepeating()) {
        repeats = HSpush( (INTSLOTS ? (*(ue++)) : (cte++)->U) );
    } else {
        repeats = 0;
    }
    const unsigned keylen = et->getKeySize(repeats);

    //
    // Loop over key portion
    //
    for (unsigned i=0; i<keylen; i++) {
        const ct_itemtype &it = et->getKeyType(i);
        if (!it.shouldBeHashed()) {
            if (INTSLOTS) {
                ue += it.intslots();
            } else {
                ++cte;
            }
            continue;
        }
        if (it.requiresTwoSlots()) {
            if (INTSLOTS) {
                HSpush(ue[0]);
                HSpush(ue[1]);
                ue += 2;
            } else {
                HSpush(cte->raw[0]);
                HSpush(cte->raw[1]);
                ++cte;
            }
        } else {
            if (INTSLOTS) {
                HSpush(*ue);
                ++ue;
            } else {
                HSpush(cte->U);
                ++cte;
            }
        }

    } // for i

    return hash32();
#endif
}

// **********************************************************************

template <bool M, bool C, bool I>
void MEDDLY::ct_tmpl<M, C, I>
    ::showEntry(output &s, unsigned long h, bool keyonly) const
{
    MEDDLY_DCASSERT(h);
    const void* hptr = MMAN->getChunkAddress(h);
    const unsigned* uptr = (unsigned*) hptr;
    const ct_entry_item* ctptr = (ct_entry_item*) hptr;

    //
    // Skip over chained portion
    //
    if (C) {
        if (I) {
            uptr += sizeof(unsigned long) / sizeof(unsigned);
        } else {
            ++ctptr;
        }
    }

    //
    // Get operation and advance
    //
    const ct_entry_type* et = M
        ?   getEntryType( I ? *uptr : ctptr->U )
        :   global_et;
    if (M) {
        if (I) {
            ++uptr;
        } else {
            ++ctptr;
        }
    }

    //
    // Get key size and advance
    //
    unsigned reps;
    if (et->isRepeating()) {
        if (I) {
            reps = *uptr;
            ++uptr;
        } else {
            reps = ctptr->U;
            ++ctptr;
        }
    } else {
        reps = 0;
    }

    //
    // Display key portion
    //

    ct_item x;
    s << "[" << et->getName() << "(";
    const unsigned klen = et->getKeySize(reps);
    for (unsigned i=0; i<klen; i++) {
        if (i) s << ", ";

        const ct_typeID tid = et->getKeyType(i).getType();
        if (I) {
            uptr = x.set(tid, uptr);
        } else {
            x.set(tid, *ctptr);
            ++ctptr;
        }
        x.show(s);

    } // for i

    //
    // Display result portion
    //
    s << "): ";

    if (keyonly) {
        s << "temp";
    } else {
        for (unsigned i=0; i<et->getResultSize(); i++) {
            if (i) s << ", ";

            const ct_typeID tid = et->getResultType(i).getType();
            if (I) {
                uptr = x.set(tid, uptr);
            } else {
                x.set(tid, *ctptr);
                ++ctptr;
            }
            x.show(s);
        }
    }
    s << "]";
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***                                                                     ***
// ***                             Front  ends                             ***
// ***                                                                     ***
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// **********************************************************************
// *                                                                    *
// *                  monolithic_chained_style methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::monolithic_chained_style::monolithic_chained_style()
    : compute_table_style(true)
{
}

MEDDLY::compute_table*
MEDDLY::monolithic_chained_style::create(const ct_settings &s) const
{
#ifdef USE_NEW_TEMPLATE
    switch (s.compression) {
        case compressionOption::None:
                return new ct_tmpl<true, true, false>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_tmpl<true, true, true>(s, 0, 0);
        default:
                return 0;
    }
#else
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<true, true>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_typebased<true, true>(s, 0, 0);
        default:
                return 0;
    }
#endif
}

// **********************************************************************
// *                                                                    *
// *                 monolithic_unchained_style methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::monolithic_unchained_style::monolithic_unchained_style()
    : compute_table_style(true)
{
}

MEDDLY::compute_table*
MEDDLY::monolithic_unchained_style::create(const ct_settings &s) const
{
#ifdef USE_NEW_TEMPLATE
    switch (s.compression) {
        case compressionOption::None:
                return new ct_tmpl<true, false, false>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_tmpl<true, false, true>(s, 0, 0);
        default:
                return nullptr;
    }
#else
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<true, false>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_typebased<true, false>(s, 0, 0);
        default:
                return nullptr;
    }
#endif
}

// **********************************************************************
// *                                                                    *
// *                  operation_chained_style  methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::operation_chained_style::operation_chained_style()
    : compute_table_style(false)
{
}

MEDDLY::compute_table*
MEDDLY::operation_chained_style::create(const ct_settings &s, operation* op, unsigned slot) const
{
#ifdef USE_NEW_TEMPLATE
    switch (s.compression) {
        case compressionOption::None:
            return new ct_tmpl<false, true, false>(s, op, slot);
        case compressionOption::TypeBased:
            return new ct_tmpl<false, true, true>(s, op, slot);
        default:
            return nullptr;
    }
#else
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<false, true>(s, op, slot);
        case compressionOption::TypeBased:
                return new ct_typebased<false, true>(s, op, slot);
        default:
                return nullptr;
    }
#endif
}

// **********************************************************************
// *                                                                    *
// *                 operation_unchained_style  methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::operation_unchained_style::operation_unchained_style()
    : compute_table_style(false)
{
}

MEDDLY::compute_table*
MEDDLY::operation_unchained_style::create(const ct_settings &s, operation* op, unsigned slot) const
{
#ifdef USE_NEW_TEMPLATE
    switch (s.compression) {
        case compressionOption::None:
                return new ct_tmpl<false, false, false>(s, op, slot);
        case compressionOption::TypeBased:
                return new ct_tmpl<false, false, false>(s, op, slot);
        default:
                return nullptr;
    }
#else
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<false, false>(s, op, slot);
        case compressionOption::TypeBased:
                return new ct_typebased<false, false>(s, op, slot);
        default:
                return nullptr;
    }
#endif
}

