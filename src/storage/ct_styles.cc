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
#include "../ct_generics.h"
#include "../ct_initializer.h"
#include "../operators.h"

// #define DEBUG_FIND
// #define DEBUG_REHASH
// #define DEBUG_CT

// #define HASHARRAY_SIMPLE

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
            @param  TTYPE           Type for the hash table entries.
                                    In practice, either 'unsigned'
                                    (limit of 4 billion entries)
                                    or 'unsigned long'.

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
            |         next pointer  (TTYPE)         |
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
    template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
    class ct_tmpl : public compute_table {
        public:
            ct_tmpl(const ct_settings &s, unsigned etid);
            virtual ~ct_tmpl();

            // Required functions

#ifdef ALLOW_DEPRECATED_0_17_6
            virtual void find(ct_entry_key* key, ct_entry_result &res);
            virtual void addEntry(ct_entry_key* key,
                const ct_entry_result& res);
            virtual void doneKey(ct_entry_key* key);
#endif

            virtual bool find(const ct_entry_type &ET, ct_vector &key,
                    ct_vector &res);

            virtual void addEntry(const ct_entry_type &ET, ct_vector &key,
                    const ct_vector &res);

            virtual void doneKey(ct_vector &key);

            // calls removeStaleEntries, then maybe shrinks table
            virtual void removeStales();
            virtual void removeAll();
            virtual void show(output &s, int verbLevel = 0);
            virtual void countNodeEntries(const forest* f,
                    std::vector <unsigned long> &counts) const;

        protected:
            // Helper functions

            void removeStaleEntries();

#ifdef ALLOW_DEPRECATED_0_17_6
            void* key2entry(const ct_entry_key& key, void* e);
            static void result2entry(const ct_entry_result& res, void* e);
#endif

            void* key2entry(const ct_entry_type &ET, const ct_vector &key,
                    void* e);
            static void result2entry(const ct_entry_type &ET,
                    const ct_vector& res, void* e);

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

            inline static TTYPE getNext(const void* raw)
            {
                ASSERT(__FILE__, __LINE__, CHAINED);

                union {
                    TTYPE tt;
                    unsigned u[2];
                    unsigned long ul;
                } thingy;

                if (INTSLOTS) {
                    const unsigned* u = (const unsigned*) raw;
                    thingy.u[0] = u[0];
                    if (sizeof(TTYPE) == sizeof(long)) {
                        thingy.u[1] = u[1];
                    }
                    return thingy.tt;
                } else {
                    const ct_entry_item* ct = (const ct_entry_item*) raw;
                    if (sizeof(TTYPE) == sizeof(unsigned)) {
                        thingy.u[0] = ct[0].U;
                        return thingy.tt;
                    }
                    thingy.ul = ct[0].UL;
                    return thingy.tt;
                }
            }
            inline static void setNext(void* raw, TTYPE n)
            {
                ASSERT(__FILE__, __LINE__, CHAINED);

                union {
                    TTYPE tt;
                    unsigned u[2];
                    unsigned long ul;
                } thingy;

                thingy.tt = n;

                if (INTSLOTS) {
                    unsigned* u = (unsigned*) raw;
                    u[0] = thingy.u[0];
                    if (sizeof(TTYPE) == sizeof(long)) {
                        u[1] = thingy.u[1];
                    }
                } else {
                    ct_entry_item* ct = (ct_entry_item*) raw;
                    if (sizeof(TTYPE) == sizeof(unsigned)) {
                        ct[0].U = thingy.u[0];
                    } else {
                        ct[0].UL = thingy.ul;
                    }
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

            /**
                Copy the result portion of an entry, and check if it is dead.
                    @param  ET      Entry type information
                    @param  sres    Pointer to result portion of entry
                    @param  dres    Where to copy the result
             */
            static bool isDead(const ct_entry_type &ET, const void* sres,
                    ct_vector &dres);


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
                Try to set table[h] to curr.
                If the slot is full, check ahead the next couple for empty.
                If none are empty, then recycle current table[h]
                and overwrite it.
            */
            inline void setTable(TTYPE h, TTYPE curr)
            {
                ASSERT(__FILE__, __LINE__, !CHAINED);
                TTYPE hfree = h;
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
                ASSERT(__FILE__, __LINE__, table[h]);
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
            void deleteEntry(TTYPE &hnd);

            /**
                Resize the table.
                    @param  newsz   New size to use
            */
            void resizeTable(TTYPE newsz);

            /**
                Compute the hash of an entry.
                    @param  entry   Raw pointer to entire entry

                    @return     Hashed value; still needs to be modded
                                by hash table size.

                TBD - switch to 64-bit hash
             */
            TTYPE hashEntry(const void* entry);


            /**
                Display an entry.
                Used for debugging, and by method(s) to display the entire CT.
                    @param  s         Stream to write to.
                    @param  h         Handle of the entry.
                    @param  keyonly     true: we don't have a result yet
                                        false: it's a proper entry with both
            */
            void showEntry(output &s, TTYPE h, bool keyonly) const;

            /// Display a chain
            inline void showChain(output &s, TTYPE L) const
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

            /// Multiply and check for overflow
            static inline TTYPE safemult(TTYPE sz, TTYPE a)
            {
                a *= sz;
                if (a < sz) {
                    // we overflowed
                    return 0;
                }
                return a;
            }


        private:
            //
            // Abstraction layer for hashing array
            //
            inline void HSreset() {
#ifdef HASHARRAY_SIMPLE
                HS.resize(0);
#else
                hsi = 0;
#endif
            }
            inline void HSreserve(unsigned sz) {
#ifdef HASHARRAY_SIMPLE
                HS.reserve(sz);
                HS.resize(0);
#else
                if (sz > HS.size()) {
                    HS.resize(sz);
                }
                hsi = 0;
#endif
            }
            inline unsigned HSpush(unsigned x) {
#ifdef HASHARRAY_SIMPLE
                HS.push_back(x);
                return x;
#else
#ifdef DEVELOPMENT_CODE
                return HS.at(hsi++) = x;
#else
                return HS[hsi++] = x;
#endif
#endif
            }
            inline unsigned hash32() {
#ifdef HASHARRAY_SIMPLE
                return hash_stream::raw_hash(HS.data(), HS.size());
#else
                return hash_stream::raw_hash(HS.data(), hsi);
#endif
            }

        private:
            /// Hashing array
            std::vector <unsigned> HS;

#ifndef HASHARRAY_SIMPLE
            /// Current index in hashing array
            unsigned hsi;
#endif

            /// List of entries to delete
            std::vector <TTYPE> toDelete;

            /// Hash table
            std::vector <TTYPE> table;

            /// When to next expand the table
            TTYPE tableExpand;

            /// When to next shrink the table
            TTYPE tableShrink;

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

template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
MEDDLY::ct_tmpl<TTYPE, MONOLITHIC, CHAINED, INTSLOTS>::ct_tmpl(
    const ct_settings &s, unsigned etid) : compute_table(s, etid)
{
    if (MONOLITHIC) {
        ASSERT(__FILE__, __LINE__, 0==etid);
    } else {
        ASSERT(__FILE__, __LINE__, etid);
    }

    /*
        Initialize memory management for entries.
    */
    ASSERT(__FILE__, __LINE__, s.MMS);
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

    mstats.incMemUsed(table.size() * sizeof(TTYPE));
    mstats.incMemAlloc(table.size() * sizeof(TTYPE));

    collisions = 0;
}

// **********************************************************************

template <class TT, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
MEDDLY::ct_tmpl<TT, MONOLITHIC, CHAINED, INTSLOTS>::~ct_tmpl()
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

template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<TTYPE,MONOLITHIC,CHAINED,INTSLOTS>
    ::find(ct_entry_key* key, ct_entry_result &res)
{
    ASSERT(__FILE__, __LINE__, key);
    const ct_entry_type* et = key->getET();
    ASSERT(__FILE__, __LINE__, et);
    if (!MONOLITHIC) {
        if (et->getID() != global_etid)
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

    const unsigned chain_slots = CHAINED
            ? ( INTSLOTS ? sizeof(TTYPE)/sizeof(unsigned) : 1 )
            : 0;

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
    ASSERT(__FILE__, __LINE__, key->my_entry);

    genentry keyentry, keyafterchain;
    keyentry.vptr = MMAN->getChunkAddress(key->my_entry);

    //
    // (2) build the entry, except for the result.
    //     If the entire key is hashable, then we
    //     hash everything all at once at the end;
    //     otherwise, we hash it as we go.
    //

    HSreserve(num_slots);
    if (INTSLOTS) {
        unsigned* e = keyentry.uptr;
        if (CHAINED) {
            // skip over NEXT pointer slot(s)
            if (sizeof(TTYPE) == sizeof(unsigned)) {
                e++;
            } else {
                e += 2;
            }
        }
        keyafterchain.uptr = e;
        if (et->isEntireKeyHashable(key->numRepeats())) {
            //
            // Copy operation and key size info
            //
            if (MONOLITHIC) {
                *e = et->getID();
                ++e;
            }
            if (et->isRepeating()) {
                *e = key->numRepeats();
                ++e;
            }
        } else {
            //
            // Copy and hash operation and key size info
            //
            if (MONOLITHIC) {
                *e = HSpush(et->getID());
                ++e;
            }
            if (et->isRepeating()) {
                *e = HSpush(key->numRepeats());
                ++e;
            }
        }
        key->result_shift = (unsigned*) key2entry(*key, e) - keyentry.uptr;
    } else {
        ct_entry_item* e = keyentry.ctptr;
        if (CHAINED) {
            // skip over NEXT pointer
            ++e;
        }
        keyafterchain.ctptr = e;
        if (et->isEntireKeyHashable(key->numRepeats())) {
            //
            // Copy operation and key size info
            //
            if (MONOLITHIC) {
                e->raw[0] = et->getID();
                e->raw[1] = 0;
                ++e;
            }
            if (et->isRepeating()) {
                e->raw[0] = key->numRepeats();
                e->raw[1] = 0;
                ++e;
            }
        } else {
            //
            // Copy and hash operation and key size info
            //
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
        }
        key->result_shift = (ct_entry_item*) key2entry(*key, e) - keyentry.ctptr;
    }

    //
    // (3) Save hash value
    //

    if (et->isEntireKeyHashable(key->numRepeats())) {
        //
        // Hash directly from the (entire) Key
        //
        if (INTSLOTS) {
            key->setHash(
                hash_stream::raw_hash(keyafterchain.uptr, entry_slots)
            );
        } else {
            key->setHash(
                hash_stream::raw_hash(keyafterchain.uptr, 2*entry_slots)
            );
        }
    } else {
        //
        // Hash from the array we have set aside
        //
        key->setHash(hash32());
    }
    ASSERT(__FILE__, __LINE__, key->getHash() == hashEntry(keyentry.vptr));
    const TTYPE hslot = key->getHash() % table.size();
    TTYPE hcurr = hslot;
    // TBD: use a 64-bit hash when TTYPE is 64 bit

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
    TTYPE currhnd = table.at(hcurr);
#else
    TTYPE currhnd = table[hcurr];
#endif
    TTYPE prevhnd = 0;
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
                TTYPE next = getNext(currentry.vptr);
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


template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<TTYPE, MONOLITHIC, CHAINED, INTSLOTS>
    ::addEntry(ct_entry_key* key, const ct_entry_result& res)
{
    ASSERT(__FILE__, __LINE__, key);
    const ct_entry_type* et = key->getET();
    ASSERT(__FILE__, __LINE__, et);
    if (!MONOLITHIC) {
        if (et->getID() != global_etid)
            throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }
    ct_entry_type::incEntries(et->getID());

    TTYPE h = key->getHash() % table.size();

    //
    // Increment cache counters for nodes in the key
    //
    key->cacheNodes();

    //
    // Copy the result portion, and increment cache counters in the result
    //
    ASSERT(__FILE__, __LINE__, key->my_entry);
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
    std::vector <bool> skipF(NF, false);
    clearForestCTBits(skipF);
    removeStaleEntries();
    sweepForestCTBits(skipF);

    //
    // Do we still need to enlarge the table?
    //
    if (CHAINED) {
        if (perf.numEntries < table.size()) return;
    } else {
        if (perf.numEntries < tableExpand / 4) return;
    }
    TTYPE newsize = safemult(table.size(), 2);
    if (!newsize) newsize = maxSize;
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
            tableExpand = std::numeric_limits<TTYPE>::max();
        } else {
            tableExpand = safemult(table.size(), 4);
            if (!tableExpand) tableExpand = std::numeric_limits<TTYPE>::max();
        }
        tableShrink = table.size() / 2;
    } else {
        if (table.size() == maxSize) {
            tableExpand = std::numeric_limits<TTYPE>::max();
        } else {
            tableExpand = table.size() / 2;
        }
        tableShrink = table.size() / 8;
    }
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <class T, bool M, bool C, bool I>
void MEDDLY::ct_tmpl<T,M,C,I>::doneKey(ct_entry_key* key)
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

template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
bool MEDDLY::ct_tmpl<TTYPE,MONOLITHIC,CHAINED,INTSLOTS>
    ::find(const ct_entry_type &ET, ct_vector &key, ct_vector &res)
{
    //
    // Go ahead and start building an entry (the key portion anyway)
    // as this makes searching for matches easier.
    //

    if (!MONOLITHIC) {
        if (ET.getID() != global_etid)
            throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }

    //
    // (1) determine entry size, and allocate
    //

    const unsigned chain_slots = CHAINED
            ? ( INTSLOTS ? sizeof(TTYPE)/sizeof(unsigned) : 1 )
            : 0;

    const unsigned op_slots = MONOLITHIC ? 1 : 0;

    ASSERT(__FILE__, __LINE__, ET.canHaveKeySize(key.size()));
    const unsigned repeats = ET.repeatsForKeySize(key.size());

    const unsigned key_slots =
        INTSLOTS  ? ET.getKeyIntslots(repeats)
                  : ET.getKeySize(repeats);

    const unsigned res_slots =
        INTSLOTS  ? ET.getResultIntslots()
                  : ET.getResultSize();

    const unsigned entry_slots =
          op_slots
        + (ET.isRepeating() ? 1 : 0)
        + key_slots
    ;

    size_t num_slots =
          chain_slots
        + entry_slots
        + res_slots
    ;

    key.entry_slots = num_slots;
    key.my_entry = MMAN->requestChunk(num_slots);
    if (num_slots < key.entry_slots) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    ASSERT(__FILE__, __LINE__, key.my_entry);

    const void* entry_vptr = MMAN->getChunkAddress(key.my_entry);
    unsigned* entry_after_chain = nullptr;

    //
    // (2) build the entry, except for the result.
    //     If the entire key is hashable, then we
    //     hash everything all at once at the end;
    //     otherwise, we hash it as we go.
    //

    HSreserve(num_slots);
    if (INTSLOTS) {
        unsigned* const entry_uptr = (unsigned*) entry_vptr;
        unsigned* e = entry_uptr;
        if (CHAINED) {
            // skip over NEXT pointer slot(s)
            if (sizeof(TTYPE) == sizeof(unsigned)) {
                e++;
            } else {
                e += 2;
            }
        }
        entry_after_chain = e;
        if (ET.isEntireKeyHashable(repeats)) {
            //
            // Copy operation and key size info
            //
            if (MONOLITHIC) {
                *e = ET.getID();
                ++e;
            }
            if (ET.isRepeating()) {
                *e = repeats;
                ++e;
            }
        } else {
            //
            // Copy and hash operation and key size info
            //
            if (MONOLITHIC) {
                *e = HSpush(ET.getID());
                ++e;
            }
            if (ET.isRepeating()) {
                *e = HSpush(repeats);
                ++e;
            }
        }
        key.result_shift = (unsigned*) key2entry(ET, key, e) - entry_uptr;
    } else {
        ct_entry_item* const entry_ctptr = (ct_entry_item*) entry_vptr;
        ct_entry_item* e = entry_ctptr;
        if (CHAINED) {
            // skip over NEXT pointer
            ++e;
        }
        entry_after_chain = (unsigned*) e;

        if (ET.isEntireKeyHashable(repeats)) {
            //
            // Copy operation and key size info
            //
            if (MONOLITHIC) {
                e->raw[0] = ET.getID();
                e->raw[1] = 0;
                ++e;
            }
            if (ET.isRepeating()) {
                e->raw[0] = repeats;
                e->raw[1] = 0;
                ++e;
            }
        } else {
            //
            // Copy and hash operation and key size info
            //
            if (MONOLITHIC) {
                e->raw[0] = HSpush(ET.getID());
                e->raw[1] = 0;
                ++e;
            }
            if (ET.isRepeating()) {
                e->raw[0] = HSpush(repeats);
                e->raw[1] = 0;
                ++e;
            }
        }
        key.result_shift = (ct_entry_item*) key2entry(ET, key, e) - entry_ctptr;
    }

    //
    // (3) Save hash value
    //

    if (ET.isEntireKeyHashable(repeats)) {
        //
        // Hash directly from the (entire) Key
        //
        if (INTSLOTS) {
            key.setHash(
                hash_stream::raw_hash(entry_after_chain, entry_slots)
            );
        } else {
            key.setHash(
                hash_stream::raw_hash(entry_after_chain, 2*entry_slots)
            );
        }
    } else {
        //
        // Hash from the array we have set aside
        //
        key.setHash(hash32());
    }
    ASSERT(__FILE__, __LINE__, key.getHash() == hashEntry(entry_vptr));
    const TTYPE hslot = key.getHash() % table.size();
    TTYPE hcurr = hslot;
    // TBD: use a 64-bit hash when TTYPE is 64 bit

#ifdef DEBUG_FIND
    ostream_output out(std::cout);
    out << "Searching for entry ";
    showEntry(out, key.my_entry, true);
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
    TTYPE currhnd = table.at(hcurr);
#else
    TTYPE currhnd = table[hcurr];
#endif
    TTYPE prevhnd = 0;
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
        void* curr_vptr = MMAN->getChunkAddress(currhnd);
        void* curr_after;
        void* curr_res;
        curr_vptr = MMAN->getChunkAddress(currhnd);
        if (INTSLOTS) {
            curr_after = ((unsigned*) curr_vptr) + chain_slots;
            curr_res = ((unsigned*) curr_after) + entry_slots;
        } else {
            curr_after = ((ct_entry_item*) curr_vptr) + chain_slots;
            curr_res = ((ct_entry_item*) curr_after) + entry_slots;
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
            if (equal_sw((unsigned*) curr_after, entry_after_chain, entry_slots))
            {
                equal = true;
                remove = isDead(ET, curr_res, res);
            } else {
                equal = false;
                if (checkStalesOnFind) {
                    remove = isStale(curr_vptr, false);
                } else {
                    remove = false;
                }
            }
        } else {
            if (equal_sw((ct_entry_item*) curr_after, (ct_entry_item*) entry_after_chain, entry_slots))
            {
                equal = true;
                remove = isDead(ET, curr_res, res);
            } else {
                equal = false;
                if (checkStalesOnFind) {
                    remove = isStale(curr_vptr, false);
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
                TTYPE next = getNext(curr_vptr);
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
                    setNext(prev, getNext(curr_vptr));
                    setNext(curr_vptr, table[hcurr]);
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
            // Clean up
            //
            sawSearch(chainlen);
            perf.hits++;
            batchDelete();
            MMAN->recycleChunk(key.my_entry, key.entry_slots);
            key.entry_slots = 0;
            return true;
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
            currhnd = getNext(curr_vptr);
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
    batchDelete();
#ifdef DEBUG_FIND
    out << "        not found\n";
#endif
    return false;
}

// **********************************************************************

template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<TTYPE,MONOLITHIC,CHAINED,INTSLOTS>
    ::addEntry(const ct_entry_type &ET, ct_vector &key, const ct_vector &res)
{
    if (!MONOLITHIC) {
        if (ET.getID() != global_etid)
            throw error(error::UNKNOWN_OPERATION, __FILE__, __LINE__);
    }
    ASSERT(__FILE__, __LINE__, res.size() == ET.getResultSize());

    ct_entry_type::incEntries(ET.getID());

    TTYPE h = key.getHash() % table.size();

    //
    // Increment cache counters for nodes in the key
    //
    for (unsigned i=0; i<key.size(); i++) {
        const ct_itemtype &item = ET.getKeyType(i);
        if (item.hasNodeType()) {
            item.cacheNode(key[i].getN());
        }
    }

    //
    // Copy the result portion, and increment cache counters in the result
    //
    ASSERT(__FILE__, __LINE__, key.my_entry);
    void* rawentry = MMAN->getChunkAddress(key.my_entry);
    if (INTSLOTS) {
        unsigned* entry = (unsigned*) rawentry;
        result2entry(ET, res, entry + key.result_shift);
    } else {
        ct_entry_item* entry = (ct_entry_item*) rawentry;
        result2entry(ET, res, entry + key.result_shift);
    }

    //
    // Add to the CT
    //
    if (CHAINED) {
        // Add this to front of chain
        setNext(rawentry, table[h]);
        table[h] = key.my_entry;
    } else {
        setTable(h, key.my_entry);
    }
    key.my_entry = 0;

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
    std::vector <bool> skipF(NF, false);
    clearForestCTBits(skipF);
    removeStaleEntries();
    sweepForestCTBits(skipF);

    //
    // Do we still need to enlarge the table?
    //
    if (CHAINED) {
        if (perf.numEntries < table.size()) return;
    } else {
        if (perf.numEntries < tableExpand / 4) return;
    }
    TTYPE newsize = safemult(table.size(), 2);
    if (!newsize) newsize = maxSize;
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
            tableExpand = std::numeric_limits<TTYPE>::max();
        } else {
            tableExpand = safemult(table.size(), 4);
            if (!tableExpand) tableExpand = std::numeric_limits<TTYPE>::max();
        }
        tableShrink = table.size() / 2;
    } else {
        if (table.size() == maxSize) {
            tableExpand = std::numeric_limits<TTYPE>::max();
        } else {
            tableExpand = table.size() / 2;
        }
        tableShrink = table.size() / 8;
    }
}

// **********************************************************************

template <class T, bool M, bool C, bool I>
void MEDDLY::ct_tmpl<T,M,C,I>::doneKey(ct_vector &key)
{
    if (key.my_entry) {
        MMAN->recycleChunk(key.my_entry, key.entry_slots);
        key.my_entry = 0;
    }
}

// **********************************************************************

template <class TTYPE, bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<TTYPE, M, CHAINED, I>::removeStales()
{
    removeStaleEntries();

    //
    // Is it time to shrink the table?
    //
    if (perf.numEntries >= tableShrink) return;
    TTYPE newsize = table.size() / 2;
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

template <class TTYPE, bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<TTYPE,M,CHAINED,I>::removeAll()
{
    for (TTYPE i=0; i<table.size(); i++) {
        while (table[i]) {
            TTYPE curr = table[i];
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

template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<TTYPE, MONOLITHIC, CHAINED, INTSLOTS>
    ::show(output &s, int verbLevel)
{
    if (verbLevel < 1) return;

    const char* chained = CHAINED ? "(chained)" : "(unchained)";

    if (MONOLITHIC) {
        s << "Monolithic " << chained << " compute table\n";
    } else {
        const ct_entry_type* et = ct_entry_type::getEntryType(global_etid);
        if (et) {
            s << "Compute table " << chained << " for " << et->getName()
              << " (index " << long(global_etid) << ")\n";
        } else {
            s << "Compute table " << chained << " for null entry"
              << " (index " << long(global_etid) << ")\n";
        }
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

    for (TTYPE i=0; i<table.size(); i++) {
        if (0==table[i]) continue;
        s << "table[";
        s.put(long(i), 9);
        s << "]: ";
        showChain(s, table[i]);
    }

    if (--verbLevel < 1) return;

    s << "\nHash table nodes:\n";

    for (TTYPE i=0; i<table.size(); i++) {
        TTYPE curr = table[i];
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

template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
void MEDDLY::ct_tmpl<TTYPE, MONOLITHIC, CHAINED, INTSLOTS>
    ::countNodeEntries(const forest* f, std::vector <unsigned long> &counts)
    const
{
    for (TTYPE h=0; h<table.size(); h++) {
        TTYPE curr = table[h];
        while (curr) {
            const void* entry = MMAN->getChunkAddress(curr);
            const unsigned* ue = (const unsigned*) entry;
            const ct_entry_item* cte = (const ct_entry_item*) entry;

            //
            // Skip over chained portion
            //
            if (CHAINED) {
                if (INTSLOTS) {
                    if (sizeof(TTYPE) == sizeof(unsigned)) {
                        ++ue;
                    } else {
                        ue += 2;
                    }
                } else {
                    ++cte;
                }
            }

            //
            // Get operation and advance
            //
            const ct_entry_type* et = ct_entry_type::getEntryType(
                MONOLITHIC
                    ?   (INTSLOTS ? *ue : cte->U)
                    :   global_etid
            );

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

template <class TTYPE, bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<TTYPE, M, CHAINED, I>::removeStaleEntries()
{
    if (CHAINED) {
        // Chained
        ostream_output out(std::cout);
        for (TTYPE i=0; i<table.size(); i++) {
            if (!table[i]) continue;

            void* preventry = nullptr;
            TTYPE curr;
            for (curr = table[i]; curr; ) {
                void* entry = MMAN->getChunkAddress(curr);
                TTYPE next = getNext(entry);

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
        for (TTYPE i=0; i<table.size(); i++) {
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

template <class T, bool M, bool C, bool I>
void* MEDDLY::ct_tmpl<T,M,C,I>::key2entry(const ct_entry_key& key, void* e)
{
    const ct_entry_type* et = key.getET();
    ASSERT(__FILE__, __LINE__, et);
    const unsigned keylen = et->getKeySize(key.numRepeats());
    const ct_entry_item* data = key.rawData();

    unsigned* ue = (unsigned*) e;
    ct_entry_item* cte = (ct_entry_item*) e;


    if (et->isEntireKeyHashable(key.numRepeats())) {
        //
        // Just copy the key over; no hashing
        // (caller will hash everything at once)
        //
        for (unsigned i=0; i<keylen; i++) {
            const ct_itemtype &it = et->getKeyType(i);
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
        } // for i
    } else {
        //
        // Copy the key over and hash as we go
        //
        for (unsigned i=0; i<keylen; i++) {
            const ct_itemtype &it = et->getKeyType(i);

            if (it.shouldBeHashed()) {
                // copy and hash
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
            } else {
                // just copy
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
        } // for i
    } // hash everything or not?

    if (I)  return ue;
    else    return cte;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <class T, bool M, bool C, bool I>
void MEDDLY::ct_tmpl<T,M,C,I>
    ::result2entry(const ct_entry_result& res, void* e)
{
    const ct_entry_type* et = res.getET();
    ASSERT(__FILE__, __LINE__, et);

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
                    ASSERT(__FILE__, __LINE__, sizeof(data[i].N) == sizeof(data[i].U));
                    ASSERT(__FILE__, __LINE__, sizeof(data[i].I) == sizeof(data[i].U));
                    ASSERT(__FILE__, __LINE__, sizeof(data[i].F) == sizeof(data[i].U));
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
                    ASSERT(__FILE__, __LINE__, sizeof(data[i].G) == sizeof(data[i].UL));
                    ASSERT(__FILE__, __LINE__, sizeof(data[i].D) == sizeof(data[i].UL));
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
                  FAIL(__FILE__, __LINE__, "Unknown ct type");
        } // switch t
    } // for i
}

#endif

// **********************************************************************

template <class T, bool M, bool C, bool I>
void* MEDDLY::ct_tmpl<T,M,C,I>::key2entry(const ct_entry_type &ET,
        const ct_vector& key, void* e)
{
    unsigned* ue = (unsigned*) e;
    ct_entry_item* cte = (ct_entry_item*) e;

    const unsigned repeats = ET.repeatsForKeySize(key.size());

    if (ET.isEntireKeyHashable(repeats)) {
        //
        // Just copy the key over; no hashing
        // (caller will hash everything at once)
        //
        for (unsigned i=0; i<key.size(); i++) {
            const ct_itemtype &it = ET.getKeyType(i);
            ASSERT(__FILE__, __LINE__, it.hasType(key[i].getType()));

            if (it.requiresTwoSlots()) {
                if (I) {
                    ue[0] = key[i].raw0();
                    ue[1] = key[i].raw1();
                    ue+=2;
                } else {
                    cte->UL = key[i].rawUL();
                    cte++;
                }
                continue;
            } else {
                if (I) {
                    *ue = key[i].raw0();
                    ue++;
                } else {
                    cte->raw[0] = key[i].raw0();
                    cte->raw[1] = 0;
                    // zero out the entire slot, makes comparisons easier
                    cte++;
                }
                continue;
            }
        } // for i
    } else {
        //
        // Copy the key over and hash as we go
        //
        for (unsigned i=0; i<key.size(); i++) {
            const ct_itemtype &it = ET.getKeyType(i);
            ASSERT(__FILE__, __LINE__, it.hasType(key[i].getType()));

            if (it.shouldBeHashed()) {
                // copy and hash
                if (it.requiresTwoSlots()) {
                    if (I) {
                        ue[0] = HSpush(key[i].raw0());
                        ue[1] = HSpush(key[i].raw1());
                        ue+=2;
                    } else {
                        cte->UL = key[i].rawUL();
                        HSpush(cte->raw[0]);
                        HSpush(cte->raw[1]);
                        cte++;
                    }
                    continue;
                } else {
                    if (I) {
                        *ue = HSpush(key[i].raw0());
                        ue++;
                    } else {
                        cte->raw[0] = HSpush(key[i].raw0());
                        cte->raw[1] = 0;
                        // zero out the entire slot, makes comparisons easier
                        cte++;
                    }
                    continue;
                }
            } else {
                // just copy
                if (it.requiresTwoSlots()) {
                    if (I) {
                        ue[0] = key[i].raw0();
                        ue[1] = key[i].raw1();
                        ue+=2;
                    } else {
                        cte->UL = key[i].rawUL();
                        cte++;
                    }
                    continue;
                } else {
                    if (I) {
                        *ue = key[i].raw0();
                        ue++;
                    } else {
                        cte->raw[0] = key[i].raw0();
                        cte->raw[1] = 0;
                        // zero out the entire slot, makes comparisons easier
                        cte++;
                    }
                    continue;
                }
            }
        } // for i
    }

    if (I)  return ue;
    else    return cte;
}

// **********************************************************************

template <class T, bool M, bool C, bool I>
void MEDDLY::ct_tmpl<T,M,C,I>
    ::result2entry(const ct_entry_type &ET, const ct_vector& res, void* e)
{
    //
    // Copy the result and increment cache counts for any nodes
    //
    unsigned* ue = (unsigned*) e;
    ct_entry_item* cte = (ct_entry_item*) e;

    for (unsigned i=0; i<res.size(); i++) {
        const ct_itemtype &it = ET.getResultType(i);
        ASSERT(__FILE__, __LINE__, it.hasType(res[i].getType()));
        if (it.hasNodeType()) {
            ASSERT(__FILE__, __LINE__, !it.requiresTwoSlots());
            it.cacheNode(res[i].getN());
            if (I) {
                *ue = res[i].raw0();
                ++ue;
            } else {
                cte->U = res[i].raw0();
                ++cte;
            }
            continue;
        }

        if (it.requiresTwoSlots()) {
            if (I) {
                ue[0] = res[i].raw0();
                ue[1] = res[i].raw1();
                ue+=2;
            } else {
                cte->UL = res[i].rawUL();
                cte++;
            }
        } else {
            if (I) {
                *ue = res[i].raw0();
                ++ue;
            } else {
                cte->U = res[i].raw0();
                ++cte;
            }
        }
    } // for i
}

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <class T, bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<T,M,C,I>::isDead(const ct_entry_type &ET,
        const void* sres, ct_entry_result &dres)
{
    ASSERT(__FILE__, __LINE__, sres);
    ct_entry_item* data = dres.rawData();
    ASSERT(__FILE__, __LINE__, dres.dataLength() == ET.getResultSize());

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

template <class T, bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<T,M,C,I>::isDead(const ct_entry_type &ET,
        const void* sres, ct_vector &dres)
{
    ASSERT(__FILE__, __LINE__, sres);
    ASSERT(__FILE__, __LINE__, dres.size() == ET.getResultSize());

    const unsigned* ures = (const unsigned*) sres;
    const ct_entry_item* ctres = (const ct_entry_item*) sres;

    //
    // Copy the entry result from sres into dres,
    // and make sure it's not a dead entry while scanning it.
    //
    for (unsigned i=0; i<ET.getResultSize(); i++) {
        const ct_itemtype &item = ET.getResultType(i);
        if (I) {
            dres[i].setType(item.getType());
            if (item.requiresTwoSlots()) {
                dres[i].raw0() = ures[0];
                dres[i].raw1() = ures[1];
                ures += 2;
            } else {
                dres[i].raw0() = ures[0];
                ures++;
            }
        } else {
            dres[i].set(item.getType(), ctres[i]);
        }
        if (item.hasNodeType()) {
            if (item.isDeadEntry(dres[i].getN())) {
                return true;
            }
        }
    }
    return false;
}

// **********************************************************************

template <class T, bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<T,M,C,I>::isStale(const void* entry, bool mark) const
{
    const unsigned* uptr = (const unsigned*) entry;
    const ct_entry_item* ctptr = (const ct_entry_item*) entry;

    //
    // Skip over next pointer
    //
    if (C) {
        if (I) {
            if (sizeof(T) == sizeof(unsigned)) {
                ++uptr;
            } else {
                uptr += 2;
            }
        } else {
            ++ctptr;
        }
    }

    //
    // Get entry type and check if those are marked
    // If monolithic, advance entry pointer
    //
    const ct_entry_type* et = ct_entry_type::getEntryType(
        M   ?   ( I ? *uptr : ctptr->U )
            :   global_etid
    );

    if (M) {
        if (I) {
            ++uptr;
        } else {
            ++ctptr;
        }
    }
    ASSERT(__FILE__, __LINE__, et);
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

template <class TTYPE, bool M, bool C, bool I>
void MEDDLY::ct_tmpl<TTYPE,M,C,I>::deleteEntry(TTYPE &h)
{
    ASSERT(__FILE__, __LINE__, h);
    const void* hptr = MMAN->getChunkAddress(h);
    const unsigned* uptr = (unsigned*) hptr;
    const ct_entry_item* ctptr = (ct_entry_item*) hptr;

    //
    // Skip over next pointer
    //
    if (C) {
        if (I) {
            if (sizeof(TTYPE) == sizeof(unsigned)) {
                ++uptr;
            } else {
                uptr += 2;
            }
        } else {
            ++ctptr;
        }
    }

    //
    // Get operation and advance
    //
    const ct_entry_type* et = ct_entry_type::getEntryType(
         M  ?   ( I ? *uptr : ctptr->U )
            :   global_etid
    );
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
            continue;
        }
        if (item.hasType(ct_typeID::GENERIC)) {
            if (I) {
                ct_item x;
                x.set(ct_typeID::GENERIC, uptr);
                delete x.getG();
            } else {
                delete ctptr->G;
            }
        }
        if (I) {
            uptr += item.intslots();
        } else {
            ++ctptr;
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
            continue;
        }
        if (item.hasType(ct_typeID::GENERIC)) {
            if (I) {
                ct_item x;
                x.set(ct_typeID::GENERIC, uptr);
                delete x.getG();
            } else {
                delete ctptr->G;
            }
        }
        if (I) {
            uptr += item.intslots();
        } else {
            ++ctptr;
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

    //
    // Decrement entry counter in entry type.
    //
    ct_entry_type::decEntries(et->getID());
}

// **********************************************************************

template <class TTYPE, bool M, bool CHAINED, bool I>
void MEDDLY::ct_tmpl<TTYPE,M,CHAINED,I>::resizeTable(TTYPE newsz)
{
#ifdef DEBUG_REHASH
    std::cout << "Resizing table, old size: " << table.size() << ", new size: " << newsz << "\n";
#endif
    //
    // Allocate the new table
    //
    std::vector <TTYPE> oldtab(newsz, 0);
    mstats.incMemUsed(oldtab.size() * sizeof(TTYPE));
    mstats.incMemAlloc(oldtab.size() * sizeof(TTYPE));
    table.swap(oldtab);

    //
    // Copy entries; requires re-hashing
    //
    if (CHAINED) {
        for (TTYPE i=0; i<oldtab.size(); i++) {
            while (oldtab[i]) {
                // Remove from old
                const TTYPE curr = oldtab[i];
                void* ptr = MMAN->getChunkAddress(curr);
                oldtab[i] = getNext(ptr);
                // Add to new
                const TTYPE h = hashEntry(ptr) % table.size();
                setNext(ptr, table[h]);
                table[h] = curr;
            }
        } // for i
    } else {
        for (TTYPE i=0; i<oldtab.size(); i++) {
            if (!oldtab[i]) continue;
            const void* ptr = MMAN->getChunkAddress(oldtab[i]);
            const TTYPE h = hashEntry(ptr) % table.size();
            setTable(h, oldtab[i]);
        } // for i
    }

    //
    // account for deletion of old table
    //
    mstats.decMemUsed(oldtab.size() * sizeof(TTYPE));
    mstats.decMemAlloc(oldtab.size() * sizeof(TTYPE));
}

// **********************************************************************

template <class TTYPE, bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
TTYPE MEDDLY::ct_tmpl<TTYPE, MONOLITHIC, CHAINED, INTSLOTS>
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
            if (sizeof(TTYPE) == sizeof(unsigned)) {
                ++ue;
            } else {
                ue += 2;
            }
        } else {
            ++cte;
        }
    }

    //
    // Determine entry type; do NOT advance pointer
    //
    const ct_entry_type* et;
    unsigned eind = 0;
    if (MONOLITHIC) {
        if (INTSLOTS) {
            et = ct_entry_type::getEntryType(ue[eind++]);
        } else {
            et = ct_entry_type::getEntryType(cte[eind++].U);
        }
    } else {
        et = ct_entry_type::getEntryType(global_etid);
    }
    ASSERT(__FILE__, __LINE__, et);

    //
    // Determine key size; do NOT advance pointer
    //
    unsigned repeats;
    if (et->isRepeating()) {
        repeats = INTSLOTS ? ue[eind++] : cte[eind++].U;
    } else {
        repeats = 0;
    }

    //
    // Hash the entire key and return, if we can!
    //
    if (et->isEntireKeyHashable(repeats)) {
        if (INTSLOTS) {
            eind += et->getKeyIntslots(repeats);
            return hash_stream::raw_hash(ue, eind);
        } else {
            eind += et->getKeySize(repeats);
            return hash_stream::raw_hash((unsigned*) cte, 2*eind);
        }
    }

    //
    // Still here? We can't hash the entire key,
    // so do everything in pieces.
    // Advance pointers to start of key.
    //
    if (INTSLOTS) {
        ue += eind;
    } else {
        cte += eind;
    }

    //
    // Hash operator and key size as needed
    //
    if (MONOLITHIC) {
        HSpush(et->getID());
    }
    if (et->isRepeating()) {
        HSpush(repeats);
    }

    //
    // Loop over key portion
    //
    const unsigned keylen = et->getKeySize(repeats);
    for (unsigned i=0; i<keylen; i++) {
        const ct_itemtype &it = et->getKeyType(i);
        if (!it.shouldBeHashed()) {
            //
            // Skip this one
            //
            if (INTSLOTS) {
                ue += it.intslots();
            } else {
                ++cte;
            }
            continue;
        }
        //
        // Hash this one
        //
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

    // TBD: use 64-bit when TTYPE is 64 bits
    return hash32();
}

// **********************************************************************

template <class TTYPE, bool M, bool C, bool I>
void MEDDLY::ct_tmpl<TTYPE, M, C, I>
    ::showEntry(output &s, TTYPE h, bool keyonly) const
{
    ASSERT(__FILE__, __LINE__, h);
    const void* hptr = MMAN->getChunkAddress(h);
    const unsigned* uptr = (unsigned*) hptr;
    const ct_entry_item* ctptr = (ct_entry_item*) hptr;

    //
    // Skip over next pointer
    //
    if (C) {
        if (I) {
            if (sizeof(TTYPE) == sizeof(unsigned)) {
                ++uptr;
            } else {
                uptr += 2;
            }
        } else {
            ++ctptr;
        }
    }

    //
    // Get operation and advance
    //
    const ct_entry_type* et = ct_entry_type::getEntryType(
         M  ?   ( I ? *uptr : ctptr->U )
            :   global_etid
    );
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
    switch (s.compression) {
        case compressionOption::None:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, true, true, false>(s, 0);
            else
                return new ct_tmpl<unsigned, true, true, false>(s, 0);

        case compressionOption::TypeBased:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, true, true, true>(s, 0);
            else
                return new ct_tmpl<unsigned, true, true, true>(s, 0);

        default:
                return 0;
    }
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
    switch (s.compression) {
        case compressionOption::None:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, true, false, false>(s, 0);
            else
                return new ct_tmpl<unsigned, true, false, false>(s, 0);

        case compressionOption::TypeBased:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, true, false, true>(s, 0);
            else
                return new ct_tmpl<unsigned, true, false, true>(s, 0);

        default:
                return nullptr;
    }
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
MEDDLY::operation_chained_style::create(const ct_settings &s, unsigned etid)
   const
{
    switch (s.compression) {
        case compressionOption::None:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, false, true, false>(s, etid);
            else
                return new ct_tmpl<unsigned, false, true, false>(s, etid);

        case compressionOption::TypeBased:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, false, true, true>(s, etid);
            else
                return new ct_tmpl<unsigned, false, true, true>(s, etid);

        default:
            return nullptr;
    }
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
MEDDLY::operation_unchained_style::create(const ct_settings &s, unsigned etid)
   const
{
    switch (s.compression) {
        case compressionOption::None:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, false, false, false>(s, etid);
            else
                return new ct_tmpl<unsigned, false, false, false>(s, etid);

        case compressionOption::TypeBased:
            if (s.allowHugeTables)
                return new ct_tmpl<unsigned long, false, false, false>(s, etid);
            else
                return new ct_tmpl<unsigned, false, false, false>(s, etid);

        default:
                return nullptr;
    }
}

