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


// #define DEBUG_SLOW

// #define DEBUG_CT_SEARCHES
// #define DEBUG_CT_SLOTS

// #define DEBUG_REHASH
// #define DEBUG_TABLE2LIST
// #define DEBUG_LIST2TABLE
// #define DEBUG_CTALLOC

// #define DEBUG_ISDEAD
// #define DEBUG_ISSTALE
// #define DEBUG_REMOVESTALES
// #define DEBUG_REMOVESTALES_DETAILS

// #define DEBUG_CT_SCAN
// #define DEBUG_CT

// #define DEBUG_VALIDATE_COUNTS

#include "../memory.h"

#include "ct_typebased.h"
#include "ct_none.h"


#include <climits>

#include "../hash_stream.h"
#include "../ct_initializer.h"

/*
#include "../node_storage.h"
#include "../compute_table.h"
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../ct_vector.h"
*/

// Type for hash table indexes.
// Use unsigned for 32-bit integers (4 billion max table size)
// or unsigned long for 64-bit integers.
typedef unsigned hslot_type;

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***                                                                     ***
// ***                    New templatized compute table                    ***
// ***                                                                     ***
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

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
#endif

            virtual bool find(const ct_entry_type &ET, ct_vector &key,
                    ct_vector &res);
            virtual void addEntry(const ct_entry_type &ET, ct_vector &key,
                    const ct_vector &res);

            virtual void removeStales();
            virtual void removeAll();
            virtual void show(output &s, int verbLevel = 0);
            virtual void countNodeEntries(const forest* f, size_t* counts)
                const;

        protected:
            // Helper functions

#ifdef ALLOW_DEPRECATED_0_17_6
            static ct_entry_item* key2entry(const ct_entry_key& key,
                    ct_entry_item* e, hash_stream &H);
            static unsigned* key2entry(const ct_entry_key& key, unsigned* e,
                    hash_stream &H);

            inline static ct_entry_item* result2entry(
                    const ct_entry_result& res, ct_entry_item* e)
            {
                MEDDLY_DCASSERT(!INTSLOTS);
                const ct_entry_type* et = res.getET();
                MEDDLY_DCASSERT(et);
                const unsigned res_slots = et->getResultSize();
                memcpy(e, res.rawData(), res_slots * sizeof(ct_entry_item));
                return e + res_slots;
            }
            static unsigned* result2entry(const ct_entry_result& res,
                    unsigned* e);
#endif
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

            //
            //  Short vector equality check
            //
            template <class VT>
            inline static bool equal_sw(const VT* a, const VT* b, unsigned N)
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
                    default:    return (0==memcmp(a, b, N*sizeof(VT)));
                };
            }

            //
            // Chain management
            //
            inline static hslot_type getNext(const void* raw)
            {
                MEDDLY_DCASSERT(CHAINED);
                if (INTSLOTS) {
                    const unsigned* u = (const unsigned*) raw;
                    if (sizeof(hslot_type) == sizeof(unsigned)) {
                        return u[0];
                    } else {
                        hslot_type n = u[0];
                        n <<= 32;
                        n |= u[1];
                        return n;
                    }
                } else {
                    const ct_entry_item* ct = (const ct_entry_item*) raw;
                    return ct[0].UL;
                }
            }
            inline static void setNext(void* raw, hslot_type n)
            {
                MEDDLY_DCASSERT(CHAINED);
                if (INTSLOTS) {
                    unsigned* u = (unsigned*) raw;
                    if (sizeof(hslot_type) == sizeof(unsigned)) {
                        u[0] = n;
                    } else {
                        u[0] = (n >> 32);
                        u[1] = n & 0xffffffff;
                        return;
                    }
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
            static bool isDead(const ct_entry_type &ET, const unsigned* sres,
                    ct_entry_result &dres);
            /**
                Copy the result portion of an entry, and check if it is dead.
                    @param  ET      Entry type information
                    @param  sres    Pointer to result portion of entry
                    @param  dres    Where to copy the result
             */
            static bool isDead(const ct_entry_type &ET, const ct_entry_item* sres,
                    ct_entry_result &dres);
#endif

            //
            // Should we discard a CT miss?
            //
#ifdef ALLOW_DEPRECATED_0_17_6
            /**
                See if an entry is stale.
                    @param  e       Pointer to the entire entry
                    @param  mark    Mark nodes in the entry
             */
            bool isStale(const unsigned* e, bool mark) const;
            /**
                See if an entry is stale.
                    @param  e       Pointer to the entire entry
                    @param  mark    Mark nodes in the entry
             */
            bool isStale(const ct_entry_item* e, bool mark) const;
#endif

        private:
            /// Hash table
            std::vector <hslot_type> table;

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

#define maxCollisionSearch 2

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
    const ct_entry_type* et = key->getET();
    MEDDLY_DCASSERT(et);

    const unsigned chain_slots =
        CHAINED ? (INTSLOTS ? sizeof(hslot_type) / sizeof(unsigned) : 1) : 0;

    const unsigned op_slots = MONOLITHIC ? 1 : 0;

    const unsigned key_slots =
        INTSLOTS  ? ( et->getKeyBytes(key->numRepeats()) / sizeof(unsigned) )
                  : et->getKeySize(key->numRepeats());

    const unsigned res_slots =
        INTSLOTS  ? ( et->getResultBytes() / sizeof(unsigned) )
                  : et->getResultSize();

    const unsigned entry_slots =
          op_slots
        + (et->isRepeating() ? 1 : 0)
        + key_slots
    ;

    const unsigned num_slots =
          chain_slots
        + entry_slots
        + res_slots
    ;

    size_t slots = num_slots;

    key->my_entry = MMAN->requestChunk(slots);
    if (slots < num_slots) {
        throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
    }
    genentry keyentry;
    keyentry.vptr = MMAN->getChunkAddress(key->my_entry);
    perf.numEntries++;

    //
    // (2) build the entry, except for the result.
    //     we hash it as we go.
    //

    hash_stream H;
    H.start();
    if (INTSLOTS) {
        unsigned* e = keyentry.uptr;
        if (CHAINED) {
            e += sizeof(hslot_type)/sizeof(unsigned);
            // skip over NEXT pointer slot(s)
        }
        if (MONOLITHIC) {
            H.push(*e = et->getID());
            ++e;
        }
        if (et->isRepeating()) {
            H.push(*e = key->numRepeats());
            ++e;
        }
        key->result_shift = key2entry(key, e, H) - keyentry.uptr;
    } else {
        ct_entry_item* e = keyentry.ctptr;
        if (CHAINED) {
            ++e;    // skip over NEXT pointer
        }
        if (MONOLITHIC) {
            H.push(e->U = et->getID());
            ++e;
        }
        if (et->isRepeating()) {
            H.push(e->U = key->numRepeats());
            ++e;
        }
        key->result_shift = key2entry(key, e, H) - keyentry.ctptr;
    }

    //
    // (3) Save hash value
    //
    setHash(key, H.finish());
    const hslot_type hslot = key->getHash() % table.size();
    hslot_type hcurr = hslot;
    // TBD: ^ unsigned long
    // TBD: use a 64-bit hash

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
    hslot_type curr = table.at(hcurr);
#else
    hslot_type curr = table[hcurr];
#endif
    genentry currentry, preventry;
    preventry.vptr = nullptr;
    for (chainlen = 0; ; )
    {
        if (!curr) {
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
            curr = table.at(hcurr);
#else
            curr = table[hcurr];
#endif
            continue;
        } else {
            if (CHAINED) chainlen++;
        }

        //
        // Get the next entry
        //
        currentry.vptr = MMAN->getChunkAddress(curr);

        //
        // Check if it is equal to the entry we built
        //
        if (INTSLOTS
            ? equal_sw(currentry.uptr + chain_slots,
                        keyentry.uptr + chain_slots, entry_slots)
            : equal_sw(currentry.etptr + chain_slots,
                        keyentry.etptr + chain_slots, entry_slots) )
        {
            //
            // Equal, that's a CT hit.
            // See if we can use the result.
            //
            if (INTSLOTS
                ?   isDead(*et, currentry.uptr+chain_slots+entry_slots,  res)
                :   isDead(*et, currentry.etptr+chain_slots+entry_slots, res)
               )
            {
                //
                // Nope.
                // Delete this entry.
                //
                if (CHAINED) {
                    hslot_type next = getNext(currentry.vptr);
                    if (preventry.vptr) {
                        setNext(preventry.vptr, next);
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

                res.setInvalid();
                sawSearch(chainlen);
                deleteEntry(currentry);
                return;
            } else {
                //
                // Yes.
                // Move to front if we're chained
                //
                if (CHAINED) {
                    if (preventry.vptr) {
                        setNext(preventry.vptr, getNext(currentry.vptr));
                        setNext(currentry.vptr, table[hcurr]);
                        table[hcurr] = curr;
                    }
                    // if preventry.vptr is null, then
                    // we are already at the front.
                }

                res.reset();
                res.setValid();
                sawSearch(chainlen);
                perf.hits++;
                return;
            }

        }  // if equal

        //
        // Not equal.
        // See if this entry is stale.
        //
        if (checkStalesOnFind && (INTSLOTS
            ?   isStale(currentry.uptr, false)
            :   isStale(currentry.etptr, false)
           ))
        {
            //
            // Stale; delete
            //
            if (CHAINED) {
                hslot_type next = getNext(currentry.vptr);
                if (preventry.vptr) {
                    setNext(preventry.vptr, next);
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

        }

        //
        // Advance
        //
        if (CHAINED) {
            preventry.vptr = currentry.vptr;
            curr = getNext(currentry.vptr);
        } else {
            ++chainlen;
            if (chainlen > maxCollisionSearch) break;
            hcurr = (1+hcurr) % table.size();
#ifdef DEVELOPMENT_CODE
            curr = table.at(hcurr);
#else
            curr = table[hcurr];
#endif
        }

    } // for chainlen

    //
    // not found
    //

    sawSearch(chainlen);
    res.setInvalid();
}

#endif

// **********************************************************************
//
//      Helper methods
//
// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
MEDDLY::ct_entry_item* MEDDLY::ct_tmpl<M,C,I>::key2entry(
        const ct_entry_key& key, ct_entry_item* e, hash_stream &H)
{
    MEDDLY_DCASSERT(!I);
    const ct_entry_type* et = key.getET();
    MEDDLY_DCASSERT(et);
    const unsigned key_slots = et->getKeySize(key.numRepeats());
    const ct_entry_item* data = key.rawData();

    for (unsigned i=0; i<key_slots; i++) {
        // zero out the entire slot, makes comparisons easier
        e->UL = 0;
        switch (et->getKeyType(i).getType())
        {
            case ct_typeID::NODE:
            case ct_typeID::INTEGER:
                    H.push(e->U = data[i].U);
                    e++;
                    continue;

            case ct_typeID::LONG:
                    e->UL = data[i].UL;
                    H.push(unsigned(e->UL >> 32));
                    H.push(unsigned(e->UL & 0xffffffff));
                    e++;
                    continue;

            case ct_typeID::FLOAT:
                    // don't hash
                    e->F = data[i].F;
                    e++;
                    continue;

            case ct_typeID::DOUBLE:
                    // don't hash
                    e->D = data[i].D;
                    e++;
                    continue;

            default:
                  MEDDLY_DCASSERT(0);
        } // switch t
    } // for i

    return e;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
unsigned* MEDDLY::ct_tmpl<M,C,I>::key2entry(const ct_entry_key& key,
        unsigned* e, hash_stream &H)
{
    MEDDLY_DCASSERT(I);
    const ct_entry_type* et = key.getET();
    MEDDLY_DCASSERT(et);

    //
    // Copy the key into temp_entry
    //
    const ct_entry_item* data = key.rawData();
    const unsigned datalen = key.dataLength();
    for (unsigned i=0; i<datalen; i++) {
        switch (et->getKeyType(i).getType())
        {
            case ct_typeID::NODE:
            case ct_typeID::INTEGER:
                    MEDDLY_DCASSERT(sizeof(data[i].N) == sizeof(data[i].U));
                    MEDDLY_DCASSERT(sizeof(data[i].I) == sizeof(data[i].U));
                    H.push(*e = data[i].U);
                    e++;
                    continue;

            case ct_typeID::GENERIC:    // probably don't hash this way!
            case ct_typeID::LONG:
                    H.push(*e = (data[i].UL >> 32));
                    e++;
                    H.push(*e = (data[i].UL & 0xffffffff));
                    e++;
                    continue;

            case ct_typeID::FLOAT:
                    // DON'T hash
                    MEDDLY_DCASSERT(sizeof(data[i].F) == sizeof(data[i].U));
                    *e = data[i].U;
                    e++;
                    continue;

            case ct_typeID::DOUBLE:
                    // DON'T hash
                    *e = (data[i].UL >> 32);
                    e++;
                    *e = (data[i].UL & 0xffffffff);
                    e++;
                    continue;

            default:
                  MEDDLY_DCASSERT(0);
        } // switch t
    } // for i
    return e;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
unsigned* MEDDLY::ct_tmpl<M,C,I>::result2entry(const ct_entry_result& res,
        unsigned* e)
{
    MEDDLY_DCASSERT(I);
    const ct_entry_type* et = res.getET();
    MEDDLY_DCASSERT(et);

    //
    // Copy the result into temp_entry
    //
    const ct_entry_item* data = res.rawData();
    const unsigned datalen = res.dataLength();
    for (unsigned i=0; i<datalen; i++) {
        switch (et->getResultType(i).getType())
        {
            case ct_typeID::NODE:
            case ct_typeID::INTEGER:
            case ct_typeID::FLOAT:
                    MEDDLY_DCASSERT(sizeof(data[i].N) == sizeof(data[i].U));
                    MEDDLY_DCASSERT(sizeof(data[i].I) == sizeof(data[i].U));
                    MEDDLY_DCASSERT(sizeof(data[i].F) == sizeof(data[i].U));
                    *e = data[i].U;
                    e++;
                    continue;

            case ct_typeID::DOUBLE:
            case ct_typeID::GENERIC:
            case ct_typeID::LONG:
                    MEDDLY_DCASSERT(sizeof(data[i].G) == sizeof(data[i].UL));
                    MEDDLY_DCASSERT(sizeof(data[i].D) == sizeof(data[i].UL));
                    *e = (data[i].UL >> 32);
                    e++;
                    *e = (data[i].UL & 0xffffffff);
                    e++;
                    continue;

            default:
                  MEDDLY_DCASSERT(0);
        } // switch t
    } // for i
    return e;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<M,C,I>::isDead(const ct_entry_type &ET,
        const unsigned* sres, ct_entry_result &dres)
{
    MEDDLY_DCASSERT(I);

    MEDDLY_DCASSERT(sres);
    ct_entry_item* data = dres.rawData();
    MEDDLY_DCASSERT(dres.dataLength() == ET.getResultSize());

    MEDDLY_DCASSERT(sizeof(data[0].N) == sizeof(data[0].U));
    MEDDLY_DCASSERT(sizeof(data[0].I) == sizeof(data[0].U));
    MEDDLY_DCASSERT(sizeof(data[0].F) == sizeof(data[0].U));

    MEDDLY_DCASSERT(sizeof(data[0].G) == sizeof(data[0].UL));
    MEDDLY_DCASSERT(sizeof(data[0].D) == sizeof(data[0].UL));

    //
    // Copy the entry result from sres into dres,
    // and make sure it's not a dead entry while scanning it.
    //
    for (unsigned i=0; i<ET.getResultSize(); i++) {
        const ct_itemtype &item = ET.getResultType(i);

        switch (item.getType())
        {
            case ct_typeID::NODE:
                    MEDDLY_DCASSERT(item.hasForest());
                    if (item.getForest()->isDeadEntry(*sres)) {
                        return true;
                    }
                    continue;

            case ct_typeID::INTEGER:
            case ct_typeID::FLOAT:
                    data[i].U = *sres;
                    sres++;
                    continue;

            case ct_typeID::DOUBLE:
            case ct_typeID::GENERIC:
            case ct_typeID::LONG:
                    data[i].UL = *sres;
                    sres++;
                    data[i].UL <<= 32;
                    data[i].UL |= *sres;
                    sres++;
                    continue;
            default:
                  MEDDLY_DCASSERT(0);
        } // switch
    }
    return false;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<M,C,I>::isDead(const ct_entry_type &ET,
        const ct_entry_item* sres, ct_entry_result &dres)
{
    MEDDLY_DCASSERT(!I);

    MEDDLY_DCASSERT(sres);
    ct_entry_item* data = dres.rawData();
    MEDDLY_DCASSERT(dres.dataLength() == ET.getResultSize());

    //
    // Copy the entry result from sres into dres,
    // and make sure it's not a dead entry while scanning it.
    //
    for (unsigned i=0; i<ET.getResultSize(); i++) {
        const ct_itemtype &item = ET.getResultType(i);
        if (item.hasType(ct_typeID::NODE)) {
            MEDDLY_DCASSERT(item.hasForest());
            if (item.getForest()->isDeadEntry(sres[i].N)) {
                return true;
            }
        }

        data[i] = sres[i];

    }
    return false;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<M,C,I>::isStale(const unsigned* entry, bool mark) const
{
    MEDDLY_DCASSERT(I);

    //
    // Ignore next pointer, if there is one
    //
    entry += C ? (sizeof(hslot_type) / sizeof(unsigned)) : 0;


    //
    // Get entry type and check if those are marked
    // If monolithic, advance entry pointer
    //
    const ct_entry_type* et = M ? getEntryType(*entry++) : global_et;
    MEDDLY_DCASSERT(et);
    if (et->isMarkedForDeletion()) return true;

    //
    // Get entry size (advancing entry pointer if needed)
    //
    const unsigned reps = (et->isRepeating()) ? (*entry++) : 0;
    const unsigned klen = et->getKeySize(reps);

    //
    // Check the key portion of the entry
    //
    for (unsigned i=0; i<klen; i++) {
        const ct_itemtype &item = et->getKeyType(i);
        if (item.hasForest()) {
            if (item.getForest()->isStaleEntry(*entry)) {
                return true;
            } else {
                // Indicate that this node is in some cache entry
                if (mark) item.getForest()->setCacheBit(*entry);
            }
            entry++;
        } else {
            MEDDLY_DCASSERT(! item.hasType(ct_typeID::NODE));
            entry += item.intslots();
        }
    } // for i

    //
    // Check result portion of the entry
    //
    for (unsigned i=0; i<et->getResultSize(); i++) {
        const ct_itemtype &item = et->getResultType(i);
        if (item.hasForest()) {
            if (item.getForest()->isStaleEntry(*entry)) {
                return true;
            } else {
                if (mark) item.getForest()->setCacheBit(*entry);
            }
            entry++;
        } else {
            MEDDLY_DCASSERT(! item.hasType(ct_typeID::NODE));
            entry += item.intslots();
        }
    } // for i

    return false;
}

#endif

// **********************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6

template <bool M, bool C, bool I>
bool MEDDLY::ct_tmpl<M,C,I>::isStale(const ct_entry_item* entry, bool mark)
    const
{
    MEDDLY_DCASSERT(I);

    //
    // Ignore next pointer, if there is one
    //
    if (C) entry++;

    //
    // Get entry type and check if those are marked
    // If monolithic, advance entry pointer
    //
    const ct_entry_type* et = M ? getEntryType((entry++)->U) : global_et;
    MEDDLY_DCASSERT(et);
    if (et->isMarkedForDeletion()) return true;

    //
    // Get entry size (advancing entry pointer if needed)
    //
    const unsigned reps = (et->isRepeating()) ? ((entry++)->U) : 0;
    const unsigned klen = et->getKeySize(reps);

    //
    // Check the key portion of the entry
    //
    for (unsigned i=0; i<klen; i++) {
        const ct_itemtype &item = et->getKeyType(i);
        if (item.hasForest()) {
            if (item.getForest()->isStaleEntry(entry->N)) {
                return true;
            } else {
                // Indicate that this node is in some cache entry
                if (mark) item.getForest()->setCacheBit(entry->N);
            }
        } else {
            MEDDLY_DCASSERT(! item.hasType(ct_typeID::NODE));
        }
        entry++;
    } // for i

    //
    // Check result portion of the entry
    //
    for (unsigned i=0; i<et->getResultSize(); i++) {
        const ct_itemtype &item = et->getResultType(i);
        if (item.hasForest()) {
            if (item.getForest()->isStaleEntry(entry->N)) {
                return true;
            } else {
                if (mark) item.getForest()->setCacheBit(entry->N);
            }
        } else {
            MEDDLY_DCASSERT(! item.hasType(ct_typeID::NODE));
        }
        entry++;
    } // for i

    return false;
}

#endif

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
                return new ct_none<true, true>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_typebased<true, true>(s, 0, 0);
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
                return new ct_none<true, false>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_typebased<true, false>(s, 0, 0);
        default:
                return 0;
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
MEDDLY::operation_chained_style::create(const ct_settings &s, operation* op, unsigned slot) const
{
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<false, true>(s, op, slot);
        case compressionOption::TypeBased:
                return new ct_typebased<false, true>(s, op, slot);
        default:
                return 0;
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
MEDDLY::operation_unchained_style::create(const ct_settings &s, operation* op, unsigned slot) const
{
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<false, false>(s, op, slot);
        case compressionOption::TypeBased:
                return new ct_typebased<false, false>(s, op, slot);
        default:
                return 0;
    }
}

