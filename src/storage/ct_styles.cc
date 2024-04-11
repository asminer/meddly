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
            virtual void updateEntry(ct_entry_key* key,
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

        private:
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

    //
    // (1) determine entry size, and allocate
    //
    const ct_entry_type* et = key->getET();
    MEDDLY_DCASSERT(et);
    const unsigned key_slots =
        INTSLOTS  ? ( et->getKeyBytes(key->numRepeats()) / sizeof(unsigned) )
                  : et->getKeySize(key->numRepeats());
    const unsigned res_slots =
        INTSLOTS  ? ( et->getResultBytes() / sizeof(unsigned) )
                  : et->getResultSize();
    size_t num_slots =
          (CHAINED ? 1 : 0)
        + (MONOLITHIC ? 1 : 0)
        + (et->isRepeating() ? 1 : 0)
        + key_slots
        + res_slots
    ;
    key->my_entry = MMAN->requestChunk(num_slots);
    perf.numEntries++;

    //
    // (2) build the entry, except for the result.
    //     we hash it as we go.
    //

    hash_stream H;
    H.start();
    if (INTSLOTS) {
        unsigned* entry = (unsigned*) MMAN->getChunkAddress(key->my_entry);
        if (CHAINED) {
            entry += 2;  // skip over NEXT pointer slots
        }
        if (MONOLITHIC) {
            H.push(*entry = et->getID());
            ++entry;
        }
        if (et->isRepeating()) {
            H.push(*entry = key->numRepeats());
            ++entry;
        }
        key->resptr = key2entry(key, entry, H);
    } else {
        ct_entry_item* entry =
            (ct_entry_item*) MMAN->getChunkAddress(key->my_entry);
        if (CHAINED) {
            ++entry;    // skip over NEXT pointer
        }
        if (MONOLITHIC) {
            H.push(entry->U = et->getID());
            ++entry;
        }
        if (et->isRepeating()) {
            H.push(entry->U = key->numRepeats());
            ++entry;
        }
        key->resptr = key2entry(key, entry, H);
    }

    //
    // (3) Save hash value, and search the table.
    //
    setHash(key, H.finish());
    const unsigned hslot = key->getHash() % table.size();
    // TBD: ^ unsigned long
    // TBD: use a 64-bit hash
}
#endif

// **********************************************************************
// Helper functions
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

