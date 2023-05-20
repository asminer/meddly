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

#ifndef MEDDLY_CT_ENTRY_KEY_H
#define MEDDLY_CT_ENTRY_KEY_H

#include "ct_entry_type.h"
#include "error.h"
#include "forest.h"

namespace MEDDLY {
    class ct_entry_key;
    class compute_table;
};

/**
    The key portion of a compute table entry.
    Internally, in the compute table, we may store
    entries differently.  This class is used to build
    keys for searching and to construct CT entries.
*/
class MEDDLY::ct_entry_key {
        friend class compute_table;
    public:
        ct_entry_key();
        ~ct_entry_key();

    protected:
        /// Start using for this operation
        void setup(const ct_entry_type* et, unsigned repeats);

    public:
        inline const ct_entry_type* getET() const { return etype; }

        /// Write a node into the next slot
        inline void writeN(node_handle nh) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, currslot, total_slots);
            MEDDLY_DCASSERT(ct_typeID::NODE == theSlotType());
            data[currslot++].N = nh;
        }

        /// Write an integer into the next slot
        inline void writeI(int i) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, currslot, total_slots);
            MEDDLY_DCASSERT(ct_typeID::INTEGER == theSlotType());
            data[currslot++].I = i;
        }

        /// Write a long into the next slot
        inline void writeL(long i) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, currslot, total_slots);
            MEDDLY_DCASSERT(ct_typeID::LONG == theSlotType());
            data[currslot++].L = i;
        }

        /// Write a float into the next slot
        inline void writeF(float f) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, currslot, total_slots);
            MEDDLY_DCASSERT(ct_typeID::FLOAT == theSlotType());
            data[currslot++].F = f;
        }

        // For templates
        inline void write_ev(long i)  { writeL(i); }
        inline void write_ev(float f) { writeF(f); }

    public:
        //
        // Interface, mostly for compute table.
        //

        inline const ct_entry_item* rawData() const { return data; }

        inline unsigned dataLength() const { return total_slots; }
        inline unsigned numRepeats() const { return num_repeats; }

        inline const void* readTempData() const { return temp_data; }
        inline unsigned numTempBytes() const { return temp_bytes; }

        inline unsigned getHash() const {
            MEDDLY_DCASSERT(has_hash);
            return hash_value;
        }

        void* allocTempData(unsigned bytes);

        /// Increase cache counters for nodes in this portion of the entry.
        inline void cacheNodes() const {
            for (unsigned i=0; i<total_slots; i++) {
                expert_forest* f = etype->getKeyForest(i);
                if (f) {
                    f->cacheNode(data[i].N);
                }
            }
        }

    protected:
        // protected interface, for compute_table.
        inline void setHash(unsigned h) {
            hash_value = h;
#ifdef DEVELOPMENT_CODE
            has_hash = true;
#endif
        }

    private:
        inline ct_typeID theSlotType() const {
            return etype->getKeyType(currslot);
        }



    private:
        const ct_entry_type* etype;
        ct_entry_item* data;
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
        ct_entry_key* next;
};


#endif // #include guard

