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
    public:
          ct_entry_key();
          ~ct_entry_key();

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
          const ct_entry_item* rawData() const;
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

    friend class compute_table;
};

// ******************************************************************
// *                                                                *
// *                  inlined ct_entry_key methods                  *
// *                                                                *
// ******************************************************************

inline void
MEDDLY::ct_entry_key::setup(const ct_entry_type* et, unsigned repeats)
{
  MEDDLY_DCASSERT(et);
  etype = et;
  num_repeats = repeats;
  MEDDLY_DCASSERT( 0==repeats || et->isRepeating() );
  total_slots = et->getKeySize(repeats);
  if (total_slots > data_alloc) {
    data_alloc = (1+(data_alloc / 8)) * 8;   // allocate in chunks of size 8
    data = (ct_entry_item*) realloc(data, data_alloc*sizeof(ct_entry_item));
    if (0==data) throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  memset(data, 0, total_slots * sizeof(ct_entry_item));
  currslot = 0;
#ifdef DEVELOPMENT_CODE
  has_hash = false;
#endif
}

inline const MEDDLY::ct_entry_type*
MEDDLY::ct_entry_key::getET() const
{
  return etype;
}

inline void MEDDLY::ct_entry_key::writeN(node_handle nh)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::NODE == theSlotType());
  data[currslot++].N = nh;
}

inline void MEDDLY::ct_entry_key::writeI(int i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::INTEGER == theSlotType());
  data[currslot++].I = i;
}

inline void MEDDLY::ct_entry_key::writeL(long i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::LONG == theSlotType());
  data[currslot++].L = i;
}

inline void MEDDLY::ct_entry_key::writeF(float f)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, currslot, total_slots);
  MEDDLY_DCASSERT(ct_typeID::FLOAT == theSlotType());
  data[currslot++].F = f;
}

inline const MEDDLY::ct_entry_item*
MEDDLY::ct_entry_key::rawData() const
{
  return data;
}

inline unsigned MEDDLY::ct_entry_key::dataLength() const
{
  return total_slots;
}

inline unsigned MEDDLY::ct_entry_key::numRepeats() const
{
  return num_repeats;
}

inline const void*
MEDDLY::ct_entry_key::readTempData() const
{
  return temp_data;
}

inline unsigned
MEDDLY::ct_entry_key::numTempBytes() const
{
  return temp_bytes;
}

inline void*
MEDDLY::ct_entry_key::allocTempData(unsigned bytes)
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
MEDDLY::ct_entry_key::cacheNodes() const
{
  for (unsigned i=0; i<total_slots; i++) {
    expert_forest* f = etype->getKeyForest(i);
    if (f) {
      f->cacheNode(data[i].N);
    }
  }
}

inline unsigned MEDDLY::ct_entry_key::getHash() const
{
  MEDDLY_DCASSERT(has_hash);
  return hash_value;
}

inline void MEDDLY::ct_entry_key::setHash(unsigned h)
{
  hash_value = h;
#ifdef DEVELOPMENT_CODE
  has_hash = true;
#endif
}

inline MEDDLY::ct_typeID MEDDLY::ct_entry_key::theSlotType() const
{
  //
  // Adjust currslot for OP entry, and number of repeats entry
  //
  // return etype->getKeyType(currslot - (etype->isRepeating() ? 2 : 1) );
  return etype->getKeyType(currslot);
}


#endif // #include guard

