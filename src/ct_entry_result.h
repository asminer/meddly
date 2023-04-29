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

#ifndef MEDDLY_CT_ENTRY_RESULT_H
#define MEDDLY_CT_ENTRY_RESULT_H

#include "ct_entry_type.h"

namespace MEDDLY {
    class ct_entry_result;
    class compute_table;
};

/**
        The result portion of an entry.
        Internally, in the compute table, we may store
        entries differently.  This class is used to return
        results from searches and to construct CT entries.
*/
class MEDDLY::ct_entry_result {
        public:
          ct_entry_result();
          ~ct_entry_result();

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
          void setValid(const ct_entry_item* d);
          void setInvalid();
          operator bool() const;
          /// Increase cache counters for nodes in this portion of the entry.
          void cacheNodes() const;

          const ct_entry_item* rawData() const;
          unsigned dataLength() const;


        private:
          const ct_entry_type* etype;
          ct_entry_item* build;
          const ct_entry_item* data;
          bool is_valid;
          unsigned currslot;
};

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
MEDDLY::compute_table::entry_result::setValid(const ct_entry_item* d)
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

inline const MEDDLY::compute_table::ct_entry_item*
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




#endif // #include guard
