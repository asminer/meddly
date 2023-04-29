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
#include "forest.h"

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
        /// Initialize after construction
        void initialize(const ct_entry_type* et);

        //
        // interface, for operations (reading).
        //

        inline node_handle readN() {
            READ_SLOT(ct_typeID::NODE);
            return data[currslot++].N;
        }

        inline int readI() {
            READ_SLOT(ct_typeID::INTEGER);
            return data[currslot++].I;
        }

        inline float readF() {
            READ_SLOT(ct_typeID::FLOAT);
            return data[currslot++].F;
        }

        inline long readL() {
            READ_SLOT(ct_typeID::LONG);
            return data[currslot++].L;
        }

        inline double readD() {
            READ_SLOT(ct_typeID::DOUBLE);
            return data[currslot++].D;
        }

        inline ct_object* readG() {
            READ_SLOT(ct_typeID::GENERIC);
            return data[currslot++].G;
        }

        // for templates
        void read_ev(long &l)   { l = readL(); }
        void read_ev(float &f)  { f = readF(); }

        //
        // interface, for operations (building).
        //

        inline void reset() { currslot = 0; }

        inline void writeN(node_handle nh) {
            WRITE_SLOT(ct_typeID::NODE);
            build[currslot++].N = nh;
        }

        inline void writeI(int i) {
            WRITE_SLOT(ct_typeID::INTEGER);
            build[currslot++].I = i;
        }

        inline void writeF(float f) {
            WRITE_SLOT(ct_typeID::FLOAT);
            build[currslot++].F = f;
        }

        inline void writeL(long L) {
            WRITE_SLOT(ct_typeID::LONG);
            build[currslot++].L = L;
        }

        inline void writeD(double D) {
            WRITE_SLOT(ct_typeID::DOUBLE);
            build[currslot++].D = D;
        }

        inline void writeG(ct_object* G) {
            WRITE_SLOT(ct_typeID::GENERIC);
            build[currslot++].G = G;
        }

        //
        // interface, for compute tables.
        //

        inline void setValid() {
            is_valid = true;
            data = build;
        }

        inline void setValid(const ct_entry_item* d) {
            is_valid = true;
            data = d;
        }

        inline void setInvalid() { is_valid = false; }

        inline operator bool() const { return is_valid; }

        inline const ct_entry_item* rawData() const { return build; }
        inline unsigned dataLength() const { return etype->getResultSize(); }

        /// Increase cache counters for nodes in this portion of the entry.
        inline void cacheNodes() const {
            for (unsigned i=0; i<etype->getResultSize(); i++) {
                expert_forest* f = etype->getResultForest(i);
                if (f) {
                    f->cacheNode(build[i].N);
                }
            }
        }

    private:
        inline void READ_SLOT(ct_typeID t)
        {
            MEDDLY_DCASSERT(data);
            MEDDLY_DCASSERT(currslot < dataLength());
            MEDDLY_DCASSERT(t == etype->getResultType(currslot));
        }
        inline void WRITE_SLOT(ct_typeID t)
        {
            MEDDLY_DCASSERT(build);
            MEDDLY_DCASSERT(currslot < dataLength());
            MEDDLY_DCASSERT(t == etype->getResultType(currslot));
        }


    private:
        const ct_entry_type* etype;
        ct_entry_item* build;
        const ct_entry_item* data;
        bool is_valid;
        unsigned currslot;
};

#endif // #include guard
