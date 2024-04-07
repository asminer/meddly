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

#ifndef MEDDLY_CT_VECTOR_H
#define MEDDLY_CT_VECTOR_H

#include "ct_entry_type.h"
#include "forest.h"

namespace MEDDLY {
    class ct_item;
    class ct_vector;
    class initializer_list;
};

/**
    A single item in a compute table entry.
*/
class MEDDLY::ct_item {
    public:
        /// Default constructor
        ct_item();

    public: // Setters

        /// Set item to a node
        inline void setN(node_handle nh) {
            type = ct_typeID::NODE;
            the_node = nh;
        }
        /// Set item from an integer
        inline void setI(int i) {
            type = ct_typeID::INTEGER;
            the_int = i;
        }
        /// Set item from a long
        inline void setL(long L) {
            type = ct_typeID::LONG;
            the_long = L;
        }
        /// Set item from a float
        inline void setF(float f) {
            type = ct_typeID::FLOAT;
            the_float = f;
        }
        /// Set item from a double
        inline void setD(double d) {
            type = ct_typeID::DOUBLE;
            the_double = d;
        }
        /// Set item from a generic
        inline void setG(ct_object* g) {
            type = ct_typeID::GENERIC;
            the_generic = g;
        }

        /// Set item from an edge value
        inline void set(const edge_value &ev) {
            switch (ev.getType()) {
                case edge_type::INT:    setI(ev.getInt());      return;
                case edge_type::LONG:   setL(ev.getLong());     return;
                case edge_type::FLOAT:  setF(ev.getFloat());    return;
                case edge_type::DOUBLE: setD(ev.getDouble());   return;
                default:    // includes VOID
                    type = ct_typeID::ERROR;
            }
        }

        // For templates
        inline void set_ev(long i)  { setL(i); }
        inline void set_ev(float f) { setF(f); }

        // For recycling and such
        inline void setNext(ct_item* n) {
            type = ct_typeID::ERROR;
            next = n;
        }

    public: // Getters
        /// Get a node
        inline node_handle getN() const {
            MEDDLY_DCASSERT(type == ct_typeID::NODE);
            return the_node;
        }
        /// Get an integer
        inline int getI() const {
            MEDDLY_DCASSERT(type == ct_typeID::INTEGER);
            return the_int;
        }
        /// Get a long
        inline long getL() const {
            MEDDLY_DCASSERT(type == ct_typeID::LONG);
            return the_long;
        }
        /// Get a float
        inline float getF() const {
            MEDDLY_DCASSERT(type == ct_typeID::FLOAT);
            return the_float;
        }
        /// Get a double
        inline double getD() const {
            MEDDLY_DCASSERT(type == ct_typeID::DOUBLE);
            return the_double;
        }
        /// Get a generic
        inline ct_object* getG() const {
            MEDDLY_DCASSERT(type == ct_typeID::GENERIC);
            return the_generic;
        }

        /// Get item into an edge value
        inline void get(edge_value &ev) const {
            switch (type) {
                case ct_typeID::INTEGER:    ev.set(the_int);        return;
                case ct_typeID::LONG:       ev.set(the_long);       return;
                case ct_typeID::FLOAT:      ev.set(the_float);      return;
                case ct_typeID::DOUBLE:     ev.set(the_double);     return;
                default:
                    MEDDLY_DCASSERT(false);
                    ev.set();
            }
        }

        // For templates
        inline void get_ev(long &i)     const   { i = getL(); }
        inline void get_ev(float &f)    const   { f = getF(); }

        // For recycling and such
        inline ct_item* getNext() const {
            MEDDLY_DCASSERT(type == ct_typeID::ERROR);
            return next;
        }

    private:
        union {
            int             the_int;
            long            the_long;
            node_handle     the_node;
            float           the_float;
            double          the_double;
            ct_object*      the_generic;
            ct_item*        next;
        };
        ct_typeID   type;
};


/**
    A vector of compute table items.
    Used to build keys and results.
    Internally, in the compute table, we may store
    entries differently.
*/
class MEDDLY::ct_vector {
    public:
        ct_vector(unsigned sz);
        ~ct_vector();

        inline const ct_item& operator[](unsigned u) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, u, size);
            MEDDLY_DCASSERT(data);
            return data[u];
        }
        inline ct_item& operator[](unsigned u) {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, u, size);
            MEDDLY_DCASSERT(data);
            return data[u];
        }

        inline unsigned getSize() const {
            return size;
        }

    private:
        const unsigned size;
        ct_item* data;

    private:
        friend class initializer_list;
        static void initStatics();
        static void doneStatics();

        static ct_item* useArray(unsigned sz);
        static void recycleArray(ct_item* v, unsigned sz);

    private:
        static ct_item* lists[16];
};

#endif

