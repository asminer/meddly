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
#include "hash_stream.h"

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

        // For recycling and such
        inline void setNext(ct_item* n) {
            next = n;
            type = ct_typeID::ERROR;
        }

        /// Set item to a node
        inline void setN(node_handle nh) {
            the_node = nh;
            type = ct_typeID::NODE;
        }
        /// Set item from an integer
        inline void setI(int i) {
            the_int = i;
            type = ct_typeID::INTEGER;
        }
        /// Set item from a long
        inline void setL(long L) {
            the_long = L;
            type = ct_typeID::LONG;
        }
        /// Set item from a float
        inline void setF(float f) {
            the_float = f;
            type = ct_typeID::FLOAT;
        }
        /// Set item from a double
        inline void setD(double d) {
            the_double = d;
            type = ct_typeID::DOUBLE;
        }
        /// Set item from a generic
        inline void setG(ct_object* g) {
            the_generic = g;
            type = ct_typeID::GENERIC;
        }
        /// Set item to a relation node
        inline void setR(node_handle nh) {
            the_node = nh;
            type = ct_typeID::RELNODE;
        }

        /// Set item from an edge value
        inline void set(const edge_value &ev) {
            switch (ev.getType()) {
                case edge_type::INT:    setI(ev.getInt());      return;
                case edge_type::LONG:   setL(ev.getLong());     return;
                case edge_type::FLOAT:  setF(ev.getFloat());    return;
                case edge_type::DOUBLE: setD(ev.getDouble());   return;
                default:                setNext(nullptr);
            }
        }

        /// Set item from a raw ct_entry_item
        inline void set(ct_typeID t, ct_entry_item ci)
        {
            switch (t) {
                case ct_typeID::ERROR:      setNext(nullptr);   return;
                case ct_typeID::NODE:       setN(ci.N);         return;
                case ct_typeID::INTEGER:    setI(ci.I);         return;
                case ct_typeID::LONG:       setL(ci.L);         return;
                case ct_typeID::FLOAT:      setF(ci.F);         return;
                case ct_typeID::DOUBLE:     setD(ci.D);         return;
                case ct_typeID::GENERIC:    setG(ci.G);         return;
                case ct_typeID::RELNODE:    setR(ci.N);         return;
                default:                    FAIL(__FILE__, __LINE__);
            }
        }

        /// Set item from a compacted CT entry, and advance the pointer.
        inline const unsigned* set(ct_typeID t, const unsigned* ptr)
        {
            ASSERT(__FILE__, __LINE__, sizeof(the_node) == sizeof(unsigned));
            ASSERT(__FILE__, __LINE__, sizeof(the_int) == sizeof(unsigned));
            ASSERT(__FILE__, __LINE__, sizeof(the_float) == sizeof(unsigned));

            ASSERT(__FILE__, __LINE__, sizeof(the_long) == 2*sizeof(unsigned));
            ASSERT(__FILE__, __LINE__, sizeof(the_double) == 2*sizeof(unsigned));
            ASSERT(__FILE__, __LINE__, sizeof(the_generic) == 2*sizeof(unsigned));

            type = t;
            switch (t) {
                case ct_typeID::ERROR:
                        setNext(nullptr);
                        return ptr;

                case ct_typeID::NODE:
                case ct_typeID::RELNODE:
                case ct_typeID::INTEGER:
                case ct_typeID::FLOAT:
                        raw[0] = ptr[0];
                        // slots = 1;
                        return ptr+1;

                case ct_typeID::LONG:
                case ct_typeID::DOUBLE:
                case ct_typeID::GENERIC:
                        raw[0] = ptr[0];
                        raw[1] = ptr[1];
                        // slots = 2;
                        return ptr+2;

                default:
                        FAIL(__FILE__, __LINE__);
            }
            // fail safe
            return nullptr;
        }

        // For templates
        inline void set_ev(long i)  { setL(i); }
        inline void set_ev(float f) { setF(f); }

    public: // Getters
        /// Get a node
        inline node_handle getN() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::NODE);
            return the_node;
        }
        /// Get an integer
        inline int getI() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::INTEGER);
            return the_int;
        }
        /// Get a long
        inline long getL() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::LONG);
            return the_long;
        }
        /// Get a float
        inline float getF() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::FLOAT);
            return the_float;
        }
        /// Get a double
        inline double getD() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::DOUBLE);
            return the_double;
        }
        /// Get a generic
        inline ct_object* getG() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::GENERIC);
            return the_generic;
        }
        /// Get a relation node
        inline node_handle getR() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::RELNODE);
            return the_node;
        }

        /// Get item into an edge value
        inline void get(edge_value &ev) const {
            switch (type) {
                case ct_typeID::INTEGER:    ev.set(the_int);        return;
                case ct_typeID::LONG:       ev.set(the_long);       return;
                case ct_typeID::FLOAT:      ev.set(the_float);      return;
                case ct_typeID::DOUBLE:     ev.set(the_double);     return;
                default:
                    FAIL(__FILE__, __LINE__, "Unknown edge type");
                    ev.set();
            }
        }

        // For templates
        inline void get_ev(long &i)     const   { i = getL(); }
        inline void get_ev(float &f)    const   { f = getF(); }

        // For recycling and such
        inline ct_item* getNext() const {
            ASSERT(__FILE__, __LINE__, type == ct_typeID::ERROR);
            return next;
        }

        // ================================================================
        // CT Entry encodings
        // ================================================================
        inline unsigned long rawUL() const {
            return UL;
        }
        inline unsigned raw0() const {
            return raw[0];
        }
        inline unsigned raw1() const {
            return raw[1];
        }

        inline unsigned long& rawUL() {
            return UL;
        }
        inline unsigned& raw0() {
            return raw[0];
        }
        inline unsigned& raw1() {
            return raw[1];
        }


    public: // Type checking

        inline ct_typeID getType() const { return type; }
        inline bool hasType(ct_typeID t) const { return type == t; }
        inline void setType(ct_typeID t) { type = t; }

    public: // miscellaneous

        void show(output &s) const;

    private:
        union {
            int             the_int;
            long            the_long;
            node_handle     the_node;
            float           the_float;
            double          the_double;
            ct_object*      the_generic;
            ct_item*        next;

            // For packing into compute table entries
            unsigned long   UL;

            // For hashing
            unsigned        raw[2];
        };
        ct_typeID   type;
};


// ===========================================================================

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
            CHECK_RANGE(__FILE__, __LINE__, 0u, u, size());
            ASSERT(__FILE__, __LINE__, data);
            return data[u];
        }
        inline ct_item& operator[](unsigned u) {
            CHECK_RANGE(__FILE__, __LINE__, 0u, u, size());
            ASSERT(__FILE__, __LINE__, data);
            return data[u];
        }

        inline unsigned size() const {
            return _size;
        }

        inline void setHash(unsigned h) {
            hashVal = h;
#ifdef DEVELOPMENT_CODE
            hasHashVal = true;
#endif
        }
        inline unsigned getHash() const {
#ifdef DEVELOPMENT_CODE
            ASSERT(__FILE__, __LINE__, hasHashVal);
#endif
            return hashVal;
        }

        void show(output &s) const;

    private:
        const unsigned _size;
        unsigned hashVal;
        ct_item* data;
#ifdef DEVELOPMENT_CODE
        bool hasHashVal;
#endif

    public:
        //
        // Stored data for CT misses,
        // to save time for CT adds later.
        //
        unsigned long my_entry;
        unsigned entry_slots;
        unsigned result_shift;

    private:
        friend class initializer_list;
        static void initStatics();
        static void doneStatics();

        static ct_item* useArray(unsigned sz);
        static void recycleArray(ct_item* v, unsigned sz);

    private:
        static ct_item* lists[64];
};

#endif

