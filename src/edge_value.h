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

#ifndef MEDDLY_EDGE_VALUE_H
#define MEDDLY_EDGE_VALUE_H

#include "defines.h"
#include "hash_stream.h"

namespace MEDDLY {

    /// Type of value on an edge.
    enum class edge_type {
        /// Nothing; for multi-terminal
        VOID,
        /// int
        INT,
        /// long
        LONG,
        /// float
        FLOAT,
        /// double
        DOUBLE
    };

    class edge_value;
};


/**
    Unified edge value object.
*/
class MEDDLY::edge_value {

    public:
        edge_value();

        //
        // Getters for the type
        //

        inline bool isVoid() const {
            return (edge_type::VOID == mytype);
        }
        inline bool isInt() const {
            return (edge_type::INT == mytype);
        }
        inline bool isLong() const {
            return (edge_type::LONG == mytype);
        }
        inline bool isFloat() const {
            return (edge_type::FLOAT == mytype);
        }
        inline bool isDouble() const {
            return (edge_type::DOUBLE == mytype);
        }
        inline bool hasType(edge_type t) const {
            return (t == mytype);
        }

        //
        // Getters for the value
        //

        inline int getInt() const {
            MEDDLY_DCASSERT(isInt());
            return ev_int;
        }
        inline long getLong() const {
            MEDDLY_DCASSERT(isLong());
            return ev_long;
        }
        inline float getFloat() const {
            MEDDLY_DCASSERT(isFloat());
            return ev_float;
        }
        inline double getDouble() const {
            MEDDLY_DCASSERT(isDouble());
            return ev_double;
        }

        //
        // Alternate getters (better for templates)
        //

        inline void get(int &v) const {
            MEDDLY_DCASSERT(isInt());
            v = ev_int;
        }
        inline void get(long &v) const {
            MEDDLY_DCASSERT(isLong());
            v = ev_long;
        }
        inline void get(float &v) const {
            MEDDLY_DCASSERT(isFloat());
            v = ev_float;
        }
        inline void get(double &v) const {
            MEDDLY_DCASSERT(isDouble());
            v = ev_double;
        }

        //
        // Getters, for low-level storage objects
        //
        inline void get_int(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((int*) p) );
        }
        inline void get_long(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((long*) p) );
        }
        inline void get_float(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((float*) p) );
        }
        inline void get_double(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((double*) p) );
        }
        inline void get(edge_type et, void* p) const {
            switch (et) {
                case edge_type::VOID:
                    return;

                case edge_type::INT:
                    get_int(p);
                    return;

                case edge_type::LONG:
                    get_long(p);
                    return;

                case edge_type::FLOAT:
                    get_float(p);
                    return;

                case edge_type::DOUBLE:
                    get_double(p);
                    return;

                default:
                   MEDDLY_DCASSERT(false);
            }
        }
        //
        // Setters
        //

        inline void set() {
            mytype = edge_type::VOID;
        }
        inline void set(int v) {
            mytype = edge_type::INT;
            ev_int = v;
        }
        inline void set(long v) {
            mytype = edge_type::LONG;
            ev_long = v;
        }
        inline void set(float v) {
            mytype = edge_type::FLOAT;
            ev_float = v;
        }
        inline void set(double v) {
            mytype = edge_type::DOUBLE;
            ev_double = v;
        }

        //
        // Setters, for low-level storage objects
        //
        inline void set_void(const void *p) {
            mytype = edge_type::VOID;
        }
        inline void set_int(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const int*) p) );
        }
        inline void set_long(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const long*) p) );
        }
        inline void set_float(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const float*) p) );
        }
        inline void set_double(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const double*) p) );
        }
        inline void set(edge_type et, const void* p) {
            switch (et) {
                case edge_type::VOID:
                    set_void(p);
                    return;

                case edge_type::INT:
                    set_int(p);
                    return;

                case edge_type::LONG:
                    set_long(p);
                    return;

                case edge_type::FLOAT:
                    set_float(p);
                    return;

                case edge_type::DOUBLE:
                    set_double(p);
                    return;

                default:
                   MEDDLY_DCASSERT(false);
            }
        }



        //
        // Equality checks
        //
        inline bool equals(int v) const {
            if (mytype != edge_type::INT) return false;
            return v == ev_int;
        }
        inline bool equals(long v) const {
            if (mytype != edge_type::LONG) return false;
            return v == ev_long;
        }
        inline bool equals(float v) const {
            if (mytype != edge_type::FLOAT) return false;
            if (v) {
                double diff = v-ev_float;
                return ABS(diff/v) < 1e-6;
            } else {
                return ABS(ev_float) < 1e-10;
            }
        }
        inline bool equals(double v) const {
            if (mytype != edge_type::DOUBLE) return false;
            if (v) {
                double diff = v-ev_double;
                return ABS(diff/v) < 1e-6;
            } else {
                return ABS(ev_double) < 1e-10;
            }
        }
        inline bool operator==(const edge_value &v) const {
            switch (mytype) {
                case edge_type::VOID:
                    return v.isVoid();

                case edge_type::INT:
                    return v.equals(ev_int);

                case edge_type::LONG:
                    return v.equals(ev_long);

                case edge_type::FLOAT:
                    return v.equals(ev_float);

                case edge_type::DOUBLE:
                    return v.equals(ev_double);

                default:
                    MEDDLY_DCASSERT(0);
                    return true;
            }
        }


        //
        // Equality check, for low-level storage objects
        //
        inline bool equals(const void* p) const {
            switch (mytype) {
                case edge_type::VOID:
                    return isVoid();

                case edge_type::INT:
                    return equals( *((const int*)p) );

                case edge_type::LONG:
                    return equals( *((const long*)p) );

                case edge_type::FLOAT:
                    return equals( *((const float*)p) );

                case edge_type::DOUBLE:
                    return equals( *((const double*)p) );

                default:
                   MEDDLY_DCASSERT(false);
            }
        }

        //
        // Hash the edge value.
        //
        inline void hash(hash_stream &h) const {
            switch (mytype) {
                case edge_type::INT:
                    h.push(&ev_int, sizeof(int));
                    return;

                case edge_type::LONG:
                    h.push(&ev_long, sizeof(long));
                    return;

                case edge_type::FLOAT:
                    h.push(&ev_float, sizeof(float));
                    return;

                case edge_type::DOUBLE:
                    h.push(&ev_double, sizeof(double));
                    return;

                case edge_type::VOID:
                default:
                    MEDDLY_DCASSERT(false);
            }
        }

    private:
        union {
            int     ev_int;
            long    ev_long;
            float   ev_float;
            double  ev_double;
        };

        edge_type mytype;
};


#endif // include guard
