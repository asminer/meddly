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
#include "error.h"

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

    class input;
    class output;
};


/**
    Unified edge value object.
*/
class MEDDLY::edge_value {

    public:
        /// Initializes to VOID
        edge_value();
        /// Initializes to INT
        edge_value(int v);
        /// Initializes to LONG
        edge_value(long v);
        /// Initializes to FLOAT
        edge_value(float v);
        /// Initializes to DOUBLE
        edge_value(double v);

        /// Conversion
        edge_value(edge_type newtype, const edge_value &v);

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
        inline edge_type getType() const {
            return mytype;
        }

        //
        // Get the value using type conversion
        //
        inline operator int() const {
            MEDDLY_DCASSERT(isInt());
            return ev_int;
        }
        inline operator long() const {
            MEDDLY_DCASSERT(isLong());
            return ev_long;
        }
        inline operator float() const {
            MEDDLY_DCASSERT(isFloat());
            return ev_float;
        }
        inline operator double() const {
            MEDDLY_DCASSERT(isDouble());
            return ev_double;
        }

        //
        // Getters for the value
        //

#ifdef ALLOW_DEPRECATED_0_18_0
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
#endif

        //
        // Alternate getters
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
        inline void getInt(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((int*) p) );
        }
        inline void getLong(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((long*) p) );
        }
        inline void getFloat(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((float*) p) );
        }
        inline void getDouble(void *p) const {
            MEDDLY_DCASSERT(p);
            get( *((double*) p) );
        }
        inline void get(edge_type et, void* p) const {
            switch (et) {
                case edge_type::VOID:
                    return;

                case edge_type::INT:
                    getInt(p);
                    return;

                case edge_type::LONG:
                    getLong(p);
                    return;

                case edge_type::FLOAT:
                    getFloat(p);
                    return;

                case edge_type::DOUBLE:
                    getDouble(p);
                    return;

                default:
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
            }
        }
        template <class T>
        inline void copyInto(T &val) const {
            switch (mytype) {
                case edge_type::INT:
                    val = T(ev_int);
                    return;

                case edge_type::LONG:
                    val = T(ev_long);
                    return;

                case edge_type::FLOAT:
                    val = T(ev_float);
                    return;

                case edge_type::DOUBLE:
                    val = T(ev_double);
                    return;

                default:
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
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

        template <class V>
        inline void setTempl(edge_type t, V v) {
            mytype = t;
            switch (mytype) {
                case edge_type::VOID:
                    return;

                case edge_type::INT:
                    ev_int = int(v);
                    return;

                case edge_type::LONG:
                    ev_long = long(v);
                    return;

                case edge_type::FLOAT:
                    ev_float = float(v);
                    return;

                case edge_type::DOUBLE:
                    ev_double = double(v);
                    return;

                default:
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
            }
        }

        //
        // Setters, for low-level storage objects
        //
        inline void setVoid(const void *p) {
            mytype = edge_type::VOID;
        }
        inline void setInt(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const int*) p) );
        }
        inline void setLong(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const long*) p) );
        }
        inline void setFloat(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const float*) p) );
        }
        inline void setDouble(const void *p) {
            MEDDLY_DCASSERT(p);
            set( *((const double*) p) );
        }
        inline void setType(edge_type et) {
            mytype = et;
        }
        inline void setRaw(const void* p) {
            switch (mytype) {
                case edge_type::VOID:
                    return;

                case edge_type::INT:
                    setInt(p);
                    return;

                case edge_type::LONG:
                    setLong(p);
                    return;

                case edge_type::FLOAT:
                    setFloat(p);
                    return;

                case edge_type::DOUBLE:
                    setDouble(p);
                    return;

                default:
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
            }
        }
        inline void set(edge_type et, const void* p) {
            setType(et);
            setRaw(p);
        }


        //
        // Add to the edge value
        //

        inline void add(int v) {
            MEDDLY_DCASSERT(isInt());
            ev_int += v;
        }
        inline void add(long v) {
            MEDDLY_DCASSERT(isLong());
            ev_long += v;
        }
        inline void add(float v) {
            MEDDLY_DCASSERT(isFloat());
            ev_float += v;
        }
        inline void add(double v) {
            MEDDLY_DCASSERT(isDouble());
            ev_double += v;
        }

        //
        // Subtract from the edge value
        // Useful for normalization
        //

        inline void subtract(int v) {
            MEDDLY_DCASSERT(isInt());
            ev_int -= v;
        }
        inline void subtract(long v) {
            MEDDLY_DCASSERT(isLong());
            ev_long -= v;
        }
        inline void subtract(float v) {
            MEDDLY_DCASSERT(isFloat());
            ev_float -= v;
        }
        inline void subtract(double v) {
            MEDDLY_DCASSERT(isDouble());
            ev_double -= v;
        }

        //
        // Multiply the edge value
        //

        inline void multiply(int v) {
            MEDDLY_DCASSERT(isInt());
            ev_int *= v;
        }
        inline void multiply(long v) {
            MEDDLY_DCASSERT(isLong());
            ev_long *= v;
        }
        inline void multiply(float v) {
            MEDDLY_DCASSERT(isFloat());
            ev_float *= v;
        }
        inline void multiply(double v) {
            MEDDLY_DCASSERT(isDouble());
            ev_double *= v;
        }

        //
        // Divide the edge value
        // Useful for normalization
        //

        inline void divide(int v) {
            MEDDLY_DCASSERT(isInt());
            MEDDLY_DCASSERT(v);
            ev_int /= v;
        }
        inline void divide(long v) {
            MEDDLY_DCASSERT(isLong());
            MEDDLY_DCASSERT(v);
            ev_long /= v;
        }
        inline void divide(float v) {
            MEDDLY_DCASSERT(isFloat());
            MEDDLY_DCASSERT(v);
            ev_float /= v;
        }
        inline void divide(double v) {
            MEDDLY_DCASSERT(isDouble());
            MEDDLY_DCASSERT(v);
            ev_double /= v;
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
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
            }
        }
        inline bool operator!=(const edge_value &v) const {
            switch (mytype) {
                case edge_type::VOID:
                    return !v.isVoid();

                case edge_type::INT:
                    return !v.equals(ev_int);

                case edge_type::LONG:
                    return !v.equals(ev_long);

                case edge_type::FLOAT:
                    return !v.equals(ev_float);

                case edge_type::DOUBLE:
                    return !v.equals(ev_double);

                default:
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
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
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
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
                    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
            }
        }

        //
        // File I/O
        //

        /**
            Read an edge value from an input file (stream).
            This will set the edge type appropriately.
        */
        void read(input &s);

        /**
            Write an edge value to an output file (stream).
            Formatted so that we can read the value back
            using read().
        */
        void write(output &s) const;

        /**
            Display an edge value for human consumption.
        */
        void show(output &s) const;

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
