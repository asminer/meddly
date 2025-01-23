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

#ifndef MEDDLY_RANGEVAL_H
#define MEDDLY_RANGEVAL_H

#include "defines.h"
#include "error.h"

namespace MEDDLY {
    class rangeval;
    class forest;
    class output;

    /// Types of return values allowed for functions.
    enum class range_type {
        /// boolean-valued functions.
        BOOLEAN,
        /// integer-valued functions.
        INTEGER,
        /// real-valued functions.
        REAL
    };

    inline const char* nameOf(range_type rt)
    {
        switch (rt) {
            case range_type::BOOLEAN:   return "boolean";
            case range_type::INTEGER:   return "integer";
            case range_type::REAL:      return "real";
            default:                    return "unknown range type";
        }
    }

    inline const char* shortNameOf(range_type rt)
    {
        switch (rt) {
            case range_type::BOOLEAN:   return "bool";
            case range_type::INTEGER:   return "int ";
            case range_type::REAL:      return "real";
            default:                    return "????";
        }
    }

    /// Special values allowed for functions.
    enum class range_special {
        /// Normal value, not special
        NORMAL,

        /// Positive infinity.
        PLUS_INFINITY

            // TBD: others? undefined?
    };

};  // namespace

/**
    Object to hold a generic return value.
    We keep track of the type and the value.
    Also we allow special values, like infinity.
*/
class MEDDLY::rangeval {
    public:
        /// Initialize as a boolean
        rangeval(bool v=true)
        {
            the_type = range_type::BOOLEAN;
            s_value = range_special::NORMAL;
            l_value = v;
        }
        /// Initialize as an integer
        rangeval(int v)
        {
            the_type = range_type::INTEGER;
            s_value = range_special::NORMAL;
            l_value = v;
        }
        /// Initialize as an integer
        rangeval(long v)
        {
            the_type = range_type::INTEGER;
            s_value = range_special::NORMAL;
            l_value = v;
        }
        /// Initialize as a real
        rangeval(float v)
        {
            the_type = range_type::REAL;
            s_value = range_special::NORMAL;
            d_value = v;
        }
        /// Initialize as a real
        rangeval(double v)
        {
            the_type = range_type::REAL;
            s_value = range_special::NORMAL;
            d_value = v;
        }
        /// Initialize as special
        rangeval(range_special v, range_type rt)
        {
            the_type = rt;
            s_value = v;
        }


        //
        // Getters for the type
        //

        inline bool hasType(range_type t) const {
            return (t == the_type);
        }
        inline bool isBoolean() const {
            return hasType(range_type::BOOLEAN);
        }
        inline bool isInteger() const {
            return hasType(range_type::INTEGER);
        }
        inline bool isReal() const {
            return hasType(range_type::REAL);
        }
        inline range_type getType() const {
            return the_type;
        }

        //
        // Getters for the value
        //

        inline bool isNormal() const {
            return range_special::NORMAL == s_value;
        }
        inline bool isPlusInfinity() const {
            return range_special::PLUS_INFINITY == s_value;
        }

        inline operator bool() const {
            MEDDLY_DCASSERT(isBoolean());
            if (!isNormal()) {
                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
            }
            return l_value;
        }
        inline operator int() const {
            MEDDLY_DCASSERT(isInteger());
            if (!isNormal()) {
                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
            }
            return l_value;
        }
        inline operator long() const {
            MEDDLY_DCASSERT(isInteger());
            if (!isNormal()) {
                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
            }
            return l_value;
        }
        inline operator float() const {
            MEDDLY_DCASSERT(isReal());
            if (!isNormal()) {
                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
            }
            return d_value;
        }
        inline operator double() const {
            MEDDLY_DCASSERT(isReal());
            if (!isNormal()) {
                throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
            }
            return d_value;
        }

        inline bool operator==(const rangeval &t) const {
            return equals(t);
        }
        inline bool operator!=(const rangeval &t) const {
            return !equals(t);
        }

        //
        // Helpful for debugging
        //

        void write(output &s) const;

    private:
        inline bool equals(const rangeval &t) const
        {
            if (the_type != t.the_type) return false;
            if (s_value != t.s_value) return false;
            if (range_special::NORMAL != s_value) {
                return true;
            }
            if (isReal()) {
                //
                // Check for a small relative difference
                // (unless we're zero)
                //
                if (0.0 == d_value) {
                    return 0.0 == t.d_value;
                }
                double delta;
                if (d_value > t.d_value) {
                    delta = d_value - t.d_value;
                } else {
                    delta = t.d_value - d_value;
                }
                delta /= d_value;
                return delta < 1e-8;
            } else {
                //
                // Integers or booleans; just check for equality
                //
                return l_value == t.l_value;
            }
        }



    private:
        union {
            long        l_value;
            double      d_value;
        };
        range_type      the_type;
        range_special   s_value;
};

#endif
