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

namespace MEDDLY {
    enum class range_type;
    enum class range_special;
    class rangeval;
};


/// Types of return values allowed for functions.
enum class MEDDLY::range_type {
    /// boolean-valued functions.
    BOOLEAN,
    /// integer-valued functions.
    INTEGER,
    /// real-valued functions.
    REAL
};


/// Special values allowed for functions.
enum class MEDDLY::range_special {
    /// Normal value, not special
    NORMAL,

    /// Positive infinity.
    PLUS_INFINITY

    // TBD: others? undefined?
};

/**
    Object to hold a generic return value.
    We keep track of the type and the value.
    Also we allow special values, like infinity.
*/
class MEDDLY::rangeval {
    public:
        /// Initialize as a boolean
        rangeval(bool v=true, range_type rt = range_type::BOOLEAN);

        /// Initialize as an integer
        rangeval(int v, range_type rt = range_type::INTEGER);

        /// Initialize as a real
        rangeval(float v, range_type rt = range_type::REAL);

        /// Initialize as special
        rangeval(range_special v, range_type rt);

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
            MEDDLY_DCASSERT(isNormal());
            MEDDLY_DCASSERT(isBoolean());
            return i_value;
        }
        inline operator int() const {
            MEDDLY_DCASSERT(isNormal());
            MEDDLY_DCASSERT(isInteger());
            return i_value;
        }
        inline operator float() const {
            MEDDLY_DCASSERT(isNormal());
            MEDDLY_DCASSERT(isReal());
            return f_value;
        }

    private:
        union {
            int         i_value;
            float       f_value;
        };
        range_type      the_type;
        range_special   s_value;
};

#endif
