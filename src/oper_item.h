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

#ifndef MEDDLY_OPER_ITEM_H
#define MEDDLY_OPER_ITEM_H

#include "defines.h"
#include "oper.h"

#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif

/*
 * This is for items other than DD edges that are either arguments,
 * or results, of operations.
 */

namespace MEDDLY {
    class oper_item;
};

/**
    Items (with types) for operations.
 */
class MEDDLY::oper_item {
    public:
        oper_item(opnd_type t = opnd_type::FOREST);
        oper_item(long v);
        oper_item(double v);
#ifdef HAVE_LIBGMP
        oper_item(mpz_ptr v);
#endif
        //
        // Delayed initialization
        //

        inline void init(long v) {
            ASSERT(__FILE__, __LINE__, mytype == opnd_type::FOREST);
            mytype = opnd_type::INTEGER;
            the_long = v;
        }

        inline void init(double v) {
            ASSERT(__FILE__, __LINE__, mytype == opnd_type::FOREST);
            mytype = opnd_type::REAL;
            the_double = v;
        }

#ifdef HAVE_LIBGMP
        inline void init(mpz_ptr v) {
            ASSERT(__FILE__, __LINE__, mytype == opnd_type::FOREST);
            mytype = opnd_type::HUGEINT;
            the_mpz = v;
        }
#endif

        //
        // getters for the type
        //

        inline opnd_type getType() const {
            return mytype;
        }
        inline bool hasType(opnd_type t) const {
            return t == mytype;
        }

        //
        // getters for the value
        //

        inline long getInteger() const {
            ASSERT(__FILE__, __LINE__, hasType(opnd_type::INTEGER));
            return the_long;
        }
        inline double getReal() const {
            ASSERT(__FILE__, __LINE__, hasType(opnd_type::REAL));
            return the_double;
        }
#ifdef HAVE_LIBGMP
        inline mpz_ptr getHugeint() const {
            ASSERT(__FILE__, __LINE__, hasType(opnd_type::HUGEINT));
            return the_mpz;
        }
#endif

        //
        // value access
        //

        inline long& integer() {
            ASSERT(__FILE__, __LINE__, hasType(opnd_type::INTEGER));
            return the_long;
        }

        inline double& real() {
            ASSERT(__FILE__, __LINE__, hasType(opnd_type::REAL));
            return the_double;
        }
#ifdef HAVE_LIBGMP
        inline mpz_ptr hugeint() {
            ASSERT(__FILE__, __LINE__, hasType(opnd_type::HUGEINT));
            return the_mpz;
        }
#endif

    //
    // TBD: should we include CT stuff here,
    // i.e., things from ct_object?
    // number of CT slots required
    // how to stuff this into a CT
    // how to extract this from a CT
    // should it be hashed if part of a CT entry?
    //

    private:
        union {
            long        the_long;
            double      the_double;
#ifdef HAVE_LIBGMP
            mpz_ptr     the_mpz;
#endif
        };

        opnd_type mytype;
};

#endif
