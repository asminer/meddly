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

    TBD: decide on philosophy for mpz_t values, maybe
        () It's always just a wrapper, we never "own" the mpz_t
           so use pointers (maybe mpz_ptr type?)
        () Use something like the old mpz_object but with reference counts
 */
class MEDDLY::oper_item {
    public:
        oper_item(opnd_type t = FOREST);
        oper_item(int v);
        oper_item(long v);
        oper_item(double v);
#ifdef HAVE_GMP
        oper_item(mpz_ptr v);
#endif

        ~oper_item();

        //
        // getters for the type
        //

        inline opnd_type getType() const {
            return type;
        }
        inline bool hasType(opnd_type t) const {
            return t == type;
        }

        //
        // getters for the value
        //

        inline void get(int &v) const {
            MEDDLY_DCASSERT(hasType(opnd_type::INTEGER));
            v = the_int;
        }
        inline void get(long &v) const {
            MEDDLY_DCASSERT(hasType(opnd_type::LONG));
            v = the_long;
        }
        inline void get(double &v) const {
            MEDDLY_DCASSERT(hasType(opnd_type::DOUBLE));
            v = the_double;
        }
#ifdef HAVE_GMP
        inline void get(mpz_t &v) const {
            MEDDLY_DCASSERT(hasType(opnd_type::HUGEINT));
            mpz_set(v, the_mpz);
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
            int         the_int;
            long        the_long;
            double      the_double;
#ifdef HAVE_GMP
            mpz_ptr     the_mpz;
#endif
        };

        opnd_type mytype;
};

#endif
