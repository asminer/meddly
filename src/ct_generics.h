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

#ifndef MEDDLY_CT_GENERICS_H
#define MEDDLY_CT_GENERICS_H

#include "oper.h" // for  opnd_type

namespace MEDDLY {
    class ct_object;
#ifdef HAVE_LIBGMP
    class cached_mpz;
#endif
};

// ******************************************************************
// *                                                                *
// *                         ct_object class                        *
// *                                                                *
// ******************************************************************

/** Generic objects in compute tables.
    Used for things other than dd_edges and simple types.
*/
class MEDDLY::ct_object {
        opnd_type type;
    public:
        ct_object(opnd_type t);
        virtual ~ct_object();
        inline opnd_type getType() const {
            return type;
        }

        // Default behavior: show 'this' pointer
        virtual void show(output &s) const;
};


// ******************************************************************
// *                                                                *
// *                        cached_mpz  class                       *
// *                                                                *
// ******************************************************************

#ifdef HAVE_LIBGMP
#include <gmp.h>

class MEDDLY::cached_mpz : public ct_object {
        mpz_t value;
    public:
        cached_mpz(const mpz_ptr v);
        virtual ~cached_mpz();

        inline mpz_ptr mpz() { return value; }
};

#endif

#endif
