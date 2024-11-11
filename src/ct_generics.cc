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

#include "ct_generics.h"

// #define DEBUG_DESTROY

// **********************************************************************
// *                                                                    *
// *                         ct_object  methods                         *
// *                                                                    *
// **********************************************************************

MEDDLY::ct_object::ct_object(opnd_type t)
{
    type = t;
}

MEDDLY::ct_object::~ct_object()
{
}

void MEDDLY::ct_object::show(output &s) const
{
    s.put_hex((unsigned long) this);
    s.put(' ');
    s.put('G');
}


// **********************************************************************
// *                                                                    *
// *                         cached_mpz methods                         *
// *                                                                    *
// **********************************************************************

#ifdef HAVE_LIBGMP

MEDDLY::cached_mpz::cached_mpz(const mpz_ptr v) : ct_object(opnd_type::HUGEINT)
{
#ifdef DEBUG_DESTROY
    std::cerr << "created cached_mpz\n";
#endif
    mpz_init_set(value, v);
}

MEDDLY::cached_mpz::~cached_mpz()
{
    mpz_clear(value);
#ifdef DEBUG_DESTROY
    std::cerr << "destroyed cached_mpz\n";
#endif
}

#endif
