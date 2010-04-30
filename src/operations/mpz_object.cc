
// $Id$

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

#include "../defines.h"
#include "mpz_object.h"

#ifdef HAVE_LIBGMP 

mpz_object::mpz_object()
{
  mpz_init(value);
}

mpz_object::mpz_object(const mpz_t &v)
{
  mpz_init_set(value, v);
}

mpz_object::mpz_object(const mpz_object &x)
{
  mpz_init_set(value, x.value);
}


mpz_object::~mpz_object()
{
  mpz_clear(value);
}

op_param::type mpz_object::getType()
{
  return op_param::HUGEINT;
}

ct_object& MEDDLY_get_mpz_wrapper()
{
  static mpz_object foo;
  return foo;
}

void MEDDLY_unwrap(const ct_object &x, mpz_t &value)
{
  const mpz_object &mx = dynamic_cast <const mpz_object &> (x);
  mx.copyInto(value);
}

#endif

