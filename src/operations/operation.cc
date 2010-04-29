
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
#include "../meddly_expert.h"

operation::operation()
{ }


operation::~operation() {}


// Defaults

compute_manager::error 
operation::compute(op_info* cc, dd_edge** operands)
{
  return compute_manager::TYPE_MISMATCH;
}

compute_manager::error 
operation::compute(op_info* cc, const dd_edge& a, dd_edge& b)
{
  return compute_manager::TYPE_MISMATCH;
}

compute_manager::error 
operation::compute(op_info* cc, const dd_edge& a, long& b)
{
  return compute_manager::TYPE_MISMATCH;
}

compute_manager::error 
operation::compute(op_info* cc, const dd_edge& a, double& b)
{
  return compute_manager::TYPE_MISMATCH;
}


compute_manager::error 
operation::compute(op_info* cc, const dd_edge& a, mpz_t &b)
{
  return compute_manager::TYPE_MISMATCH;
}


compute_manager::error 
operation::compute(op_info* cc, const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  return compute_manager::TYPE_MISMATCH;
}
