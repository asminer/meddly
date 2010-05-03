
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

// ---------------------- operation ------------------------

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
operation::compute(op_info* cc, const dd_edge& a, ct_object &b)
{
  return compute_manager::TYPE_MISMATCH;
}


compute_manager::error 
operation::compute(op_info* cc, const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  return compute_manager::TYPE_MISMATCH;
}

// ---------------------- op_info ------------------------


op_info::op_info() 
 : op(0), p(0), nParams(0), cc(0) 
{
}

op_info::op_info(operation *oper, op_param* plist, int n,
  compute_cache* cache)
: op(oper), p(0), nParams(n), cc(cache)
{
  p = (op_param *) malloc(nParams * sizeof(op_param));
  for (int i = 0; i < nParams; ++i) p[i] = plist[i];
}

op_info::op_info(const op_info& a)
: op(a.op), p(0), nParams(a.nParams), cc(a.cc)
{
  p = (op_param *) malloc(nParams * sizeof(op_param));
  for (int i = 0; i < nParams; ++i) p[i] = a.p[i];
}

op_info::~op_info()
{
  if (p != 0) free(p);
}

op_info& op_info::operator=(const op_info& a)
{
  if (this == &a) return *this;
  if (nParams != a.nParams) {
    nParams = a.nParams;
    p = (op_param *) realloc(p, nParams * sizeof(op_param));
  }
  for (int i = 0; i < nParams; ++i) p[i] = a.p[i];
  op = a.op;
  cc = a.cc;
  return *this;
}



// ---------------------- ct_object ------------------------

ct_object::ct_object()
{
}

ct_object::~ct_object()
{
}

// ---------------------- builtin_op_key ------------------------

builtin_op_key::builtin_op_key(compute_manager::op_code op,
    const op_param* const p, int n)
: opCode(op), plist(0), nParams(n)
{
  plist = (op_param*) malloc(nParams * sizeof(op_param));
  for (int i = 0; i < nParams; ++i) plist[i] = p[i];
}


builtin_op_key::builtin_op_key(const builtin_op_key& a)
: opCode(a.opCode), plist(0), nParams(a.nParams)
{
  plist = (op_param*) malloc(nParams * sizeof(op_param));
  for (int i = 0; i < nParams; ++i) plist[i] = a.plist[i];
}


builtin_op_key& builtin_op_key::operator=(const builtin_op_key& a)
{
  if (this == &a) return *this;
  if (nParams != a.nParams) {
    nParams = a.nParams;
    plist = (op_param*) realloc(plist, nParams * sizeof(op_param));
  }
  for (int i = 0; i < nParams; ++i) plist[i] = a.plist[i];
  opCode = a.opCode;
  return *this;
}

builtin_op_key::~builtin_op_key()
{
  free(plist);
}

// ---------------------- custom_op_key ------------------------

custom_op_key::custom_op_key(const operation* oper,
  const op_param* const p, int n)
: op(0), plist(0), nParams(n)
{
  op = const_cast<operation*>(oper);
  plist = (op_param*) malloc(nParams * sizeof(op_param));
  for (int i = 0; i < nParams; ++i) plist[i] = p[i];
}


custom_op_key::custom_op_key(const custom_op_key& a)
: op(a.op), plist(0), nParams(a.nParams)
{
  plist = (op_param*) malloc(nParams * sizeof(op_param));
  for (int i = 0; i < nParams; ++i) plist[i] = a.plist[i];
}


custom_op_key& custom_op_key::operator=(const custom_op_key& a)
{
  if (this == &a) return *this;
  if (nParams != a.nParams) {
    nParams = a.nParams;
    plist = (op_param*) realloc(plist, nParams * sizeof(op_param));
  }
  for (int i = 0; i < nParams; ++i) plist[i] = a.plist[i];
  op = a.op;
  return *this;
}


custom_op_key::~custom_op_key()
{
  free(plist);
}

