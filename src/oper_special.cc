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

#include "oper_special.h"
#include "error.h"

// ******************************************************************
// *                 specialized_operation  methods                 *
// ******************************************************************

MEDDLY::
specialized_operation::
specialized_operation(specialized_opname* op, unsigned et_slots)
 : operation(op, et_slots)
{
}

MEDDLY::specialized_operation::~specialized_operation()
{
}

void MEDDLY::specialized_operation::compute(const dd_edge &arg, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(const dd_edge &arg, bool &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(const dd_edge &ar1,
  const dd_edge &ar2, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(double* y, const double* x)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::specialized_operation::compute(const dd_edge &ar1,
  const dd_edge &ar2, const dd_edge &ar3, dd_edge &res)
{
  throw error(error::TYPE_MISMATCH);
}

// ******************************************************************
// *                                                                *
// *                      front-end  functions                      *
// *                                                                *
// ******************************************************************

void MEDDLY::destroyOperation(MEDDLY::specialized_operation* &op)
{
  if (!op) return;
  if (!op->isMarkedForDeletion()) {
    op->markForDeletion();
    operation::removeStalesFromMonolithic();
  }
  delete op;
  op = nullptr;
}

