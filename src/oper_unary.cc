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

#include "oper_unary.h"
#include "error.h"
#include "forest.h"
#include "dd_edge.h"

// ******************************************************************
// *                    unary_operation  methods                    *
// ******************************************************************

MEDDLY::unary_operation::unary_operation(const unary_opname* code,
  unsigned et_slots, expert_forest* arg, expert_forest* res)
: operation(code, et_slots)
{
  argF = arg;
  resultType = opnd_type::FOREST;
  resF = res;

  registerInForest(argF);
  registerInForest(resF);
}

MEDDLY::unary_operation::unary_operation(const unary_opname* code,
  unsigned et_slots, expert_forest* arg, opnd_type res)
: operation(code, et_slots)
{
  argF = arg;
  resultType = res;
  resF = 0;

  registerInForest(argF);
}

MEDDLY::unary_operation::~unary_operation()
{
  unregisterInForest(argF);
  unregisterInForest(resF);
}


inline bool
MEDDLY::unary_operation::matches(const MEDDLY::dd_edge &arg,
        const MEDDLY::dd_edge &res) const
{
  return (arg.getForest() == argF && res.getForest() == resF);
}

inline bool
MEDDLY::unary_operation::matches(const MEDDLY::dd_edge &arg, opnd_type res)
        const
{
  return (arg.getForest() == argF && resultType == res);
}



void MEDDLY::unary_operation::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, long &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, double &res)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::unary_operation::compute(const dd_edge &arg, ct_object &c)
{
  throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void
MEDDLY::unary_operation::compute(const dd_edge &arg, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(arg, res, true);
}

void
MEDDLY::unary_operation::computeTemp(const dd_edge &arg, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(arg, res, false);
}

bool
MEDDLY::unary_operation::checkForestCompatibility() const
{
  if (resultType == opnd_type::FOREST) {
    auto o1 = argF->variableOrder();
    auto o2 = resF->variableOrder();
    return o1->is_compatible_with(*o2);
  }
  else {
    return true;
  }
}


// ******************************************************************
// *                                                                *
// *                      front-end  functions                      *
// *                                                                *
// ******************************************************************

/*

MEDDLY::unary_operation*
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg,
    const dd_edge& res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg.getForest(),
      (MEDDLY::expert_forest*) res.getForest());
}

MEDDLY::unary_operation*
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg,
    opnd_type res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg.getForest(), res);
}

*/
