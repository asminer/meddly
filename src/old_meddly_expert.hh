
/*
 Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
 Copyright (C) 2009, Iowa State University Research Foundation, Inc.

 This library is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published
 by the Free Software Foundation, either version 3 of the License, orf
 (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library.  If not, see <http://www.gnu.org/licenses/>.
 */

/*! \file meddly_expert.hh

 Implementation details for interface in meddly_expert.h.
 */

#ifndef MEDDLY_EXPERT_HH
#define MEDDLY_EXPERT_HH

#include "defines.h"
#include "io.h"
#include "memstats.h"

#include "forest.h"

// ******************************************************************
// *                                                                *
// *                   inlined  operation methods                   *
// *                                                                *
// ******************************************************************

#include "compute_table.h"
#include "opname.h"

inline void
MEDDLY::operation::registerInForest(MEDDLY::forest* f)
{
  if (f)
    f->registerOperation(this);
}

inline void
MEDDLY::operation::unregisterInForest(MEDDLY::forest* f)
{
  if (f)
    f->unregisterOperation(this);
}

inline void
MEDDLY::operation::registerEntryType(unsigned slot, ct_entry_type* et)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, slot, num_etids);
  MEDDLY_DCASSERT(etype);
  MEDDLY_DCASSERT(0==etype[slot]);
  etype[slot] = et;
  compute_table::registerEntryType(first_etid + slot, et);
}

inline bool
MEDDLY::operation::isMarkedForDeletion() const
{
  return is_marked_for_deletion;
}

inline void
MEDDLY::operation::setNext(operation* n)
{
  next = n;
}

inline MEDDLY::operation*
MEDDLY::operation::getNext()
{
  return next;
}

inline bool
MEDDLY::operation::usesMonolithicComputeTable()
{
  return Monolithic_CT;
}

inline unsigned
MEDDLY::operation::getIndex() const
{
  return oplist_index;
}

inline MEDDLY::operation*
MEDDLY::operation::getOpWithIndex(unsigned i)
{
  MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, list_size);
  return op_list[i];
}

inline unsigned
MEDDLY::operation::getOpListSize()
{
  return list_size;
}

inline void
MEDDLY::operation::setFirstETid(unsigned slot)
{
  first_etid = slot;
}

inline unsigned
MEDDLY::operation::getFirstETid() const
{
  return first_etid;
}

inline unsigned
MEDDLY::operation::getNumETids() const
{
  return num_etids;
}

inline const char*
MEDDLY::operation::getName() const
{
  return theOpName->getName();
}

inline const MEDDLY::opname*
MEDDLY::operation::getOpName() const
{
  return theOpName;
}

// ******************************************************************
// *                                                                *
// *                inlined  unary_operation methods                *
// *                                                                *
// ******************************************************************

inline bool
MEDDLY::unary_operation::matches(const MEDDLY::expert_forest* arg,
    const MEDDLY::expert_forest* res) const
{
  return (arg == argF && res == resF);
}

inline bool
MEDDLY::unary_operation::matches(const MEDDLY::expert_forest* arg, opnd_type res) const
{
  return (arg == argF && resultType == res);
}

inline bool
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

inline void
MEDDLY::unary_operation::compute(const dd_edge &arg, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(arg, res, true);
}

inline void
MEDDLY::unary_operation::computeTemp(const dd_edge &arg, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(arg, res, false);
}

// ******************************************************************
// *                                                                *
// *                inlined binary_operation methods                *
// *                                                                *
// ******************************************************************


inline bool
MEDDLY::binary_operation::matches(const MEDDLY::expert_forest* arg1,
    const MEDDLY::expert_forest* arg2, const MEDDLY::expert_forest* res) const
{
  return (arg1 == arg1F && arg2 == arg2F && res == resF);
}

inline void
MEDDLY::binary_operation::operationCommutes()
{
  can_commute = (arg1F == arg2F);
}

inline bool
MEDDLY::binary_operation::checkForestCompatibility() const
{
  auto o1 = arg1F->variableOrder();
  auto o2 = arg2F->variableOrder();
  auto o3 = resF->variableOrder();
  return o1->is_compatible_with(*o2) && o1->is_compatible_with(*o3);
}

inline void
MEDDLY::binary_operation::compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(ar1, ar2, res, true);
}

inline void
MEDDLY::binary_operation::computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(ar1, ar2, res, false);
}

// ******************************************************************
// *                                                                *
// *                       inlined  functions                       *
// *                                                                *
// ******************************************************************


inline MEDDLY::unary_operation*
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg,
    const dd_edge& res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg.getForest(),
      (MEDDLY::expert_forest*) res.getForest());
}

inline MEDDLY::unary_operation*
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg,
    opnd_type res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg.getForest(), res);
}

inline MEDDLY::binary_operation*
MEDDLY::getOperation(const binary_opname* code, const dd_edge& arg1,
    const dd_edge& arg2, const dd_edge& res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg1.getForest(),
      (MEDDLY::expert_forest*) arg2.getForest(), (MEDDLY::expert_forest*) res.getForest());
}

#endif



