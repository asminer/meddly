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
#include "minus.h"
#include "apply_base.h"

namespace MEDDLY {
  class minus_mdd;
  class minus_mxd;
  class minus_evplus;
  class minus_evtimes;

  class minus_opname;
};

// ******************************************************************
// *                                                                *
// *                        minus_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_mdd : public generic_binary_mdd {
  public:
    minus_mdd(binary_opname* opcode, forest* arg1,
      forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::minus_mdd::minus_mdd(binary_opname* opcode,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
}

bool MEDDLY::minus_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    }
    return true;
  }
  if (0==b) {
    if (arg1F == resF) {
      c = arg1F->linkNode(a);
      return true;
    }
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                        minus_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_mxd : public generic_binary_mxd {
  public:
    minus_mxd(binary_opname* opcode, forest* arg1,
      forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::minus_mxd::minus_mxd(binary_opname* opcode,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(opcode, arg1, arg2, res)
{
}

bool MEDDLY::minus_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       minus_evplus class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_evplus : public generic_binary_evplus {
  public:
    minus_evplus(binary_opname* opcode, forest* arg1,
      forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long& cev, node_handle& c);
};

MEDDLY::minus_evplus::minus_evplus(binary_opname* opcode,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_evplus(opcode, arg1, arg2, res)
{
}

bool MEDDLY::minus_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1; cev = aev - bev;
    return true;
  }
  if (0 == a || 0 == b) {
    c = 0;
    cev = 0;
    return true;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                      minus_evtimes  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_evtimes : public generic_binary_evtimes {
  public:
    minus_evtimes(binary_opname* opcode, forest* arg1,
      forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(float aev, node_handle a, float bev, node_handle b,
      float& cev, node_handle& c);
};

MEDDLY::minus_evtimes::minus_evtimes(binary_opname* opcode,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_evtimes(opcode, arg1, arg2, res)
{
}

bool MEDDLY::minus_evtimes::checkTerminals(float aev, node_handle a,
  float bev, node_handle b, float& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1; cev = aev - bev;
    return true;
  }
  if (0 == a && 0 == b) {
    c = 0;
    cev = 0;
    return true;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                       minus_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_opname : public binary_opname {
  public:
    minus_opname();
    virtual binary_operation* buildOperation(forest* a1,
      forest* a2, forest* r);
};

MEDDLY::minus_opname::minus_opname()
 : binary_opname("Minus")
{
}

MEDDLY::binary_operation*
MEDDLY::minus_opname::buildOperation(forest* a1, forest* a2,
  forest* r)
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (
    (a1->getDomain() != r->getDomain()) ||
    (a2->getDomain() != r->getDomain())
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    (a1->isForRelations() != r->isForRelations()) ||
    (a2->isForRelations() != r->isForRelations()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (r->getRangeType() == range_type::BOOLEAN)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (r->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    if (r->isForRelations())
      return new minus_mxd(this, a1, a2, r);
    else
      return new minus_mdd(this, a1, a2, r);
  }

  if (
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType())
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (r->getEdgeLabeling() == edge_labeling::EVPLUS)
    return new minus_evplus(this, a1, a2, r);

  if (r->getEdgeLabeling() == edge_labeling::EVTIMES)
    return new minus_evtimes(this, a1, a2, r);

  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeMinus()
{
  return new minus_opname;
}

