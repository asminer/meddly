
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "plus.h"
#include "apply_base.h"

namespace MEDDLY {
  class plus_mdd;
  class plus_mxd;
  class plus_evplus;
  class plus_evtimes;

  class plus_opname;
};

// ******************************************************************
// *                                                                *
// *                         plus_mdd class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_mdd : public generic_binary_mdd {
  public:
    plus_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::plus_mdd::plus_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::plus_mdd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == forest::INTEGER) {
      c = resF->getTerminalNode(arg1F->getInteger(a) + arg2F->getInteger(b));
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == forest::REAL);
      c = resF->getTerminalNode(arg1F->getReal(a) + arg2F->getReal(b));
    }
    return true;
  }
  if (0==a) {
    if (arg2F == resF) {
      c = arg2F->linkNode(b);
      return true;
    }
    return false;
  }
  if (0==b) {
    if (arg1F == resF) {
      c = arg1F->linkNode(a);
      return true;
    }
    return false;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                         plus_mxd class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_mxd : public generic_binary_mxd {
  public:
    plus_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::plus_mxd::plus_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mxd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::plus_mxd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == forest::INTEGER) {
      c = resF->getTerminalNode(arg1F->getInteger(a) + arg2F->getInteger(b));
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == forest::REAL);
      c = resF->getTerminalNode(arg1F->getReal(a) + arg2F->getReal(b));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       plus_evplus  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_evplus : public generic_binary_evplus {
  public:
    plus_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int aev, int a, int bev, int b, 
      int& cev, int& c);
};

MEDDLY::plus_evplus::plus_evplus(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_evplus(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::plus_evplus::checkTerminals(int aev, int a, int bev, int b,
  int& cev, int& c)
{
  if (a == -1 && b == -1) {
    c = -1; cev = aev + bev;
    return true;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                       plus_evtimes class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_evtimes : public generic_binary_evtimes {
  public:
    plus_evtimes(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(float aev, int a, float bev, int b, 
      float& cev, int& c);
};

MEDDLY::plus_evtimes::plus_evtimes(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_evtimes(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::plus_evtimes::checkTerminals(float aev, int a, 
  float bev, int b, float& cev, int& c)
{
  if (a == -1 && b == -1) {
    c = -1; cev = aev + bev;
    return true;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                       plus_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_opname : public binary_opname {
  public:
    plus_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::plus_opname::plus_opname()
 : binary_opname("Plus")
{
}

MEDDLY::binary_operation* 
MEDDLY::plus_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    (a1->isForRelations() != r->isForRelations()) ||
    (a2->isForRelations() != r->isForRelations()) ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (r->getRangeType() == forest::BOOLEAN)
  )
    throw error(error::TYPE_MISMATCH);

  if (r->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    if (r->isForRelations())
      return new plus_mxd(this, a1, a2, r);
    else
      return new plus_mdd(this, a1, a2, r);
  }

  if (r->getEdgeLabeling() == forest::EVPLUS)
    return new plus_evplus(this, a1, a2, r);

  if (r->getEdgeLabeling() == forest::EVTIMES)
    return new plus_evtimes(this, a1, a2, r);

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializePlus(const settings &s)
{
  return new plus_opname;
}

