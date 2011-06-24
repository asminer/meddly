
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
#include "comp_eq.h"
#include "../compute_table.h"
#include "apply_base.h"

namespace MEDDLY {
  class equal_mdd;
  class equal_mxd;
  class equal_evtimes;

  class equal_opname;
};


// ******************************************************************
// *                                                                *
// *                        equal_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::equal_mdd : public generic_binary_mdd {
  public:
    equal_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::equal_mdd::equal_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::equal_mdd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    int term = (a==b) ? 1 : 0;
    if (resF->getRangeType() == forest::INTEGER) {
      c = resF->getTerminalNode(term);
    } else {
      DCASSERT(resF->getRangeType() == forest::REAL);
      c = resF->getTerminalNode(float(term));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                        equal_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::equal_mxd : public generic_binbylevel_mxd {
  public:
    equal_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::equal_mxd::equal_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binbylevel_mxd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::equal_mxd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    int term = (a==b) ? 1 : 0;
    if (resF->getRangeType() == forest::INTEGER) {
      c = resF->getTerminalNode(term);
    } else {
      DCASSERT(resF->getRangeType() == forest::REAL);
      c = resF->getTerminalNode(float(term));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                      equal_evtimes  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::equal_evtimes : public generic_binary_evtimes {
  public:
    equal_evtimes(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(float av, int a, float bv, int b, 
      float &cv, int& c);
};

MEDDLY::equal_evtimes::equal_evtimes(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_evtimes(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::equal_evtimes
::checkTerminals(float aev, int a, float bev, int b, float &cev, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    if (a == b && ((aev == bev) || (isNan(aev) && isNan(bev)))) {
      c = -1;
      cev = 1.0;
    } else {
      c = 0;
      cev = NAN;
    }
    return true;
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                       equal_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::equal_opname : public binary_opname {
  public:
    equal_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::equal_opname::equal_opname()
 : binary_opname("Equal")
{
}

MEDDLY::binary_operation* 
MEDDLY::equal_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
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
      return new equal_mxd(this, a1, a2, r);
    else
      return new equal_mdd(this, a1, a2, r);
  }

  if (r->getEdgeLabeling() == forest::EVTIMES)
    return new equal_evtimes(this, a1, a2, r);

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeEQ(const settings &s)
{
  return new equal_opname;
}

