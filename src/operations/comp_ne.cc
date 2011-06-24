
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
#include "comp_ne.h"
#include "../compute_table.h"
#include "apply_base.h"

namespace MEDDLY {
  class unequal_mdd;
  class unequal_mxd;

  class unequal_opname;
};


// ******************************************************************
// *                                                                *
// *                       unequal_mdd  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::unequal_mdd : public generic_binary_mdd {
  public:
    unequal_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::unequal_mdd::unequal_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::unequal_mdd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    int term = (a!=b) ? 1 : 0;
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
// *                       unequal_mxd  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::unequal_mxd : public generic_binbylevel_mxd {
  public:
    unequal_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::unequal_mxd::unequal_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binbylevel_mxd(opcode, arg1, arg2, res)
{
  operationCommutes();
}

bool MEDDLY::unequal_mxd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    int term = (a!=b) ? 1 : 0;
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
// *                      unequal_opname class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::unequal_opname : public binary_opname {
  public:
    unequal_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::unequal_opname::unequal_opname()
 : binary_opname("Unequal")
{
}

MEDDLY::binary_operation* 
MEDDLY::unequal_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
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
      return new unequal_mxd(this, a1, a2, r);
    else
      return new unequal_mdd(this, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeNE(const settings &s)
{
  return new unequal_opname;
}

