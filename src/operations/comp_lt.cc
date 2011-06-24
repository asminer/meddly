
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
#include "comp_lt.h"
#include "../compute_table.h"
#include "apply_base.h"

namespace MEDDLY {
  class lessthan_mdd;
  class lessthan_mxd;

  class lessthan_opname;
};


// ******************************************************************
// *                                                                *
// *                       lessthan_mdd class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::lessthan_mdd : public generic_binary_mdd {
  public:
    lessthan_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::lessthan_mdd::lessthan_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
}

bool MEDDLY::lessthan_mdd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == forest::INTEGER) {
      bool lt = arg1F->getInteger(a) < arg2F->getInteger(b);
      c = resF->getTerminalNode(lt ? 1 : 0);
    } else {
      DCASSERT(resF->getRangeType() == forest::REAL);
      bool lt = arg1F->getReal(a) < arg2F->getReal(b);
      c = resF->getTerminalNode(float(lt ? 1 : 0));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       lessthan_mxd class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::lessthan_mxd : public generic_binbylevel_mxd {
  public:
    lessthan_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::lessthan_mxd::lessthan_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binbylevel_mxd(opcode, arg1, arg2, res)
{
}

bool MEDDLY::lessthan_mxd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == forest::INTEGER) {
      bool lt = arg1F->getInteger(a) < arg2F->getInteger(b);
      c = resF->getTerminalNode(lt ? 1 : 0);
    } else {
      DCASSERT(resF->getRangeType() == forest::REAL);
      bool lt = arg1F->getReal(a) < arg2F->getReal(b);
      c = resF->getTerminalNode(float(lt ? 1 : 0));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                     lessthan_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::lessthan_opname : public binary_opname {
  public:
    lessthan_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::lessthan_opname::lessthan_opname()
 : binary_opname("LessThan")
{
}

MEDDLY::binary_operation* 
MEDDLY::lessthan_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
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
      return new lessthan_mxd(this, a1, a2, r);
    else
      return new lessthan_mdd(this, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeLT(const settings &s)
{
  return new lessthan_opname;
}

