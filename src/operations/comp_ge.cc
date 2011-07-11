
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
#include "comp_ge.h"
#include "apply_base.h"

namespace MEDDLY {
  class moreequal_mdd;
  class moreequal_mxd;

  class moreequal_opname;
};


// ******************************************************************
// *                                                                *
// *                      moreequal_mdd  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::moreequal_mdd : public generic_binary_mdd {
  public:
    moreequal_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::moreequal_mdd::moreequal_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
}

bool MEDDLY::moreequal_mdd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == forest::INTEGER) {
      bool lt = arg1F->getInteger(a) >= arg2F->getInteger(b);
      c = resF->getTerminalNode(lt ? 1 : 0);
    } else {
      DCASSERT(resF->getRangeType() == forest::REAL);
      bool lt = arg1F->getReal(a) >= arg2F->getReal(b);
      c = resF->getTerminalNode(float(lt ? 1 : 0));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                      moreequal_mxd  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::moreequal_mxd : public generic_binbylevel_mxd {
  public:
    moreequal_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(int a, int b, int& c);
};

MEDDLY::moreequal_mxd::moreequal_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binbylevel_mxd(opcode, arg1, arg2, res)
{
}

bool MEDDLY::moreequal_mxd::checkTerminals(int a, int b, int& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == forest::INTEGER) {
      bool lt = arg1F->getInteger(a) >= arg2F->getInteger(b);
      c = resF->getTerminalNode(lt ? 1 : 0);
    } else {
      DCASSERT(resF->getRangeType() == forest::REAL);
      bool lt = arg1F->getReal(a) >= arg2F->getReal(b);
      c = resF->getTerminalNode(float(lt ? 1 : 0));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                     moreequal_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::moreequal_opname : public binary_opname {
  public:
    moreequal_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::moreequal_opname::moreequal_opname()
 : binary_opname("Unequal")
{
}

MEDDLY::binary_operation* 
MEDDLY::moreequal_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
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
      return new moreequal_mxd(this, a1, a2, r);
    else
      return new moreequal_mdd(this, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeGE(const settings &s)
{
  return new moreequal_opname;
}

