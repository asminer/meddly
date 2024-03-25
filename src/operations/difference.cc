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
#include "difference.h"
#include "apply_base.h"

namespace MEDDLY {
  class diffr_mdd;
  class diffr_mxd;

  class diffr_opname;
};

// ******************************************************************
// *                                                                *
// *                        diffr_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mdd : public generic_binary_mdd {
  public:
    diffr_mdd(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::diffr_mdd::diffr_mdd(binary_list& opcode,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
  //  difference does NOT commute
}

bool MEDDLY::diffr_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0 || b == -1) {
    c = 0;
    return true;
  }
  if (a == -1 && b == 0) {
    c = -1;
    return true;
  }
  if (a == b) {
    if (arg1F == arg2F) {
      c = 0;
      return true;
    } else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } else {
      return false;
    }
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                        diffr_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mxd : public generic_binary_mxd {
  public:
    diffr_mxd(binary_list& opcode, forest* arg1,
      forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::diffr_mxd::diffr_mxd(binary_list& opcode,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(opcode, arg1, arg2, res)
{
  //  difference does NOT commute
}

bool MEDDLY::diffr_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0) {
    c = 0;
    return true;
  }
  if (a == -1 && b == 0) {
    c = -1;
    return true;
  }
  if (a == b) {
    if (arg1F == arg2F) {
      c = 0;
      return true;
    } else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } else {
      return false;
    }
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       diffr_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_opname : public binary_opname {
  public:
    diffr_opname();
    virtual binary_operation* buildOperation(binary_list &c, forest* a1,
      forest* a2, forest* r);
};

MEDDLY::diffr_opname::diffr_opname()
 : binary_opname("Difference")
{
}

MEDDLY::binary_operation*
MEDDLY::diffr_opname::buildOperation(binary_list &c, forest* a1, forest* a2,
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
    (a1->getRangeType() != r->getRangeType()) ||
    (a2->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (r->getRangeType() != range_type::BOOLEAN)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (r->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    if (r->isForRelations())
      return new diffr_mxd(c, a1, a2, r);
    else
      return new diffr_mdd(c, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeDifference()
{
  return new diffr_opname;
}

