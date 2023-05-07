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
#include "comp_ge.h"
#include "apply_base.h"

namespace MEDDLY {
  class moreequal_opname;
};


// ******************************************************************
// *                                                                *
// *                      moreequal_mdd  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class moreequal_mdd : public generic_binary_mdd {
  public:
    moreequal_mdd(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res)
      : generic_binary_mdd(opcode, arg1, arg2, res) { }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool moreequal_mdd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av >= bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                      moreequal_mxd  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class moreequal_mxd : public generic_binbylevel_mxd {
  public:
    moreequal_mxd(binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res)
      : generic_binbylevel_mxd(opcode, arg1, arg2, res) { }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool moreequal_mxd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av >= bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                     moreequal_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::moreequal_opname : public binary_opname {
  public:
    moreequal_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r);
};

MEDDLY::moreequal_opname::moreequal_opname()
 : binary_opname("Unequal")
{
}

MEDDLY::binary_operation*
MEDDLY::moreequal_opname::buildOperation(expert_forest* a1, expert_forest* a2,
  expert_forest* r)
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
    (a2->getEdgeLabeling() != r->getEdgeLabeling())
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  bool use_reals = (
    a1->getRangeType() == range_type::REAL || a2->getRangeType() == range_type::REAL
  );
  if (r->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    if (use_reals) {
      if (r->isForRelations())
        return new moreequal_mxd<float>(this, a1, a2, r);
      else
        return new moreequal_mdd<float>(this, a1, a2, r);
    } else {
      if (r->isForRelations())
        return new moreequal_mxd<int>(this, a1, a2, r);
      else
        return new moreequal_mdd<int>(this, a1, a2, r);
    }
  }

  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeGE()
{
  return new moreequal_opname;
}

