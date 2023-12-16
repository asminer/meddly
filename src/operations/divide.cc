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
#include "divide.h"
#include "apply_base.h"

namespace MEDDLY {
  class divide_opname;
};


// ******************************************************************
// *                                                                *
// *                        divide_mdd class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename REAL>
class divide_mdd : public generic_binary_mdd {
  public:
    divide_mdd(binary_opname* opcode, forest* arg1,
      forest* arg2, forest* res)
      : generic_binary_mdd(opcode, arg1, arg2, res) { }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename REAL>
bool divide_mdd<REAL>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b))
  {
    REAL av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    if (0 == bv) throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
    c = resF->handleForValue( av / bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                        divide_mxd class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename REAL>
class divide_mxd : public generic_binbylevel_mxd {
  public:
    divide_mxd(binary_opname* opcode, forest* arg1,
      forest* arg2, forest* res)
      : generic_binbylevel_mxd(opcode, arg1, arg2, res) { }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename REAL>
bool divide_mxd<REAL>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b))
  {
    REAL av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    if (0 == bv) throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
    c = resF->handleForValue( av / bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                      divide_opname  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::divide_opname : public binary_opname {
  public:
    divide_opname();
    virtual binary_operation* buildOperation(forest* a1,
      forest* a2, forest* r);
};

MEDDLY::divide_opname::divide_opname()
 : binary_opname("Divide")
{
}

MEDDLY::binary_operation*
MEDDLY::divide_opname::buildOperation(forest* a1, forest* a2,
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
    (a2->getEdgeLabeling() != r->getEdgeLabeling())
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (r->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    switch (r->getRangeType()) {

      case range_type::INTEGER:
          if (r->isForRelations())
            return new divide_mxd<int>(this, a1, a2, r);
          else
            return new divide_mdd<int>(this, a1, a2, r);

      case range_type::REAL:
          if (r->isForRelations())
            return new divide_mxd<float>(this, a1, a2, r);
          else
            return new divide_mdd<float>(this, a1, a2, r);

      default:
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

  }

  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeDivide()
{
  return new divide_opname;
}

