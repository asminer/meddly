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
    binary_list DIVIDE_cache;
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
    divide_mdd(forest* arg1, forest* arg2, forest* res)
      : generic_binary_mdd(DIVIDE_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, SET);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

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
    divide_mxd(forest* arg1, forest* arg2, forest* res)
      : generic_binbylevel_mxd(DIVIDE_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, RELATION);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::DIVIDE(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;

    binary_operation* bop =  DIVIDE_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        switch (c->getRangeType()) {

            case range_type::INTEGER:
                if (c->isForRelations())
                    return DIVIDE_cache.add(new divide_mxd<long>(a, b, c));
                else
                    return DIVIDE_cache.add(new divide_mdd<long>(a, b, c));

            case range_type::REAL:
                if (c->isForRelations())
                    return DIVIDE_cache.add(new divide_mxd<float>(a, b, c));
                else
                    return DIVIDE_cache.add(new divide_mdd<float>(a, b, c));

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::DIVIDE_init()
{
    DIVIDE_cache.reset("Minus");
}

void MEDDLY::DIVIDE_done()
{
    MEDDLY_DCASSERT(DIVIDE_cache.isEmpty());
}


