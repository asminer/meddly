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
#include "comp_ne.h"
#include "apply_base.h"

namespace MEDDLY {
    binary_list NEQ_cache;
};


// ******************************************************************
// *                                                                *
// *                       unequal_mdd  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class unequal_mdd : public generic_binary_mdd {
  public:
    unequal_mdd(forest* arg1, forest* arg2, forest* res)
      : generic_binary_mdd(NEQ_cache, arg1, arg2, res)
      {
        operationCommutes();
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, SET);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
      }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool unequal_mdd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av != bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                       unequal_mxd  class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class unequal_mxd : public generic_binbylevel_mxd {
  public:
    unequal_mxd(forest* arg1, forest* arg2, forest* res)
      : generic_binbylevel_mxd(NEQ_cache, arg1, arg2, res)
      {
        operationCommutes();
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, RELATION);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
      }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool unequal_mxd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av != bv );
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

MEDDLY::binary_operation* MEDDLY::NOT_EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  NEQ_cache.find(a, b, c);
    if (bop) {
        return bop;
    }
    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        bool use_reals = (
            a->getRangeType() == range_type::REAL ||
            b->getRangeType() == range_type::REAL
        );
        if (use_reals) {
            if (c->isForRelations())
                return NEQ_cache.add(new unequal_mxd<float>(a, b, c));
            else
                return NEQ_cache.add(new unequal_mdd<float>(a, b, c));
        } else {
            if (c->isForRelations())
                return NEQ_cache.add(new unequal_mxd<long>(a, b, c));
            else
                return NEQ_cache.add(new unequal_mdd<long>(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::NOT_EQUAL_init()
{
    NEQ_cache.reset("Unequal");
}

void MEDDLY::NOT_EQUAL_done()
{
    MEDDLY_DCASSERT(NEQ_cache.isEmpty());
}

