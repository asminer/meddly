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
#include "comp_eq.h"
#include "apply_base.h"

namespace MEDDLY {
    // TBD: this operation is strange
    class equal_evtimes;

    binary_list EQUAL_cache;
};


// ******************************************************************
// *                                                                *
// *                        equal_mdd  class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class equal_mdd : public generic_binary_mdd {
  public:
    equal_mdd(forest* arg1, forest* arg2, forest* res)
      : generic_binary_mdd(EQUAL_cache, arg1, arg2, res)
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
bool equal_mdd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == b && arg1F == arg2F) {
    c = resF->handleForValue(true);
    return true;
  }
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av == bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                        equal_mxd  class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class equal_mxd : public generic_binbylevel_mxd {
  public:
    equal_mxd(forest* arg1, forest* arg2, forest* res)
      : generic_binbylevel_mxd(EQUAL_cache, arg1, arg2, res)
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
bool equal_mxd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av == bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                      equal_evtimes  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::equal_evtimes : public generic_binary_evtimes {
  public:
    equal_evtimes(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(float av, node_handle a, float bv, node_handle b,
      float &cv, node_handle& c);
};

MEDDLY::equal_evtimes::equal_evtimes(forest* arg1, forest* arg2, forest* res)
  : generic_binary_evtimes(EQUAL_cache, arg1, arg2, res)
{
  operationCommutes();
  checkDomains(__FILE__, __LINE__);
  checkAllRelations(__FILE__, __LINE__, RELATION);
  checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVTIMES);
}

bool MEDDLY::equal_evtimes
::checkTerminals(float aev, node_handle a, float bev, node_handle b, float &cev, node_handle& c)
{
  if (arg1F->isTerminalNode(a) &&
      arg2F->isTerminalNode(b)) {
    if (a == b && ((aev == bev) || (isNan(aev) && isNan(bev)))) {
      c = -1;
      cev = 1.0;
    } else {
      c = 0;
      cev = Nan();
    }
    return true;
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  EQUAL_cache.find(a, b, c);
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
                return EQUAL_cache.add(new equal_mxd<float>(a, b, c));
            else
                return EQUAL_cache.add(new equal_mdd<float>(a, b, c));
        } else {
            if (c->isForRelations())
                return EQUAL_cache.add(new equal_mxd<long>(a, b, c));
            else
                return EQUAL_cache.add(new equal_mdd<long>(a, b, c));
        }
    }

    if (c->getEdgeLabeling() == edge_labeling::EVTIMES) {
        return EQUAL_cache.add(new equal_evtimes(a, b, c));
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::EQUAL_init()
{
    EQUAL_cache.reset("Equal");
}

void MEDDLY::EQUAL_done()
{
    MEDDLY_DCASSERT(EQUAL_cache.isEmpty());
}

