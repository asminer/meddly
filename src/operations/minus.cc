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
#include "minus.h"
#include "apply_base.h"

namespace MEDDLY {
    class minus_mdd;
    class minus_mxd;
    class minus_evplus;
    class minus_evtimes;

    binary_list MINUS_cache;
};

// ******************************************************************
// *                                                                *
// *                        minus_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_mdd : public generic_binary_mdd {
  public:
    minus_mdd(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::minus_mdd::minus_mdd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mdd(MINUS_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::minus_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    }
    return true;
  }
  if (0==b) {
    if (arg1F == resF) {
      c = arg1F->linkNode(a);
      return true;
    }
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                        minus_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_mxd : public generic_binary_mxd {
  public:
    minus_mxd(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::minus_mxd::minus_mxd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(MINUS_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::minus_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av - bv);
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       minus_evplus class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_evplus : public generic_binary_evplus {
  public:
    minus_evplus(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long& cev, node_handle& c);
};

MEDDLY::minus_evplus::minus_evplus(forest* arg1, forest* arg2, forest* res)
  : generic_binary_evplus(MINUS_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

bool MEDDLY::minus_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1; cev = aev - bev;
    return true;
  }
  if (0 == a || 0 == b) {
    c = 0;
    cev = 0;
    return true;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                      minus_evtimes  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::minus_evtimes : public generic_binary_evtimes {
  public:
    minus_evtimes(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(float aev, node_handle a, float bev, node_handle b,
      float& cev, node_handle& c);
};

MEDDLY::minus_evtimes::minus_evtimes(forest* arg1, forest* arg2, forest* res)
  : generic_binary_evtimes(MINUS_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVTIMES);
}

bool MEDDLY::minus_evtimes::checkTerminals(float aev, node_handle a,
  float bev, node_handle b, float& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1; cev = aev - bev;
    return true;
  }
  if (0 == a && 0 == b) {
    c = 0;
    cev = 0;
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::MINUS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;

    binary_operation* bop =  MINUS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations())
            return MINUS_cache.add(new minus_mxd(a, b, c));
        else
            return MINUS_cache.add(new minus_mdd(a, b, c));
    }

    if  (
            (a->getRangeType() != b->getRangeType()) ||
            (a->getRangeType() != c->getRangeType())
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS)
        return MINUS_cache.add(new minus_evplus(a, b, c));

    if (c->getEdgeLabeling() == edge_labeling::EVTIMES)
        return MINUS_cache.add(new minus_evtimes(a, b, c));

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::MINUS_init()
{
    MINUS_cache.reset("Minus");
}

void MEDDLY::MINUS_done()
{
    MEDDLY_DCASSERT(MINUS_cache.isEmpty());
}

