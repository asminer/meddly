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
#include "multiply.h"
#include "apply_base.h"

namespace MEDDLY {
    class multiply_mdd;
    class multiply_mxd;
    class multiply_evplus;
    class multiply_evtimes;

    binary_list MULTIPLY_cache;
};

// ******************************************************************
// *                                                                *
// *                       multiply_mdd class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::multiply_mdd : public generic_binary_mdd {
  public:
    multiply_mdd(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::multiply_mdd::multiply_mdd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mdd(MULTIPLY_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::multiply_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0 || b == 0) {
    c = 0;
    return true;
  }
  if (arg1F->isTerminalNode(a)) {
    if (arg2F->isTerminalNode(b)) {
      if (resF->getRangeType() == range_type::INTEGER) {
        int av, bv;
        arg1F->getValueFromHandle(a, av);
        arg2F->getValueFromHandle(b, bv);
        c = resF->handleForValue(av * bv);
      } else {
        MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
        float av, bv;
        arg1F->getValueFromHandle(a, av);
        arg2F->getValueFromHandle(b, bv);
        c = resF->handleForValue(av * bv);
      }
      return true;
    }
    if (arg2F != resF) return false;
    if (resF->getRangeType() == range_type::INTEGER) {
      if (1==arg1F->getIntegerFromHandle(a)) {
        c = arg2F->linkNode(b);
        return true;
      }
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      if (1.0==arg1F->getRealFromHandle(a)) {
        c = arg2F->linkNode(b);
        return true;
      }
    }
  } // a is terminal
  if (arg2F->isTerminalNode(b)) {
    if (arg1F != resF) return false;
    if (resF->getRangeType() == range_type::INTEGER) {
      if (1==arg2F->getIntegerFromHandle(b)) {
        c = arg1F->linkNode(a);
        return true;
      }
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      if (1.0==arg2F->getRealFromHandle(b)) {
        c = arg1F->linkNode(a);
        return true;
      }
    }
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                       multiply_mxd class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::multiply_mxd : public generic_binary_mxd {
  public:
    multiply_mxd(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::multiply_mxd::multiply_mxd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mxd(MULTIPLY_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::multiply_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0 || b == 0) {
    c = 0;
    return true;
  }

  if (arg1F->isTerminalNode(a)) {
    if (arg2F->isTerminalNode(b)) {
      if (resF->getRangeType() == range_type::INTEGER) {
        int av, bv;
        arg1F->getValueFromHandle(a, av);
        arg2F->getValueFromHandle(b, bv);
        c = resF->handleForValue(av * bv);
      } else {
        MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
        float av, bv;
        arg1F->getValueFromHandle(a, av);
        arg2F->getValueFromHandle(b, bv);
        c = resF->handleForValue(av * bv);
      }
      return true;
    }
    if (arg2F != resF) return false;
    if (resF->getRangeType() == range_type::INTEGER) {
      if (1==arg1F->getIntegerFromHandle(a)) {
        c = arg2F->linkNode(b);
        return true;
      }
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      if (1.0==arg1F->getRealFromHandle(a)) {
        c = arg2F->linkNode(b);
        return true;
      }
    }
  } // a is terminal
  if (arg2F->isTerminalNode(b)) {
    if (arg1F != resF) return false;
    if (resF->getRangeType() == range_type::INTEGER) {
      if (1==arg2F->getIntegerFromHandle(b)) {
        c = arg1F->linkNode(a);
        return true;
      }
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      if (1.0==arg2F->getRealFromHandle(b)) {
        c = arg1F->linkNode(a);
        return true;
      }
    }
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                     multiply_evplus  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::multiply_evplus : public generic_binary_evplus {
  public:
    multiply_evplus(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long& cev, node_handle& c);
};

MEDDLY::multiply_evplus::multiply_evplus(forest* arg1, forest* arg2,
    forest* res) : generic_binary_evplus(MULTIPLY_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, arg1->isForRelations());
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

bool MEDDLY::multiply_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == 0 || b == 0) {
    c = 0; cev = 0;
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev * bev;
    return true;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                     multiply_evtimes class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::multiply_evtimes : public generic_binary_evtimes {
  public:
    multiply_evtimes(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(float aev, node_handle a, float bev, node_handle b,
      float& cev, node_handle& c);
};

MEDDLY::multiply_evtimes::multiply_evtimes(forest* arg1, forest* arg2,
    forest* res) : generic_binary_evtimes(MULTIPLY_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, arg1->isForRelations());
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVTIMES);
}

bool MEDDLY::multiply_evtimes::checkTerminals(float aev, node_handle a,
  float bev, node_handle b, float& cev, node_handle& c)
{
  if (a == 0 || b == 0) {
    c = 0; cev = Nan();
    return true;
  }
  if (a == -1 && b == -1) {
    c = -1; cev = aev * bev;
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::MULTIPLY(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;

    binary_operation* bop =  MULTIPLY_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations())
            return MULTIPLY_cache.add(new multiply_mxd(a, b, c));
        else
            return MULTIPLY_cache.add(new multiply_mdd(a, b, c));
    }

    if  (
            (a->getRangeType() != b->getRangeType()) ||
            (a->getRangeType() != c->getRangeType())
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS)
        return MULTIPLY_cache.add(new multiply_evplus(a, b, c));

    if (c->getEdgeLabeling() == edge_labeling::EVTIMES)
        return MULTIPLY_cache.add(new multiply_evtimes(a, b, c));

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::MULTIPLY_init()
{
    MULTIPLY_cache.reset("Multiply");
}

void MEDDLY::MULTIPLY_done()
{
    MEDDLY_DCASSERT(MULTIPLY_cache.isEmpty());
}


