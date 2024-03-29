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
#include "plus.h"
#include "apply_base.h"

namespace MEDDLY {
    class plus_mdd;
    class plus_mxd;
    class plus_evplus;
    class plus_evtimes;

    binary_list PLUS_cache;
};

// ******************************************************************
// *                                                                *
// *                         plus_mdd class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_mdd : public generic_binary_mdd {
    public:
        plus_mdd(forest* arg1, forest* arg2, forest* res);

    protected:
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};

MEDDLY::plus_mdd::plus_mdd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mdd(PLUS_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::plus_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av + bv);
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av + bv);
    }
    return true;
  }
  if (0==a) {
    if (arg2F == resF) {
      c = arg2F->linkNode(b);
      return true;
    }
    return false;
  }
  if (0==b) {
    if (arg1F == resF) {
      c = arg1F->linkNode(a);
      return true;
    }
    return false;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                         plus_mxd class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_mxd : public generic_binary_mxd {
    public:
        plus_mxd(forest* arg1, forest* arg2, forest* res);

    protected:
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};

MEDDLY::plus_mxd::plus_mxd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(PLUS_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::plus_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av + bv);
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(av + bv);
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       plus_evplus  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_evplus : public generic_binary_evplus {
  public:
    plus_evplus(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual ct_entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c);
    virtual void saveResult(ct_entry_key* key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c);

    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long& cev, node_handle& c);
};

MEDDLY::plus_evplus::plus_evplus(forest* arg1, forest* arg2, forest* res)
  : generic_binary_evplus(PLUS_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::plus_evplus::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  if (can_commute && a > b) {
    CTsrch->writeL(0);
    CTsrch->writeN(b);
    CTsrch->writeL(0);
    CTsrch->writeN(a);
  } else {
    CTsrch->writeL(0);
    CTsrch->writeN(a);
    CTsrch->writeL(0);
    CTsrch->writeN(b);
  }
  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  cev = CTresult[0].readL();
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev += aev + bev;
  }
  else {
    MEDDLY_DCASSERT(cev == 0);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::plus_evplus::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  CTresult[0].reset();
  CTresult[0].writeL(c == 0 ? 0L : cev - aev - bev);
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::plus_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1;
    cev = aev + bev;
    MEDDLY_DCASSERT(cev >= 0);
    return true;
  }
  if (a == 0 || b == 0) {
    c = 0;
    cev = 0;
    return true;
  }
  if (a == b && arg1F == arg2F && arg2F == resF) {
    c = resF->linkNode(a);
    cev = aev + bev;
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       plus_evtimes class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::plus_evtimes : public generic_binary_evtimes {
  public:
    plus_evtimes(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(float aev, node_handle a, float bev, node_handle b,
      float& cev, node_handle& c);
};

MEDDLY::plus_evtimes::plus_evtimes(forest* arg1, forest* arg2, forest* res)
  : generic_binary_evtimes(PLUS_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVTIMES);
}

bool MEDDLY::plus_evtimes::checkTerminals(float aev, node_handle a,
  float bev, node_handle b, float& cev, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (0 == a) {
      c = b;
      cev = bev;
    } else {
      c = a;
      if (0 == b) cev = aev; else cev = aev + bev;
    }
    return true;
  }
  if (0 == a) {
    if (arg2F == resF) {
      c = resF->linkNode(b);
      cev = bev;
      return true;
    }
    return false;
  }
  if (0 == b) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      cev = aev;
      return true;
    }
    return false;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::PLUS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;

    binary_operation* bop =  PLUS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations())
            return PLUS_cache.add(new plus_mxd(a, b, c));
        else
            return PLUS_cache.add(new plus_mdd(a, b, c));
    }

    if  (
            (a->getRangeType() != b->getRangeType()) ||
            (a->getRangeType() != c->getRangeType())
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS)
        return PLUS_cache.add(new plus_evplus(a, b, c));

    if (c->getEdgeLabeling() == edge_labeling::EVTIMES)
        return PLUS_cache.add(new plus_evtimes(a, b, c));

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::PLUS_init()
{
    PLUS_cache.reset("Plus");
}

void MEDDLY::PLUS_done()
{
    MEDDLY_DCASSERT(PLUS_cache.isEmpty());
}

