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
#include "maxmin.h"
#include "apply_base.h"

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"

namespace MEDDLY {
    class maximum_mdd;
    class maximum_mxd;
    binary_list MAXIMUM_cache;

    class minimum_mdd;
    class minimum_mxd;
    binary_list MINIMUM_cache;
};

// ******************************************************************
// *                                                                *
// *                       maximum_mdd  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::maximum_mdd : public generic_binary_mdd {
    public:
        maximum_mdd(forest* arg1, forest* arg2, forest* res);

    protected:
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};

MEDDLY::maximum_mdd::maximum_mdd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mdd(MAXIMUM_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::maximum_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MAX(av, bv));
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MAX(av, bv));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       maximum_mxd  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::maximum_mxd : public generic_binary_mxd {
    public:
        maximum_mxd(forest* arg1, forest* arg2, forest* res);

    protected:
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};

MEDDLY::maximum_mxd::maximum_mxd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mxd(MAXIMUM_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::maximum_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MAX(av, bv));
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MAX(av, bv));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       minimum_mdd  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::minimum_mdd : public generic_binary_mdd {
    public:
        minimum_mdd(forest* arg1, forest* arg2, forest* res);

    protected:
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};

MEDDLY::minimum_mdd::minimum_mdd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mdd(MINIMUM_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::minimum_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MIN(av, bv));
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MIN(av, bv));
    }
    return true;
  }
  return false;
}


// ******************************************************************
// *                                                                *
// *                       minimum_mxd  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::minimum_mxd : public generic_binary_mxd {
    public:
        minimum_mxd(forest* arg1, forest* arg2, forest* res);

    protected:
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};

MEDDLY::minimum_mxd::minimum_mxd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mxd(MINIMUM_cache, arg1, arg2, res)
{
    operationCommutes();
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::minimum_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    if (resF->getRangeType() == range_type::INTEGER) {
      int av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MIN(av, bv));
    } else {
      MEDDLY_DCASSERT(resF->getRangeType() == range_type::REAL);
      float av, bv;
      arg1F->getValueFromHandle(a, av);
      arg2F->getValueFromHandle(b, bv);
      c = resF->handleForValue(MIN(av, bv));
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

MEDDLY::binary_operation* MEDDLY::MAXIMUM(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  MAXIMUM_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations())
            return MAXIMUM_cache.add(new maximum_mxd(a, b, c));
        else
            return MAXIMUM_cache.add(new maximum_mdd(a, b, c));
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::MAXIMUM_init()
{
    MAXIMUM_cache.reset("Maximum");
}

void MEDDLY::MAXIMUM_done()
{
    MEDDLY_DCASSERT(MAXIMUM_cache.isEmpty());
}

MEDDLY::binary_operation* MEDDLY::MINIMUM(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  MINIMUM_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations())
            return MINIMUM_cache.add(new minimum_mxd(a, b, c));
        else
            return MINIMUM_cache.add(new minimum_mdd(a, b, c));
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::MINIMUM_init()
{
    MINIMUM_cache.reset("Minimum");
}

void MEDDLY::MINIMUM_done()
{
    MEDDLY_DCASSERT(MINIMUM_cache.isEmpty());
}

