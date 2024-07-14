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
#include "compare.h"

// TBD

// ******************************************************************
// ******************************************************************
// ******************************************************************
//
// OLD IMPLEMENTATION HERE
//
// ******************************************************************
// ******************************************************************
// ******************************************************************

#include "apply_base.h"

// ******************************************************************
// ******************************************************************
// **                            EQUAL                             **
// ******************************************************************
// ******************************************************************

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
// ******************************************************************
// **                          NOT_EQUAL                           **
// ******************************************************************
// ******************************************************************

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
// ******************************************************************
// **                        GREATER  THAN                         **
// ******************************************************************
// ******************************************************************

namespace MEDDLY {
    binary_list GT_cache;
};

// ******************************************************************
// *                                                                *
// *                       morethan_mdd class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class morethan_mdd : public generic_binary_mdd {
  public:
    morethan_mdd(forest* arg1, forest* arg2, forest* res)
      : generic_binary_mdd(GT_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, SET);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool morethan_mdd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av > bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                       morethan_mxd class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class morethan_mxd : public generic_binbylevel_mxd {
  public:
    morethan_mxd(forest* arg1, forest* arg2, forest* res)
      : generic_binbylevel_mxd(GT_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, RELATION);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool morethan_mxd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av > bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// ******************************************************************
// **                      GREATER  OR EQUAL                       **
// ******************************************************************
// ******************************************************************

namespace MEDDLY {
    binary_list GE_cache;
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
    moreequal_mdd(forest* arg1, forest* arg2, forest* res)
      : generic_binary_mdd(GE_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, SET);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

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
    moreequal_mxd(forest* arg1, forest* arg2, forest* res)
      : generic_binbylevel_mxd(GE_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, RELATION);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

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
// ******************************************************************
// **                          LESS THAN                           **
// ******************************************************************
// ******************************************************************

namespace MEDDLY {
    binary_list LT_cache;
};

// ******************************************************************
// *                                                                *
// *                       lessthan_mdd class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class lessthan_mdd : public generic_binary_mdd {
  public:
    lessthan_mdd(forest* arg1, forest* arg2, forest* res)
      : generic_binary_mdd(LT_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, SET);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool lessthan_mdd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av <  bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                       lessthan_mxd class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class lessthan_mxd : public generic_binbylevel_mxd {
  public:
    lessthan_mxd(forest* arg1, forest* arg2, forest* res)
      : generic_binbylevel_mxd(LT_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, RELATION);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool lessthan_mxd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av <  bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// ******************************************************************
// **                        LESS OR EQUAL                         **
// ******************************************************************
// ******************************************************************

namespace MEDDLY {
    binary_list LE_cache;
};

// ******************************************************************
// *                                                                *
// *                      lessequal_mdd  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class lessequal_mdd : public generic_binary_mdd {
  public:
    lessequal_mdd(forest* arg1, forest* arg2, forest* res)
      : generic_binary_mdd(LE_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, SET);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool lessequal_mdd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == b) {
    if (arg1F == arg2F) {
      c = resF->handleForValue(true);
      return true;
    }
  }
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av <= bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY

// ******************************************************************
// *                                                                *
// *                      lessequal_mxd  class                      *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

template <typename T>
class lessequal_mxd : public generic_binbylevel_mxd {
  public:
    lessequal_mxd(forest* arg1, forest* arg2, forest* res)
      : generic_binbylevel_mxd(LE_cache, arg1, arg2, res)
    {
        checkDomains(__FILE__, __LINE__);
        checkAllRelations(__FILE__, __LINE__, RELATION);
        checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    }

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

template <typename T>
bool lessequal_mxd<T>
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b)) {
    T av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    c = resF->handleForValue( av <= bv );
    return true;
  }
  return false;
}

};  // namespace MEDDLY


// ******************************************************************
// *                                                                *
// *                                                                *
// *                           Front ends                           *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                             EQUAL                              *
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

// ******************************************************************
// *                            NOT_EQUAL                           *
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

// ******************************************************************
// *                         GREATER  THAN                          *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::GREATER_THAN(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  GT_cache.find(a, b, c);
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
                return GT_cache.add(new morethan_mxd<float>(a, b, c));
            else
                return GT_cache.add(new morethan_mdd<float>(a, b, c));
        } else {
            if (c->isForRelations())
                return GT_cache.add(new morethan_mxd<long>(a, b, c));
            else
                return GT_cache.add(new morethan_mdd<long>(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::GREATER_THAN_init()
{
    GT_cache.reset("MoreThan");
}

void MEDDLY::GREATER_THAN_done()
{
    MEDDLY_DCASSERT(GT_cache.isEmpty());
}

// ******************************************************************
// *                       GREATER  OR EQUAL                        *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::GREATER_THAN_EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  GE_cache.find(a, b, c);
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
                return GE_cache.add(new moreequal_mxd<float>(a, b, c));
            else
                return GE_cache.add(new moreequal_mdd<float>(a, b, c));
        } else {
            if (c->isForRelations())
                return GE_cache.add(new moreequal_mxd<long>(a, b, c));
            else
                return GE_cache.add(new moreequal_mdd<long>(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::GREATER_THAN_EQUAL_init()
{
    GE_cache.reset("MoreEqual");
}

void MEDDLY::GREATER_THAN_EQUAL_done()
{
    MEDDLY_DCASSERT(GE_cache.isEmpty());
}

// ******************************************************************
// *                           LESS THAN                            *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::LESS_THAN(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  LT_cache.find(a, b, c);
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
                return LT_cache.add(new lessthan_mxd<float>(a, b, c));
            else
                return LT_cache.add(new lessthan_mdd<float>(a, b, c));
        } else {
            if (c->isForRelations())
                return LT_cache.add(new lessthan_mxd<long>(a, b, c));
            else
                return LT_cache.add(new lessthan_mdd<long>(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::LESS_THAN_init()
{
    LT_cache.reset("LessThan");
}

void MEDDLY::LESS_THAN_done()
{
    MEDDLY_DCASSERT(LT_cache.isEmpty());
}

// ******************************************************************
// *                         LESS OR EQUAL                          *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::LESS_THAN_EQUAL(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  LE_cache.find(a, b, c);
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
                return LE_cache.add(new lessequal_mxd<float>(a, b, c));
            else
                return LE_cache.add(new lessequal_mdd<float>(a, b, c));
        } else {
            if (c->isForRelations())
                return LE_cache.add(new lessequal_mxd<long>(a, b, c));
            else
                return LE_cache.add(new lessequal_mdd<long>(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::LESS_THAN_EQUAL_init()
{
    LE_cache.reset("LessEqual");
}

void MEDDLY::LESS_THAN_EQUAL_done()
{
    MEDDLY_DCASSERT(LE_cache.isEmpty());
}

