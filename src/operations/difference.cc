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
#include "difference.h"
#include "apply_base.h"

namespace MEDDLY {
    class diffr_mdd;
    class diffr_mxd;
};

// ******************************************************************
// *                                                                *
// *                        diffr_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mdd : public generic_binary_mdd {
    protected:
        diffr_mdd(forest* arg1, forest* arg2, forest* res);
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
    public:
        static binary_list cache;

        inline static binary_operation* build(forest* a, forest* b, forest* c)
        {
            binary_operation* bop =  cache.findOperation(a, b, c);
            if (bop) {
                return bop;
            }
            return cache.addOperation(new diffr_mdd(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::diffr_mdd::cache;

MEDDLY::diffr_mdd::diffr_mdd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mdd(cache, arg1, arg2, res)
{
    //  difference does NOT commute

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::diffr_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0 || b == -1) {
    c = 0;
    return true;
  }
  if (a == -1 && b == 0) {
    c = -1;
    return true;
  }
  if (a == b) {
    if (arg1F == arg2F) {
      c = 0;
      return true;
    } else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } else {
      return false;
    }
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                        diffr_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mxd : public generic_binary_mxd {
    protected:
        diffr_mxd(forest* arg1, forest* arg2, forest* res);
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);

    public:
        static binary_list cache;

        inline static binary_operation* build(forest* a, forest* b, forest* c)
        {
            binary_operation* bop =  cache.findOperation(a, b, c);
            if (bop) {
                return bop;
            }
            return cache.addOperation(new diffr_mxd(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::diffr_mxd::cache;

MEDDLY::diffr_mxd::diffr_mxd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(cache, arg1, arg2, res)
{
    //  difference does NOT commute

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::diffr_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0) {
    c = 0;
    return true;
  }
  if (a == -1 && b == 0) {
    c = -1;
    return true;
  }
  if (a == b) {
    if (arg1F == arg2F) {
      c = 0;
      return true;
    } else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::DIFFERENCE(MEDDLY::forest* a, MEDDLY::forest* b, MEDDLY::forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return diffr_mxd::build(a, b, c);
        } else {
            return diffr_mdd::build(a, b, c);
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::DIFFERENCE_init()
{
    diffr_mdd::cache.reset("Difference");
    diffr_mxd::cache.reset("Difference");
}

void MEDDLY::DIFFERENCE_done()
{
    MEDDLY_DCASSERT(diffr_mdd::cache.isEmpty());
    MEDDLY_DCASSERT(diffr_mxd::cache.isEmpty());
}
