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
#include "modulo.h"
#include "apply_base.h"

namespace MEDDLY {
    class modulo_mdd;
    class modulo_mxd;

    binary_list MODULO_cache;
};


// ******************************************************************
// *                                                                *
// *                        modulo_mdd class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::modulo_mdd : public generic_binary_mdd {
  public:
    modulo_mdd(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::modulo_mdd::modulo_mdd(forest* arg1, forest* arg2, forest* res)
    : generic_binary_mdd(MODULO_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    checkAllRanges(__FILE__, __LINE__, range_type::INTEGER);
}

bool MEDDLY::modulo_mdd::
checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b))
  {
    long av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    if (0 == bv) throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
    c = resF->handleForValue( av % bv );
    return true;
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                        modulo_mxd class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::modulo_mxd : public generic_binbylevel_mxd {
  public:
    modulo_mxd(forest* arg1, forest* arg2, forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::modulo_mxd::modulo_mxd(forest* arg1, forest* arg2, forest* res)
    : generic_binbylevel_mxd(MODULO_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    checkAllRanges(__FILE__, __LINE__, range_type::INTEGER);
}


bool MEDDLY::modulo_mxd
::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (arg1F->isTerminalNode(a) && arg2F->isTerminalNode(b))
  {
    long av, bv;
    arg1F->getValueFromHandle(a, av);
    arg2F->getValueFromHandle(b, bv);
    if (0 == bv) throw error(error::DIVIDE_BY_ZERO, __FILE__, __LINE__);
    c = resF->handleForValue( av / bv );
    return true;
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::MODULO(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;

    binary_operation* bop =  MODULO_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->isForRelations())
        return MODULO_cache.add(new modulo_mxd(a, b, c));
    else
        return MODULO_cache.add(new modulo_mdd(a, b, c));
}

void MEDDLY::MODULO_init()
{
    MODULO_cache.reset("Modulo");
}

void MEDDLY::MODULO_done()
{
    MEDDLY_DCASSERT(MODULO_cache.isEmpty());
}

