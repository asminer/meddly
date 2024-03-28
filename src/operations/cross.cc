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
#include "cross.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"

// #define TRACE_ALL_OPS
// #define DEBUG_CROSS

namespace MEDDLY {
    class cross_bool;
};

// ******************************************************************
// *                                                                *
// *                        cross_bool class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::cross_bool : public binary_operation {
    protected:
        cross_bool(forest* a1, forest* a2, forest* res);

        virtual void computeDDEdge(const dd_edge& a, const dd_edge& b,
                dd_edge &c, bool userFlag);

        node_handle compute_pr(unsigned in, int ht, node_handle a,
                node_handle b);

        node_handle compute_un(int ht, node_handle a, node_handle b);

    public:
        static binary_list cache;

        inline static binary_operation* build(forest* a, forest* b, forest* c)
        {
            binary_operation* bop =  cache.find(a, b, c);
            if (bop) {
                return bop;
            }
            return cache.add(new cross_bool(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::cross_bool::cache;

MEDDLY::cross_bool::cross_bool(forest* a1, forest* a2, forest* res)
: binary_operation(cache, 1, a1, a2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
    checkRelations(__FILE__, __LINE__, SET, SET, RELATION);

    ct_entry_type* et = new ct_entry_type(cache.getName(), "INN:N");
    et->setForestForSlot(1, a1);
    et->setForestForSlot(2, a2);
    et->setForestForSlot(4, res);
    registerEntryType(0, et);
    buildCTs();
}

void
MEDDLY::cross_bool::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c, bool userFlag)
{
  int L = arg1F->getMaxLevelIndex();
  node_handle cnode = compute_un(L, a.getNode(), b.getNode());
  c.set(cnode);
}

MEDDLY::node_handle MEDDLY::cross_bool::compute_un(int k, node_handle a, node_handle b)
{
#ifdef DEBUG_CROSS
  printf("calling compute_un(%d, %d, %d)\n", k, a, b);
#endif
  MEDDLY_DCASSERT(k>=0);
  if (0==a || 0==b) return 0;
  if (0==k) {
    return a;
  }

  // check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeI(k);
  CTsrch->writeN(a);
  CTsrch->writeN(b);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    CT0->recycle(CTsrch);
    return resF->linkNode(CTresult[0].readN());
  }

  // Initialize unpacked node
  unpacked_node *A = unpacked_node::New(arg1F);
  if (arg1F->getNodeLevel(a) < k) {
    A->initRedundant(arg1F, k, a, FULL_ONLY);
  } else {
    arg1F->unpackNode(A, a, FULL_ONLY);
  }

  const unsigned resultSize = unsigned(resF->getLevelSize(k));
  unpacked_node *C = unpacked_node::newFull(resF, k, resultSize);

  // recurse
  for (unsigned i=0; i<resultSize; i++) {
    C->setFull(i, compute_pr(i, -k, A->down(i), b));
  }

  // reduce, save in compute table
  unpacked_node::Recycle(A);
  node_handle c = resF->createReducedNode(-1, C);

  CTresult[0].reset();
  CTresult[0].writeN(c);
  CT0->addEntry(CTsrch, CTresult[0]);

#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d, %d) = %d\n", getName(), k, a, b, c);
#endif

  return c;
}

MEDDLY::node_handle MEDDLY::cross_bool::compute_pr(unsigned in, int k, node_handle a, node_handle b)
{
#ifdef DEBUG_CROSS
  printf("calling compute_pr(%d, %d, %d, %d)\n", in, k, a, b);
#endif
  MEDDLY_DCASSERT(k<0);
  if (0==a || 0==b) return 0;

  // DON'T check compute table

  // Initialize unpacked node
  unpacked_node *B = unpacked_node::New(arg2F);
  if (arg2F->getNodeLevel(b) < -k) {
    B->initRedundant(arg2F, -k, b, FULL_ONLY);
  } else {
    arg2F->unpackNode(B, b, FULL_ONLY);
  }

  const unsigned resultSize = unsigned(resF->getLevelSize(k));
  unpacked_node *C = unpacked_node::newFull(resF, k, resultSize);

  // recurse
  for (unsigned i=0; i<resultSize; i++) {
    C->setFull(i, compute_un(-(k+1), a, B->down(i)));
  }

  // reduce
  unpacked_node::Recycle(B);
  node_handle c = resF->createReducedNode(int(in), C);

  // DON'T save in compute table

#ifdef TRACE_ALL_OPS
  printf("computed %s((%d), %d, %d, %d) = %d\n", getName(), in, k, a, b, c);
#endif

  return c;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::CROSS(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }

    return cross_bool::build(a, b, c);
}

void MEDDLY::CROSS_init()
{
    cross_bool::cache.reset("Cross");
}

void MEDDLY::CROSS_done()
{
    MEDDLY_DCASSERT(cross_bool::cache.isEmpty());
}

