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
#include "union.h"
#include "apply_base.h"

namespace MEDDLY {
    class union_mdd;
    class union_mxd;

    class union_min_evplus;
    class union_min_evplus_mxd;
};

// ******************************************************************
// *                                                                *
// *                        union_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mdd : public generic_binary_mdd {
    protected:
        union_mdd(forest* arg1, forest* arg2, forest* res);

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
            return cache.addOperation(new union_mdd(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::union_mdd::cache;

MEDDLY::union_mdd::union_mdd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mdd(cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::union_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a < 0 || b < 0) {
    c = resF->handleForValue(true);
    return true;
  }
  if (a == 0) {
    if (b==0) {
      c = 0;
      return true;
    }
    if (arg2F == resF) {
      c = resF->linkNode(b);
      return true;
    }
    return false;
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    }
    return false;
  }
  if (a == b) {
    if (arg1F == arg2F && arg1F == resF) {
      c = resF->linkNode(b);
      return true;
    }
    return false;
  }
  return false;
}



// ******************************************************************
// *                                                                *
// *                        union_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mxd : public generic_binary_mxd {
    protected:
        union_mxd(forest* arg1, forest* arg2, forest* res);

        virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
        virtual MEDDLY::node_handle compute_ext(node_handle a, node_handle b);

    public:
        static binary_list cache;

        inline static binary_operation* build(forest* a, forest* b, forest* c)
        {
            binary_operation* bop =  cache.findOperation(a, b, c);
            if (bop) {
                return bop;
            }
            return cache.addOperation(new union_mxd(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::union_mxd::cache;

MEDDLY::union_mxd::union_mxd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::union_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a < 0 && b < 0) {
    c = resF->handleForValue(true);
    return true;
  }
  if (0 == a) {
    if (0 == b) {
      c = 0;
      return true;
    }
    if (arg2F == resF) {
      c = resF->linkNode(b);
      return true;
    } else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg1F == resF) {
      c = resF->linkNode(b);
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

MEDDLY::node_handle
MEDDLY::union_mxd::compute_ext(node_handle a, node_handle b)
{
  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));

  MEDDLY_DCASSERT(resF->isExtensibleLevel(resultLevel));

  const int dwnLevel = resF->downLevel(resultLevel);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, SPARSE_ONLY)
    : arg1F->newUnpacked(a, SPARSE_ONLY)
    ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, SPARSE_ONLY)
    : arg2F->newUnpacked(b, SPARSE_ONLY)
    ;

  // Initialize result writer
  unsigned resultSize = A->getSize() + B->getSize() + 1 + 1;
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = A->getSize()-1;
  for ( ; last_nz >= 0 && A->down(last_nz) == 0; last_nz--);
  const unsigned int A_nnzs = last_nz + 1;
  const int A_last_index = last_nz >= 0? A->index(last_nz): -1;

  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getSize()-1;
  for ( ; last_nz >= 0 && B->down(last_nz) == 0; last_nz--);
  const unsigned int B_nnzs = last_nz + 1;
  const int B_last_index = last_nz >= 0? B->index(last_nz): -1;

  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  unsigned nnz = 0;
  unsigned A_curr_index = 0;
  unsigned B_curr_index = 0;
  for ( ; A_curr_index < A_nnzs && B_curr_index < B_nnzs; ) {
    // get a_i, a_d, b_i, b_d
    unsigned a_i, b_i;
    node_handle a_d, b_d;
    a_i = A->index(A_curr_index);
    b_i = B->index(B_curr_index);
    if (a_i <= b_i) {
      a_d = A->down(A_curr_index);
      A_curr_index++;
    } else {
      a_d = 0;
    }
    if (a_i >= b_i) {
      b_d = B->down(B_curr_index);
      B_curr_index++;
    } else {
      b_d = 0;
    }

    MEDDLY_DCASSERT(a_d != 0 || b_d != 0);

    // compute union(a_d, b_d)
    unsigned index = (a_d? a_i: b_i);
    node_handle down = compute_r(int(index), dwnLevel, a_d, b_d);

    // if union is non-zero, add it to the new node
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  }
  for ( ; A_curr_index < A_nnzs; A_curr_index++) {
    // do union(a_i, b_ext_i)
    unsigned index = A->index(A_curr_index);
    node_handle down = compute_r(int(index), dwnLevel, A->down(A_curr_index), B_ext_d);
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  }
  for ( ; B_curr_index < B_nnzs; B_curr_index++) {
    // do union(a_ext_i, b_i)
    unsigned index = B->index(B_curr_index);
    node_handle down = compute_r(int(index), dwnLevel, A_ext_d, B->down(B_curr_index));
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  }
  if (A->isExtensible() || B->isExtensible()) {
    const unsigned index = max_a_b_last_index+1;
    node_handle down = compute_r(index, dwnLevel, A_ext_d, B_ext_d);
    if (down) {
      MEDDLY_DCASSERT(index >= 0);
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = unsigned(index);
      // C->d_ref(nnz) = down;
      C->markAsExtensible();
      nnz++;
    } else {
      C->markAsNotExtensible();
    }
  } else {
    C->markAsNotExtensible();
  }
  C->shrink(nnz);

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce result
  node_handle result = resF->createReducedNode(-1, C);
  return result;
}

// ******************************************************************
// *                                                                *
// *                    union_min_evplus  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::union_min_evplus : public generic_binary_evplus {
    protected:
        union_min_evplus(forest* arg1, forest* arg2, forest* res);

        virtual ct_entry_key* findResult(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle &c);
        virtual void saveResult(ct_entry_key* key, long aev, node_handle a,
            long bev, node_handle b, long cev, node_handle c);

        virtual bool checkTerminals(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle& c);

    public:
        static binary_list cache;

        inline static binary_operation* build(forest* a, forest* b, forest* c)
        {
            binary_operation* bop =  cache.findOperation(a, b, c);
            if (bop) {
                return bop;
            }
            return cache.addOperation(new union_min_evplus(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::union_min_evplus::cache;

MEDDLY::union_min_evplus::union_min_evplus(forest* a, forest* b, forest* c)
  : generic_binary_evplus(cache, a, b, c)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::union_min_evplus::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  if (can_commute && a > b) {
    CTsrch->writeL(0);
    CTsrch->writeN(b);
    CTsrch->writeL(aev - bev);
    CTsrch->writeN(a);
  } else {
    CTsrch->writeL(0);
    CTsrch->writeN(a);
    CTsrch->writeL(bev - aev);
    CTsrch->writeN(b);
  }
  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  cev = CTresult[0].readL();
  MEDDLY_DCASSERT(cev == 0);
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev = MIN(aev, bev);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::union_min_evplus::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  MEDDLY_DCASSERT(c == 0 || cev == MIN(aev, bev));
  CTresult[0].reset();
  CTresult[0].writeL(0);   //   Why always 0?
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::union_min_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == 0) {
    if (b == 0) {
      cev = 0;
      c = 0;
      return true;
    }
    else if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg1F->isTerminalNode(a) && aev <= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && bev <= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MIN(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                  union_min_evplus_mxd  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::union_min_evplus_mxd : public generic_binary_evplus_mxd {
    protected:
        union_min_evplus_mxd(forest* arg1, forest* arg2, forest* res);

        virtual ct_entry_key* findResult(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle &c);
        virtual void saveResult(ct_entry_key* key, long aev, node_handle a,
            long bev, node_handle b, long cev, node_handle c);

        virtual bool checkTerminals(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle& c);

    public:
        static binary_list cache;

        inline static binary_operation* build(forest* a, forest* b, forest* c)
        {
            binary_operation* bop =  cache.findOperation(a, b, c);
            if (bop) {
                return bop;
            }
            return cache.addOperation(new union_min_evplus_mxd(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::union_min_evplus_mxd::cache;

MEDDLY::union_min_evplus_mxd::union_min_evplus_mxd(
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_evplus_mxd(cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::union_min_evplus_mxd::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  if (can_commute && a > b) {
    CTsrch->writeL(0);
    CTsrch->writeN(b);
    CTsrch->writeL(aev - bev);
    CTsrch->writeN(a);
  } else {
    CTsrch->writeL(0);
    CTsrch->writeN(a);
    CTsrch->writeL(bev - aev);
    CTsrch->writeN(b);
  }
  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  cev = CTresult[0].readL();
  MEDDLY_DCASSERT(cev == 0);
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev = MIN(aev, bev);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::union_min_evplus_mxd::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  MEDDLY_DCASSERT(c == 0 || cev == MIN(aev, bev));
  CTresult[0].reset();
  CTresult[0].writeL(0);   // why always 0?
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::union_min_evplus_mxd::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == 0) {
    if (b == 0) {
      cev = 0;
      b = 0;
      return true;
    }
    else if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg1F->isTerminalNode(a) && aev <= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && bev <= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MIN(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
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
MEDDLY::UNION(MEDDLY::forest* a, MEDDLY::forest* b, MEDDLY::forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return union_mxd::build(a, b, c);
        } else {
            return union_mdd::build(a, b, c);
        }
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS) {
        if (c->isForRelations()) {
            return union_min_evplus_mxd::build(a, b, c);
        } else {
            return union_min_evplus::build(a, b, c);
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::UNION_init()
{
    union_mdd::cache.reset("Union");
    union_mxd::cache.reset("Union");
    union_min_evplus::cache.reset("Union");
    union_min_evplus_mxd::cache.reset("Union");
}

void MEDDLY::UNION_done()
{
    MEDDLY_DCASSERT(union_mdd::cache.isEmpty());
    MEDDLY_DCASSERT(union_mxd::cache.isEmpty());
    MEDDLY_DCASSERT(union_min_evplus::cache.isEmpty());
    MEDDLY_DCASSERT(union_min_evplus_mxd::cache.isEmpty());
}

