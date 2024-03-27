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
#include "intersection.h"
#include "apply_base.h"

namespace MEDDLY {
    class inter_mdd;
    class inter_mxd;
    class inter_max_evplus;
};


// ******************************************************************
// *                                                                *
// *                        inter_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_mdd : public generic_binary_mdd {
    protected:
        inter_mdd(forest* arg1, forest* arg2, forest* res);

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
            return cache.addOperation(new inter_mdd(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::inter_mdd::cache;

MEDDLY::inter_mdd::inter_mdd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mdd(cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::inter_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0 || b == 0) {
    c = 0;
    return true;
  }
  if (a==-1 && b==-1) {
    c = -1;
    return true;
  }
  if (a == -1) {
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
  if (b == -1) {
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
// *                        inter_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_mxd : public generic_binary_mxd {
    protected:
        inter_mxd(forest* arg1, forest* arg2, forest* res);

        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
        virtual MEDDLY::node_handle compute_ext(node_handle a, node_handle b);
    public:
        static binary_list cache;

        inline static binary_operation* build(forest* a, forest* b, forest* c)
        {
            binary_operation* bop =  cache.findOperation(a, b, c);
            if (bop) {
                return bop;
            }
            return cache.addOperation(new inter_mxd(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::inter_mxd::cache;

MEDDLY::inter_mxd::inter_mxd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::inter_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0 || b == 0) {
    c = 0;
    return true;
  }
  if (a==-1 && b==-1) {
    c = -1;
    return true;
  }
  if (a == b) {
    if (arg1F == arg2F && arg1F == resF) {
      c = resF->linkNode(b);
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                     inter_max_evplus  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_max_evplus : public generic_binary_evplus {
    protected:
        inter_max_evplus(forest* arg1, forest* arg2, forest* res);

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
            return cache.addOperation(new inter_max_evplus(a, b, c));
        }
};

MEDDLY::binary_list MEDDLY::inter_max_evplus::cache;

MEDDLY::inter_max_evplus::inter_max_evplus(forest* arg1, forest* arg2,
        forest* res) : generic_binary_evplus(cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::inter_max_evplus::findResult(long aev, node_handle a,
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
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev += (a > b ? bev : aev);
  }
  else {
    MEDDLY_DCASSERT(cev == 0);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::inter_max_evplus::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  if (c == 0) {
    CTresult[0].writeL(0);
  }
  else {
    CTresult[0].writeL(cev - (a > b ? bev : aev));
  }
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::inter_max_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
    long& cev, node_handle& c)
{
  if (a == 0 || b == 0) {
    cev = 0;
    c = 0;
    return true;
  }
  if (arg1F->isTerminalNode(a) && bev >= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && aev >= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MAX(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}

MEDDLY::node_handle
MEDDLY::inter_mxd::compute_ext(node_handle a, node_handle b)
{
  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));
  const int dwnLevel = resF->downLevel(resultLevel);

  MEDDLY_DCASSERT(resF->isExtensibleLevel(resultLevel));

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, SPARSE_ONLY)
    : arg1F->newUnpacked(a, SPARSE_ONLY)
    ;
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = int(A->getSize())-1;
  for ( ; last_nz >= 0 && A->down(unsigned(last_nz)) == 0; last_nz--);
  const unsigned int A_nnzs = last_nz + 1;
  const int A_last_index = last_nz >= 0? int(A->index(unsigned(last_nz))): -1;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, SPARSE_ONLY)
    : arg2F->newUnpacked(b, SPARSE_ONLY)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = int(B->getSize())-1;
  for ( ; last_nz >= 0 && B->down(unsigned(last_nz)) == 0; last_nz--);
  const unsigned int B_nnzs = last_nz + 1;
  const int B_last_index = last_nz >= 0? int(B->index(unsigned(last_nz))): -1;

  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  unsigned resultSize = A->getSize() + B->getSize() + 1 + 1;
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  unsigned nnz = 0;
  unsigned A_curr_index = 0;
  unsigned B_curr_index = 0;
  for ( ; A_curr_index < A_nnzs && B_curr_index < B_nnzs; ) {
    // get a_i, a_d, b_i, b_d
    node_handle a_d, b_d;
    const unsigned a_i = A->index(A_curr_index);
    const unsigned b_i = B->index(B_curr_index);
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

    if (a_d == 0 || b_d == 0) continue;

    // compute inter(a_d, b_d)
    const unsigned index = (a_d? a_i: b_i);
    const node_handle down = compute_r(int(index), dwnLevel, a_d, b_d);

    // if inter is non-zero, add it to the new node
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  } // for loop
  if (B_ext_d != 0) {
    for ( ; A_curr_index < A_nnzs; A_curr_index++) {
      // do inter(a_i, b_ext_i)
      const unsigned index = A->index(A_curr_index);
      const node_handle down = compute_r(int(index), dwnLevel, A->down(A_curr_index), B_ext_d);
      if (down) {
        C->setSparse(nnz, index, down);
        // C->i_ref(nnz) = index;
        // C->d_ref(nnz) = down;
        nnz++;
      }
    }
  }
  if (A_ext_d != 0) {
    for ( ; B_curr_index < B_nnzs; B_curr_index++) {
      // do inter(a_ext_i, b_i)
      const unsigned index = B->index(B_curr_index);
      node_handle down = compute_r(int(index), dwnLevel, A_ext_d, B->down(B_curr_index));
      if (down) {
        C->setSparse(nnz, index, down);
        // C->i_ref(nnz) = index;
        // C->d_ref(nnz) = down;
        nnz++;
      }
    }
  }
  if (A_ext_d != 0 && B_ext_d != 0) {
    const unsigned index = max_a_b_last_index+1;
    MEDDLY_DCASSERT(index >= 0);
    const node_handle down = compute_r(index, dwnLevel, A_ext_d, B_ext_d);
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = unsigned(index);
      // C->d_ref(nnz) = down;
      C->markAsExtensible();
      nnz++;
    } else {
      C->markAsNotExtensible();
    }
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
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::INTERSECTION(MEDDLY::forest* a, MEDDLY::forest* b, MEDDLY::forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return inter_mxd::build(a, b, c);
        } else {
            return inter_mdd::build(a, b, c);
        }
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS) {
        if (! c->isForRelations()) {
            return inter_max_evplus::build(a, b, c);
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::INTERSECTION_init()
{
    inter_mdd::cache.reset("Intersection");
    inter_mxd::cache.reset("Intersection");
    inter_max_evplus::cache.reset("Intersection");
}

void MEDDLY::INTERSECTION_done()
{
    MEDDLY_DCASSERT(inter_mdd::cache.isEmpty());
    MEDDLY_DCASSERT(inter_mxd::cache.isEmpty());
    MEDDLY_DCASSERT(inter_max_evplus::cache.isEmpty());
}

