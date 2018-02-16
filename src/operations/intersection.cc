
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "intersection.h"
#include "apply_base.h"

namespace MEDDLY {
  class inter_mdd;
  class inter_mxd;
  class inter_max_evplus;

  class inter_opname;
};


// ******************************************************************
// *                                                                *
// *                        inter_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_mdd : public generic_binary_mdd {
  public:
    inter_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::inter_mdd::inter_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
  operationCommutes();
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
  public:
    inter_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
#ifdef USE_XDDS
    virtual MEDDLY::node_handle compute(node_handle a, node_handle b);
#endif
};

MEDDLY::inter_mxd::inter_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mxd(opcode, arg1, arg2, res)
{
  operationCommutes();
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
  public:
    inter_max_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual compute_table::search_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c);
    virtual void saveResult(compute_table::search_key* key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c);

    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
        long& cev, node_handle& c);
};

MEDDLY::inter_max_evplus::inter_max_evplus(const binary_opname* opcode,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_evplus(opcode, arg1, arg2, res)
{
  operationCommutes();
}

MEDDLY::compute_table::search_key* MEDDLY::inter_max_evplus::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  compute_table::search_key* CTsrch = useCTkey();
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->reset();
  if (can_commute && a > b) {
    CTsrch->write(0L);
    CTsrch->writeNH(b);
    CTsrch->write(aev - bev);
    CTsrch->writeNH(a);
  } else {
    CTsrch->write(0L);
    CTsrch->writeNH(a);
    CTsrch->write(bev - aev);
    CTsrch->writeNH(b);
  }
  compute_table::search_result &cacheFind = CT->find(CTsrch);
  if (!cacheFind) return CTsrch;
  cacheFind.read(cev);
  c = resF->linkNode(cacheFind.readNH());
  if (c != 0) {
    cev += (a > b ? bev : aev);
  }
  else {
    MEDDLY_DCASSERT(cev == 0);
  }
  doneCTkey(CTsrch);
  return 0;
}

void MEDDLY::inter_max_evplus::saveResult(compute_table::search_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  arg1F->cacheNode(a);
  arg2F->cacheNode(b);
  compute_table::entry_builder &entry = CT->startNewEntry(key);
  if (c == 0) {
    entry.writeResult(0L);
  }
  else {
    entry.writeResult(cev - (a > b ? bev : aev));
  }
  entry.writeResultNH(resF->cacheNode(c));
  CT->addEntry();
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

#ifdef USE_XDDS
MEDDLY::node_handle 
MEDDLY::inter_mxd::compute(node_handle a, node_handle b) 
{
  //  Compute for the unprimed levels.
  //
  node_handle result = 0;
  if (checkTerminals(a, b, result))
    return result;

  compute_table::search_key* Key = findResult(a, b, result);
  if (0==Key) return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));
  const int dwnLevel = resF->downLevel(resultLevel);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel) 
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, false)
    : unpacked_node::newFromNode(arg1F, a, false)
    ;
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = A->getNNZs()-1;
  for ( ; last_nz >= 0 && A->d(last_nz) == 0; last_nz--);
  const int A_nnzs = last_nz + 1;
  const int A_last_index = last_nz >= 0? A->i(last_nz): -1;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, false)
    : unpacked_node::newFromNode(arg2F, b, false)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getNNZs()-1;
  for ( ; last_nz >= 0 && B->d(last_nz) == 0; last_nz--);
  const int B_nnzs = last_nz + 1;
  const int B_last_index = last_nz >= 0? B->i(last_nz): -1;

  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  int resultSize = A->getNNZs() + B->getNNZs() + 1 + 1;
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  int nnz = 0;
  int A_curr_index = 0;
  int B_curr_index = 0;
  for ( ; A_curr_index < A_nnzs && B_curr_index < B_nnzs; ) {
    // get a_i, a_d, b_i, b_d
    int a_i, a_d, b_i, b_d;
    a_i = A->i(A_curr_index);
    b_i = B->i(B_curr_index);
    if (a_i <= b_i) {
      a_d = A->d(A_curr_index);
      A_curr_index++;
    } else {
      a_d = 0;
    }
    if (a_i >= b_i) {
      b_d = B->d(B_curr_index);
      B_curr_index++;
    } else {
      b_d = 0;
    }

    if (a_d == 0 || b_d == 0) continue;

    // compute inter(a_d, b_d)
    int index = (a_d? a_i: b_i);
    node_handle down = compute_r(index, dwnLevel, a_d, b_d);

    // if inter is non-zero, add it to the new node
    if (down) {
      C->i_ref(nnz) = index;
      C->d_ref(nnz) = down;
      nnz++;
    }
  }
  if (B_ext_d != 0) {
    for ( ; A_curr_index < A_nnzs; A_curr_index++) {
      // do inter(a_i, b_ext_i)
      int index = A->i(A_curr_index);
      node_handle down = compute_r(index, dwnLevel, A->d(A_curr_index), B_ext_d);
      if (down) {
        C->i_ref(nnz) = index;
        C->d_ref(nnz) = down;
        nnz++;
      }
    }
  }
  if (A_ext_d != 0) {
    for ( ; B_curr_index < B_nnzs; B_curr_index++) {
      // do inter(a_ext_i, b_i)
      int index = B->i(B_curr_index);
      node_handle down = compute_r(index, dwnLevel, A_ext_d, B->d(B_curr_index));
      if (down) {
        C->i_ref(nnz) = index;
        C->d_ref(nnz) = down;
        nnz++;
      }
    }
  }
  if (A_ext_d != 0 && B_ext_d != 0) {
    int index = max_a_b_last_index+1;
    int down = compute_r(index, dwnLevel, A_ext_d, B_ext_d);
    if (down) {
      C->i_ref(nnz) = index;
      C->d_ref(nnz) = down;
      C->markAsExtensible();
      nnz++;
    } else {
      C->markAsNotExtensible();
    }
  }
  C->shrinkSparse(nnz);

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // reduce and save result
  result = resF->createReducedNode(-1, C);
  saveResult(Key, a, b, result);

#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
#endif

  return result;
}
#endif


// ******************************************************************
// *                                                                *
// *                       inter_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_opname : public binary_opname {
  public:
    inter_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::inter_opname::inter_opname()
 : binary_opname("Intersection")
{
}

MEDDLY::binary_operation* 
MEDDLY::inter_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    (a1->isForRelations() != r->isForRelations()) ||
    (a2->isForRelations() != r->isForRelations()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != r->getEdgeLabeling())
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (r->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    if (r->isForRelations())
      return new inter_mxd(this, a1, a2, r);
    else
      return new inter_mdd(this, a1, a2, r);
  }

  if (r->getEdgeLabeling() == forest::EVPLUS) {
    if (r->isForRelations()) {
      throw error(error::NOT_IMPLEMENTED);
    }
    else {
      return new inter_max_evplus(this, a1, a2, r);
    }
  }

  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeIntersection()
{
  return new inter_opname;
}

