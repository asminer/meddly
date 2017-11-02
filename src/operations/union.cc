
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
#include "union.h"
#include "apply_base.h"

namespace MEDDLY {
  class union_mdd;
  class union_mxd;

  class union_opname;
};

// ******************************************************************
// *                                                                *
// *                        union_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mdd : public generic_binary_mdd {
  public:
    union_mdd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
};

MEDDLY::union_mdd::union_mdd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mdd(opcode, arg1, arg2, res)
{
  operationCommutes();
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
  public:
    union_mxd(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
#ifdef USE_XDDS
  virtual MEDDLY::node_handle compute(node_handle a, node_handle b);
#endif
};

MEDDLY::union_mxd::union_mxd(const binary_opname* opcode, 
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_mxd(opcode, arg1, arg2, res)
{
  operationCommutes();
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

#ifdef USE_XDDS
MEDDLY::node_handle 
MEDDLY::union_mxd::compute(node_handle a, node_handle b) 
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
  const int A_last_index = last_nz >= 0? A->i(last_nz): -1;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, false)
    : unpacked_node::newFromNode(arg2F, b, false)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getNNZs()-1;
  for ( ; last_nz >= 0 && B->d(last_nz) == 0; last_nz--);
  const int B_last_index = last_nz >= 0? B->i(last_nz): -1;

  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  int resultSize = A->getNNZs() + B->getNNZs() + 1 + 1;
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  int nnz = 0;
  int A_curr_index = 0;
  int B_curr_index = 0;
  for ( ; A_curr_index < A->getNNZs() && B_curr_index < B->getNNZs(); ) {
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

    MEDDLY_DCASSERT(a_d != 0 || b_d != 0);

    // compute union(a_d, b_d)
    int index = (a_d? a_i: b_i);
    node_handle down = compute_r(index, dwnLevel, a_d, b_d);

    // if union is non-zero, add it to the new node
    if (down) {
      C->i_ref(nnz) = index;
      C->d_ref(nnz) = down;
      nnz++;
    }
  }
  for ( ; A_curr_index < A->getNNZs(); A_curr_index++) {
    // do union(a_i, b_ext_i)
    int index = A->i(A_curr_index);
    node_handle down = compute_r(index, dwnLevel, A->d(A_curr_index), B_ext_d);
    if (down) {
      C->i_ref(nnz) = index;
      C->d_ref(nnz) = down;
      nnz++;
    }
  }
  for ( ; B_curr_index < B->getNNZs(); B_curr_index++) {
    // do union(a_ext_i, b_i)
    int index = B->i(B_curr_index);
    node_handle down = compute_r(index, dwnLevel, A_ext_d, B->d(B_curr_index));
    if (down) {
      C->i_ref(nnz) = index;
      C->d_ref(nnz) = down;
      nnz++;
    }
  }
  if (A->isExtensible() || B->isExtensible()) {
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
// *                       union_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::union_opname : public binary_opname {
  public:
    union_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::union_opname::union_opname()
 : binary_opname("Union")
{
}

MEDDLY::binary_operation* 
MEDDLY::union_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
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
      return new union_mxd(this, a1, a2, r);
    else
      return new union_mdd(this, a1, a2, r);
  }

  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeUnion()
{
  return new union_opname;
}

