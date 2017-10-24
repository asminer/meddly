
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
#include "apply_base.h"

// #define TRACE_ALL_OPS
// #define DISABLE_CACHE

// ******************************************************************
// *                                                                *
// *                   generic_binary_mdd methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_mdd::generic_binary_mdd(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : binary_operation(code, 2, 1, arg1, arg2, res)
{
}

MEDDLY::generic_binary_mdd::~generic_binary_mdd()
{
}

#ifndef USE_NODE_STATUS
bool MEDDLY::generic_binary_mdd::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}
#else
MEDDLY::forest::node_status
MEDDLY::generic_binary_mdd::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = arg1F->getNodeStatus(data[0]);
  MEDDLY::forest::node_status b = arg2F->getNodeStatus(data[1]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[2]);

  if (a == MEDDLY::forest::DEAD ||
      b == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      b == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif


void MEDDLY::generic_binary_mdd::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::generic_binary_mdd ::showEntry(output &strm, const node_handle *data) const
{
  strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) 
       << "): " << long(data[2]) << "]";
}

void MEDDLY::generic_binary_mdd::computeDDEdge(const dd_edge &a, const dd_edge &b, 
  dd_edge &c)
{
  node_handle cnode = compute(a.getNode(), b.getNode());
  c.set(cnode);
#ifdef TRACE_ALL_OPS
  printf("completed %s(%d, %d) = %d\n", 
    getName(), a.getNode(), b.getNode(), cnode);
#endif
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}


MEDDLY::node_handle 
MEDDLY::generic_binary_mdd::compute(node_handle a, node_handle b)
{
  node_handle result = 0;
  if (checkTerminals(a, b, result))
    return result;

  compute_table::search_key* Key = findResult(a, b, result);
  if (0==Key) return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const int resultLevel = MAX(aLevel, bLevel);
  const int resultSize = resF->getLevelSize(resultLevel);

  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel) 
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, true)
    : unpacked_node::newFromNode(arg1F, a, true)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, true)
    : unpacked_node::newFromNode(arg2F, b, true)
  ;

  MEDDLY_DCASSERT(A->isFull() && resultSize == A->getSize());
  MEDDLY_DCASSERT(B->isFull() && resultSize == B->getSize());
  MEDDLY_DCASSERT(C->isFull() && resultSize == C->getSize());

  // do computation
  for (int i=0; i<resultSize; i++) {
    C->d_ref(i) = compute(A->d(i), B->d(i));
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // reduce and save result
  result = resF->createReducedNode(-1, C);
  saveResult(Key, a, b, result);

#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
  fflush(stdout);
#endif
  return result;
}



// ******************************************************************
// *                                                                *
// *                   generic_binary_mxd methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_mxd::generic_binary_mxd(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : binary_operation(code, 2, 1, arg1, arg2, res)
{
  // data[0] : arg1
  // data[1] : arg2
  // data[2] : result
}

MEDDLY::generic_binary_mxd::~generic_binary_mxd()
{
}

#ifndef USE_NODE_STATUS
bool MEDDLY::generic_binary_mxd::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}
#else
MEDDLY::forest::node_status
MEDDLY::generic_binary_mxd::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = arg1F->getNodeStatus(data[0]);
  MEDDLY::forest::node_status b = arg2F->getNodeStatus(data[1]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[2]);

  if (a == MEDDLY::forest::DEAD ||
      b == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      b == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif

void MEDDLY::generic_binary_mxd::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::generic_binary_mxd ::showEntry(output &strm, const node_handle *data) const
{
  strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) 
       << "): " << long(data[2]) << "]";
}

void MEDDLY::generic_binary_mxd::computeDDEdge(const dd_edge &a, const dd_edge &b, 
  dd_edge &c)
{
  node_handle cnode = compute(a.getNode(), b.getNode());
  c.set(cnode);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

#ifndef USE_XDDS

MEDDLY::node_handle 
MEDDLY::generic_binary_mxd::compute(node_handle a, node_handle b) 
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

  int resultSize = resF->getLevelSize(resultLevel);

  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel) 
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, true)
    : unpacked_node::newFromNode(arg1F, a, true)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, true)
    : unpacked_node::newFromNode(arg2F, b, true)
  ;

  // Do computation
  for (int j=0; j<resultSize; j++) {
    C->d_ref(j) = compute_r(j, resF->downLevel(resultLevel), A->d(j), B->d(j));
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // reduce and save result
  result = resF->createReducedNode(-1, C);
  saveResult(Key, a, b, result);

  /*
  printf("Saving %s un (%d, %d) = %d\n", getName(), a, b, result);
  printf("\tLevels %d and %d\n", aLevel, bLevel);
  */
#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
  fflush(stdout);
#endif

  return result;
}

MEDDLY::node_handle 
MEDDLY::generic_binary_mxd::compute_r(int in, int k, node_handle a, node_handle b)
{
  //  Compute for the primed levels.
  //
  MEDDLY_DCASSERT(k<0);

  node_handle result;
  /*
  //
  // Note - we cache the primed levels, but only when "safe"
  //
  compute_table::search_key* Key = findResult(a, b, result);
  if (0==Key) {
    printf("Found %s pr (%d, %d) = %d\n", getName(), a, b, result);
    printf("\tat level %d\n", k);
    return result;
  }
  */

  // bool canSaveResult = true;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const int resultSize = resF->getLevelSize(k);

  unpacked_node* C = unpacked_node::newFull(resF, k, resultSize);

  // Initialize readers
  unpacked_node *A = unpacked_node::useUnpackedNode();
  unpacked_node *B = unpacked_node::useUnpackedNode();

  if (aLevel == k) {
    A->initFromNode(arg1F, a, true);
  } else if (arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, k, a, true);
  } else {
    A->initIdentity(arg1F, k, in, a, true);
  }
  MEDDLY_DCASSERT(A->getSize() == C->getSize());

  if (bLevel == k) {
    B->initFromNode(arg2F, b, true);
  } else if (arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, k, b, true);
  } else {
    B->initIdentity(arg2F, k, in, b, true);
  }
  MEDDLY_DCASSERT(B->getSize() == C->getSize());

  // Do computation
  for (int j=0; j<resultSize; j++) {
    C->d_ref(j) = compute(A->d(j), B->d(j));
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // reduce 
  result = resF->createReducedNode(in, C);

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, %d, %d) = %d\n", getName(), in, a, b, result);
#endif

  return result;
}

#else

MEDDLY::node_handle 
MEDDLY::generic_binary_mxd::compute(node_handle a, node_handle b) 
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

  const int min_a_b_last_index = MIN(A_last_index, B_last_index);
  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  const node_handle C_ext_d =
    (A->isExtensible() || B->isExtensible())
    ? compute_r(max_a_b_last_index+1, dwnLevel, A_ext_d, B_ext_d)
    : 0;
  const bool C_is_extensible = (C_ext_d != 0);
  int resultSize = max_a_b_last_index + 1 + (C_is_extensible? 1: 0);
  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);

  //
  // Three loops to reduce the amount of checking needed:
  //
  // Loop 1: deal with indices in
  //         [0,min(a_last_index,b_last_index)]
  // Loop 2: deal with indices in
  //         [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  // 
  // Last index: result of A_ext_d and B_ext_d

  //
  // Loop 1: [0, min(a_last_index,b_last_index)]
  //
  int A_curr_index = 0;
  int B_curr_index = 0;
  int j = 0;
  for ( ; j <= min_a_b_last_index; j++) {
    const int a_d = ((j == A->i(A_curr_index))? A->d(A_curr_index++): 0);
    const int b_d = ((j == B->i(B_curr_index))? B->d(B_curr_index++): 0);
    C->d_ref(j) = compute_r(j, dwnLevel, a_d, b_d);
  }

  //
  // At most one of the next two loops will execute
  // Loop 2: [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  //
  for ( ; j <= A_last_index; j++) {
    const int a_d = ((j == A->i(A_curr_index))? A->d(A_curr_index++): 0);
    C->d_ref(j) = compute_r(j, dwnLevel, a_d, B_ext_d);
  }
  for ( ; j <= B_last_index; j++) {
    const int b_d = ((j == B->i(B_curr_index))? B->d(B_curr_index++): 0);
    C->d_ref(j) = compute_r(j, dwnLevel, A_ext_d, b_d);
  }
  MEDDLY_DCASSERT(j == max_a_b_last_index+1);

  //
  // Last index
  //
  MEDDLY_DCASSERT(!C->isExtensible());  // default: not extensible
  if (C_is_extensible) {
    C->d_ref(j) = C_ext_d;
    C->markAsExtensible();
  }

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

MEDDLY::node_handle 
MEDDLY::generic_binary_mxd::compute_r(int in, int k, node_handle a, node_handle b)
{
  //  Compute for the primed levels.
  //
  MEDDLY_DCASSERT(k<0);

  node_handle result = 0;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  // Initialize readers
  unpacked_node *A =
  (aLevel == k) 
    ? unpacked_node::newFromNode(arg1F, a, false)
    : arg1F->isFullyReduced()
    ? unpacked_node::newRedundant(arg1F, k, a, false)
    : unpacked_node::newIdentity(arg1F, k, in, a, false)
    ;
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = A->getNNZs()-1;
  for ( ; last_nz >= 0 && A->d(last_nz) == 0; last_nz--);
  const int A_last_index = last_nz >= 0? A->i(last_nz): -1;

  unpacked_node *B =
  (bLevel == k) 
    ? unpacked_node::newFromNode(arg2F, b, false)
    : arg2F->isFullyReduced()
    ? unpacked_node::newRedundant(arg2F, k, b, false)
    : unpacked_node::newIdentity(arg2F, k, in, b, false)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getNNZs()-1;
  for ( ; last_nz >= 0 && B->d(last_nz) == 0; last_nz--);
  const int B_last_index = last_nz >= 0? B->i(last_nz): -1;

  const int min_a_b_last_index = MIN(A_last_index, B_last_index);
  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  const node_handle C_ext_d =
    (A->isExtensible() || B->isExtensible())
    ? compute(A_ext_d, B_ext_d)
    : 0;
  const bool C_is_extensible = (C_ext_d != 0);
  int resultSize = max_a_b_last_index + 1 + (C_is_extensible? 1: 0);

  if (resultSize > 0) {
    unpacked_node* C = unpacked_node::newFull(resF, k, resultSize);

    // Loop 1: [0, min(a_last_index,b_last_index)]
    //
    int A_curr_index = 0;
    int B_curr_index = 0;
    int j = 0;
    for ( ; j <= min_a_b_last_index; j++) {
      const int a_d = ((j == A->i(A_curr_index))? A->d(A_curr_index++): 0);
      const int b_d = ((j == B->i(B_curr_index))? B->d(B_curr_index++): 0);
      C->d_ref(j) = compute(a_d, b_d);
    }

    // At most one of the next two loops will execute
    // Loop 2: [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
    //
    for ( ; j <= A_last_index; j++) {
      const int a_d = ((j == A->i(A_curr_index))? A->d(A_curr_index++): 0);
      C->d_ref(j) = compute(a_d, B_ext_d);
    }
    for ( ; j <= B_last_index; j++) {
      const int b_d = ((j == B->i(B_curr_index))? B->d(B_curr_index++): 0);
      C->d_ref(j) = compute(A_ext_d, b_d);
    }
    MEDDLY_DCASSERT(j == max_a_b_last_index+1);

    // Last index
    //
    MEDDLY_DCASSERT(!C->isExtensible());  // default: not extensible
    if (C_is_extensible) {
      C->d_ref(j) = C_ext_d;
      C->markAsExtensible();
    }

    // reduce 
    result = resF->createReducedNode(in, C);
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, %d, %d) = %d\n", getName(), in, a, b, result);
  fflush(stdout);
#endif

  return result;
}

#endif

// ******************************************************************
// *                                                                *
// *                 generic_binbylevel_mxd methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binbylevel_mxd
::generic_binbylevel_mxd(const binary_opname* code, expert_forest* arg1, 
  expert_forest* arg2, expert_forest* res)
 : binary_operation(code, 3, 1, arg1, arg2, res)
{
  can_commute = false;
}

MEDDLY::generic_binbylevel_mxd::~generic_binbylevel_mxd()
{
}

#ifndef USE_NODE_STATUS
bool MEDDLY::generic_binbylevel_mxd::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[2]) ||
         resF->isStale(data[3]);
}
#else
MEDDLY::forest::node_status
MEDDLY::generic_binbylevel_mxd::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = arg1F->getNodeStatus(data[1]);
  MEDDLY::forest::node_status b = arg2F->getNodeStatus(data[2]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[3]);

  if (a == MEDDLY::forest::DEAD ||
      b == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      b == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif

void MEDDLY::generic_binbylevel_mxd::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[2]);
  resF->uncacheNode(data[3]);
}

void
MEDDLY::generic_binbylevel_mxd
::showEntry(output &strm, const node_handle *data) const
{
  strm << "[" << getName() << "(" << long(data[0]) << ", " << long(data[1]) 
       << ", " << long(data[2]) << "): " << long(data[3]) << "]";
}

void MEDDLY::generic_binbylevel_mxd
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  node_handle result = compute(
    resF->getDomain()->getNumVariables(), a.getNode(), b.getNode()
  );
  c.set(result);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

MEDDLY::node_handle 
MEDDLY::generic_binbylevel_mxd::compute(int level, node_handle a, node_handle b) 
{
  return compute_r(-1, level, a, b);
}

// #ifndef USE_XDDS
#if 1
MEDDLY::node_handle 
MEDDLY::generic_binbylevel_mxd
::compute_r(int in, int resultLevel, node_handle a, node_handle b)
{
  node_handle result = 0;
  if (0==resultLevel) {
    if (checkTerminals(a, b, result))
      return result;
  }

  //
  // Note - we cache the primed levels, but only when "safe"
  //
  compute_table::search_key* Key = findResult(resultLevel, a, b, result);
  if (0==Key) return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultSize = resF->getLevelSize(resultLevel);

  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);

  bool canSaveResult = true;

  // Initialize readers
  unpacked_node* A = unpacked_node::useUnpackedNode();
  unpacked_node* B = unpacked_node::useUnpackedNode();

  if (aLevel == resultLevel) {
    A->initFromNode(arg1F, a, true);
  } else if (resultLevel>0 || arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, resultLevel, a, true);
  } else {
    A->initIdentity(arg1F, resultLevel, in, a, true);
    canSaveResult = false;
  }

  MEDDLY_DCASSERT(!A->isExtensible());

  if (bLevel == resultLevel) {
    B->initFromNode(arg2F, b, true);
  } else if (resultLevel>0 || arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, resultLevel, b, true);
  } else {
    B->initIdentity(arg2F, resultLevel, in, b, true);
    canSaveResult = false;
  }

  MEDDLY_DCASSERT(!B->isExtensible());

  // Do computation
  int nextLevel = resF->downLevel(resultLevel);
  int nnz = 0;
  for (int j=0; j<resultSize; j++) {
    node_handle d = compute_r(j, nextLevel, A->d(j), B->d(j));
    C->d_ref(j) = d;
    if (d) nnz++;
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // reduce
  result = resF->createReducedNode(in, C);

  // save result in compute table, when we can
  if (resultLevel<0 && 1==nnz) canSaveResult = false;
  if (canSaveResult)  saveResult(Key, resultLevel, a, b, result);
  else                doneCTkey(Key);

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, %d, %d) = %d\n", getName(), in, a, b, result);
  fflush(stdout);
#endif

  return result;
}

#else

MEDDLY::node_handle 
MEDDLY::generic_binbylevel_mxd
::compute_r(int in, int resultLevel, node_handle a, node_handle b)
{
  node_handle result = 0;
  // TODO: can be made more efficient (for example, a==a)
  if (0==resultLevel) {
    if (checkTerminals(a, b, result))
      return result;
  }

  //
  // Note - we cache the primed levels, but only when "safe"
  //
  compute_table::search_key* Key = findResult(resultLevel, a, b, result);
  if (0==Key) return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  const int dwnLevel = resF->downLevel(resultLevel);

  // Initialize readers
  unpacked_node* A = unpacked_node::useUnpackedNode();
  unpacked_node* B = unpacked_node::useUnpackedNode();

  bool canSaveResult = true;

  if (aLevel == resultLevel) {
    A->initFromNode(arg1F, a, false);
  } else if (resultLevel>0 || arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, resultLevel, a, false);
  } else {
    A->initIdentity(arg1F, resultLevel, in, a, false);
    canSaveResult = false;
  }
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;

  if (bLevel == resultLevel) {
    B->initFromNode(arg2F, b, false);
  } else if (resultLevel>0 || arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, resultLevel, b, false);
  } else {
    B->initIdentity(arg2F, resultLevel, in, b, false);
    canSaveResult = false;
  }
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;

  const int A_last_index = A->i(A->getNNZs()-1);
  const int B_last_index = B->i(B->getNNZs()-1);
  const int min_a_b_last_index = MIN(A_last_index, B_last_index);
  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  const node_handle C_ext_d =
    (A->isExtensible() || B->isExtensible())
    ? compute_r(max_a_b_last_index+1, dwnLevel, A_ext_d, B_ext_d)
    : 0;
  const bool C_is_extensible = (C_ext_d != 0);
  int resultSize = max_a_b_last_index + 1 + (C_is_extensible? 1: 0);
  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);

  //
  // Three loops to reduce the amount of checking needed:
  //
  // Loop 1: deal with indices in
  //         [0,min(a_last_index,b_last_index)]
  // Loop 2: deal with indices in
  //         [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  // 
  // Last index: result of A_ext_d and B_ext_d

  //
  // Loop 1: [0, min(a_last_index,b_last_index)]
  //
  int A_curr_index = 0;
  int B_curr_index = 0;
  int j = 0;
  int nnz = 0;
  for ( ; j <= min_a_b_last_index; j++) {
    const node_handle a_d = ((j == A->i(A_curr_index))? A->d(A_curr_index++): 0);
    const node_handle b_d = ((j == B->i(B_curr_index))? B->d(B_curr_index++): 0);
    node_handle d = compute_r(j, dwnLevel, a_d, b_d);
    C->d_ref(j) = d;
    if (d) nnz++;
  }

  //
  // At most one of the next two loops will execute
  // Loop 2: [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  //
  for ( ; j <= A_last_index; j++) {
    const node_handle a_d = ((j == A->i(A_curr_index))? A->d(A_curr_index++): 0);
    node_handle d = compute_r(j, dwnLevel, a_d, B_ext_d);
    C->d_ref(j) = d;
    if (d) nnz++;
  }
  for ( ; j <= B_last_index; j++) {
    const node_handle b_d = ((j == B->i(B_curr_index))? B->d(B_curr_index++): 0);
    node_handle d = compute_r(j, dwnLevel, A_ext_d, b_d);
    C->d_ref(j) = d;
    if (d) nnz++;
  }
  MEDDLY_DCASSERT(j == max_a_b_last_index+1);

  //
  // Last index
  //
  MEDDLY_DCASSERT(!C->isExtensible());  // default: not extensible
  if (C_is_extensible) {
    C->d_ref(j) = C_ext_d;
    C->markAsExtensible();
    nnz++;
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // reduce
  result = resF->createReducedNode(in, C);

  // save result in compute table, when we can
  if (resultLevel<0 && 1==nnz) canSaveResult = false;
  if (canSaveResult)  saveResult(Key, resultLevel, a, b, result);
  else                doneCTkey(Key);

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, %d, %d) = %d\n", getName(), in, a, b, result);
  fflush(stdout);
#endif

  return result;
}

#endif


// ******************************************************************
// *                                                                *
// *                   generic_binary_ev  methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_ev::generic_binary_ev(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : binary_operation(code, 4, 2, arg1, arg2, res)
{
  can_commute = false;
}

MEDDLY::generic_binary_ev::~generic_binary_ev()
{
}

#ifndef USE_NODE_STATUS
bool MEDDLY::generic_binary_ev::isStaleEntry(const node_handle* data)
{
  bool a = arg1F->isStale(data[1]);
  bool b = arg2F->isStale(data[3]);
  bool c = resF->isStale(data[5]);

  return (a | b | c);
}

#else

MEDDLY::forest::node_status
MEDDLY::generic_binary_ev::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = arg1F->getNodeStatus(data[1]);
  MEDDLY::forest::node_status b = arg2F->getNodeStatus(data[3]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[5]);

  if (a == MEDDLY::forest::DEAD ||
      b == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      b == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif

void MEDDLY::generic_binary_ev::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[3]);
  resF->uncacheNode(data[5]);
}

// ******************************************************************
// *                                                                *
// *                 generic_binary_evplus  methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evplus::generic_binary_evplus(const binary_opname* code,
  expert_forest* arg1, expert_forest* arg2, expert_forest* res)
  : generic_binary_ev(code, arg1, arg2, res)
{
}

MEDDLY::generic_binary_evplus::~generic_binary_evplus()
{
}

void MEDDLY::generic_binary_evplus
::showEntry(output &strm, const node_handle *data) const
{
  strm << "[" << getName() << "(<" << long(data[0]) << ":" << long(data[1]) 
       << ">, <" << long(data[2]) << ":" << long(data[3]) << ">): <"
       << long(data[4]) << ":" << long(data[5]) << ">]";
}

void MEDDLY::generic_binary_evplus
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  node_handle result;
  int ev, aev, bev;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(result, ev);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

void MEDDLY::generic_binary_evplus
::compute(int aev, node_handle a, int bev, node_handle b, 
  int& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  compute_table::search_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const int resultLevel = aLevel > bLevel? aLevel: bLevel;
  const int resultSize = resF->getLevelSize(resultLevel);

  // Initialize result
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel) 
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0, a, true)
    : unpacked_node::newFromNode(arg1F, a, true)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0, b, true)
    : unpacked_node::newFromNode(arg2F, b, true)
  ;

  MEDDLY_DCASSERT(!A->isExtensible() && !B->isExtensible());


  // do computation
  for (int i=0; i<resultSize; i++) {
    int ev;
    node_handle ed;
    compute(aev + A->ei(i), A->d(i), 
            bev + B->ei(i), B->d(i), 
            ev, ed);
    nb->d_ref(i) = ed;
    nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(-1, nb, cev, cl);
  c = cl;

  // Add to CT
  saveResult(Key, aev, a, bev, b, cev, c);
}


// ******************************************************************
// *                                                                *
// *                 generic_binary_evtimes methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evtimes
::generic_binary_evtimes(const binary_opname* code, expert_forest* arg1, 
  expert_forest* arg2, expert_forest* res)
: generic_binary_ev(code, arg1, arg2, res)
{
}

MEDDLY::generic_binary_evtimes::~generic_binary_evtimes()
{
}

void MEDDLY::generic_binary_evtimes
::showEntry(output &strm, const node_handle *data) const
{
  float ev0;
  float ev2;
  float ev4;
  compute_table::readEV(data+0, ev0);
  compute_table::readEV(data+2, ev2);
  compute_table::readEV(data+4, ev4);
  strm << "[" << getName() << "(<" << ev0 << ":" << long(data[1]) 
       << ">, <" << ev2 << ":" << long(data[3]) << ">): <"
       << ev4 << ":" << long(data[5]) << ">]";
}

void MEDDLY::generic_binary_evtimes
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  node_handle result; 
  float ev, aev, bev;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(result, ev);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

void MEDDLY::generic_binary_evtimes
::compute(float aev, node_handle a, float bev, node_handle b, 
  float& cev, node_handle& c)
{
  // Compute for the unprimed levels.
  //

  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

#ifndef DISABLE_CACHE
  compute_table::search_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;
#endif

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  // Initialize result
  const int resultLevel = ABS(topLevel(aLevel, bLevel));
  MEDDLY_DCASSERT(resultLevel>0);
  const int resultSize = resF->getLevelSize(resultLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel) 
    ? unpacked_node::newRedundant(arg1F, resultLevel, 1.0f, a, true)
    : unpacked_node::newFromNode(arg1F, a, true)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 1.0f, b, true)
    : unpacked_node::newFromNode(arg2F, b, true)
  ;

  MEDDLY_DCASSERT(!A->isExtensible() && !B->isExtensible());

  // do computation
  for (int i=0; i<resultSize; i++) {
    float ev;
    node_handle ed;
    compute_k(
        i, -resultLevel,
        aev * A->ef(i), A->d(i), 
        bev * B->ef(i), B->d(i), 
        ev, ed);
    nb->d_ref(i) = ed;
    nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(-1, nb, cev, cl);
  c = cl;

#ifndef DISABLE_CACHE
  // Add to CT
  saveResult(Key, aev, a, bev, b, cev, c);
#endif
}

void MEDDLY::generic_binary_evtimes
::compute_k(
  int in, int resultLevel,
  float aev, node_handle a,
  float bev, node_handle b,
  float& cev, node_handle& c)
{
  // Compute for the primed levels.
  //
  MEDDLY_DCASSERT(resultLevel < 0);

  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  // Initialize result
  const int resultSize = resF->getLevelSize(resultLevel);
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = unpacked_node::useUnpackedNode();
  unpacked_node *B = unpacked_node::useUnpackedNode();

  if (aLevel == resultLevel) {
    A->initFromNode(arg1F, a, true);
  } else if (arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, resultLevel, 1.0f, a, true);
  } else {
    A->initIdentity(arg1F, resultLevel, in, 1.0f, a, true);
  }

  if (bLevel == resultLevel) {
    B->initFromNode(arg2F, b, true);
  } else if (arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, resultLevel, 1.0f, b, true);
  } else {
    B->initIdentity(arg2F, resultLevel, in, 1.0f, b, true);
  }

  MEDDLY_DCASSERT(!A->isExtensible() && !B->isExtensible());

  // do computation
  for (int i=0; i<resultSize; i++) {
    float ev;
    node_handle ed;
    compute(
        aev * A->ef(i), A->d(i),
        bev * B->ef(i), B->d(i), 
        ev, ed);
    nb->d_ref(i) = ed;
    nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::recycle(B);
  unpacked_node::recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(in, nb, cev, cl);
  c = cl;
}
