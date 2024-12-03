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
#include "../forests/mt.h"
#include "apply_base.h"

// #define TRACE_ALL_OPS
// #define DISABLE_CACHE


// ******************************************************************
// *                                                                *
// *                   generic_binary_mdd methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_mdd::generic_binary_mdd(binary_list& code,
  forest* arg1, forest* arg2, forest* res)
  : binary_operation(code, 1, arg1, arg2, res)
{
    ct_entry_type* et = new ct_entry_type(code.getName(), "NN:N");
    et->setForestForSlot(0, arg1);
    et->setForestForSlot(1, arg2);
    et->setForestForSlot(3, res);
    registerEntryType(0, et);
    buildCTs();
}

MEDDLY::generic_binary_mdd::~generic_binary_mdd()
{
}

void MEDDLY::generic_binary_mdd::computeDDEdge(const dd_edge &a, const dd_edge &b,
  dd_edge &c, bool userFlag)
{
#ifdef TRACE_ALL_OPS
    printf("computing Top %s(%d, %d)\n",
        getName(), a.getNode(), b.getNode()
    );
#endif
    node_handle cnode = compute(a.getNode(), b.getNode());
    const int num_levels = resF->getMaxLevelIndex();
    if ( userFlag && resF->isQuasiReduced() && cnode != resF->getTransparentNode()
        && resF->getNodeLevel(cnode) < num_levels)
    {
        node_handle temp = ((mt_forest*)resF)->makeNodeAtLevel(num_levels, cnode);
        resF->unlinkNode(cnode);
        cnode = temp;
    }
#ifdef TRACE_ALL_OPS
    printf("completed Top %s(%d, %d) = %d\n",
        getName(), a.getNode(), b.getNode(), cnode);
#endif
    c.set(cnode);
#ifdef DEVELOPMENT_CODE
    resF->validateIncounts(true);
#endif
}


MEDDLY::node_handle
MEDDLY::generic_binary_mdd::compute(node_handle a, node_handle b)
{
    node_handle result = 0;
    if (checkTerminals(a, b, result)) {
        return result;
    }

    ct_entry_key* Key = findResult(a, b, result);
    if (!Key) {
#ifdef TRACE_ALL_OPS
        printf("computing %s(%d, %d), got %d from cache\n", getName(), a, b, result);
        fflush(stdout);
#endif

        return result;
    }

#ifdef TRACE_ALL_OPS
    printf("computing %s(%d, %d)\n", getName(), a, b);
    fflush(stdout);
#endif

#ifdef ALLOW_EXTENSIBLE
    // Get level information
    const int aLevel = arg1F->getNodeLevel(a);
    const int bLevel = arg2F->getNodeLevel(b);
    const int resultLevel = MAX(aLevel, bLevel);

    result = resF->isExtensibleLevel(resultLevel)
        ? compute_ext(a, b)
        : compute_normal(a, b);
#else
    result = compute_normal(a, b);
#endif

    // save result
    saveResult(Key, a, b, result);

#ifdef TRACE_ALL_OPS
    printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
    fflush(stdout);
#endif

    return result;
}

MEDDLY::node_handle
MEDDLY::generic_binary_mdd::compute_normal(node_handle a, node_handle b)
{
  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  const int resultLevel = MAX(aLevel, bLevel);
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  MEDDLY_DCASSERT(!resF->isExtensibleLevel(resultLevel));

  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);
  MEDDLY_DCASSERT(!C->isExtensible());

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY)
  ;
  MEDDLY_DCASSERT(!A->isExtensible());

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY)
  ;
  MEDDLY_DCASSERT(!B->isExtensible());

  MEDDLY_DCASSERT(A->isFull() && resultSize == A->getSize());
  MEDDLY_DCASSERT(B->isFull() && resultSize == B->getSize());
  MEDDLY_DCASSERT(C->isFull() && resultSize == C->getSize());

  // do computation
  for (unsigned i=0; i<resultSize; i++) {
      C->setFull(i, compute(A->down(i), B->down(i)));
    // C->d_ref(i) = compute(A->down(i), B->down(i));
  }

  if (resF->isQuasiReduced()) {
    int nextLevel = resultLevel - 1;
    for (unsigned i = 0; i < C->getSize(); i++) {
      if (resF->getNodeLevel(C->down(i)) < nextLevel) {
        node_handle temp = ((mt_forest*)resF)->makeNodeAtLevel(nextLevel, C->down(i));
        resF->unlinkNode(C->down(i));
        C->setFull(i, temp);
        // C->d_ref(i) = temp;
      }
    }
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce and return result
  edge_value ev;
  node_handle result;
  resF->createReducedNode(C, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
  return result;
}


#ifdef ALLOW_EXTENSIBLE
MEDDLY::node_handle
MEDDLY::generic_binary_mdd::compute_ext(node_handle a, node_handle b)
{
  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = topLevel(aLevel, bLevel);

  MEDDLY_DCASSERT(resF->isExtensibleLevel(resultLevel));

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, SPARSE_ONLY)
    : arg1F->newUnpacked(a, SPARSE_ONLY)
    ;
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = A->getSize()-1;
  for ( ; last_nz >= 0 && A->down(last_nz) == 0; last_nz--);
  const int A_last_index = last_nz >= 0? A->index(last_nz): -1;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, SPARSE_ONLY)
    : arg2F->newUnpacked(b, SPARSE_ONLY)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getSize()-1;
  for ( ; last_nz >= 0 && B->down(last_nz) == 0; last_nz--);
  const int B_last_index = last_nz >= 0? B->index(last_nz): -1;

  const int min_a_b_last_index = MIN(A_last_index, B_last_index);
  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  const node_handle C_ext_d =
    (A->isExtensible() || B->isExtensible())
    ? compute(A_ext_d, B_ext_d)
    : 0;
  const bool C_is_extensible = (C_ext_d != 0);

  //
  // Three loops to reduce the amount of checking needed:
  //
  // Loop 1: deal with indices in
  //         [0,min(a_last_index,b_last_index)]
  // Loop 2: deal with indices in
  //         [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  //
  // Last index: result of A_ext_d and B_ext_d

  int resultSize = max_a_b_last_index + 1 + (C_is_extensible? 1: 0);
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  //
  // Loop 1: [0, min(a_last_index,b_last_index)]
  //
  unsigned A_curr_index = 0;
  unsigned B_curr_index = 0;
  int j = 0;
  unsigned nnz = 0;
  for ( ; j <= min_a_b_last_index; j++) {
    const node_handle a_d = ((j == A->index(A_curr_index))? A->down(A_curr_index++): 0);
    const node_handle b_d = ((j == B->index(B_curr_index))? B->down(B_curr_index++): 0);
    node_handle down = compute(a_d, b_d);
    if (down) {
      C->setSparse(nnz, j, down);
      // C->d_ref(nnz) = down;
      // C->i_ref(nnz) = j;
      nnz++;
    }
  }

  //
  // At most one of the next two loops will execute
  // Loop 2: [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  //
  for ( ; j <= A_last_index; j++) {
    const node_handle a_d = ((j == A->index(A_curr_index))? A->down(A_curr_index++): 0);
    node_handle down = compute(a_d, B_ext_d);
    if (down) {
      C->setSparse(nnz, j, down);
      // C->d_ref(nnz) = down;
      // C->i_ref(nnz) = j;
      nnz++;
    }
  }
  for ( ; j <= B_last_index; j++) {
    const node_handle b_d = ((j == B->index(B_curr_index))? B->down(B_curr_index++): 0);
    node_handle down = compute(A_ext_d, b_d);
    if (down) {
      C->setSparse(nnz, j, down);
      // C->d_ref(nnz) = down;
      // C->i_ref(nnz) = j;
      nnz++;
    }
  }
  MEDDLY_DCASSERT(j == max_a_b_last_index+1);

  //
  // Last index
  //
  MEDDLY_DCASSERT(!C->isExtensible());  // default: not extensible
  if (C_is_extensible) {
    MEDDLY_DCASSERT(C_ext_d);
    C->setSparse(nnz, j, C_ext_d);
    // C->d_ref(nnz) = C_ext_d;
    // C->i_ref(nnz) = j;
    nnz++;
    C->markAsExtensible();
  }
  C->shrink(nnz);

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  if (resF->isQuasiReduced()) {
    int nextLevel = resultLevel - 1;
    for (unsigned i = 0; i < C->getSize(); i++) {
      if (resF->getNodeLevel(C->down(i)) < nextLevel) {
        node_handle temp = ((mt_forest*)resF)->makeNodeAtLevel(nextLevel, C->down(i));
        resF->unlinkNode(C->down(i));
        // C->d_ref(i) = temp;
        C->setSparse(i, C->index(i), temp);
      }
    }
  }

  // reduce and return result
  edge_value ev;
  node_handle result;
  resF->createReducedNode(C, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
  return result;
}

#endif // ALLOW_EXTENSIBLE


// ******************************************************************
// *                                                                *
// *                   generic_binary_mxd methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_mxd::generic_binary_mxd(binary_list& code,
  forest* arg1, forest* arg2, forest* res)
  : binary_operation(code, 1, arg1, arg2, res)
{
  ct_entry_type* et = new ct_entry_type(code.getName(), "NN:N");
  et->setForestForSlot(0, arg1);
  et->setForestForSlot(1, arg2);
  et->setForestForSlot(3, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::generic_binary_mxd::~generic_binary_mxd()
{
}

void MEDDLY::generic_binary_mxd::computeDDEdge(const dd_edge &a, const dd_edge &b,
  dd_edge &c, bool userFlag)
{
  node_handle cnode = compute(a.getNode(), b.getNode());
  c.set(cnode);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

MEDDLY::node_handle
MEDDLY::generic_binary_mxd::compute(node_handle a, node_handle b)
{
  //  Compute for the unprimed levels.
  //
  node_handle result = 0;
  if (checkTerminals(a, b, result))
    return result;

  ct_entry_key* Key = findResult(a, b, result);
  if (0==Key) return result;

#ifdef ALLOW_EXTENSIBLE
  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  const int resultLevel = ABS(topLevel(aLevel, bLevel));
  result =
    resF->isExtensibleLevel(resultLevel)
    ? compute_ext(a, b)
    : compute_normal(a, b);
#else
  result = compute_normal(a, b);
#endif

  // save result
  saveResult(Key, a, b, result);

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
  ct_entry_key* Key = findResult(a, b, result);
  if (0==Key) {
  printf("Found %s pr (%d, %d) = %d\n", getName(), a, b, result);
  printf("\tat level %d\n", k);
  return result;
  }
  */

#ifdef ALLOW_EXTENSIBLE
  result =
    resF->isExtensibleLevel(k)
    ? compute_r_ext(in, k, a, b)
    : compute_r_normal(in, k, a, b);
#else
  result = compute_r_normal(in, k, a, b);
#endif

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, level %d, %d, %d) = %d\n", getName(), in, k, a, b, result);
  fflush(stdout);
#endif


  return result;
}

MEDDLY::node_handle
MEDDLY::generic_binary_mxd::compute_normal(node_handle a, node_handle b)
{
  node_handle result = 0;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));
  unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  MEDDLY_DCASSERT(!resF->isExtensibleLevel(resultLevel));

  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);
  MEDDLY_DCASSERT(!C->isExtensible());

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY)
  ;
  MEDDLY_DCASSERT(!A->isExtensible());

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY)
  ;
  MEDDLY_DCASSERT(!B->isExtensible());

  // Do computation
  for (unsigned j=0; j<resultSize; j++) {
    C->setFull(j,
        compute_r(j, MXD_levels::downLevel(resultLevel), A->down(j), B->down(j))
    );
    // C->d_ref(j) = compute_r(j, MXD_levels::downLevel(resultLevel), A->down(j), B->down(j));
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce and return result
  edge_value ev;
  resF->createReducedNode(C, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
  return result;
}

#ifdef ALLOW_EXTENSIBLE
MEDDLY::node_handle
MEDDLY::generic_binary_mxd::compute_ext(node_handle a, node_handle b)
{
  node_handle result = 0;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));
  const int dwnLevel = MXD_levels::downLevel(resultLevel);

  MEDDLY_DCASSERT(resF->isExtensibleLevel(resultLevel));

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, SPARSE_ONLY)
    : arg1F->newUnpacked(a, SPARSE_ONLY)
    ;
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = A->getSize()-1;
  for ( ; last_nz >= 0 && A->down(last_nz) == 0; last_nz--);
  const int A_last_index = last_nz >= 0? A->index(last_nz): -1;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, SPARSE_ONLY)
    : arg2F->newUnpacked(b, SPARSE_ONLY)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getSize()-1;
  for ( ; last_nz >= 0 && B->down(last_nz) == 0; last_nz--);
  const int B_last_index = last_nz >= 0? B->index(last_nz): -1;

  const int min_a_b_last_index = MIN(A_last_index, B_last_index);
  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  const node_handle C_ext_d =
    (A->isExtensible() || B->isExtensible())
    ? compute_r(max_a_b_last_index+1, dwnLevel, A_ext_d, B_ext_d)
    : 0;
  const bool C_is_extensible = (C_ext_d != 0);

  //
  // Three loops to reduce the amount of checking needed:
  //
  // Loop 1: deal with indices in
  //         [0,min(a_last_index,b_last_index)]
  // Loop 2: deal with indices in
  //         [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  //
  // Last index: result of A_ext_d and B_ext_d

  int resultSize = max_a_b_last_index + 1 + (C_is_extensible? 1: 0);
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  //
  // Loop 1: [0, min(a_last_index,b_last_index)]
  //
  unsigned A_curr_index = 0;
  unsigned B_curr_index = 0;
  int j = 0;
  unsigned nnz = 0;
  for ( ; j <= min_a_b_last_index; j++) {
    const node_handle a_d = ((j == A->index(A_curr_index))? A->down(A_curr_index++): 0);
    const node_handle b_d = ((j == B->index(B_curr_index))? B->down(B_curr_index++): 0);
    node_handle down = compute_r(j, dwnLevel, a_d, b_d);
    if (down) {
      C->setSparse(nnz, j, down);
      // C->d_ref(nnz) = down;
      // C->i_ref(nnz) = j;
      nnz++;
    }
  }

  //
  // At most one of the next two loops will execute
  // Loop 2: [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
  //
  for ( ; j <= A_last_index; j++) {
    const node_handle a_d = ((j == A->index(A_curr_index))? A->down(A_curr_index++): 0);
    node_handle down = compute_r(j, dwnLevel, a_d, B_ext_d);
    if (down) {
      C->setSparse(nnz, j, down);
      // C->d_ref(nnz) = down;
      // C->i_ref(nnz) = j;
      nnz++;
    }
  }
  for ( ; j <= B_last_index; j++) {
    const node_handle b_d = ((j == B->index(B_curr_index))? B->down(B_curr_index++): 0);
    node_handle down = compute_r(j, dwnLevel, A_ext_d, b_d);
    if (down) {
      C->setSparse(nnz, j, down);
      // C->d_ref(nnz) = down;
      // C->i_ref(nnz) = j;
      nnz++;
    }
  }
  MEDDLY_DCASSERT(j == max_a_b_last_index+1);

  //
  // Last index
  //
  MEDDLY_DCASSERT(!C->isExtensible());  // default: not extensible
  if (C_is_extensible) {
    MEDDLY_DCASSERT(C_ext_d);
    C->setSparse(nnz, j, C_ext_d);
    // C->d_ref(nnz) = C_ext_d;
    // C->i_ref(nnz) = j;
    nnz++;
    C->markAsExtensible();
  }
  C->shrink(nnz);

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce and return result
  edge_value ev;
  node_handle result;
  resF->createReducedNode(C, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
  return result;
}
#endif // ALLOW_EXTENSIBLE


MEDDLY::node_handle
MEDDLY::generic_binary_mxd::compute_r_normal(int in, int k, node_handle a, node_handle b)
{
  MEDDLY_DCASSERT(!resF->isExtensibleLevel(k));

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const unsigned resultSize = unsigned(resF->getLevelSize(k));

  unpacked_node* C = unpacked_node::newFull(resF, k, resultSize);

  // Initialize readers
  unpacked_node *A = unpacked_node::New(arg1F);
  unpacked_node *B = unpacked_node::New(arg2F);

  if (aLevel == k) {
    arg1F->unpackNode(A, a, FULL_ONLY);
  } else if (arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, k, a, FULL_ONLY);
  } else {
    A->initIdentity(arg1F, k, in, a, FULL_ONLY);
  }
  MEDDLY_DCASSERT(A->getSize() == C->getSize());

  if (bLevel == k) {
    arg2F->unpackNode(B, b, FULL_ONLY);
  } else if (arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, k, b, FULL_ONLY);
  } else {
    B->initIdentity(arg2F, k, in, b, FULL_ONLY);
  }
  MEDDLY_DCASSERT(B->getSize() == C->getSize());

  // Do computation
  for (unsigned j=0; j<resultSize; j++) {
    C->setFull(j, compute(A->down(j), B->down(j)));
    // C->d_ref(j) = compute(A->down(j), B->down(j));
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce
  edge_value ev;
  node_handle result;
  resF->createReducedNode(C, ev, result, in);
  MEDDLY_DCASSERT(ev.isVoid());
  return result;
}

#ifdef ALLOW_EXTENSIBLE
MEDDLY::node_handle
MEDDLY::generic_binary_mxd::compute_r_ext(int in, int k, node_handle a, node_handle b)
{
  MEDDLY_DCASSERT(resF->isExtensibleLevel(k));

  node_handle result = 0;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  // Initialize readers
  unpacked_node *A =
  (aLevel == k)
    ? arg1F->newUnpacked(a, SPARSE_ONLY)
    : arg1F->isFullyReduced()
    ? unpacked_node::newRedundant(arg1F, k, a, SPARSE_ONLY)
    : unpacked_node::newIdentity(arg1F, k, in, a, SPARSE_ONLY)
    ;
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = A->getSize()-1;
  for ( ; last_nz >= 0 && A->down(last_nz) == 0; last_nz--);
  const int A_last_index = last_nz >= 0? A->index(last_nz): -1;

  unpacked_node *B =
  (bLevel == k)
    ? arg2F->newUnpacked(b, SPARSE_ONLY)
    : arg2F->isFullyReduced()
    ? unpacked_node::newRedundant(arg2F, k, b, SPARSE_ONLY)
    : unpacked_node::newIdentity(arg2F, k, in, b, SPARSE_ONLY)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getSize()-1;
  for ( ; last_nz >= 0 && B->down(last_nz) == 0; last_nz--);
  const int B_last_index = last_nz >= 0? B->index(last_nz): -1;

  const int min_a_b_last_index = MIN(A_last_index, B_last_index);
  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  const node_handle C_ext_d =
    (A->isExtensible() || B->isExtensible())
    ? compute(A_ext_d, B_ext_d)
    : 0;
  const bool C_is_extensible = (C_ext_d != 0);
  int resultSize = max_a_b_last_index + 1 + (C_is_extensible? 1: 0);

  if (resultSize > 0) {
    unpacked_node* C = unpacked_node::newSparse(resF, k, resultSize);

    // Loop 1: [0, min(a_last_index,b_last_index)]
    //
    unsigned A_curr_index = 0;
    unsigned B_curr_index = 0;
    int j = 0;
    unsigned nnz = 0;
    for ( ; j <= min_a_b_last_index; j++) {
      const node_handle a_d = ((j == A->index(A_curr_index))? A->down(A_curr_index++): 0);
      const node_handle b_d = ((j == B->index(B_curr_index))? B->down(B_curr_index++): 0);
      node_handle down = compute(a_d, b_d);
      if (down) {
        C->setSparse(nnz, j, down);
        // C->d_ref(nnz) = down;
        // C->i_ref(nnz) = j;
        nnz++;
      }
    }

    // At most one of the next two loops will execute
    // Loop 2: [min(a_last_index,b_last_index)+1, max(a_last_index,b_last_index)]
    //
    for ( ; j <= A_last_index; j++) {
      const node_handle a_d = ((j == A->index(A_curr_index))? A->down(A_curr_index++): 0);
      node_handle down = compute(a_d, B_ext_d);
      if (down) {
        C->setSparse(nnz, j, down);
        // C->d_ref(nnz) = down;
        // C->i_ref(nnz) = j;
        nnz++;
      }
    }
    for ( ; j <= B_last_index; j++) {
      const node_handle b_d = ((j == B->index(B_curr_index))? B->down(B_curr_index++): 0);
      node_handle down = compute(A_ext_d, b_d);
      if (down) {
        C->setSparse(nnz, j, down);
        // C->d_ref(nnz) = down;
        // C->i_ref(nnz) = j;
        nnz++;
      }
    }
    MEDDLY_DCASSERT(j == max_a_b_last_index+1);

    // Last index
    //
    MEDDLY_DCASSERT(!C->isExtensible());  // default: not extensible
    if (C_is_extensible) {
      MEDDLY_DCASSERT(C_ext_d);
      C->setSparse(nnz, j, C_ext_d);
      // C->d_ref(nnz) = C_ext_d;
      // C->i_ref(nnz) = j;
      nnz++;
      C->markAsExtensible();
    }
    C->shrink(nnz);

    // reduce
    edge_value ev;
    node_handle result;
    resF->createReducedNode(C, ev, result, in);
    MEDDLY_DCASSERT(ev.isVoid());
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  return result;
}
#endif // ALLOW_EXTENSIBLE

// ******************************************************************
// *                                                                *
// *                 generic_binbylevel_mxd methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binbylevel_mxd
::generic_binbylevel_mxd(binary_list& code, forest* arg1,
  forest* arg2, forest* res)
 : binary_operation(code, 1, arg1, arg2, res)
{
  ct_entry_type* et = new ct_entry_type(code.getName(), "INN:N");
  et->setForestForSlot(1, arg1);
  et->setForestForSlot(2, arg2);
  et->setForestForSlot(4, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::generic_binbylevel_mxd::~generic_binbylevel_mxd()
{
}

void MEDDLY::generic_binbylevel_mxd
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c, bool userFlag)
{
  node_handle result = compute(
    resF->getMaxLevelIndex(), a.getNode(), b.getNode()
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
  ct_entry_key* Key = findResult(resultLevel, a, b, result);
  if (0==Key) return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  unpacked_node* C = unpacked_node::newFull(resF, resultLevel, resultSize);

  bool canSaveResult = true;

  // Initialize readers
  unpacked_node* A = unpacked_node::New(arg1F);
  unpacked_node* B = unpacked_node::New(arg2F);

  if (aLevel == resultLevel) {
    arg1F->unpackNode(A, a, FULL_ONLY);
  } else if (resultLevel>0 || arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, resultLevel, a, FULL_ONLY);
  } else {
    A->initIdentity(arg1F, resultLevel, in, a, FULL_ONLY);
    canSaveResult = false;
  }

  MEDDLY_DCASSERT(!A->isExtensible());

  if (bLevel == resultLevel) {
    arg2F->unpackNode(B, b, FULL_ONLY);
  } else if (resultLevel>0 || arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, resultLevel, b, FULL_ONLY);
  } else {
    B->initIdentity(arg2F, resultLevel, in, b, FULL_ONLY);
    canSaveResult = false;
  }

  MEDDLY_DCASSERT(!B->isExtensible());

  // Do computation
  int nextLevel = MXD_levels::downLevel(resultLevel);
  unsigned nnz = 0;
  for (unsigned j=0; j<resultSize; j++) {
    node_handle d = compute_r(j, nextLevel, A->down(j), B->down(j));
    C->setFull(j, d);
    // C->d_ref(j) = d;
    if (d) nnz++;
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce
  edge_value ev;
  resF->createReducedNode(C, ev, result, in);
  MEDDLY_DCASSERT(ev.isVoid());

  // save result in compute table, when we can
  if (resultLevel<0 && 1==nnz) canSaveResult = false;
  if (canSaveResult)  saveResult(Key, resultLevel, a, b, result);
  else                CT0->recycle(Key);

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, level %d, %d, %d) = %d\n", getName(), in, resultLevel, a, b, result);
  fflush(stdout);
#endif

  return result;
}


// ******************************************************************
// *                                                                *
// *                   generic_binary_ev  methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_ev::generic_binary_ev(binary_list& code,
  forest* arg1, forest* arg2, forest* res)
  : binary_operation(code, 1, arg1, arg2, res)
{
}

MEDDLY::generic_binary_ev::~generic_binary_ev()
{
}


// ******************************************************************
// *                                                                *
// *                 generic_binary_evplus  methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evplus::generic_binary_evplus(binary_list& code,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_ev(code, arg1, arg2, res)
{
  ct_entry_type* et = new ct_entry_type(code.getName(), "LNLN:LN");
  et->setForestForSlot(1, arg1);
  et->setForestForSlot(3, arg2);
  et->setForestForSlot(6, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::generic_binary_evplus::~generic_binary_evplus()
{
}

void MEDDLY::generic_binary_evplus
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c, bool userFlag)
{
  node_handle result;
  long ev = Inf<long>(), aev = Inf<long>(), bev = Inf<long>();
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(ev, result);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

void MEDDLY::generic_binary_evplus
::compute(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

  ct_entry_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const int resultLevel = MAX(aLevel, bLevel);
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  // Initialize result
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0L, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0L, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY)
  ;

  MEDDLY_DCASSERT(!A->isExtensible() && !B->isExtensible());


  // do computation
  for (unsigned i=0; i<resultSize; i++) {
    long ev = 0;
    node_handle ed = 0;
    compute(aev + A->edge_long(i), A->down(i),
            bev + B->edge_long(i), B->down(i),
            ev, ed);
    MEDDLY_DCASSERT(ed != 0 || ev == 0);
    nb->setFull(i, edge_value(ev), ed);
    // nb->d_ref(i) = ed;
    // nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // Reduce
  edge_value ev;
  resF->createReducedNode(nb, ev, c);
  cev = ev.getLong();

  // Add to CT
  saveResult(Key, aev, a, bev, b, cev, c);
}

// ******************************************************************
// *                                                                *
// *               generic_binary_evplus_mxd methods                *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evplus_mxd::generic_binary_evplus_mxd(binary_list& code,
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_ev(code, arg1, arg2, res)
{
  if (!arg1->isForRelations() || !arg2->isForRelations() || !res->isForRelations()) {
    throw error::TYPE_MISMATCH;
  }
  ct_entry_type* et = new ct_entry_type(code.getName(), "LNLN:LN");
  et->setForestForSlot(1, arg1);
  et->setForestForSlot(3, arg2);
  et->setForestForSlot(6, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::generic_binary_evplus_mxd::~generic_binary_evplus_mxd()
{
}

void MEDDLY::generic_binary_evplus_mxd
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c, bool userFlag)
{
  node_handle result;
  long ev = Inf<long>(), aev = Inf<long>(), bev = Inf<long>();
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(ev, result);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

void MEDDLY::generic_binary_evplus_mxd
::compute(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c)) {
    return;
  }

  ct_entry_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0 == Key) {
    return;
  }

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const int resultLevel = ABS(isLevelAbove(aLevel, bLevel) ? aLevel: bLevel);
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));

  // Initialize result
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = isLevelAbove(resultLevel, aLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 0L, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY);

  unpacked_node *B = isLevelAbove(resultLevel, bLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 0L, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY);

  // do computation
  for (unsigned i = 0; i < resultSize; i++) {
//    long ev = Inf<long>();
    long ev = 0;
    node_handle ed = 0;
    compute_r(i, MXD_levels::downLevel(resultLevel),
      aev + A->edge_long(i), A->down(i),
      bev + B->edge_long(i), B->down(i),
      ev, ed);
    MEDDLY_DCASSERT(ed != 0 || ev == 0);
    nb->setFull(i, edge_value(ev), ed);
    // nb->d_ref(i) = ed;
    // nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // Reduce
  edge_value ev;
  resF->createReducedNode(nb, ev, c);
  cev = ev.getLong();

  // Add to CT
  saveResult(Key, aev, a, bev, b, cev, c);
}

void MEDDLY::generic_binary_evplus_mxd
::compute_r(int in, int level, long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c)
{
  //  Compute for the primed levels.
  //
  MEDDLY_DCASSERT(level < 0);

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  const unsigned resultSize = unsigned(resF->getLevelSize(level));

  unpacked_node* C = unpacked_node::newFull(resF, level, resultSize);

  // Initialize readers
  unpacked_node *A = unpacked_node::New(arg1F);
  unpacked_node *B = unpacked_node::New(arg2F);

  if (aLevel == level) {
    arg1F->unpackNode(A, a, FULL_ONLY);
  } else if (arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, level, 0L, a, FULL_ONLY);
  } else {
    A->initIdentity(arg1F, level, in, 0L, a, FULL_ONLY);
  }

  if (bLevel == level) {
    arg2F->unpackNode(B, b, FULL_ONLY);
  } else if (arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, level, 0L, b, FULL_ONLY);
  } else {
    B->initIdentity(arg2F, level, in, 0L, b, FULL_ONLY);
  }

  // Do computation
  for (unsigned i = 0; i < resultSize; i++) {
    long ev = 0;
    node_handle e = 0;
    compute(aev + A->edge_long(i), A->down(i), bev + B->edge_long(i), B->down(i), ev, e);
    MEDDLY_DCASSERT(e != 0 || ev == 0);
    C->setFull(i, edge_value(ev), e);
    // C->setEdge(i, ev);
    // C->d_ref(i) = e;
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce
  edge_value ev;
  resF->createReducedNode(C, ev, c, in);
  cev = ev.getLong();

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, %d, %d) = %d\n", getName(), in, a, b, c);
#endif
}

// ******************************************************************
// *                                                                *
// *                 generic_binary_evtimes methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::generic_binary_evtimes
::generic_binary_evtimes(binary_list& code, forest* arg1,
  forest* arg2, forest* res)
: generic_binary_ev(code, arg1, arg2, res)
{
  ct_entry_type* et = new ct_entry_type(code.getName(), "FNFN:FN");
  et->setForestForSlot(1, arg1);
  et->setForestForSlot(3, arg2);
  et->setForestForSlot(6, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::generic_binary_evtimes::~generic_binary_evtimes()
{
}

void MEDDLY::generic_binary_evtimes
::computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge& c, bool userFlag)
{
  node_handle result;
  float ev, aev, bev;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(ev, result);
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

void MEDDLY::generic_binary_evtimes
::compute(float aev, node_handle a, float bev, node_handle b,
  float& cev, node_handle& c)
{
#ifdef TRACE_ALL_OPS
   printf("computing %s(<%f, %d>, <%f, %d>)\n", getName(), aev, a, bev, b);
   fflush(stdout);
#endif

  // Compute for the unprimed levels.
  //

  if (checkTerminals(aev, a, bev, b, cev, c))
    return;

#ifndef DISABLE_CACHE
  ct_entry_key* Key = findResult(aev, a, bev, b, cev, c);
  if (0==Key) return;
#endif

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  // Initialize result
  const int resultLevel = ABS(topLevel(aLevel, bLevel));
  MEDDLY_DCASSERT(resultLevel>0);
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, 1.0f, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY)
  ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, 1.0f, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY)
  ;

  MEDDLY_DCASSERT(!A->isExtensible() && !B->isExtensible());

  // do computation
  for (unsigned i=0; i<resultSize; i++) {
    float ev;
    node_handle ed;
    compute_k(
        i, -resultLevel,
        aev * A->edge_float(i), A->down(i),
        bev * B->edge_float(i), B->down(i),
        ev, ed);
    nb->setFull(i, edge_value(ev), ed);
    // nb->d_ref(i) = ed;
    // nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // Reduce
  edge_value ev;
  resF->createReducedNode(nb, ev, c);
  cev = ev.getFloat();

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
  const unsigned resultSize = unsigned(resF->getLevelSize(resultLevel));
  unpacked_node* nb = unpacked_node::newFull(resF, resultLevel, resultSize);

  // Initialize readers
  unpacked_node *A = unpacked_node::New(arg1F);
  unpacked_node *B = unpacked_node::New(arg2F);

  if (aLevel == resultLevel) {
    arg1F->unpackNode(A, a, FULL_ONLY);
  } else if (arg1F->isFullyReduced()) {
    A->initRedundant(arg1F, resultLevel, 1.0f, a, FULL_ONLY);
  } else {
    A->initIdentity(arg1F, resultLevel, in, 1.0f, a, FULL_ONLY);
  }

  if (bLevel == resultLevel) {
    arg2F->unpackNode(B, b, FULL_ONLY);
  } else if (arg2F->isFullyReduced()) {
    B->initRedundant(arg2F, resultLevel, 1.0f, b, FULL_ONLY);
  } else {
    B->initIdentity(arg2F, resultLevel, in, 1.0f, b, FULL_ONLY);
  }

  MEDDLY_DCASSERT(!A->isExtensible() && !B->isExtensible());

  // do computation
  for (unsigned i=0; i<resultSize; i++) {
    float ev;
    node_handle ed;
    compute(
        aev * A->edge_float(i), A->down(i),
        bev * B->edge_float(i), B->down(i),
        ev, ed);
    nb->setFull(i, edge_value(ev), ed);
    // nb->d_ref(i) = ed;
    // nb->setEdge(i, ev);
  }

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // Reduce
  edge_value ev;
  resF->createReducedNode(nb, ev, c, in);
  cev = ev.getFloat();
}
