
// $Id$

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
// #define DISABLE_CACHE 1

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

bool MEDDLY::generic_binary_mdd::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::generic_binary_mdd::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::generic_binary_mdd ::showEntry(FILE* strm, const node_handle *data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::generic_binary_mdd::compute(const dd_edge &a, const dd_edge &b, 
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
  if (findResult(a, b, result))
    return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = MAX(aLevel, bLevel);
  int resultSize = resF->getLevelSize(resultLevel);
  node_builder& nb = resF->useNodeBuilder(resultLevel, resultSize);

  // Initialize readers
  node_reader* A = (aLevel < resultLevel) 
    ? arg1F->initRedundantReader(resultLevel, a, true)
    : arg1F->initNodeReader(a, true);

  node_reader* B = (bLevel < resultLevel) 
    ? arg2F->initRedundantReader(resultLevel, b, true)
    : arg2F->initNodeReader(b, true);

  // do computation
  for (int i=0; i<resultSize; i++) {
    nb.d(i) = compute(A->d(i), B->d(i));
  }

  // cleanup
  node_reader::recycle(B);
  node_reader::recycle(A);

  // reduce and save result
  result = resF->createReducedNode(-1, nb);
  saveResult(a, b, result);

#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
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
  : binary_operation(code, 3, 1, arg1, arg2, res)
{
  // data[0] : incoming index
  // data[1] : arg1
  // data[2] : arg2
  // data[3] : result
}

MEDDLY::generic_binary_mxd::~generic_binary_mxd()
{
}

bool MEDDLY::generic_binary_mxd::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[2]) ||
         resF->isStale(data[3]);
}

void MEDDLY::generic_binary_mxd::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[2]);
  resF->uncacheNode(data[3]);
}

void
MEDDLY::generic_binary_mxd ::showEntry(FILE* strm, const node_handle *data) const
{
  fprintf(strm, "[%s(in: %d, %d, %d): %d]", getName(), 
    data[0], data[1], data[2], data[3]);
}

void MEDDLY::generic_binary_mxd::compute(const dd_edge &a, const dd_edge &b, 
  dd_edge &c)
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

  if (findResult(-1, a, b, result))
    return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));

  int resultSize = resF->getLevelSize(resultLevel);
  node_builder& nb = resF->useNodeBuilder(resultLevel, resultSize);

  // Initialize readers
  node_reader* A;
  if (aLevel == resultLevel) {
    A = arg1F->initNodeReader(a, true);
  } else {
    A = arg1F->initRedundantReader(resultLevel, a, true);
  } 

  node_reader* B;
  if (bLevel == resultLevel) {
    B = arg2F->initNodeReader(b, true);
  } else {
    B = arg2F->initRedundantReader(resultLevel, b, true);
  } 

  for (int j=0; j<resultSize; j++) {
    nb.d(j) = compute(j, -resultLevel, A->d(j), B->d(j));
  }

  // cleanup
  node_reader::recycle(B);
  node_reader::recycle(A);

  // reduce and save result
  result = resF->createReducedNode(-1, nb);
  saveResult(-1, a, b, result);

#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d) = %d\n", getName(), a, b, result);
#endif

  return result;
}

MEDDLY::node_handle 
MEDDLY::generic_binary_mxd::compute(int in, int k, node_handle a, node_handle b)
{
  //  Compute for the primed levels.
  //
  MEDDLY_DCASSERT(k<0);

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultSize = resF->getLevelSize(k);
  node_builder& nb = resF->useNodeBuilder(k, resultSize);

  // Initialize readers
  node_reader* A;
  if (aLevel == k) {
    A = arg1F->initNodeReader(a, true);
  } else if (arg1F->isFullyReduced()) {
    A = arg1F->initRedundantReader(k, a, true);
  } else {
    A = arg1F->initIdentityReader(k, in, a, true);
  }

  node_reader* B;
  if (bLevel == k) {
    B = arg2F->initNodeReader(b, true);
  } else if (arg2F->isFullyReduced()) {
    B = arg2F->initRedundantReader(k, b, true);
  } else {
    B = arg2F->initIdentityReader(k, in, b, true);
  }

  for (int j=0; j<resultSize; j++) {
    nb.d(j) = compute(A->d(j), B->d(j));
  }

  // cleanup
  node_reader::recycle(B);
  node_reader::recycle(A);

  // reduce 
  node_handle result = resF->createReducedNode(in, nb);

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, %d, %d) = %d\n", getName(), in, a, b, result);
#endif

  return result;
}

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

bool MEDDLY::generic_binbylevel_mxd::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[2]) ||
         resF->isStale(data[3]);
}

void MEDDLY::generic_binbylevel_mxd::discardEntry(const node_handle* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[2]);
  resF->uncacheNode(data[3]);
}

void
MEDDLY::generic_binbylevel_mxd
::showEntry(FILE* strm, const node_handle *data) const
{
  fprintf(strm, "[%s(%d, %d, %d): %d]", getName(), 
    data[0], data[1], data[2], data[3]
  );
}

void MEDDLY::generic_binbylevel_mxd
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
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
  return compute(-1, level, a, b);
}

MEDDLY::node_handle 
MEDDLY::generic_binbylevel_mxd
::compute(int in, int resultLevel, node_handle a, node_handle b)
{
  node_handle result = 0;
  if (0==resultLevel) {
    if (checkTerminals(a, b, result))
      return result;
  }

  if (findResult(resultLevel, a, b, result))
    return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultSize = resF->getLevelSize(resultLevel);
  node_builder& nb = resF->useNodeBuilder(resultLevel, resultSize);

  bool canSaveResult = true;

  // Initialize readers
  node_reader* A;
  if (aLevel == resultLevel) {
    A = arg1F->initNodeReader(a, true);
  } else if (resultLevel>0 || arg1F->isFullyReduced()) {
    A = arg1F->initRedundantReader(resultLevel, a, true);
  } else {
    A = arg1F->initIdentityReader(resultLevel, in, a, true);
    canSaveResult = false;
  }

  node_reader* B;
  if (bLevel == resultLevel) {
    B = arg2F->initNodeReader(b, true);
  } else if (resultLevel>0 || arg2F->isFullyReduced()) {
    B = arg2F->initRedundantReader(resultLevel, b, true);
  } else {
    B = arg2F->initIdentityReader(resultLevel, in, b, true);
    canSaveResult = false;
  }

  int nextLevel = (resultLevel > 0)? -resultLevel: -resultLevel-1;
  int nnz = 0;
  for (int j=0; j<resultSize; j++) {
    int d = compute(j, nextLevel, A->d(j), B->d(j));
    nb.d(j) = d;
    if (d) nnz++;
  }

  // cleanup
  node_reader::recycle(B);
  node_reader::recycle(A);

  // reduce and save result
  result = resF->createReducedNode(in, nb);
  if (resultLevel<0 && 1==nnz) canSaveResult = false;
  if (canSaveResult) saveResult(resultLevel, a, b, result);

#ifdef TRACE_ALL_OPS
  printf("computed %s(in %d, %d, %d) = %d\n", getName(), in, a, b, result);
#endif

  return result;
}


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

bool MEDDLY::generic_binary_ev::isStaleEntry(const node_handle* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[3]) ||
         resF->isStale(data[5]);
}

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
::showEntry(FILE* strm, const node_handle *data) const
{
  fprintf(strm, "[%s(<%d:%d>, <%d:%d>): <%d:%d>]",
      getName(),
      data[0], data[1], data[2], data[3], data[4], data[5]
  );
}

void MEDDLY::generic_binary_evplus
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
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
  if (findResult(aev, a, bev, b, cev, c))
    return;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = resF->getLevelSize(resultLevel);

  // Initialize result
  node_builder& nb = resF->useNodeBuilder(resultLevel, resultSize);

  // Initialize readers
  node_reader* A = (aLevel < resultLevel)
    ? arg1F->initRedundantReader(resultLevel, 0, a, true)
    : arg1F->initNodeReader(a, true);

  node_reader* B = (bLevel < resultLevel)
    ? arg2F->initRedundantReader(resultLevel, 0, b, true)
    : arg2F->initNodeReader(b, true);

  // do computation
  for (int i=0; i<resultSize; i++) {
    int ev;
    node_handle ed;
    compute(aev + A->ei(i), A->d(i), 
            bev + B->ei(i), B->d(i), 
            ev, ed);
    nb.d(i) = ed;
    nb.setEdge(i, ev);
  }

  // cleanup
  node_reader::recycle(B);
  node_reader::recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(-1, nb, cev, cl);
  c = cl;

  // Add to CT
  saveResult(aev, a, bev, b, cev, c);
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
::showEntry(FILE* strm, const node_handle *data) const
{
  float ev0;
  float ev2;
  float ev4;
  compute_table::readEV(data+0, ev0);
  compute_table::readEV(data+2, ev2);
  compute_table::readEV(data+4, ev4);
  fprintf(strm, "[%s(<%f:%d>, <%f:%d>): <%f:%d>]",
      getName(), ev0, data[1], ev2, data[3], ev4, data[5]
  );
}

void MEDDLY::generic_binary_evtimes
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
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
  if (findResult(aev, a, bev, b, cev, c))
    return;
#endif

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  // Initialize result
  int resultLevel = ABS(topLevel(aLevel, bLevel));
  int resultSize = resF->getLevelSize(resultLevel);
  node_builder& nb = resF->useNodeBuilder(resultLevel, resultSize);

  // Initialize readers
  node_reader* A =
    (aLevel == resultLevel)
    ? arg1F->initNodeReader(a, true)
    : arg1F->initRedundantReader(resultLevel, 1.0f, a, true);

  node_reader* B =
    (bLevel == resultLevel)
    ? arg2F->initNodeReader(b, true)
    : arg2F->initRedundantReader(resultLevel, 1.0f, b, true);

  // do computation
  for (int i=0; i<resultSize; i++) {
    float ev;
    node_handle ed;
    compute(
        i, -resultLevel,
        aev * A->ef(i), A->d(i), 
        bev * B->ef(i), B->d(i), 
        ev, ed);
    nb.d(i) = ed;
    nb.setEdge(i, ev);
    MEDDLY_DCASSERT(ed == 0 || ev > 0.0f);
  }

  // cleanup
  node_reader::recycle(B);
  node_reader::recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(-1, nb, cev, cl);
  c = cl;

#ifndef DISABLE_CACHE
  // Add to CT
  saveResult(aev, a, bev, b, cev, c);
#endif
}

void MEDDLY::generic_binary_evtimes
::compute(
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
  int resultSize = resF->getLevelSize(resultLevel);
  node_builder& nb = resF->useNodeBuilder(resultLevel, resultSize);

  // Initialize readers
  node_reader* A =
    (aLevel == resultLevel)
    ? arg1F->initNodeReader(a, true)
    : arg1F->isFullyReduced()
    ? arg1F->initRedundantReader(resultLevel, 1.0f, a, true)
    : arg1F->initIdentityReader(resultLevel, in, 1.0f, a, true);

  node_reader* B =
    (bLevel == resultLevel)
    ? arg2F->initNodeReader(b, true)
    : arg2F->isFullyReduced()
    ? arg2F->initRedundantReader(resultLevel, 1.0f, b, true)
    : arg2F->initIdentityReader(resultLevel, in, 1.0f, b, true);

  // do computation
  for (int i=0; i<resultSize; i++) {
    float ev;
    node_handle ed;
    compute(
        aev * A->ef(i), A->d(i),
        bev * B->ef(i), B->d(i), 
        ev, ed);
    nb.d(i) = ed;
    nb.setEdge(i, ev);
    MEDDLY_DCASSERT(ed == 0 || ev > 0.0f);
  }

  // cleanup
  node_reader::recycle(B);
  node_reader::recycle(A);

  // Reduce
  node_handle cl;
  resF->createReducedNode(in, nb, cev, cl);
  c = cl;
}


