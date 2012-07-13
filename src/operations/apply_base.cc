
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

bool MEDDLY::generic_binary_mdd::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[0]) ||
         arg2F->isStale(data[1]) ||
         resF->isStale(data[2]);
}

void MEDDLY::generic_binary_mdd::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[0]);
  arg2F->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void
MEDDLY::generic_binary_mdd ::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::generic_binary_mdd::compute(const dd_edge &a, const dd_edge &b, 
  dd_edge &c)
{
  int cnode = compute(a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
#ifdef TRACE_ALL_OPS
  printf("completed %s(%d, %d) = %d\n", 
    getName(), a.getNode(), b.getNode(), cnode);
#endif
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}


int MEDDLY::generic_binary_mdd::compute(int a, int b)
{
  int result = 0;
  if (checkTerminals(a, b, result))
    return result;
  if (findResult(a, b, result))
    return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = MAX(aLevel, bLevel);
  int resultSize = resF->getLevelSize(resultLevel);
  expert_forest::nodeBuilder& 
    nb = resF->useNodeBuilder(resultLevel, resultSize);

  // Initialize readers
  expert_forest::nodeReader* A = (aLevel < resultLevel) 
    ? arg1F->initRedundantReader(resultLevel, a, true)
    : arg1F->initNodeReader(a, true);

  expert_forest::nodeReader* B = (bLevel < resultLevel) 
    ? arg2F->initRedundantReader(resultLevel, b, true)
    : arg2F->initNodeReader(b, true);

  // do computation
  for (int i=0; i<resultSize; i++) {
    nb.d(i) = compute(A->d(i), B->d(i));
  }

  // cleanup
  arg2F->recycle(B);
  arg1F->recycle(A);

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

bool MEDDLY::generic_binary_mxd::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[2]) ||
         resF->isStale(data[3]);
}

void MEDDLY::generic_binary_mxd::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[2]);
  resF->uncacheNode(data[3]);
}

void
MEDDLY::generic_binary_mxd ::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(in: %d, %d, %d): %d]", getName(), 
    data[0], data[1], data[2], data[3]);
}

void MEDDLY::generic_binary_mxd::compute(const dd_edge &a, const dd_edge &b, 
  dd_edge &c)
{
  int cnode = compute(a.getNode(), b.getNode());
  c.set(cnode, 0, resF->getNodeLevel(cnode));
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

int MEDDLY::generic_binary_mxd::compute(int a, int b) 
{
  return compute(-1, a, b);
}

int MEDDLY::generic_binary_mxd::compute(int in, int a, int b)
{
  int result = 0;
  if (checkTerminals(a, b, result))
    return result;

  if (findResult(-1, a, b, result))
    return result;

  if (in>=0) if (findResult(in, a, b, result))
    return result;

  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = topLevel(aLevel, bLevel);
  // trick: make sure we "start" at an unprimed level
  if (in<0) resultLevel = ABS(resultLevel);

  int resultSize = resF->getLevelSize(resultLevel);
  expert_forest::nodeBuilder& 
    nb = resF->useNodeBuilder(resultLevel, resultSize);

  bool needsIndex = false;
 
  // Initialize readers
  expert_forest::nodeReader* A;
  if (aLevel == resultLevel) {
    A = arg1F->initNodeReader(a, true);
  } else if (resultLevel>0 || arg1F->isFullyReduced()) {
    A = arg1F->initRedundantReader(resultLevel, a, true);
  } else {
    A = arg1F->initIdentityReader(resultLevel, in, a, true);
    needsIndex = true;
  }

  expert_forest::nodeReader* B;
  if (bLevel == resultLevel) {
    B = arg2F->initNodeReader(b, true);
  } else if (resultLevel>0 || arg2F->isFullyReduced()) {
    B = arg2F->initRedundantReader(resultLevel, b, true);
  } else {
    B = arg2F->initIdentityReader(resultLevel, in, b, true);
    needsIndex = true;
  }

  int nnz = 0;
  for (int j=0; j<resultSize; j++) {
    int d = compute(j, A->d(j), B->d(j));
    nb.d(j) = d;
    if (d) nnz++;
  }

  // cleanup
  arg2F->recycle(B);
  arg1F->recycle(A);

  // reduce and save result
  result = resF->createReducedNode(in, nb);

  if (resultLevel < 0 && 1==nnz) needsIndex = true;

  if (needsIndex) {
    saveResult(in, a, b, result);
  } else {
    saveResult(-1, a, b, result);
  }

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

bool MEDDLY::generic_binbylevel_mxd::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[2]) ||
         resF->isStale(data[3]);
}

void MEDDLY::generic_binbylevel_mxd::discardEntry(const int* data)
{
  arg1F->uncacheNode(data[1]);
  arg2F->uncacheNode(data[2]);
  resF->uncacheNode(data[3]);
}

void
MEDDLY::generic_binbylevel_mxd::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d, %d, %d): %d]", getName(), 
    data[0], data[1], data[2], data[3]
  );
}

void MEDDLY::generic_binbylevel_mxd
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  int result = compute(
    resF->getDomain()->getNumVariables(), a.getNode(), b.getNode()
  );
  c.set(result, 0, resF->getNodeLevel(result));
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

int MEDDLY::generic_binbylevel_mxd::compute(int level, int a, int b) 
{
  return compute(-1, level, a, b);
}

int MEDDLY::generic_binbylevel_mxd
::compute(int in, int resultLevel, int a, int b)
{
  int result = 0;
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
  expert_forest::nodeBuilder& 
    nb = resF->useNodeBuilder(resultLevel, resultSize);

  bool canSaveResult = true;

  // Initialize readers
  expert_forest::nodeReader* A;
  if (aLevel == resultLevel) {
    A = arg1F->initNodeReader(a, true);
  } else if (resultLevel>0 || arg1F->isFullyReduced()) {
    A = arg1F->initRedundantReader(resultLevel, a, true);
  } else {
    A = arg1F->initIdentityReader(resultLevel, in, a, true);
    canSaveResult = false;
  }

  expert_forest::nodeReader* B;
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
  arg2F->recycle(B);
  arg1F->recycle(A);

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

bool MEDDLY::generic_binary_ev::isStaleEntry(const int* data)
{
  return arg1F->isStale(data[1]) ||
         arg2F->isStale(data[3]) ||
         resF->isStale(data[5]);
}

void MEDDLY::generic_binary_ev::discardEntry(const int* data)
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

void MEDDLY::generic_binary_evplus::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(<%d:%d>, <%d:%d>): <%d:%d>]",
      getName(),
      data[0], data[1], data[2], data[3], data[4], data[5]
  );
}

void MEDDLY::generic_binary_evplus
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  int result, ev, aev, bev;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(result, ev, resF->getNodeLevel(result));
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}


void MEDDLY::generic_binary_evplus
::compute(int aev, int a, int bev, int b, int& cev, int& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;
  if (findResult(aev, a, bev, b, cev, c))
    return;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = resF->getLevelSize(resultLevel);

  // Three vectors: operands a and b, and result c

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Bev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int aZero, aZeroev;
    compute(aev, a, INF, 0, aZeroev, aZero);

    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        resF->linkNode(aZero);
        *iterCev = aZeroev;
      } else {
        compute(aev, a, *iterBev + bev, *iterB, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(aZero);
    c = resF->createTempNode(resultLevel, C, Cev);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    std::vector<int> A(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Aev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int zeroB, zeroBev;
    compute(INF, 0, bev, b, zeroBev, zeroB);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterAev = Aev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        resF->linkNode(zeroB);
        *iterCev = zeroBev;
      } else {
        compute(*iterAev + aev, *iterA, bev, b, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(zeroB);
    c = resF->createTempNode(resultLevel, C, Cev);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    std::vector<int> A(resultSize, 0);
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<int> Aev(resultSize, INF);
    std::vector<int> Bev(resultSize, INF);
    std::vector<int> Cev(resultSize, INF);

    int zeroZero, zeroZeroEv;
    compute(INF, 0, INF, 0, zeroZeroEv, zeroZero);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<int>::iterator iterAev = Aev.begin();
    std::vector<int>::iterator iterBev = Bev.begin();
    std::vector<int>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end();
        iterA++, iterB++, iterC++, iterAev++, iterBev++, iterCev++)
    {
      if (*iterA == 0) {
        if (*iterB == 0) {
          *iterC = zeroZero;
          resF->linkNode(zeroZero);
          *iterCev = zeroZeroEv;
        } else {
          compute(INF, 0, *iterBev + bev, *iterB, *iterCev, *iterC);
        }
      } else if (*iterB == 0) {
        compute(*iterAev + aev, *iterA, INF, 0, *iterCev, *iterC);
      } else {
        compute(*iterAev + aev, *iterA, *iterBev + bev, *iterB, 
          *iterCev, *iterC);
      }
    }

    resF->unlinkNode(zeroZero);
    c = resF->createTempNode(resultLevel, C, Cev);
  }

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = resF->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  resF->normalizeAndReduceNode(c, cev);
#endif

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
::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(<%f:%d>, <%f:%d>): <%f:%d>]",
      getName(),
      toFloat(data[0]), data[1], 
      toFloat(data[2]), data[3], 
      toFloat(data[4]), data[5]
  );
}

void MEDDLY::generic_binary_evtimes
::compute(const dd_edge& a, const dd_edge& b, dd_edge& c)
{
  int result; 
  float ev, aev, bev;
  a.getEdgeValue(aev);
  b.getEdgeValue(bev);
  compute(aev, a.getNode(), bev, b.getNode(), ev, result);
  c.set(result, ev, resF->getNodeLevel(result));
#ifdef DEVELOPMENT_CODE
  resF->validateIncounts(true);
#endif
}

void MEDDLY::generic_binary_evtimes
::compute(float aev, int a, float bev, int b, float& cev, int& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c))
    return;
  if (findResult(aev, a, bev, b, cev, c))
    return;

  // 0. initialize result
  // 1. if a is at a lower level than b, expand b
  //    else if b is at a lower level than a, expand a
  //    else expand both

  // 0. initialize result
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);

  int resultLevel = aLevel > bLevel? aLevel: bLevel;
  int resultSize = resF->getLevelSize(resultLevel);

  // Three vectors: operands a and b, and result c

  if (aLevel < resultLevel) {
    // expand b
    // result[i] = a op b[i]
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<float> Bev(resultSize, NAN);
    std::vector<float> Cev(resultSize, NAN);

    int aZero; 
    float aZeroev;
    compute(aev, a, NAN, 0, aZeroev, aZero);

    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterBev = Bev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterB != B.end(); iterB++, iterC++, iterBev++, iterCev++)
    {
      if (*iterB == 0) {
        *iterC = aZero;
        resF->linkNode(aZero);
        *iterCev = aZeroev;
      } else {
        compute(aev, a, *iterBev + bev, *iterB, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(aZero);
    c = resF->createTempNode(resultLevel, C, Cev);
  }
  else if (bLevel < resultLevel) {
    // expand a
    // result[i] = a[i] op b
    std::vector<int> A(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<float> Aev(resultSize, NAN);
    std::vector<float> Cev(resultSize, NAN);

    int zeroB; 
    float zeroBev;
    compute(NAN, 0, bev, b, zeroBev, zeroB);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterAev = Aev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end(); iterA++, iterC++, iterAev++, iterCev++)
    {
      if (*iterA == 0) {
        *iterC = zeroB;
        resF->linkNode(zeroB);
        *iterCev = zeroBev;
      } else {
        compute(*iterAev + aev, *iterA, bev, b, *iterCev, *iterC);
      }
    }
    resF->unlinkNode(zeroB);
    c = resF->createTempNode(resultLevel, C, Cev);
  }
  else {
    // expand both a and b
    // result[i] = a[i] op b[i]
    std::vector<int> A(resultSize, 0);
    std::vector<int> B(resultSize, 0);
    std::vector<int> C(resultSize, 0);
    std::vector<float> Aev(resultSize, NAN);
    std::vector<float> Bev(resultSize, NAN);
    std::vector<float> Cev(resultSize, NAN);

    int zeroZero; 
    float zeroZeroEv;
    compute(NAN, 0, NAN, 0, zeroZeroEv, zeroZero);

    arg1F->getDownPtrsAndEdgeValues(a, A, Aev);
    arg2F->getDownPtrsAndEdgeValues(b, B, Bev);
    std::vector<int>::iterator iterA = A.begin();
    std::vector<int>::iterator iterB = B.begin();
    std::vector<int>::iterator iterC = C.begin();
    std::vector<float>::iterator iterAev = Aev.begin();
    std::vector<float>::iterator iterBev = Bev.begin();
    std::vector<float>::iterator iterCev = Cev.begin();
    for ( ; iterA != A.end();
        iterA++, iterB++, iterC++, iterAev++, iterBev++, iterCev++)
    {
      if (*iterA == 0) {
        if (*iterB == 0) {
          *iterC = zeroZero;
          resF->linkNode(zeroZero);
          *iterCev = zeroZeroEv;
        } else {
          compute(NAN, 0, *iterBev + bev, *iterB, *iterCev, *iterC);
        }
      } else if (*iterB == 0) {
        compute(*iterAev + aev, *iterA, NAN, 0, *iterCev, *iterC);
      } else {
        compute(*iterAev + aev, *iterA, *iterBev + bev, *iterB, 
          *iterCev, *iterC);
      }
    }

    resF->unlinkNode(zeroZero);
    c = resF->createTempNode(resultLevel, C, Cev);
  }

  // save result in compute cache and return it

#if 0
  printf("reduce(%d): ", result);
  result = resF->reduceNode(result);
  printf("%d  [", result);
  for (unsigned i = 0; i < C.size(); i++ )
  {
    printf("%d ", C[i]);
  }
  printf("]\n");
#else
  resF->normalizeAndReduceNode(c, cev);
#endif

  saveResult(aev, a, bev, b, cev, c);
}



