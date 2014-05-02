
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
#include "../forests/mt.h"
#include "complement.h"

namespace MEDDLY {
  class compl_mdd;
  class compl_mxd;
  class compl_opname;
};

// #define DEBUG_MXD_COMPL

// ******************************************************************
// *                                                                *
// *                        compl_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mdd : public unary_operation {
  public:
    compl_mdd(const unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(FILE* strm, const node_handle* entryData) const;
    virtual void compute(const dd_edge& a, dd_edge& b);

  protected:
    node_handle compute(node_handle a);

    inline bool findResult(node_handle a, node_handle &b) {
      CTsrch.key(0) = a;
      const node_handle* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      b = resF->linkNode(cacheFind[1]);
      return true;
    }
    inline node_handle saveResult(node_handle a, node_handle b) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      entry.key(0) = argF->cacheNode(a);
      entry.result(0) = resF->cacheNode(b);
      CT->addEntry();
      return b;
    }
};

MEDDLY::compl_mdd
::compl_mdd(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 1, 1, arg, res)
{
  // ct entry 0: input node
  // ct entry 1: output node
}

bool MEDDLY::compl_mdd::isStaleEntry(const node_handle* data)
{
  return argF->isStale(data[0]) || resF->isStale(data[1]);
}

void MEDDLY::compl_mdd::discardEntry(const node_handle* data)
{
  argF->uncacheNode(data[0]);
  resF->uncacheNode(data[1]);
}

void MEDDLY::compl_mdd::showEntry(FILE* strm, const node_handle* data) const
{
  fprintf(strm, "[%s(%d): %d]", getName(), data[0], data[1]);
}

void MEDDLY::compl_mdd::compute(const dd_edge& a, dd_edge& b) 
{
  int result = compute(a.getNode());
  b.set(result);
}

MEDDLY::node_handle MEDDLY::compl_mdd::compute(node_handle a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    return expert_forest::bool_Tencoder::value2handle(
      !expert_forest::bool_Tencoder::handle2value(a)
    );
  }

  // Check compute table
  node_handle b;
  if (findResult(a, b)) return b;

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  const int size = resF->getLevelSize(level);
  node_builder& nb = resF->useNodeBuilder(level, size);

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, true);

  bool addRedundentNode=(resF->isQuasiReduced() && level>1);

  // recurse
  for (int i=0; i<size; i++) {
    nb.d(i) = compute(A->d(i));

    if(addRedundentNode && resF->isTerminalNode(nb.d(i)) && nb.d(i)!=resF->getTransparentNode()){
    	nb.d(i)=((mt_forest*)resF)->makeNodeAtLevel(level-1, nb.d(i));
    }
  }

  // Cleanup
  node_reader::recycle(A);

  // Reduce
  b = resF->createReducedNode(-1, nb);

  // Add to compute table
  return saveResult(a, b);
}


// ******************************************************************
// *                                                                *
// *                        compl_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mxd : public unary_operation {
  public:
    compl_mxd(const unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(FILE* strm, const node_handle* entryData) const;
    virtual void compute(const dd_edge& a, dd_edge& b);

    node_handle compute(int in, int k, node_handle a);
};

MEDDLY::compl_mxd
::compl_mxd(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 2, 1, arg, res)
{
  // ct entry 0: level
  // ct entry 1: input node
  // ct entry 2: output node
}

bool MEDDLY::compl_mxd::isStaleEntry(const node_handle* data)
{
  return argF->isStale(data[1]) || resF->isStale(data[2]);
}

void MEDDLY::compl_mxd::discardEntry(const node_handle* data)
{
  argF->uncacheNode(data[1]);
  resF->uncacheNode(data[2]);
}

void MEDDLY::compl_mxd::showEntry(FILE* strm, const node_handle* data) const
{
  fprintf(strm, "[%s(%d, %d): %d]", getName(), data[0], data[1], data[2]);
}

void MEDDLY::compl_mxd::compute(const dd_edge& a, dd_edge& b) 
{
  node_handle result = compute(-1, argF->getDomain()->getNumVariables(), a.getNode());
  b.set(result);
}

MEDDLY::node_handle MEDDLY::compl_mxd::compute(int in, int k, node_handle a)
{
  if (0==k) {
    return expert_forest::bool_Tencoder::value2handle(
      !expert_forest::bool_Tencoder::handle2value(a)
    );
  }
  if (argF->isTerminalNode(a) &&
      resF->isFullyReduced())
  {
    return expert_forest::bool_Tencoder::value2handle(
      !expert_forest::bool_Tencoder::handle2value(a)
    );
  }
  // Check compute table
  CTsrch.key(0) = k;
  CTsrch.key(1) = a;
  const node_handle* cacheFind = CT->find(CTsrch);
  if (cacheFind) {
#ifdef DEBUG_MXD_COMPL
    fprintf(stderr, "\tin CT:   compl_mxd(%d, %d) : %d\n", ht, a, cacheFind[2]);
#endif
    return resF->linkNode(cacheFind[2]);
  }

#ifdef DEBUG_MXD_COMPL
  fprintf(stderr, "\tstarting compl_mxd(%d, %d)\n", ht, a);
#endif

  // Initialize node builder
  const int size = resF->getLevelSize(k);
  node_builder& nb = resF->useNodeBuilder(k, size);

  // Initialize node reader
  const int aLevel = argF->getNodeLevel(a);
  MEDDLY_DCASSERT(!isLevelAbove(aLevel, k));
  node_reader* A;
  bool canSave = true;
  if (aLevel == k) {
    A = argF->initNodeReader(a, true);
  } else if (k>0 || argF->isFullyReduced()) {
    A = argF->initRedundantReader(k, a, true);
  } else {
    MEDDLY_DCASSERT(in>=0);
    A = argF->initIdentityReader(k, in, a, true);
    canSave = false;
  }

  // recurse
  int nextLevel = (k>0) ? -k : -k-1;
  int nnz = 0;
  bool addRedundentNode=(resF->isQuasiReduced() && (k>0 || k<-1));

  // recurse
  for (int i=0; i<size; i++) {
    int d = compute(i, nextLevel, A->d(i));
    nb.d(i) = d;
    if (d!=resF->getTransparentNode()) nnz++;

    if(addRedundentNode && resF->isTerminalNode(nb.d(i)) && nb.d(i)!=resF->getTransparentNode()){
      nb.d(i)=((mt_forest*)resF)->makeNodeAtLevel(nextLevel, nb.d(i));
    }
  }

  // cleanup
  node_reader::recycle(A);

  // reduce, save in CT
  node_handle result = resF->createReducedNode(in, nb);
  if (k<0 && 1==nnz) canSave = false;
  if (canSave) {
    compute_table::temp_entry &entry = CT->startNewEntry(this);
    entry.key(0) = k;
    entry.key(1) = argF->cacheNode(a);
    entry.result(0) = resF->cacheNode(result);
    CT->addEntry();
  }
  return result;
}


// ******************************************************************
// *                                                                *
// *                       compl_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_opname : public unary_opname {
  public:
    compl_opname();
    virtual unary_operation*
      buildOperation(expert_forest* ar, expert_forest* res) const;
};

MEDDLY::compl_opname::compl_opname()
 : unary_opname("Complement")
{
}

MEDDLY::unary_operation*
MEDDLY::compl_opname
::buildOperation(expert_forest* arg, expert_forest* res) const
{
  if (0==arg || 0==res) return 0;

  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH);

  if (arg->getRangeType() != forest::BOOLEAN ||
      arg->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      res->getRangeType() != forest::BOOLEAN ||
      res->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      arg->isForRelations() != res->isForRelations()
  ) throw error(error::TYPE_MISMATCH);

  if (arg->isForRelations())
    return new compl_mxd(this,  arg,  res);
  else
    return new compl_mdd(this,  arg,  res);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeComplement(const settings &s)
{
  return new compl_opname;
}

