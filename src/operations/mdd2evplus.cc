
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
#include "mdd2evplus.h"
#include "../compute_table.h"

namespace MEDDLY {
  class mdd2evplus_operation;
  class mdd2evplus_opname;
};

// ******************************************************************
// *                                                                *
// *                   mdd2evplus_operation class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd2evplus_operation : public unary_operation {
  public:
    mdd2evplus_operation(const unary_opname* oc, expert_forest* arg, 
      expert_forest* res);

    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;

    virtual void compute(const dd_edge &arg, dd_edge &res);

    void compute(int ht, int a, int &bdn, int &bcard);
};

MEDDLY::mdd2evplus_operation::mdd2evplus_operation(const unary_opname* oc, 
  expert_forest* arg, expert_forest* res) : unary_operation(oc, arg, res)
{
  key_length = 1;
  ans_length = 2; // pointer, cardinality
}

bool 
MEDDLY::mdd2evplus_operation
::isEntryStale(const int* entryData)
{
  return 
    argF->isStale(entryData[0]) ||
    resF->isStale(entryData[1]);
}

void 
MEDDLY::mdd2evplus_operation
::discardEntry(const int* entryData)
{
  argF->uncacheNode(entryData[0]);
  resF->uncacheNode(entryData[1]);
}

void 
MEDDLY::mdd2evplus_operation
::showEntry(FILE* strm, const int *entryData) const
{
  fprintf(strm, "[%s %d %d (card %d)]", getName(), entryData[0], 
    entryData[1], entryData[2]);
}

void
MEDDLY::mdd2evplus_operation
::compute(const dd_edge &arg, dd_edge &res)
{
  DCASSERT(arg.getForest() == argF);
  DCASSERT(res.getForest() == resF);
  DCASSERT(resF->getTerminalNode(false) == 0);
  DCASSERT(argF->getTerminalNode(true) == -1);
  DCASSERT(argF->getTerminalNode(false) == 0);
  int down, card;
  int nVars = argF->getDomain()->getNumVariables();
  compute(nVars, arg.getNode(), down, card);
  res.set(down, 0, resF->getNodeLevel(down));
}

void
MEDDLY::mdd2evplus_operation
::compute(int height, int a, int &bdn, int &bcard)
{
  // Deal with terminals
  if (a == 0) {
    bdn = 0;
    bcard = 0;
    return;
  }
  if (height == 0) {
    DCASSERT(argF->getTerminalNode(true) == a);
    bdn = resF->getTerminalNode(true);
    bcard = 1;
    return;
  }

  int aHeight = argF->getNodeHeight(a);
  DCASSERT(aHeight <= height);

  // Check compute table
  if (aHeight == height) {
    const int* cacheEntry = CT->find(this, &a);
    if (cacheEntry) {
      bdn = cacheEntry[1];
      bcard = cacheEntry[2];
      resF->linkNode(bdn);
      return;
    }
  }

  // Create node at appropriate height 
  const int resultSize = resF->getDomain()->getVariableBound(height);
  bdn = resF->createTempNode(height, resultSize, true);
  bcard = 0;

  if (aHeight < height) {
    int ddn, dcard;
    compute(height-1, a, ddn, dcard);
    if (ddn) {
      for (int i=0; i<resultSize; i++) {
        resF->setDownPtrWoUnlink(bdn, i, ddn);
        resF->setEdgeValue(bdn, i, bcard);
        bcard += dcard;
      } // for i
      resF->unlinkNode(ddn);
    } // if ddn
  } // aHeight < height
  else {
    if (argF->isFullNode(a)) {
      int aSize = argF->getFullNodeSize(a);
      for (int i=0; i<aSize; i++) {
        int ddn, dcard;
        compute(height-1, argF->getFullNodeDownPtr(a, i), ddn, dcard);
        if (0==ddn) continue;
        resF->setDownPtrWoUnlink(bdn, i, ddn);
        resF->setEdgeValue(bdn, i, bcard);
        resF->unlinkNode(ddn);
        bcard += dcard;
      } // for i
    } // a is full
    else {
      int aNnz = argF->getSparseNodeSize(a);
      for (int i=0; i<aNnz; i++) {
        int ddn, dcard;
        compute(height-1, argF->getSparseNodeDownPtr(a, i), ddn, dcard);
        if (0==ddn) continue;
        resF->setDownPtrWoUnlink(bdn, i, ddn);
        resF->setEdgeValue(bdn, i, bcard);
        resF->unlinkNode(ddn);
        bcard += dcard;
      } // for i
    } // a is sparse
  } // aHeight == height

  // reduce, save in compute table
  int dummy;
  resF->normalizeAndReduceNode(bdn, dummy);
  DCASSERT((bdn==0 && dummy==INF) || (bdn!= 0 && dummy==0));
  static int cacheEntry[3];
  cacheEntry[0] = a;
  cacheEntry[1] = bdn;
  cacheEntry[2] = bcard;
  argF->cacheNode(a);
  resF->cacheNode(bdn);
  CT->add(this, cacheEntry);
  resF->setIndexSetCardinality(bdn, bcard);
}

// ******************************************************************
// *                                                                *
// *                    mdd2evplus_opname  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd2evplus_opname : public unary_opname {
  public:
    mdd2evplus_opname();
    virtual unary_operation* 
      buildOperation(const forest* ar, const forest* res) const;
};

MEDDLY::mdd2evplus_opname::mdd2evplus_opname()
 : unary_opname("ConvertToIndexSet")
{
}

MEDDLY::unary_operation*
MEDDLY::mdd2evplus_opname
::buildOperation(const forest* arg, const forest* res) const
{
  if (0==arg || 0==res) return 0;
  if (arg->isForRelations() || 
      arg->getRangeType() != forest::BOOLEAN ||
      arg->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      res->isForRelations() ||
      res->getRangeType() != forest::INTEGER ||
      res->getEdgeLabeling() != forest::EVPLUS
  ) throw error(error::TYPE_MISMATCH);

  return new mdd2evplus_operation(
    this, (expert_forest*) arg, (expert_forest*) res
  );
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeMDD2EVPLUS(const settings &s)
{
  return new mdd2evplus_opname;
}

// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************

#if 0

// Old implementation, yanked and stuck here 

/// Convert MTMDDs to EVMDDs
class mtmdd_to_evmdd : public old_operation {
  public:
    static mtmdd_to_evmdd* getInstance();

    virtual const char* getName() const { return "Convert MtMdd to EvMdd"; }

    virtual int getKeyLength() const { return 1; }
    virtual int getAnsLength() const { return 2; }
    virtual int getCacheEntryLength() const { return 3; }

    virtual int getKeyLengthInBytes() const { return 4; }
    virtual int getAnsLengthInBytes() const { return 8; }
    virtual int getCacheEntryLengthInBytes() const { return 12; }

    virtual void typeCheck(const op_info* owner);
    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;
    // Calls compute(op_info*, dd_edge, dd_edge)
    virtual void compute(op_info* owner, dd_edge** operands);
    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual void compute(op_info* owner, const dd_edge& a,
        dd_edge& b);
    // Returns an error
    virtual void compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

  protected:
    mtmdd_to_evmdd();
    mtmdd_to_evmdd(const mtmdd_to_evmdd& copy);
    mtmdd_to_evmdd& operator=(const mtmdd_to_evmdd& copy);
    ~mtmdd_to_evmdd();

    virtual void compute(op_info* owner, int a, int& b, int &ev);
    virtual bool checkTerminals(op_info* op, int a, int& b, int &ev);
    virtual bool findResult(op_info* owner, int a, int& b, int &ev);
    virtual void saveResult(op_info* owner, int a, int b, int ev);

    // for later: once evmdds can use real edge-values
    virtual bool checkTerminals(op_info* op, int a, int& b, float &ev);
    virtual bool findResult(op_info* owner, int a, int& b, float &ev);
    virtual void saveResult(op_info* owner, int a, int b, float ev);
};


// -------------------- Convert MTMDD to EVMDD ------------------


mtmdd_to_evmdd* mtmdd_to_evmdd::getInstance()
{
  static mtmdd_to_evmdd instance;
  return &instance;
}


mtmdd_to_evmdd::mtmdd_to_evmdd()
{ }


mtmdd_to_evmdd::~mtmdd_to_evmdd() {}


void
mtmdd_to_evmdd::typeCheck(const op_info* owner)
{
  if (owner == 0)
    throw error(error::UNKNOWN_OPERATION);
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    throw error(error::TYPE_MISMATCH);
  if (!owner->areAllForests())
    throw error(error::TYPE_MISMATCH);
  if (owner->nParams != 2)
    throw error(error::WRONG_NUMBER);
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    throw error(error::TYPE_MISMATCH);
  if (!getExpertForest(owner, 0)->isMtMdd())
    throw error(error::TYPE_MISMATCH);
  if (owner->p[1].isMT())
    // !MULTI_TERMINAL is same as EV+ or EV*
    throw error(error::TYPE_MISMATCH);
}


bool mtmdd_to_evmdd::isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call isStale for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  return 
    getExpertForest(owner, 0)->isStale(data[0]) ||
    getExpertForest(owner, 1)->isStale(data[1]);
}


void
mtmdd_to_evmdd::discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call uncacheNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  getExpertForest(owner, 0)->uncacheNode(data[0]);
  getExpertForest(owner, 1)->uncacheNode(data[1]);
}


void
mtmdd_to_evmdd::showEntry(const op_info* owner, FILE* strm,
  const int* data) const
{
  // data[] is of size owner.nParams + 1 (last int is edge value)
  // data[i] <--> forest[i]
  // call showNode for each forest[i] and data[i]
  DCASSERT(owner->nParams == 2);
  if (getExpertForest(owner, 1)->getRangeType() == forest::REAL) {
    fprintf(strm, "[%s %d %d %f]",
        owner->op->getName(), data[0], data[1], toFloat(data[2]));
  } else {
    fprintf(strm, "[%s %d %d %d]", owner->op->getName(),
        data[0], data[1], data[2]);
  }
}


bool
mtmdd_to_evmdd::findResult(op_info* owner, int a, int& b, int& ev)
{
  static int key[1];
  // create cache entry
  key[0] = a;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  b = cacheEntry[1];
  ev = cacheEntry[2];
  getExpertForest(owner, 1)->linkNode(b);
  return true;
}


bool
mtmdd_to_evmdd::findResult(op_info* owner, int a, int& b, float& ev)
{
  static int key[1];
  // create cache entry
  key[0] = a;

  const int* cacheEntry = owner->cc->find(owner, const_cast<const int*>(key));

  if (cacheEntry == 0) return false;
  b = cacheEntry[1];
  ev = toFloat(cacheEntry[2]);
  getExpertForest(owner, 1)->linkNode(b);
  return true;
}


void
mtmdd_to_evmdd::saveResult(op_info* owner, int a, int b, int c)
{
  static int cacheEntry[3];

#ifdef IGNORE_INCOUNT
  if (!getExpertForest(owner, 0)->isTerminalNode(a) &&
      getExpertForest(owner, 0)->getInCount(a) < IGNORE_INCOUNT)
    return;
#endif

  // create cache entry
  cacheEntry[0] = a; cacheEntry[1] = b; cacheEntry[2] = c;

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


void
mtmdd_to_evmdd::saveResult(op_info* owner, int a, int b, float c)
{
  static int cacheEntry[3];

#ifdef IGNORE_INCOUNT
  if (!getExpertForest(owner, 0)->isTerminalNode(a) &&
      getExpertForest(owner, 0)->getInCount(a) < IGNORE_INCOUNT)
    return;
#endif

  // create cache entry
  cacheEntry[0] = a; cacheEntry[1] = b; cacheEntry[2] = toInt(c);

  getExpertForest(owner, 0)->cacheNode(a);
  getExpertForest(owner, 1)->cacheNode(b);

  owner->cc->add(owner, const_cast<const int*>(cacheEntry));
}


bool mtmdd_to_evmdd::checkTerminals(op_info* op, int a, int& b, int& ev)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    ev = getExpertForest(op, 0)->getInteger(a);
    if (ev == 0) {
      // replace with non-reachable node
      b = getExpertForest(op, 1)->getTerminalNode(false);
      ev = INF;
    } else {
      b = getExpertForest(op, 1)->getTerminalNode(true);
    }
    return true;
  }

  return false;
}


bool mtmdd_to_evmdd::checkTerminals(op_info* op, int a, int& b, float& ev)
{
  if (getExpertForest(op, 0)->isTerminalNode(a)) {
    b = getExpertForest(op, 1)->getTerminalNode(true);
    ev = getExpertForest(op, 0)->getReal(a);
    return true;
  }

  return false;
}


void
mtmdd_to_evmdd::compute(op_info* owner, dd_edge** operands)
{
  if (operands == 0) 
    throw error(error::TYPE_MISMATCH);
  // compute(owner, dd_edge, dd_edge, dd_edge) checks for owner == 0
  compute(owner, *operands[0], *operands[1]);
}


void
mtmdd_to_evmdd::compute(op_info* owner, const dd_edge& a,
    const dd_edge& b, dd_edge& c)
{
  throw error(error::TYPE_MISMATCH);
}


void
mtmdd_to_evmdd::compute(op_info* owner, const dd_edge& a, dd_edge& b)
{
  if (owner == 0) 
    throw error(error::TYPE_MISMATCH);
  int result = 0;
  int ev = 0;
  compute(owner, a.getNode(), result, ev);
  b.set(result, ev, getExpertForest(owner, 1)->getNodeLevel(result));
}


// result[i] = a[i]
void mtmdd_to_evmdd::compute(op_info* owner, int a, int& b, int &bev)
{
  b = 0; bev = 0;
  if (checkTerminals(owner, a, b, bev)) return;
  if (findResult(owner, a, b, bev)) return;

  // expand nodes
  // 0. initialize result
  // 1. copy node a to node result
  // 2. do operation between contents of result and node b

  // 0. initialize result
  // a and result belong to the same domain so the levels must have the
  // same handles and the same sizes
  const int resultLevel = getExpertForest(owner, 0)->getNodeLevel(a);
  int resultSize = getExpertForest(owner, 1)->getLevelSize(resultLevel);
  int result =
    getExpertForest(owner, 1)->createTempNode(resultLevel, resultSize, false);

  // 1. copy node a to node result
  if (getExpertForest(owner, 0)->isFullNode(a)) {
    // a is a full-node
    const int aSize = getExpertForest(owner, 0)->getFullNodeSize(a);
    DCASSERT(aSize <= resultSize);
    int i = 0;
    int node = 0;
    int ev = 0;
    for ( ; i < aSize; ++i) {
      node = ev = 0;
      compute(owner, getExpertForest(owner, 0)->getFullNodeDownPtr(a, i),
          node, ev);
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
      getExpertForest(owner, 1)->unlinkNode(node);
    }
    node = ev = 0;
    compute(owner, 0, node, ev);
    for ( ; i < resultSize; ++i) {
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
    }
    getExpertForest(owner, 1)->unlinkNode(node);
  }
  else {
    // a is a sparse-node
    const int aSize =
      getExpertForest(owner, 0)->getSparseNodeSize(a);

    DCASSERT(getExpertForest(owner, 0)->getSparseNodeIndex(a, aSize - 1)
        < resultSize);

    int i = 0;        // goes from 0..resultSize
    int aIndex = 0;   // j is the sparse index; aIndex is the equiv full index
    int node = 0;
    int ev = 0;

    for (int j = 0; j < aSize; ++i, ++j) {
      aIndex = getExpertForest(owner, 0)->getSparseNodeIndex(a, j);
      for ( ; i < aIndex; ++i) {
        // a[i] is 0
        node = ev = 0;
        compute(owner, 0, node, ev);
        getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
        getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
        getExpertForest(owner, 1)->unlinkNode(node);
      }
      DCASSERT(i == aIndex && i < resultSize);
      node = ev = 0;
      compute(owner, getExpertForest(owner, 0)->getSparseNodeDownPtr(a, j),
          node, ev);
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
      getExpertForest(owner, 1)->unlinkNode(node);
    }

    node = ev = 0;
    compute(owner, 0, node, ev);
    for ( ; i < resultSize; ++i) {
      getExpertForest(owner, 1)->setDownPtrWoUnlink(result, i, node);
      getExpertForest(owner, 1)->setEdgeValue(result, i, ev);
    }
    getExpertForest(owner, 1)->unlinkNode(node);
  }

  // save result in compute cache and return it
  b = result;
  getExpertForest(owner, 1)->normalizeAndReduceNode(b, bev);
  saveResult(owner, a, b, bev);
}




class mdd_to_evplusmdd_index_set : public mtmdd_to_evmdd {
  public:
    static mdd_to_evplusmdd_index_set* getInstance();
    virtual const char* getName() const { return "Convert Mdd to Index Set"; }
    virtual void typeCheck(const op_info* owner);
    virtual void compute(op_info* owner,
        const dd_edge& a, dd_edge& b);

  protected:
    mdd_to_evplusmdd_index_set();
    mdd_to_evplusmdd_index_set(const mdd_to_evplusmdd_index_set& copy);
    mdd_to_evplusmdd_index_set& operator=(const
        mdd_to_evplusmdd_index_set& copy);
    ~mdd_to_evplusmdd_index_set();

    virtual void compute(op_info* owner, int a, int height, int& b, int& bev);
};



// -------------------- Convert MDD to EV+MDD ------------------


mdd_to_evplusmdd_index_set* mdd_to_evplusmdd_index_set::getInstance()
{
  static mdd_to_evplusmdd_index_set instance;
  return &instance;
}


mdd_to_evplusmdd_index_set::mdd_to_evplusmdd_index_set()
{ }


mdd_to_evplusmdd_index_set::~mdd_to_evplusmdd_index_set() {}


void
mdd_to_evplusmdd_index_set::typeCheck(const op_info* owner)
{
  if (owner == 0)
    throw error(error::UNKNOWN_OPERATION);
  if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
    throw error(error::TYPE_MISMATCH);
  if (!owner->areAllForests())
    throw error(error::TYPE_MISMATCH);
  if (owner->nParams != 2)
    throw error(error::WRONG_NUMBER);
  if (owner->p[0].getDomain() != owner->p[1].getDomain())
    throw error(error::TYPE_MISMATCH);
  if (!getExpertForest(owner, 0)->isMdd())
    throw error(error::TYPE_MISMATCH);
  if (!getExpertForest(owner, 1)->isEvplusMdd())
    throw error(error::TYPE_MISMATCH);
}


void
mdd_to_evplusmdd_index_set::compute(op_info* owner,
    const dd_edge& a, dd_edge& b)
{
  if (owner == 0) 
    throw error(error::TYPE_MISMATCH);
  int result = 0;
  int ev = 0;
  int nVars = owner->p[0].getDomain()->getNumVariables();
  compute(owner, a.getNode(), nVars, result, ev);
  // note that ev will be equal to the cardinality of the ev+mdd;
  // but we do not use that number here.
  b.set(result, 0, getExpertForest(owner, 1)->getNodeLevel(result));
}


void mdd_to_evplusmdd_index_set::compute(op_info* owner,
    int a, int height, int& b, int& bev)
{
  DCASSERT(getExpertForest(owner, 1)->getTerminalNode(false) == 0);
  DCASSERT(getExpertForest(owner, 0)->getTerminalNode(true) == -1);
  DCASSERT(getExpertForest(owner, 0)->getTerminalNode(false) == 0);

  b = 0; bev = INF;

  // Deal with terminals
  if (a == 0) return;
  if (height == 0) {
    DCASSERT(getExpertForest(owner, 0)->getTerminalNode(true) == a);
    b = getExpertForest(owner, 1)->getTerminalNode(true);
    bev = 1;
    return;
  }

  expert_forest* mddf = getExpertForest(owner, 0);
  expert_forest* evmddf = getExpertForest(owner, 1);
  DCASSERT(mddf->getDomain() == evmddf->getDomain());
  expert_domain* d = smart_cast<expert_domain*>(mddf->useDomain());

  int aHeight = mddf->getNodeHeight(a);
  DCASSERT(aHeight <= height);

  // Search compute table
  if (aHeight == height) {
    if (findResult(owner, a, b, bev)) return;
  }

  // Create node at appropriate height
  const int resultSize = d->getVariableBound(height);
  int result = evmddf->createTempNode(height, resultSize, true);

  // Total counts all the paths leading from result
  int total = 0;

  if (aHeight < height) {
    int node = 0;
    int ev = INF;
    compute(owner, a, height - 1, node, ev);
    // ev gives the cardinality of node

    // Set the down-pointers and edge-values
    // Down-pointers are the same for all indexes
    // Edge-value at index i is ev*i
    if (node != 0) {
      for (int i = 0; i < resultSize; i++)
      {
        evmddf->setDownPtrWoUnlink(result, i, node);
        // edge-value indicates the number of paths from all the previous
        // down-pointers
        evmddf->setEdgeValue(result, i, total);
        // update the number of paths
        total += ev;
      }
      evmddf->unlinkNode(node);
    }
  } // aHeight < height
  else {
    // create node at level corresponding to variable 'aHeight'
    // note that aHeight == height
    
    if (mddf->isFullNode(a)) {
      // a is full (or truncated-full)
      int aSize = mddf->getFullNodeSize(a);
      for (int i = 0; i < aSize; i++)
      {
        int node = 0;
        int ev = INF;
        compute(owner, mddf->getFullNodeDownPtr(a, i), height - 1, node, ev);
        if (node == 0) continue;

        evmddf->setDownPtrWoUnlink(result, i, node);
        // edge-value indicates the number of paths from all the previous
        // down-pointers
        evmddf->setEdgeValue(result, i, total);
        evmddf->unlinkNode(node);
        // update the number of paths
        total += ev;
      }
    } // a is full
    else {
      // a is sparse
      int aNnz = mddf->getSparseNodeSize(a);
      for (int i = 0; i < aNnz; i++)
      {
        int node = 0;
        int ev = INF;
        compute(owner, mddf->getSparseNodeDownPtr(a, i), height - 1, node, ev);
        if (node == 0) continue;

        int index = mddf->getSparseNodeIndex(a, i);
        evmddf->setDownPtrWoUnlink(result, index, node);
        // edge-value indicates the number of paths from all the previous
        // down-pointers
        evmddf->setEdgeValue(result, index, total);
        evmddf->unlinkNode(node);
        // update the number of paths
        total += ev;
      }
    } // a is sparse
  } // aHeight == height

  // Save result in compute cache and return it
  b = result;
  evmddf->normalizeAndReduceNode(b, bev);
  DCASSERT((b == 0 && bev == INF) || (b != 0 && bev == 0));
  bev = total;
  saveResult(owner, a, b, bev);
  // Store the cardinality along with node
  evmddf->setIndexSetCardinality(b, bev);
}

#endif
