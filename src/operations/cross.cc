
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
#include "cross.h"

// #define TRACE_ALL_OPS
// #define DEBUG_CROSS

namespace MEDDLY {
  class cross_bool;
  class cross_opname;
};

// ******************************************************************
// *                                                                *
// *                        cross_bool class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::cross_bool : public binary_operation {
#ifdef OLD_OP_CT
    static const int LEVEL_INDEX = 0;
    static const int OPNDA_INDEX = 1;
    static const int OPNDB_INDEX = 2;
    static const int RESLT_INDEX = 3;
#endif
  public:
    cross_bool(const binary_opname* oc, expert_forest* a1,
      expert_forest* a2, expert_forest* res);

#ifdef OLD_OP_CT
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle*);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle*);
#endif
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;
#endif // OLD_OP_CT
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    node_handle compute_pr(int in, int ht, node_handle a, node_handle b);
    node_handle compute_un(int ht, node_handle a, node_handle b);
};

#ifdef OLD_OP_CT
MEDDLY::cross_bool::cross_bool(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) 
: binary_operation(oc, 3, 1, a1, a2, res)
{
  // data[LEVEL_INDEX] : level
  // data[OPNDA_INDEX] : a
  // data[OPNDB_INDEX] : b
  // data[RESLT_INDEX] : c
}
#else
MEDDLY::cross_bool::cross_bool(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) 
: binary_operation(oc, 1, a1, a2, res)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "INN:N");
  et->setForestForSlot(1, a1);
  et->setForestForSlot(2, a2);
  et->setForestForSlot(4, res);
  registerEntryType(0, et);
  buildCTs();
}
#endif

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
bool MEDDLY::cross_bool::isStaleEntry(const node_handle* data)
{
  // data[0] is the level number
  return arg1F->isStale(data[OPNDA_INDEX]) ||
         arg2F->isStale(data[OPNDB_INDEX]) ||
         resF->isStale(data[RESLT_INDEX]);
}
#else
MEDDLY::forest::node_status
MEDDLY::cross_bool::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = arg1F->getNodeStatus(data[OPNDA_INDEX]);
  MEDDLY::forest::node_status b = arg2F->getNodeStatus(data[OPNDB_INDEX]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[RESLT_INDEX]);

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

void MEDDLY::cross_bool::discardEntry(const node_handle* data)
{
  // data[0] is the level number
  arg1F->uncacheNode(data[OPNDA_INDEX]);
  arg2F->uncacheNode(data[OPNDB_INDEX]);
  resF->uncacheNode(data[RESLT_INDEX]);
}

void
MEDDLY::cross_bool ::showEntry(output &strm, const node_handle* data) const
{
  strm  << "[" << getName() << "(level: " << long(data[1]) << ", " 
        << long(data[2]) << ", " << long(data[3]) << "): " << long(data[4])
        << "]";
}

#endif // OLD_OP_CT

void
MEDDLY::cross_bool::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  int L = arg1F->getDomain()->getNumVariables();
  node_handle cnode = compute_un(L, a.getNode(), b.getNode());
  c.set(cnode);
}

MEDDLY::node_handle MEDDLY::cross_bool::compute_un(int k, node_handle a, node_handle b)
{
#ifdef DEBUG_CROSS
  printf("calling compute_un(%d, %d, %d)\n", k, a, b);
#endif
  MEDDLY_DCASSERT(k>=0);
  if (0==a || 0==b) return 0;
  if (0==k) {
    return a; 
  }

  // check compute table
#ifdef OLD_OP_CT
  compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeI(k);
  CTsrch->writeN(a);
  CTsrch->writeN(b);
#ifdef OLD_OP_CT
  compute_table::entry_result& cacheFind = CT0->find(CTsrch);
#else
  static compute_table::entry_result cacheFind(etype[0]);
  CT0->find(CTsrch, cacheFind);
#endif
  if (cacheFind) {
    CT0->recycle(CTsrch);
    return resF->linkNode(cacheFind.readN());
  }

  // Initialize unpacked node
  unpacked_node *A = unpacked_node::useUnpackedNode();
  if (arg1F->getNodeLevel(a) < k) {
    A->initRedundant(arg1F, k, a, true);
  } else {
    A->initFromNode(arg1F, a, true);
  }

  const int resultSize = resF->getLevelSize(k);
  unpacked_node *C = unpacked_node::newFull(resF, k, resultSize);

  // recurse
  for (int i=0; i<resultSize; i++) {
    C->d_ref(i) = compute_pr(i, -k, A->d(i), b);
  }

  // reduce, save in compute table
  unpacked_node::recycle(A);
  node_handle c = resF->createReducedNode(-1, C);

  arg1F->cacheNode(a);
  arg2F->cacheNode(b);
#ifdef OLD_OP_CT
  static compute_table::entry_result result(1);
#else
  static compute_table::entry_result result(etype[0]);
#endif
  result.reset();
  result.writeN(resF->cacheNode(c));
  CT0->addEntry(CTsrch, result);

#ifdef TRACE_ALL_OPS
  printf("computed %s(%d, %d, %d) = %d\n", getName(), k, a, b, c);
#endif

  return c;
}

MEDDLY::node_handle MEDDLY::cross_bool::compute_pr(int in, int k, node_handle a, node_handle b)
{
#ifdef DEBUG_CROSS
  printf("calling compute_pr(%d, %d, %d, %d)\n", in, k, a, b);
#endif
  MEDDLY_DCASSERT(k<0);
  if (0==a || 0==b) return 0;

  // DON'T check compute table 

  // Initialize unpacked node
  unpacked_node *B = unpacked_node::useUnpackedNode();
  if (arg2F->getNodeLevel(b) < -k) {
    B->initRedundant(arg2F, -k, b, true);
  } else {
    B->initFromNode(arg2F, b, true);
  }

  const int resultSize = resF->getLevelSize(k);
  unpacked_node *C = unpacked_node::newFull(resF, k, resultSize);

  // recurse
  for (int i=0; i<resultSize; i++) {
    C->d_ref(i) = compute_un(-(k+1), a, B->d(i));
  }

  // reduce
  unpacked_node::recycle(B);
  node_handle c = resF->createReducedNode(in, C);

  // DON'T save in compute table

#ifdef TRACE_ALL_OPS
  printf("computed %s((%d), %d, %d, %d) = %d\n", getName(), in, k, a, b, c);
#endif

  return c;
}


// ******************************************************************
// *                                                                *
// *                       cross_opname class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::cross_opname : public binary_opname {
  public:
    cross_opname();
    virtual binary_operation* buildOperation(expert_forest* a1, 
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::cross_opname::cross_opname()
 : binary_opname("Cross")
{
}

MEDDLY::binary_operation* 
MEDDLY::cross_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()  ||
    (a1->getRangeType() != forest::BOOLEAN) ||
    (a1->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    a2->isForRelations()  ||
    (a2->getRangeType() != forest::BOOLEAN) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) ||
    (!r->isForRelations())  ||
    (r->getRangeType() != forest::BOOLEAN) ||
    (r->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  return new cross_bool(this, a1, a2, r);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeCross()
{
  return new cross_opname;
}

