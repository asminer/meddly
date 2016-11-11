
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
#include "reach_bfs.h"

// #define DEBUG_BFS
// #define VERBOSE_BFS

namespace MEDDLY {
  class common_bfs_mt;
  class forwd_bfs_mt;
  class bckwd_bfs_mt;

  class common_bfs_evplus;
  class forwd_bfs_evplus;
  class bckwd_bfs_evplus;

  class forwd_bfs_opname;
  class bckwd_bfs_opname;
};

// ******************************************************************
// *                                                                *
// *                      common_bfs_mt  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs_mt : public binary_operation {
  public:
    common_bfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual node_handle compute(node_handle a, node_handle b) = 0;
  protected:
    binary_operation* unionOp;
    binary_operation* imageOp;

    inline node_handle iterate(node_handle init, node_handle R) {
      node_handle reachableStates = arg1F->linkNode(init);
      node_handle prevReachable = 0;
#ifdef DEBUG_BFS
      fprintf(stderr, "Relation: %d\n", R);
      arg2F->showNodeGraph(stderr, R);
      fprintf(stderr, "Initial states: %d\n", init);
      arg1F->showNodeGraph(stderr, init);
      long iters = 0;
#endif
#ifdef VERBOSE_BFS
      long iters = 0;
#endif
      while (prevReachable != reachableStates) {
#ifdef VERBOSE_BFS
        iters++;
        fprintf(stderr, "Iteration %d:\n", iters);
#endif
        resF->unlinkNode(prevReachable);
        prevReachable = reachableStates;
        node_handle front = imageOp->compute(reachableStates, R);
#ifdef VERBOSE_BFS
        fprintf(stderr, "\timage done %d\n", front);
#endif
#ifdef DEBUG_BFS
        iters++;
        fprintf(stderr, "Iteration %d\npseudo-frontier: %d\n", iters, front);
        arg1F->showNodeGraph(stderr, front);
#endif
        reachableStates = unionOp->compute(reachableStates, front);
#ifdef VERBOSE_BFS
        fprintf(stderr, "\tunion done %d\n", reachableStates);
#endif
#ifdef DEBUG_BFS
        fprintf(stderr, "Reachable so far: %d\n", reachableStates);
        arg1F->showNodeGraph(stderr, reachableStates);
#endif
        resF->unlinkNode(front);
      }
      resF->unlinkNode(prevReachable);
      return reachableStates;
    }
};

MEDDLY::common_bfs_mt::common_bfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 0, 0, a1, a2, res)
{
  unionOp = 0;
  imageOp = 0;
}

bool MEDDLY::common_bfs_mt::isStaleEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::common_bfs_mt::discardEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::common_bfs_mt::showEntry(output &strm, const node_handle* entryData) const
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::common_bfs_mt
::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  node_handle cnode = compute(a.getNode(), b.getNode());
  c.set(cnode);
}

// ******************************************************************
// *                                                                *
// *                       forwd_bfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_mt : public common_bfs_mt {
  public:
    forwd_bfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual node_handle compute(node_handle a, node_handle b);
};

MEDDLY::forwd_bfs_mt::forwd_bfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs_mt(oc, a1, a2, res)
{
}

MEDDLY::node_handle MEDDLY::forwd_bfs_mt::compute(node_handle a, node_handle b)
{
  if (resF->getRangeType() == forest::BOOLEAN) {
    unionOp = getOperation(UNION, resF, resF, resF);
  } else {
    unionOp = getOperation(MAXIMUM, resF, resF, resF);
  }
  imageOp = getOperation(POST_IMAGE, arg1F, arg2F, resF);

  return iterate(a, b);
}


// ******************************************************************
// *                                                                *
// *                       bckwd_bfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_mt : public common_bfs_mt {
  public:
    bckwd_bfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual node_handle compute(node_handle a, node_handle b);
};

MEDDLY::bckwd_bfs_mt::bckwd_bfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs_mt(oc, a1, a2, res)
{
}

MEDDLY::node_handle MEDDLY::bckwd_bfs_mt::compute(node_handle a, node_handle b)
{
  if (resF->getRangeType() == forest::BOOLEAN) {
    unionOp = getOperation(UNION, resF, resF, resF);
  } else {
    unionOp = getOperation(MAXIMUM, resF, resF, resF);
  }
  imageOp = getOperation(PRE_IMAGE, arg1F, arg2F, resF);

  return iterate(a, b);
}

// ******************************************************************
// *                                                                *
// *                    common_bfs_evplus  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs_evplus : public binary_operation {
  public:
  common_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);
    virtual void compute(long aev, node_handle a, node_handle b, long& resEv, node_handle& resEvmdd) = 0;
  protected:
    binary_operation* unionMinOp;
    binary_operation* imageOp;

    inline void iterate(long ev, node_handle init, node_handle mxd, long& resEv, node_handle& resEvmdd) {
      resEv = ev;
      resEvmdd = arg1F->linkNode(init);

      node_handle prevReachable = 0;
#ifdef DEBUG_BFS
      fprintf(stderr, "Relation: %d\n", mxd);
      arg2F->showNodeGraph(stderr, mxd);
      fprintf(stderr, "Initial states: <%ld, %d>\n", ev, init);
      arg1F->showNodeGraph(stderr, init);
      long iters = 0;
#endif
#ifdef VERBOSE_BFS
      long iters = 0;
#endif
      while (prevReachable != resEvmdd) {
#ifdef VERBOSE_BFS
        iters++;
        fprintf(stimageOpderr, "Iteration %d:\n", iters);
#endif
        resF->unlinkNode(prevReachable);
        prevReachable = resEvmdd;
        long front_ev = Inf<long>();
        node_handle front = 0;
        imageOp->compute(resEv, resEvmdd, mxd, front_ev, front);
#ifdef VERBOSE_BFS
        fprintf(stderr, "\timage done <%ld, %d>\n", front_ev, front);
#endif
#ifdef DEBUG_BFS
        iters++;
        fprintf(stderr, "Iteration %d\npseudo-frontier: <%ld, %d>\n", iters, front_ev, front);
        arg1F->showNodeGraph(stderr, front);
#endif
        unionMinOp->compute(resEv, resEvmdd, front_ev, front, resEv, resEvmdd);
#ifdef VERBOSE_BFS
        fprintf(stderr, "\tunion done <%ld, %d>\n", resEv, resEvmdd);
#endif
#ifdef DEBUG_BFS
        fprintf(stderr, "Reachable so far: <%ld, %d>\n", resEv, resEvmdd);
        arg1F->showNodeGraph(stderr, resEvmdd);
#endif
        resF->unlinkNode(front);
      }
      resF->unlinkNode(prevReachable);
    }
};

MEDDLY::common_bfs_evplus::common_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 0, 0, a1, a2, res)
{
  unionMinOp = 0;
  imageOp = 0;
}

bool MEDDLY::common_bfs_evplus::isStaleEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::common_bfs_evplus::discardEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::common_bfs_evplus::showEntry(output &strm, const node_handle* entryData) const
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::common_bfs_evplus::computeDDEdge(const dd_edge &a, const dd_edge &b, dd_edge &c)
{
  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  compute(aev, a.getNode(), b.getNode(), cev, cnode);
  c.set(cnode, cev);
}

// ******************************************************************
// *                                                                *
// *                     forwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_evplus : public common_bfs_evplus {
  public:
  forwd_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::forwd_bfs_evplus::forwd_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
{
}

void MEDDLY::forwd_bfs_evplus::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  if (resF->getRangeType() == forest::INTEGER) {
    unionMinOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(POST_IMAGE, arg1F, arg2F, resF);

  iterate(ev, evmdd, mxd, resEv, resEvmdd);
}

// ******************************************************************
// *                                                                *
// *                     bckwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_evplus : public common_bfs_evplus {
  public:
    bckwd_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::bckwd_bfs_evplus::bckwd_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
{
}

void MEDDLY::bckwd_bfs_evplus::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  if (resF->getRangeType() == forest::INTEGER) {
    unionMinOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(PRE_IMAGE, arg1F, arg2F, resF);

  iterate(ev, evmdd, mxd, resEv, resEvmdd);
}


// ******************************************************************
// *                                                                *
// *                     forwd_bfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_opname : public binary_opname {
  public:
    forwd_bfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::forwd_bfs_opname::forwd_bfs_opname()
 : binary_opname("ReachableBFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::forwd_bfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH);

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new forwd_bfs_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
    return new forwd_bfs_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}

// ******************************************************************
// *                                                                *
// *                     bckwd_bfs_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_opname : public binary_opname {
  public:
    bckwd_bfs_opname();
    virtual binary_operation* buildOperation(expert_forest* a1,
      expert_forest* a2, expert_forest* r) const;
};

MEDDLY::bckwd_bfs_opname::bckwd_bfs_opname()
 : binary_opname("ReverseReachableBFS")
{
}

MEDDLY::binary_operation* 
MEDDLY::bckwd_bfs_opname::buildOperation(expert_forest* a1, expert_forest* a2, 
  expert_forest* r) const
{
  if (0==a1 || 0==a2 || 0==r) return 0;

  if (  
    (a1->getDomain() != r->getDomain()) || 
    (a2->getDomain() != r->getDomain()) 
  )
    throw error(error::DOMAIN_MISMATCH);

  if (a1 != r)
    throw error(error::FOREST_MISMATCH);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != forest::MULTI_TERMINAL) 
  )
    throw error(error::TYPE_MISMATCH);

  if (a1->getEdgeLabeling() == forest::MULTI_TERMINAL) {
    return new bckwd_bfs_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == forest::EVPLUS) {
    return new bckwd_bfs_evplus(this, a1, a2, r);
  }
  else {
    throw error(error::TYPE_MISMATCH);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname* MEDDLY::initializeForwardBFS()
{
  return new forwd_bfs_opname;
}

MEDDLY::binary_opname* MEDDLY::initializeBackwardBFS()
{
  return new bckwd_bfs_opname;
}

