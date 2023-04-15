
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
#include "old_meddly.h"
#include "old_meddly.hh"
#include "old_meddly_expert.h"
#include "old_meddly_expert.hh"
#include "reach_bfs.h"

// #define DEBUG_BFS
// #define VERBOSE_BFS

namespace MEDDLY {
  class common_bfs;
  // class common_bfs_mt;
  class forwd_bfs_mt;
  class bckwd_bfs_mt;

  // class common_bfs_evplus;
  class forwd_bfs_evplus;
  class bckwd_bfs_evplus;

  class forwd_bfs_opname;
  class bckwd_bfs_opname;
};

// ******************************************************************
// *                                                                *
// *                        common_bfs class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs : public binary_operation {
  public:
    common_bfs(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

  protected:
    inline void setUnionOp(binary_operation* uop)
    {
      MEDDLY_DCASSERT(uop);
      MEDDLY_DCASSERT(0==unionOp);
      unionOp = uop;
    }

    inline void setImageOp(binary_operation* iop)
    {
      MEDDLY_DCASSERT(iop);
      MEDDLY_DCASSERT(0==imageOp);
      imageOp = iop;
    }

  private:
    binary_operation* unionOp;
    binary_operation* imageOp;

};


MEDDLY::common_bfs::common_bfs(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res)
: binary_operation(oc, 0, a1, a2, res)
{
  unionOp = 0;
  imageOp = 0;
}

void MEDDLY::common_bfs::computeDDEdge(const dd_edge &init, const dd_edge &R, dd_edge &reachableStates, bool userFlag)
{
  MEDDLY_DCASSERT(unionOp);
  MEDDLY_DCASSERT(imageOp);

  reachableStates = init;
  dd_edge prevReachable(resF);
  dd_edge front(resF);
#ifdef DEBUG_BFS
  FILE_output debug(stderr);
  debug << "Relation: ";
  R.show(debug, 2);
  debug << "Initial states: ";
  init.show(debug, 2);
  long iters = 0;
#endif
#ifdef VERBOSE_BFS
  long iters = 0;
  FILE_OUTPUT verbose(stderr);
#endif
  while (prevReachable != reachableStates) {
#ifdef VERBOSE_BFS
    iters++;
    verbose << "Iteration " << iters << ":\n";
#endif
    prevReachable = reachableStates;
    imageOp->computeDDEdge(reachableStates, R, front, userFlag);
#ifdef VERBOSE_BFS
    verbose << "\timage done ";
    front.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    iters++;
    debug << "Iteration " << iters << "\npseudo-frontier: ";
    front.show(debug, 2);
#endif
    unionOp->computeDDEdge(reachableStates, front, reachableStates, userFlag);
#ifdef VERBOSE_BFS
    verbose << "\tunion done ";
    reachableStates.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    debug << "Reachable so far: ";
    reachable.show(debug, 2);
#endif
  }

}

// ******************************************************************
// *                                                                *
// *                       forwd_bfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_mt : public common_bfs {
  public:
    forwd_bfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);
};

MEDDLY::forwd_bfs_mt::forwd_bfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == range_type::BOOLEAN) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    setUnionOp( getOperation(MAXIMUM, res, res, res) );
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                       bckwd_bfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_mt : public common_bfs {
  public:
    bckwd_bfs_mt(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

};

MEDDLY::bckwd_bfs_mt::bckwd_bfs_mt(const binary_opname* oc, expert_forest* a1,
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == range_type::BOOLEAN) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    setUnionOp( getOperation(MAXIMUM, res, res, res) );
  }
  setImageOp( getOperation(PRE_IMAGE, a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                    common_bfs_evplus  class                    *
// *                                                                *
// ******************************************************************

/*

class MEDDLY::common_bfs_evplus : public binary_operation {
  public:
  common_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

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
        imageOp->computeTemp(resEv, resEvmdd, mxd, front_ev, front);
#ifdef VERBOSE_BFS
        fprintf(stderr, "\timage done <%ld, %d>\n", front_ev, front);
#endif
#ifdef DEBUG_BFS
        iters++;
        fprintf(stderr, "Iteration %d\npseudo-frontier: <%ld, %d>\n", iters, front_ev, front);
        arg1F->showNodeGraph(stderr, front);
#endif
        unionMinOp->computeTemp(resEv, resEvmdd, front_ev, front, resEv, resEvmdd);
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
: binary_operation(oc, 0, a1, a2, res)
{
  unionMinOp = 0;
  imageOp = 0;
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

*/

// ******************************************************************
// *                                                                *
// *                     forwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_evplus : public common_bfs {
  public:
  forwd_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

//     virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::forwd_bfs_evplus::forwd_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  // expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == range_type::INTEGER) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( getOperation(POST_IMAGE, a1, a2, res) );
}

/*
void MEDDLY::forwd_bfs_evplus::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  if (resF->getRangeType() == range_type::INTEGER) {
    unionMinOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(POST_IMAGE, arg1F, arg2F, resF);

  iterate(ev, evmdd, mxd, resEv, resEvmdd);
}
*/

// ******************************************************************
// *                                                                *
// *                     bckwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_evplus : public common_bfs {
  public:
    bckwd_bfs_evplus(const binary_opname* opcode, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

    // virtual void compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd);
};

MEDDLY::bckwd_bfs_evplus::bckwd_bfs_evplus(const binary_opname* oc, expert_forest* a1,
  // expert_forest* a2, expert_forest* res) : common_bfs_evplus(oc, a1, a2, res)
  expert_forest* a2, expert_forest* res) : common_bfs(oc, a1, a2, res)
{
  if (res->getRangeType() == range_type::INTEGER) {
    setUnionOp( getOperation(UNION, res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( getOperation(PRE_IMAGE, a1, a2, res) );
}

/*
void MEDDLY::bckwd_bfs_evplus::compute(long ev, node_handle evmdd, node_handle mxd, long& resEv, node_handle& resEvmdd)
{
  if (resF->getRangeType() == range_type::INTEGER) {
    unionMinOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(PRE_IMAGE, arg1F, arg2F, resF);

  iterate(ev, evmdd, mxd, resEv, resEvmdd);
}
*/


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
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (a1->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    return new forwd_bfs_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == edge_labeling::EVPLUS) {
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
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (a1 != r)
    throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);

  if (
    a1->isForRelations()    ||
    !a2->isForRelations()   ||
    r->isForRelations()     ||
    (a1->getRangeType() != r->getRangeType()) ||
    (a1->getEdgeLabeling() != r->getEdgeLabeling()) ||
    (a2->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
  )
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (a1->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
    return new bckwd_bfs_mt(this, a1, a2, r);
  }
  else if (a1->getEdgeLabeling() == edge_labeling::EVPLUS) {
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

