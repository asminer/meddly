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
#include "reach_bfs.h"

#include "../forest.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"

// #define DEBUG_BFS
// #define VERBOSE_BFS

namespace MEDDLY {
    class common_bfs;
    class forwd_bfs_mt;
    class bckwd_bfs_mt;

    class forwd_bfs_evplus;
    class bckwd_bfs_evplus;

    binary_list FWD_BFS_cache;
    binary_list REV_BFS_cache;
};

// ******************************************************************
// *                                                                *
// *                        common_bfs class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::common_bfs : public binary_operation {
  public:
    common_bfs(binary_list& opcode, forest* arg1, forest* arg2, forest* res);

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


MEDDLY::common_bfs::common_bfs(binary_list& oc, forest* a1,
  forest* a2, forest* res)
: binary_operation(oc, 0, a1, a2, res)
{
    unionOp = 0;
    imageOp = 0;
    checkDomains(__FILE__, __LINE__);
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
    if (a1 != res) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }
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
  R.showGraph(debug);
  debug << "Initial states: ";
  init.showGraph(debug);
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
    front.showGraph(debug);
#endif
    unionOp->computeDDEdge(reachableStates, front, reachableStates, userFlag);
#ifdef VERBOSE_BFS
    verbose << "\tunion done ";
    reachableStates.show(verbose, 0);
    verbose << "\n";
#endif
#ifdef DEBUG_BFS
    debug << "Reachable so far: ";
    reachable.showGraph(debug);
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
    forwd_bfs_mt(forest* arg1, forest* arg2, forest* res);
};

MEDDLY::forwd_bfs_mt::forwd_bfs_mt(forest* a1, forest* a2, forest* res)
    : common_bfs(FWD_BFS_cache, a1, a2, res)
{
  if (res->getRangeType() == range_type::BOOLEAN) {
    setUnionOp( UNION(res, res, res) );
  } else {
    setUnionOp( MAXIMUM(res, res, res) );
  }
  setImageOp( POST_IMAGE(a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                       bckwd_bfs_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_mt : public common_bfs {
  public:
    bckwd_bfs_mt(forest* arg1, forest* arg2, forest* res);

};

MEDDLY::bckwd_bfs_mt::bckwd_bfs_mt(forest* a1, forest* a2, forest* res)
    : common_bfs(REV_BFS_cache, a1, a2, res)
{
  if (res->getRangeType() == range_type::BOOLEAN) {
    setUnionOp( UNION(res, res, res) );
  } else {
    setUnionOp( MAXIMUM(res, res, res) );
  }
  setImageOp( PRE_IMAGE(a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                     forwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::forwd_bfs_evplus : public common_bfs {
  public:
  forwd_bfs_evplus(forest* arg1, forest* arg2, forest* res);
};

MEDDLY::forwd_bfs_evplus::forwd_bfs_evplus(forest* a1, forest* a2, forest* res)
    : common_bfs(FWD_BFS_cache, a1, a2, res)
{
  if (res->getRangeType() == range_type::INTEGER) {
    setUnionOp( UNION(res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( POST_IMAGE(a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                     bckwd_bfs_evplus class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::bckwd_bfs_evplus : public common_bfs {
  public:
    bckwd_bfs_evplus(forest* arg1, forest* arg2, forest* res);

};

MEDDLY::bckwd_bfs_evplus::bckwd_bfs_evplus(forest* a1, forest* a2, forest* res)
    : common_bfs(REV_BFS_cache, a1, a2, res)
{
  if (res->getRangeType() == range_type::INTEGER) {
    setUnionOp( UNION(res, res, res) );
  } else {
    throw error(error::INVALID_OPERATION);
  }
  setImageOp( PRE_IMAGE(a1, a2, res) );
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation* MEDDLY::REACHABLE_STATES_BFS(forest* a,
        forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  FWD_BFS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if  (
            (a->getRangeType() != c->getRangeType()) ||
            (a->getEdgeLabeling() != c->getEdgeLabeling()) ||
            (b->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return FWD_BFS_cache.add( new forwd_bfs_mt(a, b, c) );
    }
    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {
        return FWD_BFS_cache.add(new forwd_bfs_evplus(a, b, c) );
    }
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::REACHABLE_STATES_BFS_init()
{
    FWD_BFS_cache.reset("ReachableBFS");
}

void MEDDLY::REACHABLE_STATES_BFS_done()
{
    MEDDLY_DCASSERT(FWD_BFS_cache.isEmpty());
}

// ******************************************************************

MEDDLY::binary_operation* MEDDLY::REVERSE_REACHABLE_BFS(forest* a,
        forest* b, forest* c)
{
    if (!a || !b || !c) return nullptr;
    binary_operation* bop =  REV_BFS_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if  (
            (a->getRangeType() != c->getRangeType()) ||
            (a->getEdgeLabeling() != c->getEdgeLabeling()) ||
            (b->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        return REV_BFS_cache.add( new bckwd_bfs_mt(a, b, c) );
    }
    if (a->getEdgeLabeling() == edge_labeling::EVPLUS) {
        return REV_BFS_cache.add(new bckwd_bfs_evplus(a, b, c) );
    }
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

void MEDDLY::REVERSE_REACHABLE_BFS_init()
{
    REV_BFS_cache.reset("ReverseReachableBFS");
}

void MEDDLY::REVERSE_REACHABLE_BFS_done()
{
    MEDDLY_DCASSERT(REV_BFS_cache.isEmpty());
}

