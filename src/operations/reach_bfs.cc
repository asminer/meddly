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
#include "../forest_levels.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"

// #define VERBOSE_BFS

namespace MEDDLY {
    class reach_bfs;

    binary_list FWD_BFS_cache;
    binary_list REV_BFS_cache;
};

// ******************************************************************
// *                                                                *
// *                        reach_bfs  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::reach_bfs : public binary_operation {
    public:
        reach_bfs(binary_operation* ImageOp, binary_operation* UnionOp);
        virtual ~reach_bfs();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    protected:
        binary_operation* ImageOp;
        binary_operation* UnionOp;
        unary_operation*  CopyOp;
};

// ******************************************************************

MEDDLY::reach_bfs::reach_bfs(binary_operation* Img,
        binary_operation* Un)
    : binary_operation(Img->getOp1F(), Img->getOp2F(), Img->getResF())
{
    ImageOp = Img;
    UnionOp = Un;

    CopyOp = build(COPY, arg1F, resF);

    MEDDLY_DCASSERT(ImageOp);
    MEDDLY_DCASSERT(UnionOp);
    MEDDLY_DCASSERT(CopyOp);
}

MEDDLY::reach_bfs::~reach_bfs()
{
}

void MEDDLY::reach_bfs::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
    //
    // Previous iter, and frontier set.
    //
    dd_edge prevReachable(resF);

    CopyOp->compute(L, in, av, ap, cv, cp);

#ifdef VERBOSE_BFS
    long iters = 0;
    std::cerr << "Traditional reachability\n";
#endif
    do {
#ifdef VERBOSE_BFS
        std::cerr << "    Iteration " << ++iters << "\n";
#endif
        prevReachable.set(cv, cp);
        edge_value fv;
        node_handle fp;
        ImageOp->compute(L, in, cv, cp, bv, bp, fv, fp);
        dd_edge front(resF);
        front.set(fv, fp);
#ifdef VERBOSE_BFS
        std::cerr << "        image\n";
#endif
        UnionOp->compute(L, in, cv, cp, fv, fp, cv, cp);
#ifdef VERBOSE_BFS
        std::cerr << "        union\n";
#endif
    } while ((cp != prevReachable.getNode()) || !(cv == prevReachable.getEdgeValue()));

#ifdef VERBOSE_BFS
    std::cerr << "Traditional reachability took " << iters << " iterations\n";
#endif
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

    binary_operation *img = build(POST_IMAGE, a, b, c);
    binary_operation *acc = nullptr;

    // if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (a->getRangeType() == range_type::BOOLEAN) {
            acc = build(UNION, c, c, c);
        } else {
            acc = build(MINIMUM, c, c, c);
        }

        return FWD_BFS_cache.add( new reach_bfs(img, acc) );
    // }

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

    binary_operation *img = build(PRE_IMAGE, a, b, c);
    binary_operation *acc = nullptr;

    // if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (a->getRangeType() == range_type::BOOLEAN) {
            acc = build(UNION, c, c, c);
        } else {
            acc = build(MINIMUM, c, c, c);
        }

        return REV_BFS_cache.add( new reach_bfs(img, acc) );
    // }

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

