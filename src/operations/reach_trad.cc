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
#include "../forest.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"
#include "reach_trad.h"

namespace MEDDLY {
    class reachset_frontier;
    class reachset_no_frontier;
};

// ******************************************************************
// *                                                                *
// *                    reachset_frontier  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::reachset_frontier : public binary_operation {
    public:
        reachset_frontier(binary_operation* image, binary_operation* acc,
                binary_operation* diff);


        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    private:
        unary_operation*  CopyOp;
        binary_operation* imageOp;
        binary_operation* accumulateOp;
        binary_operation* differenceOp;
};

// ************************************************************************

MEDDLY::reachset_frontier::reachset_frontier(binary_operation* image,
        binary_operation* acc, binary_operation* diff)

    : binary_operation(image->getOp1F(), image->getOp2F(), image->getResF()),
      imageOp(image), accumulateOp(acc), differenceOp(diff)
{
    checkDomains(__FILE__, __LINE__);
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
    checkLabelings(__FILE__, __LINE__,
        resF->getEdgeLabeling(),
        edge_labeling::MULTI_TERMINAL,
        resF->getEdgeLabeling()
    );

    if (arg1F->getRangeType() != resF->getRangeType()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    CopyOp = build(COPY, arg1F, resF);

    if ((!CopyOp) || (!imageOp) || (!accumulateOp) || (!differenceOp))
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    usesProgress();
}

// ************************************************************************

void MEDDLY::reachset_frontier::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
    //
    // Set reachable to initial states
    //
    CopyOp->compute(L, in, av, ap, cv, cp);
    dd_edge saveC(resF);
    saveC.set(cv, cp);

    //
    // Set frontier set to initial states
    //

    edge_value fv;
    node_handle fp;

    CopyOp->compute(L, in, av, ap, fv, fp);

    dd_edge frnt(resF);
    frnt.set(fv, fp);

    //
    // Iterate:
    //

    unsigned iters = 0;
    do {
        ++iters;
        notifyProgress(iters, ' ');

        edge_value nv;
        node_handle np;
        dd_edge next(resF);

        //
        // Compute next states
        //
        imageOp->compute(L, in, fv, fp, bv, bp, nv, np);
        next.set(nv, np);
        notifyProgress(iters, 'N');

        //
        // Subtract reachable to get new frontier
        //
        differenceOp->compute(L, in, nv, np, cv, cp, fv, fp);
        frnt.set(fv, fp);
        notifyProgress(iters, 'F');

        if (0 == fp) break; // frontier is empty (or all infinity?), stop

        //
        // Add frontier to reachable
        //
        accumulateOp->compute(L, in, cv, cp, fv, fp, cv, cp);
        saveC.set(cv, cp);
        notifyProgress(iters, ';');

    } while (true);

    //
    // saveC will unlink cp, so "make a copy"
    //
    cp = resF->linkNode(saveC.getNode());
    notifyProgress(iters, ';');
}

// ******************************************************************
// *                                                                *
// *                   reachset_no_frontier class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::reachset_no_frontier : public binary_operation {
    public:
        reachset_no_frontier(binary_operation* image, binary_operation* acc);


        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    private:
        unary_operation*  CopyOp;
        binary_operation* imageOp;
        binary_operation* accumulateOp;
};

// ************************************************************************

MEDDLY::reachset_no_frontier::reachset_no_frontier(binary_operation* image,
        binary_operation* acc)

    : binary_operation(image->getOp1F(), image->getOp2F(), image->getResF()),
      imageOp(image), accumulateOp(acc)
{
    checkDomains(__FILE__, __LINE__);
    checkRelations(__FILE__, __LINE__, SET, RELATION, SET);
    checkLabelings(__FILE__, __LINE__,
        resF->getEdgeLabeling(),
        edge_labeling::MULTI_TERMINAL,
        resF->getEdgeLabeling()
    );

    if (arg1F->getRangeType() != resF->getRangeType()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    CopyOp = build(COPY, arg1F, resF);

    if ((!CopyOp) || (!imageOp) || (!accumulateOp))
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    usesProgress();
}

// ************************************************************************

void MEDDLY::reachset_no_frontier::compute(int L, unsigned in,
        const edge_value &av, node_handle ap,
        const edge_value &bv, node_handle bp,
        edge_value &cv, node_handle &cp)
{
    //
    // Previous reachable states
    //
    edge_value  oldv;
    node_handle oldp;

    //
    // Set reachable to initial states
    //
    CopyOp->compute(L, in, av, ap, cv, cp);
    dd_edge saveC(resF);
    saveC.set(cv, cp);

    //
    // Iterate:
    //

    unsigned iters = 0;
    do {
        ++iters;
        notifyProgress(iters, ' ');

        oldv = cv;
        oldp = cp;

        edge_value nv;
        node_handle np;
        dd_edge next(resF);

        //
        // Compute next states
        //
        imageOp->compute(L, in, cv, cp, bv, bp, nv, np);
        next.set(nv, np);
        notifyProgress(iters, 'N');

        //
        // Union with reachable
        //
        accumulateOp->compute(L, in, cv, cp, nv, np, cv, cp);
        saveC.set(cv, cp);
        notifyProgress(iters, ';');

    } while ((cp != oldp) || (cv != oldv));

    //
    // saveC will unlink cp, so "make a copy"
    //
    cp = resF->linkNode(saveC.getNode());
}

// ******************************************************************
// *                                                                *
// *                  reachset_tradf_factory class                  *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <bool FWD>
    class reachset_tradf_factory : public binary_factory {
        public:
            virtual void setup();

            virtual binary_operation*
                build_new(forest* a, forest* b, forest* c);
    };
};

// ******************************************************************

template <bool FWD>
void MEDDLY::reachset_tradf_factory <FWD>::setup()
{
    if (FWD) {
        _setup(__FILE__, "REACHABLE_TRAD_FS(true)", "Build forward reachability set using a traditional breadth-first iteration, with a frontier set. The iterations stop once the frontier set is empty.");
    } else {
        _setup(__FILE__, "REACHABLE_TRAD_FS(false)", "Build backward reachability set using a traditional breadth-first iteration, with a frontier set. The iterations stop once the frontier set is empty.");
    }
}

template <bool FWD>
MEDDLY::binary_operation*
MEDDLY::reachset_tradf_factory <FWD>::build_new(forest* a, forest* b, forest* c)
{
    binary_operation* imageOp = nullptr;
    binary_operation* unionOp = nullptr;
    binary_operation* diffrOp = nullptr;

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {

        switch (c->getRangeType()) {
            case range_type::BOOLEAN:
                imageOp = FWD ? MEDDLY::build(POST_IMAGE, a, b, c)
                              : MEDDLY::build(PRE_IMAGE,  a, b, c);

                unionOp = MEDDLY::build(UNION, c, c, c);
                diffrOp = MEDDLY::build(DIFFERENCE, c, c, c);

                return new reachset_frontier(imageOp, unionOp, diffrOp);

            default:
                return nullptr;
        }
    }

    return nullptr;
}


// ******************************************************************
// *                                                                *
// *                 reachset_tradnof_factory class                 *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    template <bool FWD>
    class reachset_tradnof_factory : public binary_factory {
        public:
            virtual void setup();

            virtual binary_operation*
                build_new(forest* a, forest* b, forest* c);
    };
};

// ******************************************************************

template <bool FWD>
void MEDDLY::reachset_tradnof_factory <FWD>::setup()
{
    if (FWD) {
        _setup(__FILE__, "REACHABLE_TRAD_NOFS(true)", "Build forward reachability set using a traditional breadth-first iteration, without a frontier set. The transition relation is applied to all reachable states, to obtain the new reachable states. The iterations stop once the frontier set is empty.");
    } else {
        _setup(__FILE__, "REACHABLE_TRAD_NOFS(false)", "Build backward reachability set using a traditional breadth-first iteration, without a frontier set. The transition relation is applied to all reachable states, to obtain the new reachable states. The iterations stop once the frontier set is empty.");
    }
}

template <bool FWD>
MEDDLY::binary_operation*
MEDDLY::reachset_tradnof_factory <FWD>::build_new(forest* a, forest* b, forest* c)
{
    binary_operation* imageOp = nullptr;
    binary_operation* unionOp = nullptr;

    if (a->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {

        switch (c->getRangeType()) {
            case range_type::BOOLEAN:
                imageOp = FWD ? MEDDLY::build(POST_IMAGE, a, b, c)
                              : MEDDLY::build(PRE_IMAGE,  a, b, c);

                unionOp = MEDDLY::build(UNION, c, c, c);

                return new reachset_no_frontier(imageOp, unionOp);

            default:
                return nullptr;
        }
    }

    return nullptr;
}



// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_factory& MEDDLY::REACHABLE_TRAD_FS(bool fwd)
{
    static reachset_tradf_factory<true>  forwd;
    static reachset_tradf_factory<false> bckwd;

    if (fwd) return forwd;
    else     return bckwd;
}

MEDDLY::binary_factory& MEDDLY::REACHABLE_TRAD_NOFS(bool fwd)
{
    static reachset_tradnof_factory<true>  forwd;
    static reachset_tradnof_factory<false> bckwd;

    if (fwd) return forwd;
    else     return bckwd;
}

