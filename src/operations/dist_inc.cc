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
#include "dist_inc.h"

#include "../ct_vector.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../forest_levels.h"

namespace MEDDLY {
    class dist_inc_mt;
    class DIST_INC_factory;
};

// #define TRACE

// #define TRACE_IC

#ifdef TRACE
#include "../operators.h"
#endif
#ifdef TRACE_IC
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                       dist_inc_mt  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::dist_inc_mt : public unary_operation {
    public:
        dist_inc_mt(forest* arg, forest* res);
        virtual ~dist_inc_mt();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

    protected:
        /// Implement compute(), recursively
        void _compute(int L, unsigned in, node_handle A, node_handle &C);

    private:
        ct_entry_type* ct;

#ifdef TRACE
        ostream_output out;
        unsigned top_count;
#endif
};

// ******************************************************************

MEDDLY::dist_inc_mt::dist_inc_mt(forest* arg, forest* res)
    : unary_operation(arg, res)
#ifdef TRACE
      , out(std::cout), top_count(0)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    checkAllRanges(__FILE__, __LINE__, range_type::INTEGER);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("dist_inc");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::dist_inc_mt::~dist_inc_mt()
{
    ct->markForDestroy();
}

void MEDDLY::dist_inc_mt::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    ++top_count;
    out.indentation(0);
    out << "Starting top-level dist_inc_mt::compute #" << top_count << "\n";
#endif
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    _compute(L, in, ap, cp);
#ifdef TRACE
    out.indentation(0);
    out << "Finishing top-level dist_inc_mt::compute #" << top_count << "\n";
#endif
}

void MEDDLY::dist_inc_mt::_compute(int L, unsigned in,
        node_handle A, node_handle &cp)
{
#ifdef TRACE
    out << "dist_inc_mt::_compute(" << L << ", " << A << ")\n";
#endif
    //
    // Terminal cases.
    //
    if (argF->isTerminalNode(A)) {
        long ta;
        argF->getValueFromHandle(A, ta);
        const long tc = (ta < 0) ? ta : (ta+1);
        cp = resF->handleForValue(tc);

        //
        // Add nodes in case result forest has a different
        // reduction rule than the input forest.
        //
        if (argF->isIdentityReduced()) {
            cp = resF->makeIdentitiesTo(cp, 0, L, in);
        } else {
            cp = resF->makeRedundantsTo(cp, 0, L);
        }
#ifdef TRACE
        out << "terminal case, value " << ta << " becomes " << tc << "\n";
        out << "result node is " << cp << "\n";
        out << "dist_inc_mt::_compute(" << L << ", " << A << ") returning\n";
#endif
        return;
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);

    //
    // Check compute table
    //
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        //
        // compute table 'hit'
        //
        cp = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "  CT hit " << cp << "\n";
        out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
        //
        // done: compute table 'hit'
        //
    } else {
        //
        // compute table 'miss'; do computation
        //

        //
        // Initialize unpacked nodes
        //
        unpacked_node* Au = unpacked_node::newFromNode(argF, A, FULL_ONLY);
        unpacked_node* Cu = unpacked_node::newWritable(resF, Alevel, FULL_ONLY);
        MEDDLY_DCASSERT(Au->getSize() == Cu->getSize());
#ifdef TRACE
        out << "A: ";
        Au->show(out, true);
        out.indent_more();
        out.put('\n');
#endif
        const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

        //
        // Recurse over child edges
        //
        for (unsigned i=0; i<Cu->getSize(); i++) {
            node_handle d;
            _compute(Cnextlevel, i, Au->down(i), d);
            Cu->setFull(i, d);
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "dist_inc_mt::_compute(" << L << ", " << A << ") done\n";
        out << "A: ";
        Au->show(out, true);
        out << "\nC: ";
        Cu->show(out, true);
        out << "\n";
#endif

        //
        // Reduce
        //
        edge_value cv;
        resF->createReducedNode(Cu, cv, cp);
        MEDDLY_DCASSERT(cv.isVoid());
#ifdef TRACE
        out << "reduced to " << cp << ": ";
        resF->showNode(out, cp, SHOW_DETAILS);
        out << "\n";
#endif

        //
        // Add to CT
        //
        res[0].setN(cp);
        ct->addCT(key, res);

        //
        // Cleanup
        //
        unpacked_node::Recycle(Au);
        // Cu is recycled when we reduce it :)

        //
        // done: compute table 'miss'
        //
    }

    //
    // Adjust result for singletons, added identities/redundants.
    // Do this for both CT hits and misses.
    //
    const int Clevel = resF->getNodeLevel(cp);
    if (Clevel == L) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        cp = resF->redirectSingleton(in, cp);
#ifdef TRACE
        out << "after singleton redir: " << cp << "\n";
        out << "dist_inc_mt::_compute(" << L << ", " << A << ") returning\n";
#endif
        return;
    }

    //
    // Add nodes in case result forest has a different
    // reduction rule than the input forest.
    //
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, 0, L, in);
    } else {
        cp = resF->makeRedundantsTo(cp, 0, L);
    }

#ifdef TRACE
    out << "dist_inc_mt::_compute(" << L << ", " << A << ") returning\n";
#endif
}

// ******************************************************************
// *                                                                *
// *                     DIST_INC_factory class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::DIST_INC_factory : public unary_factory {
    public:
        virtual void setup();
        virtual unary_operation* build_new(forest* arg, forest* res);
};

// ******************************************************************

void MEDDLY::DIST_INC_factory::setup()
{
    _setup(__FILE__, "DIST_INC", "Distance increment. Builds a new function g(x) from input function f(x) such that g(x) = { f(x)+1, if f(x) >= 0; g(x), otherwise. Used for adding one to a distance function, encoded as a multi-terminal MDD/MXD. The input and output forests must be over the same domain, and must have integer range.");
}

MEDDLY::unary_operation*
MEDDLY::DIST_INC_factory::build_new(forest* arg, forest* res)
{
    return new dist_inc_mt(arg, res);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_factory& MEDDLY::DIST_INC()
{
    static DIST_INC_factory F;
    return F;
}

