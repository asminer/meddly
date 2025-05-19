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
#include "mdd2index.h"

#include "../oper_unary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class mdd2index_operation;

    unary_list MDD2INDEX_cache;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                   mdd2index_operation  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd2index_operation : public unary_operation {
    public:
        mdd2index_operation(forest* arg, forest* res);
        virtual ~mdd2index_operation();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

    protected:
        void _compute(int L, node_handle A, long &cv, node_handle &cp);

    private:
        ct_entry_type* ct;

#ifdef TRACE
        ostream_output out;
#endif
};

MEDDLY::mdd2index_operation::mdd2index_operation(forest* arg, forest* res)
    : unary_operation(arg, res)
#ifdef TRACE
      , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkRanges(__FILE__, __LINE__,
        range_type::BOOLEAN, range_type::INTEGER
    );
    checkLabelings(__FILE__, __LINE__,
        edge_labeling::MULTI_TERMINAL, edge_labeling::INDEX_SET
    );

    ct = new ct_entry_type("convert2index");
    ct->setFixed(arg);
    // answer[0] : node
    // answer[1] : cardinality
    ct->setResult(res, 'L');

    ct->doneBuilding();
}

MEDDLY::mdd2index_operation::~mdd2index_operation()
{
    ct->markForDestroy();
}

void MEDDLY::mdd2index_operation::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    MEDDLY_DCASSERT(av.isVoid());
    long lcv;
    _compute(L, ap, lcv, cp);
    cv.set(0L);  // NOT lcv!
}

void MEDDLY::mdd2index_operation::_compute(int L,
        node_handle A, long &cv, node_handle &cp)
{
    //
    // Terminal cases
    //
    if (0 == A) {
        cv = 0;
        cp = OMEGA_INFINITY;
        return;
    }
    if (0 == L) {
        cv = 1;
        cp = OMEGA_NORMAL;
        return;
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
    MEDDLY_DCASSERT(Alevel <= L);
#ifdef TRACE
    out << "mdd2index::_compute(" << L << ", " << A << ")\n";
#endif
    if (Alevel == L) {
    }

    //
    // Check compute table, but only if we're at the right level
    //
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    if (Alevel == L) {
        key[0].setN(A);
        if (ct->findCT(key, res)) {
            //
            // compute table 'hit'
            //
            cp = resF->linkNode(res[0].getN());
            cv = res[1].getL();
#ifdef TRACE
            out << "  CT hit (" << cp << ", " << cv << ")\n";
            out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
            //
            // done: compute table 'hit'
            //
            return;
        }
    }

    //
    // compute table 'miss'; do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel < L)
        ? unpacked_node::newRedundant(argF, L, A, FULL_ONLY)
        : unpacked_node::newFromNode(argF, A, FULL_ONLY);

    unpacked_node* Cu = unpacked_node::newWritable(resF, L, FULL_ONLY);
    MEDDLY_DCASSERT(Au->getSize() == Cu->getSize());
#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
    out.indent_more();
    out.put('\n');
#endif

    //
    // Recurse over child edges
    //
    cv = 0;
    for (unsigned i=0; i<Cu->getSize(); i++) {
        node_handle ddn;
        long dcard = 0;
        _compute(L-1, Au->down(i), dcard, ddn);
        if (dcard) {
            Cu->setFull(i, edge_value(cv), ddn);
            cv += dcard;
        } else {
            Cu->setFull(i, edge_value(0L), ddn);
        }
    }
    Cu->setUHdata(&cv);

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "mdd2index::_compute(" << L << ", " << A << ") done\n";
    out << "A: ";
    Au->show(out, true);
    out << "\nC: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    edge_value ev;
    resF->createReducedNode(Cu, ev, cp);
    MEDDLY_DCASSERT(0 == ev.getLong());

#ifdef TRACE
    out << "reduced to (" << cp << ", " << cv << "):";
    resF->showNode(out, cp, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Add to CT
    //
    if (Alevel == L) {
        res[0].setN(cp);
        res[1].set(cv);
        ct->addCT(key, res);
    }

    //
    // Cleanup
    //
    unpacked_node::Recycle(Au);
    // Cu is recycled when we reduce it :)
}

/*
MEDDLY::mdd2index_operation::mdd2index_operation(forest* arg, forest* res)
    : unary_operation(MDD2INDEX_cache, 1, arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkRanges(__FILE__, __LINE__,
        range_type::BOOLEAN, range_type::INTEGER
    );
    checkLabelings(__FILE__, __LINE__,
        edge_labeling::MULTI_TERMINAL, edge_labeling::INDEX_SET
    );


    // answer[0] : node
    // answer[1] : cardinality
    ct_entry_type* et = new ct_entry_type(MDD2INDEX_cache.getName(), "N:NL");
    et->setForestForSlot(0, arg);
    et->setForestForSlot(2, res);
    registerEntryType(0, et);
    buildCTs();
}

void
MEDDLY::mdd2index_operation
::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  MEDDLY_DCASSERT(arg.getForest() == argF);
  MEDDLY_DCASSERT(res.getForest() == resF);
  MEDDLY_DCASSERT(argF->handleForValue(true) < 0);
  MEDDLY_DCASSERT(argF->handleForValue(false) == 0);
  node_handle down;
  long card = 0;
  int nVars = argF->getMaxLevelIndex();
  compute_r(nVars, arg.getNode(), down, card);
  res.set(0L, down);
}

void
MEDDLY::mdd2index_operation
::compute_r(int k, node_handle a, node_handle &bdn, long &bcard)
{
  // Deal with terminals
  if (0 == a) {
    bdn = 0;
    bcard = 0;
    return;
  }
  if (0 == k) {
    terminal t(true);
    bdn = t.getHandle();
    bcard = 1;
    return;
  }

  int aLevel = argF->getNodeLevel(a);
  MEDDLY_DCASSERT(aLevel <= k);

  // Check compute table
  ct_entry_key* CTsrch = 0;
  if (aLevel == k) {
    CTsrch = CT0->useEntryKey(etype[0], 0);
    MEDDLY_DCASSERT(CTsrch);
    CTsrch->writeN(a);
    CT0->find(CTsrch, CTresult[0]);
    if (CTresult[0]) {
      bdn = resF->linkNode(CTresult[0].readN());
      bcard = CTresult[0].readL();
      CT0->recycle(CTsrch);
#ifdef TRACE_ALL_OPS
      printf("CT hit  mdd2index::compute(%d, %d) = %ld\n", k, a, bcard);
#endif
      return;
    }
  }

#ifdef TRACE_ALL_OPS
  printf("calling mdd2index::compute(%d, %d)\n", k, a);
#endif

  // Initialize node builder
  const unsigned size = unsigned(resF->getLevelSize(k));
  unpacked_node* nb = unpacked_node::newWritable(resF, k, size, FULL_ONLY);

  // Initialize node reader
  unpacked_node *A = unpacked_node::New(argF, FULL_ONLY);
  if (aLevel < k) {
    A->initRedundant(k, a);
  } else {
    A->initFromNode(a);
  }

  // recurse
  bcard = 0;
  for (unsigned i=0; i<size; i++) {
    node_handle ddn;
    long dcard = 0;
    compute_r(k-1, A->down(i), ddn, dcard);
    // nb->d_ref(i) = ddn;
    if (ddn) {
      nb->setFull(i, edge_value(bcard), ddn);
      // nb->setEdge(i, bcard);
      bcard += dcard;
    } else {
      MEDDLY_DCASSERT(0 == dcard);
      nb->setFull(i, edge_value(0L), ddn);
      // nb->setEdge(i, 0L);
    }
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  nb->setUHdata(&bcard);
  // memcpy(nb->UHdata(), &bcard, sizeof(bcard));

  edge_value ev;
  resF->createReducedNode(nb, ev, bdn);
  MEDDLY_DCASSERT(0 == ev.getLong());

  // Add to compute table
  if (CTsrch) {
    CTresult[0].reset();
    CTresult[0].writeN(bdn);
    CTresult[0].writeL(bcard);
    CT0->addEntry(CTsrch, CTresult[0]);
  }

#ifdef TRACE_ALL_OPS
  printf("finish  mdd2index::compute(%d, %d) = %ld\n", k, a, bcard);
#endif
}

*/

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_operation*
MEDDLY::CONVERT_TO_INDEX_SET(forest* arg, forest* res)
{
    if (!arg || !res) {
        return nullptr;
    }

    unary_operation* uop =  MDD2INDEX_cache.find(arg, res);
    if (uop) {
        return uop;
    }

    return MDD2INDEX_cache.add(new mdd2index_operation(arg, res));
}

void MEDDLY::CONVERT_TO_INDEX_SET_init()
{
    MDD2INDEX_cache.reset("ConvertToIndexSet");
}

void MEDDLY::CONVERT_TO_INDEX_SET_done()
{
    MEDDLY_DCASSERT(MDD2INDEX_cache.isEmpty());
}

