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
#include "difference.h"
#include "apply_base.h" // remove this when we can

#include "../oper_binary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class diffr_mdd;
    class diffr_mxd;

    binary_list DIFFR_cache;
};

// #define OLD_DIFF

// ******************************************************************
// *                                                                *
// *                        diffr_mdd  class                        *
// *                                                                *
// ******************************************************************

#ifdef OLD_DIFF

class MEDDLY::diffr_mdd : public generic_binary_mdd {
    public:
        diffr_mdd(forest* arg1, forest* arg2, forest* res);
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};

MEDDLY::diffr_mdd::diffr_mdd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mdd(DIFFR_cache, arg1, arg2, res)
{
    //  difference does NOT commute

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::diffr_mdd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0 || b == -1) {
    c = 0;
    return true;
  }
  if (a == -1 && b == 0) {
    c = -1;
    return true;
  }
  if (a == b) {
    if (arg1F == arg2F) {
      c = 0;
      return true;
    } else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } else {
      return false;
    }
  }
  return false;
}

#else

class MEDDLY::diffr_mdd : public binary_operation {
    public:
        diffr_mdd(forest* arg1, forest* arg2, forest* res);
        virtual ~diffr_mdd();

        virtual void compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);

    private:
        ct_entry_type* ct;
};

// ******************************************************************

MEDDLY::diffr_mdd::diffr_mdd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(DIFFR_cache, arg1, arg2, res)
{
    ct = new ct_entry_type("difference");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

MEDDLY::diffr_mdd::~diffr_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::diffr_mdd::compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    cv.set();
    cp = _compute(ap, bp, L);
}

MEDDLY::node_handle
MEDDLY::diffr_mdd::_compute(node_handle A, node_handle B, int L)
{
    //
    // Check terminals
    //
    if (A==0) {
        // 0 - B = 0
        return 0;
    }
    if (arg2F->isTerminalNode(B) && B != 0) {
        // A - 1 = 0
        return 0;
    }
    if (arg1F->isTerminalNode(A) && (B==0)) {
        // 1 - 0 = 1
        terminal tt(true);
        return resF->makeRedundantsTo(tt.getHandle(), L);
    }

    if (B==0) {
        // A - 0 = A
        // return A if we can (same forest as result)
        if (arg1F == resF) {
            return resF->makeRedundantsTo(resF->linkNode(A), L);
        }
    }

    if (A == B) {
        if (arg1F == arg2F) {
            return 0;
        }
    }

    //
    // Check compute table
    //
    ct_vector key(2);
    ct_vector res(1);
    key[0].setN(A);
    key[1].setN(B);
    if (ct->findCT(key, res)) {
        return resF->makeRedundantsTo(resF->linkNode(res[0].getN()), L);
    }

    //
    // Do computation
    //

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = MAX(Alevel, Blevel);

    //
    // Initialize unpacked nodes.
    // Sparse for A, full for B
    //
    unpacked_node* Au = (Alevel < Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, SPARSE_ONLY)
        :   arg1F->newUnpacked(A, SPARSE_ONLY);

    unpacked_node* Bu = (Blevel <Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel, Au->getSize());

    //
    // Build result node
    // Scan through (sparse) entries of A,
    // and subtract the corresponding entries in B.
    //
    unsigned zc = 0;
    for (unsigned z=0; z<Au->getSize(); z++) {
        const unsigned i = Au->index(z);
        const node_handle cd = _compute(Au->down(z), Bu->down(i),
                Cu->getLevel()-1);
        if (cd) {
            Cu->setSparse(zc++, i, cd);
        }
    }
    Cu->resize(zc);

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);

    return resF->makeRedundantsTo(C, L);
}


#endif

// ******************************************************************
// *                                                                *
// *                        diffr_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::diffr_mxd : public generic_binary_mxd {
    public:
        diffr_mxd(forest* arg1, forest* arg2, forest* res);
        virtual bool checkTerminals(node_handle a, node_handle b,
                node_handle& c);
};


MEDDLY::diffr_mxd::diffr_mxd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(DIFFR_cache, arg1, arg2, res)
{
    //  difference does NOT commute

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::diffr_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a == 0) {
    c = 0;
    return true;
  }
  if (a == -1 && b == 0) {
    c = -1;
    return true;
  }
  if (a == b) {
    if (arg1F == arg2F) {
      c = 0;
      return true;
    } else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      c = resF->linkNode(a);
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::DIFFERENCE(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }
    binary_operation* bop =  DIFFR_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return DIFFR_cache.add(new diffr_mxd(a, b, c));
        } else {
            return DIFFR_cache.add(new diffr_mdd(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::DIFFERENCE_init()
{
    DIFFR_cache.reset("Difference");
}

void MEDDLY::DIFFERENCE_done()
{
    MEDDLY_DCASSERT(DIFFR_cache.isEmpty());
}
