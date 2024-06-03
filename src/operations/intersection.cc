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
#include "intersection.h"

#include "apply_base.h" // remove this when we can

#include "../oper_binary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class inter_mdd;
    class inter_mxd;
    class inter_max_evplus;

    binary_list INTER_cache;
};

// #define TRACE

// ******************************************************************
// *                                                                *
// *                        inter_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_mdd : public binary_operation {
    public:
        inter_mdd(forest* arg1, forest* arg2, forest* res);
        virtual ~inter_mdd();

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

MEDDLY::inter_mdd::inter_mdd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(INTER_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("inter_mdd");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::inter_mdd::~inter_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::inter_mdd::compute(const edge_value &av, node_handle ap,
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
MEDDLY::inter_mdd::_compute(node_handle A, node_handle B, int L)
{
    //
    // Check terminals
    //
    if (A==0 || B==0) {
        return 0;
    }

    if (arg1F->isTerminalNode(A)) {
        //
        // A is terminal 1 (0 case taken care of already)
        //

        if (arg2F->isTerminalNode(B)) {
            terminal tt(true);
            return resF->makeRedundantsTo(tt.getHandle(), L);
        }

        //
        // Return B if we can
        //
        if (arg2F == resF) {
            return resF->makeRedundantsTo(resF->linkNode(B), L);
        }
    }

    if (arg2F->isTerminalNode(B)) {
        //
        // B is terminal 1; return A if we can
        //
        MEDDLY_DCASSERT(!arg2F->isTerminalNode(A));
        if (arg1F == resF) {
            return resF->makeRedundantsTo(resF->linkNode(A), L);
        }
    }

    if (A == B) {
        if ((arg1F == arg2F) && (arg1F == resF)) {
            return resF->makeRedundantsTo(resF->linkNode(A), L);
        }
    }

    //
    // Reorder A and B if they commute (same forest)
    //

    if (arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
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
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel < Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, SPARSE_ONLY)
        :   arg1F->newUnpacked(A, SPARSE_ONLY);

    unpacked_node* Bu = (Blevel <Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, SPARSE_ONLY)
        :   arg2F->newUnpacked(B, SPARSE_ONLY);

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel,
            MIN(Au->getSize(), Bu->getSize()));

    //
    // Build result node
    //
    unsigned za=0, zb=0, zc=0;
    while ( (za < Au->getSize()) && (zb < Bu->getSize()) ) {
        if (Au->index(za) < Bu->index(zb)) {
            ++za;
            continue;
        }
        if (Au->index(za) > Bu->index(zb)) {
            ++zb;
            continue;
        }
        MEDDLY_DCASSERT(Au->index(za) == Bu->index(zb));

        const node_handle cd = _compute(Au->down(za), Bu->down(zb),
                Cu->getLevel()-1);

        if (cd) {
            Cu->setSparse(zc++, Au->index(za), cd);
        }
        za++;
        zb++;
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

// ******************************************************************
// *                                                                *
// *                        inter_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_mxd : public binary_operation {
    public:
        inter_mxd(forest* arg1, forest* arg2, forest* res);
        virtual ~inter_mxd();

        virtual void compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);
        node_handle _compute_primed(int in, node_handle A, node_handle B, int L);

        inline node_handle makeChainTo(node_handle p, int L)
        {
            return identity_chains
                    ? resF->makeIdentitiesTo(p, L)
                    : resF->makeRedundantsTo(p, L);
        }

        inline unpacked_node* patternNode(forest* F, int L, int in,
                node_handle A, node_storage_flags fs)
        {
            MEDDLY_DCASSERT(L<0);
            if (F->isFullyReduced()) {
                return unpacked_node::newRedundant(F, L, A, fs);
            }
            if (F->isIdentityReduced()) {
                return unpacked_node::newIdentity(F, L, in, A, fs);
            }
            MEDDLY_DCASSERT(false);
        }

    private:
        ct_entry_type* ct;

        bool identity_chains;
};

// ******************************************************************

MEDDLY::inter_mxd::inter_mxd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(INTER_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("inter_mxd");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();

    // If either argument forest is identity reduced,
    // then create identity chains when skipping levels;
    // otherwise create redundant chains.
    identity_chains = arg1->isIdentityReduced() || arg2->isIdentityReduced();
}

MEDDLY::inter_mxd::~inter_mxd()
{
    ct->markForDestroy();
}

void MEDDLY::inter_mxd::compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
    MEDDLY_DCASSERT(L>=0);
    cv.set();
    cp = _compute(ap, bp, L);
}

MEDDLY::node_handle
MEDDLY::inter_mxd::_compute(node_handle A, node_handle B, int L)
{
    MEDDLY_DCASSERT(L>=0);

    //
    // Check terminals
    //
    if (A==0 || B==0) {
        return 0;
    }

    if (arg1F->isTerminalNode(A)) {
        //
        // A is terminal 1 (0 case taken care of already)
        //

        if (arg2F->isTerminalNode(B)) {
            terminal tt(true);
            return makeChainTo(tt.getHandle(), L);
        }

        //
        // Return B if we can
        //
        if (arg1F->isFullyReduced() && (arg2F == resF)) {
            return makeChainTo(resF->linkNode(B), L);
        }
    }

    if (arg2F->isTerminalNode(B)) {
        //
        // B is terminal 1; return A if we can
        //
        MEDDLY_DCASSERT(!arg2F->isTerminalNode(A));
        if (arg2F->isFullyReduced() && arg1F == resF) {
            return makeChainTo(resF->linkNode(A), L);
        }
    }

    if (A == B) {
        if ((arg1F == arg2F) && (arg1F == resF)) {
            return makeChainTo(resF->linkNode(A), L);
        }
    }

#ifdef TRACE
    std::cout << "inter_mxd::_compute(" << A << ", " << B << ", " << L << ")\n";
#endif

    //
    // Reorder A and B if they commute (same forest)
    //

    if (arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
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
        node_handle C = makeChainTo(resF->linkNode(res[0].getN()), L);
#ifdef TRACE
        std::cout << "\tCT hit " << res[0].getN() << "\n";
        std::cout << "\tafter chain " << C << "\n";
#endif
        return C;
    }

    //
    // Do computation
    //

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = MAX(ABS(Alevel), ABS(Blevel));

#ifdef TRACE
    std::cout << "\t" << A << " level " << Alevel << "\n";
    std::cout << "\t" << B << " level " << Blevel << "\n";
    std::cout << "\tresult level " << Clevel << " before chain\n";
#endif

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, SPARSE_ONLY)
        :   arg1F->newUnpacked(A, SPARSE_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, SPARSE_ONLY)
        :   arg2F->newUnpacked(B, SPARSE_ONLY);

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel,
            MIN(Au->getSize(), Bu->getSize()));

    //
    // Build result node
    //
    unsigned za=0, zb=0, zc=0;
    while ( (za < Au->getSize()) && (zb < Bu->getSize()) ) {
        if (Au->index(za) < Bu->index(zb)) {
            ++za;
            continue;
        }
        if (Au->index(za) > Bu->index(zb)) {
            ++zb;
            continue;
        }
        MEDDLY_DCASSERT(Au->index(za) == Bu->index(zb));

        const node_handle cd = _compute_primed(Au->index(za),
                Au->down(za), Bu->down(zb),
                forest::downLevel(Cu->getLevel())
        );

        if (cd) {
            Cu->setSparse(zc++, Au->index(za), cd);
        }
        za++;
        zb++;
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

    C = makeChainTo(C, L);
#ifdef TRACE
    std::cout << "inter_mxd::_compute(" << A << ", " << B << ", " << L << ") = " << C << "\n\t";
    ostream_output out(std::cout);
    resF->showNode(out, C, SHOW_DETAILS);
    out.put('\n');
#endif
    return C;
}

MEDDLY::node_handle
MEDDLY::inter_mxd::_compute_primed(int in, node_handle A, node_handle B,
        const int Clevel)
{
    MEDDLY_DCASSERT(Clevel<0);

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
#ifdef TRACE
    std::cout << "inter_mxd::_compute_primed(" << in << ", " << A << ", " << B << ", " << Clevel << ")\n";
    std::cout << "\t" << A << " level " << Alevel << "\n";
    std::cout << "\t" << B << " level " << Blevel << "\n";
#endif

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   patternNode(arg1F, Clevel, in, A, SPARSE_ONLY)
        :   arg1F->newUnpacked(A, SPARSE_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   patternNode(arg2F, Clevel, in, B, SPARSE_ONLY)
        :   arg2F->newUnpacked(B, SPARSE_ONLY);

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel,
            MIN(Au->getSize(), Bu->getSize()));

    //
    // Build result node
    //
    unsigned za=0, zb=0, zc=0;
    while ( (za < Au->getSize()) && (zb < Bu->getSize()) ) {
        if (Au->index(za) < Bu->index(zb)) {
            ++za;
            continue;
        }
        if (Au->index(za) > Bu->index(zb)) {
            ++zb;
            continue;
        }
        MEDDLY_DCASSERT(Au->index(za) == Bu->index(zb));

        const node_handle cd = _compute(
                Au->down(za), Bu->down(zb),
                forest::downLevel(Cu->getLevel())
        );

        if (cd) {
            Cu->setSparse(zc++, Au->index(za), cd);
        }
        za++;
        zb++;
    }
    Cu->resize(zc);

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C, in);
    MEDDLY_DCASSERT(dummy.isVoid());

#ifdef TRACE
    std::cout << "inter_mxd::_compute_primed(" << in << ", " << A << ", "
              << B << ", " << Clevel << ") = " << C << "\n\t";
    ostream_output out(std::cout);
    resF->showNode(out, C, SHOW_DETAILS);
    out.put('\n');
#endif
    return C;
}

// ******************************************************************
// *                                                                *
// *                     inter_max_evplus  class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_max_evplus : public generic_binary_evplus {
    public:
        inter_max_evplus(forest* arg1, forest* arg2, forest* res);

        virtual ct_entry_key* findResult(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle &c);
        virtual void saveResult(ct_entry_key* key, long aev, node_handle a,
            long bev, node_handle b, long cev, node_handle c);

        virtual bool checkTerminals(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle& c);
};

MEDDLY::inter_max_evplus::inter_max_evplus(forest* arg1, forest* arg2,
        forest* res) : generic_binary_evplus(INTER_cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::inter_max_evplus::findResult(long aev, node_handle a,
  long bev, node_handle b, long& cev, node_handle &c)
{
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  if (canCommute() && a > b) {
    CTsrch->writeL(0);
    CTsrch->writeN(b);
    CTsrch->writeL(aev - bev);
    CTsrch->writeN(a);
  } else {
    CTsrch->writeL(0);
    CTsrch->writeN(a);
    CTsrch->writeL(bev - aev);
    CTsrch->writeN(b);
  }
  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  cev = CTresult[0].readL();
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev += (a > b ? bev : aev);
  }
  else {
    MEDDLY_DCASSERT(cev == 0);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::inter_max_evplus::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  if (c == 0) {
    CTresult[0].writeL(0);
  }
  else {
    CTresult[0].writeL(cev - (a > b ? bev : aev));
  }
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::inter_max_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
    long& cev, node_handle& c)
{
  if (a == 0 || b == 0) {
    cev = 0;
    c = 0;
    return true;
  }
  if (arg1F->isTerminalNode(a) && bev >= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && aev >= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MAX(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}

#ifdef ALLOW_EXTENSIBLE
MEDDLY::node_handle
MEDDLY::inter_mxd::compute_ext(node_handle a, node_handle b)
{
  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));
  const int dwnLevel = resF->downLevel(resultLevel);

  MEDDLY_DCASSERT(resF->isExtensibleLevel(resultLevel));

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, SPARSE_ONLY)
    : arg1F->newUnpacked(a, SPARSE_ONLY)
    ;
  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = int(A->getSize())-1;
  for ( ; last_nz >= 0 && A->down(unsigned(last_nz)) == 0; last_nz--);
  const unsigned int A_nnzs = last_nz + 1;
  const int A_last_index = last_nz >= 0? int(A->index(unsigned(last_nz))): -1;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, SPARSE_ONLY)
    : arg2F->newUnpacked(b, SPARSE_ONLY)
    ;
  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = int(B->getSize())-1;
  for ( ; last_nz >= 0 && B->down(unsigned(last_nz)) == 0; last_nz--);
  const unsigned int B_nnzs = last_nz + 1;
  const int B_last_index = last_nz >= 0? int(B->index(unsigned(last_nz))): -1;

  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  unsigned resultSize = A->getSize() + B->getSize() + 1 + 1;
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  unsigned nnz = 0;
  unsigned A_curr_index = 0;
  unsigned B_curr_index = 0;
  for ( ; A_curr_index < A_nnzs && B_curr_index < B_nnzs; ) {
    // get a_i, a_d, b_i, b_d
    node_handle a_d, b_d;
    const unsigned a_i = A->index(A_curr_index);
    const unsigned b_i = B->index(B_curr_index);
    if (a_i <= b_i) {
      a_d = A->down(A_curr_index);
      A_curr_index++;
    } else {
      a_d = 0;
    }
    if (a_i >= b_i) {
      b_d = B->down(B_curr_index);
      B_curr_index++;
    } else {
      b_d = 0;
    }

    if (a_d == 0 || b_d == 0) continue;

    // compute inter(a_d, b_d)
    const unsigned index = (a_d? a_i: b_i);
    const node_handle down = compute_r(int(index), dwnLevel, a_d, b_d);

    // if inter is non-zero, add it to the new node
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  } // for loop
  if (B_ext_d != 0) {
    for ( ; A_curr_index < A_nnzs; A_curr_index++) {
      // do inter(a_i, b_ext_i)
      const unsigned index = A->index(A_curr_index);
      const node_handle down = compute_r(int(index), dwnLevel, A->down(A_curr_index), B_ext_d);
      if (down) {
        C->setSparse(nnz, index, down);
        // C->i_ref(nnz) = index;
        // C->d_ref(nnz) = down;
        nnz++;
      }
    }
  }
  if (A_ext_d != 0) {
    for ( ; B_curr_index < B_nnzs; B_curr_index++) {
      // do inter(a_ext_i, b_i)
      const unsigned index = B->index(B_curr_index);
      node_handle down = compute_r(int(index), dwnLevel, A_ext_d, B->down(B_curr_index));
      if (down) {
        C->setSparse(nnz, index, down);
        // C->i_ref(nnz) = index;
        // C->d_ref(nnz) = down;
        nnz++;
      }
    }
  }
  if (A_ext_d != 0 && B_ext_d != 0) {
    const unsigned index = max_a_b_last_index+1;
    MEDDLY_DCASSERT(index >= 0);
    const node_handle down = compute_r(index, dwnLevel, A_ext_d, B_ext_d);
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = unsigned(index);
      // C->d_ref(nnz) = down;
      C->markAsExtensible();
      nnz++;
    } else {
      C->markAsNotExtensible();
    }
  }
  C->shrink(nnz);

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce result
  node_handle result = resF->createReducedNode(-1, C);
  return result;
}
#endif // ALLOW_EXTENSIBLE

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::binary_operation*
MEDDLY::INTERSECTION(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }

    binary_operation* bop =  INTER_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return INTER_cache.add(new inter_mxd(a, b, c));
        } else {
            return INTER_cache.add(new inter_mdd(a, b, c));
        }
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS) {
        if (! c->isForRelations()) {
            return INTER_cache.add(new inter_max_evplus(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::INTERSECTION_init()
{
    INTER_cache.reset("Intersection");
}

void MEDDLY::INTERSECTION_done()
{
    MEDDLY_DCASSERT(INTER_cache.isEmpty());
}

