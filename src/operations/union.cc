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
#include "union.h"
#include "apply_base.h" // remove this when we can

#include "../oper_binary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class union_mdd;
    class union_mxd;

    class union_min_evplus;
    class union_min_evplus_mxd;

    binary_list UNION_cache;
};

#define OLD_OPER

// ******************************************************************
// *                                                                *
// *                        union_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mdd : public binary_operation {
    public:
        union_mdd(forest* arg1, forest* arg2, forest* res);
        virtual ~union_mdd();

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

MEDDLY::union_mdd::union_mdd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(UNION_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("union");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::union_mdd::~union_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::union_mdd::compute(const edge_value &av, node_handle ap,
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
MEDDLY::union_mdd::_compute(node_handle A, node_handle B, int L)
{
    //
    // Check terminals
    //
    if (A < 0 || B < 0) {
        terminal tt(true);
        return resF->makeRedundantsTo(tt.getHandle(), 0, L);
    }

    if (A == 0) {
        if (B == 0) {
            return 0;
        }

        //
        // Return B if we can
        //
        if (arg2F == resF) {
            return resF->linkNode(B);
            // Don't need to make redundant chain b/c same forest
        }
    } // zero A

    if (B == 0) {
        //
        // Return A if we can
        //
        MEDDLY_DCASSERT(A);
        if (arg1F == resF) {
            return resF->linkNode(A);
        }
    } // zero B

    if (A == B) {
        if ((arg1F == arg2F) && (arg1F == resF)) {
            return resF->linkNode(A);
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
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = MAX(Alevel, Blevel);

    //
    // Check compute table
    //
    ct_vector key(2);
    ct_vector res(1);
    key[0].setN(A);
    key[1].setN(B);
    if (ct->findCT(key, res)) {
        return resF->makeRedundantsTo(resF->linkNode(res[0].getN()), Clevel, L);
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel < Clevel)
        ?   unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY)
        :   arg1F->newUnpacked(A, FULL_ONLY);

    unpacked_node* Bu = (Blevel <Clevel)
        ?   unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i, _compute(Au->down(i), Bu->down(i), Cu->getLevel()-1));
    }

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

    return resF->makeRedundantsTo(C, Clevel, L);
}

// ******************************************************************
// *                                                                *
// *                        union_mxd  class                        *
// *                                                                *
// ******************************************************************

#ifdef OLD_OPER

class MEDDLY::union_mxd : public generic_binary_mxd {
    public:
        union_mxd(forest* arg1, forest* arg2, forest* res);

        virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c);
#ifdef ALLOW_EXTENSIBLE
        virtual MEDDLY::node_handle compute_ext(node_handle a, node_handle b);
#endif
};

MEDDLY::union_mxd::union_mxd(forest* arg1, forest* arg2, forest* res)
  : generic_binary_mxd(UNION_cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);
}

bool MEDDLY::union_mxd::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (a < 0 && b < 0) {
    c = resF->handleForValue(true);
    return true;
  }
  if (0 == a) {
    if (0 == b) {
      c = 0;
      return true;
    }
    if (arg2F == resF) {
      c = resF->linkNode(b);
      return true;
    } else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg1F == resF) {
      c = resF->linkNode(b);
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

#ifdef ALLOW_EXTENSIBLE
MEDDLY::node_handle
MEDDLY::union_mxd::compute_ext(node_handle a, node_handle b)
{
  // Get level information
  const int aLevel = arg1F->getNodeLevel(a);
  const int bLevel = arg2F->getNodeLevel(b);
  int resultLevel = ABS(topLevel(aLevel, bLevel));

  MEDDLY_DCASSERT(resF->isExtensibleLevel(resultLevel));

  const int dwnLevel = resF->downLevel(resultLevel);

  // Initialize readers
  unpacked_node *A = (aLevel < resultLevel)
    ? unpacked_node::newRedundant(arg1F, resultLevel, a, SPARSE_ONLY)
    : arg1F->newUnpacked(a, SPARSE_ONLY)
    ;

  unpacked_node *B = (bLevel < resultLevel)
    ? unpacked_node::newRedundant(arg2F, resultLevel, b, SPARSE_ONLY)
    : arg2F->newUnpacked(b, SPARSE_ONLY)
    ;

  // Initialize result writer
  unsigned resultSize = A->getSize() + B->getSize() + 1 + 1;
  unpacked_node* C = unpacked_node::newSparse(resF, resultLevel, resultSize);

  const node_handle A_ext_d = A->isExtensible()? A->ext_d(): 0;
  int last_nz = A->getSize()-1;
  for ( ; last_nz >= 0 && A->down(last_nz) == 0; last_nz--);
  const unsigned int A_nnzs = last_nz + 1;
  const int A_last_index = last_nz >= 0? A->index(last_nz): -1;

  const node_handle B_ext_d = B->isExtensible()? B->ext_d(): 0;
  last_nz = B->getSize()-1;
  for ( ; last_nz >= 0 && B->down(last_nz) == 0; last_nz--);
  const unsigned int B_nnzs = last_nz + 1;
  const int B_last_index = last_nz >= 0? B->index(last_nz): -1;

  const int max_a_b_last_index = MAX(A_last_index, B_last_index);

  unsigned nnz = 0;
  unsigned A_curr_index = 0;
  unsigned B_curr_index = 0;
  for ( ; A_curr_index < A_nnzs && B_curr_index < B_nnzs; ) {
    // get a_i, a_d, b_i, b_d
    unsigned a_i, b_i;
    node_handle a_d, b_d;
    a_i = A->index(A_curr_index);
    b_i = B->index(B_curr_index);
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

    MEDDLY_DCASSERT(a_d != 0 || b_d != 0);

    // compute union(a_d, b_d)
    unsigned index = (a_d? a_i: b_i);
    node_handle down = compute_r(int(index), dwnLevel, a_d, b_d);

    // if union is non-zero, add it to the new node
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  }
  for ( ; A_curr_index < A_nnzs; A_curr_index++) {
    // do union(a_i, b_ext_i)
    unsigned index = A->index(A_curr_index);
    node_handle down = compute_r(int(index), dwnLevel, A->down(A_curr_index), B_ext_d);
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  }
  for ( ; B_curr_index < B_nnzs; B_curr_index++) {
    // do union(a_ext_i, b_i)
    unsigned index = B->index(B_curr_index);
    node_handle down = compute_r(int(index), dwnLevel, A_ext_d, B->down(B_curr_index));
    if (down) {
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = index;
      // C->d_ref(nnz) = down;
      nnz++;
    }
  }
  if (A->isExtensible() || B->isExtensible()) {
    const unsigned index = max_a_b_last_index+1;
    node_handle down = compute_r(index, dwnLevel, A_ext_d, B_ext_d);
    if (down) {
      MEDDLY_DCASSERT(index >= 0);
      C->setSparse(nnz, index, down);
      // C->i_ref(nnz) = unsigned(index);
      // C->d_ref(nnz) = down;
      C->markAsExtensible();
      nnz++;
    } else {
      C->markAsNotExtensible();
    }
  } else {
    C->markAsNotExtensible();
  }
  C->shrink(nnz);

  // cleanup
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(A);

  // reduce result
  edge_value ev;
  node_handle result;
  resF->createReducedNode(C, ev, result);
  MEDDLY_DCASSERT(ev.isVoid());
  return result;
}
#endif

#else

// ******************************************************************
// *                                                                *
// *                        union_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::union_mxd : public binary_operation {
    public:
        union_mxd(forest* arg1, forest* arg2, forest* res);
        virtual union_mxd();

        virtual void compute(const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        inline node_handle makeChainTo(node_handle p, int K, int L)
        {
            if (go_by_levels) return p;
            if (identity_chains) {
                MEDDLY_DCASSERT(!redundant_chains);
                return resF->makeIdentitiesTo(p, K, L);
            }
            if (redundant_chains) {
                return resF->makeRedundantsTo(p, K, L);
            }
            return p;
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

    protected:
        node_handle _compute(node_handle A, node_handle B, int L);
        node_handle _compute_primed(int in, node_handle A, node_handle B,
                int L);


    private:
        ct_entry_type* ct;

        bool identity_chains;
        bool redundant_chains;

        bool go_by_levels;
};

// ******************************************************************

MEDDLY::union_mxd::union_mxd(forest* arg1, forest* arg2, forest* res)
  : binary_operation(INTER_cache, arg1, arg2, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    //
    // Determine if we can 'skip levels', and how to add chains
    // of nodes when needed.
    //
    identity_chains = arg1->isIdentityReduced() || arg2->isIdentityReduced();
    redundant_chains = arg1->isFullyReduced() || arg2->isFullyReduced();
    go_by_levels = identity_chains && redundant_chains;

    //
    // Possible input forests:
    //
    //      quasi, quasi    :   no level skipping will occur, moot
    //      quasi, fully    :   (same)
    //      quasi, ident    :   (same)
    //      fully, quasi    :   (same)
    //      ident, quasi    :   (same)
    //
    //      fully, fully    :   can skip levels, build redundant chains
    //      ident, ident    :   can skip levels, build identity chains
    //
    //      fully, ident    :
    //      ident, fully    :   must go by levels to handle unioning
    //                          of skipped redundant/identity nodes.
    //                          compute table entries add level info.
    //                          go_by_levels will be true only for
    //                          these two cases.
    //

    ct = new ct_entry_type("union_mxd");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    if (go_by_levels) {
        ct->setFixed(INTEGER, arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::union_mdd::~union_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::union_mdd::compute(const edge_value &av, node_handle ap,
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
MEDDLY::union_mxd::_compute(node_handle A, node_handle B, int L)
{
    MEDDLY_DCASSERT(L>=0);
    // ...
    // BIG TBD HERE
    //
}

MEDDLY::node_handle
MEDDLY::union_mxd::_compute_primed(int in, node_handle A, node_handle B,
        const int Clevel)
{
    MEDDLY_DCASSERT(Clevel<0);

    //
    // Determine level information
    //
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
#ifdef TRACE
    std::cout << "union_mxd::_compute_primed(" << in << ", " << A << ", " << B << ", " << Clevel << ")\n";
    std::cout << "\t" << A << " level " << Alevel << "\n";
    std::cout << "\t" << B << " level " << Blevel << "\n";
#endif

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   patternNode(arg1F, Clevel, in, A, FULL_ONLY)
        :   arg1F->newUnpacked(A, FULL_ONLY);

    unpacked_node* Bu = (Blevel != Clevel)
        ?   patternNode(arg2F, Clevel, in, B, FULL_ONLY)
        :   arg2F->newUnpacked(B, FULL_ONLY);

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i,
            _compute(Au->down(i), Bu->down(i),
                forest::downLevel(Cu->getLevel()))
        );
    }

#ifdef TRACE
    ostream_output out(std::cout);
    std::cout << "\t";
    Cu->show(out, true);
    std::cout << "\n";
#endif

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
    std::cout << "union_mxd::_compute_primed(" << in << ", " << A << ", "
              << B << ", " << Clevel << ") = " << C << "\n\t";
    resF->showNode(out, C, SHOW_DETAILS);
    out.put('\n');
#endif
    return C;
}


#endif

// ******************************************************************
// *                                                                *
// *                    union_min_evplus  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::union_min_evplus : public generic_binary_evplus {
    public:
        union_min_evplus(forest* arg1, forest* arg2, forest* res);

        virtual ct_entry_key* findResult(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle &c);
        virtual void saveResult(ct_entry_key* key, long aev, node_handle a,
            long bev, node_handle b, long cev, node_handle c);

        virtual bool checkTerminals(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle& c);
};

MEDDLY::union_min_evplus::union_min_evplus(forest* a, forest* b, forest* c)
  : generic_binary_evplus(UNION_cache, a, b, c)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::union_min_evplus::findResult(long aev, node_handle a,
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
  MEDDLY_DCASSERT(cev == 0);
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev = MIN(aev, bev);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::union_min_evplus::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  MEDDLY_DCASSERT(c == 0 || cev == MIN(aev, bev));
  CTresult[0].reset();
  CTresult[0].writeL(0);   //   Why always 0?
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::union_min_evplus::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == 0) {
    if (b == 0) {
      cev = 0;
      c = 0;
      return true;
    }
    else if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg1F->isTerminalNode(a) && aev <= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && bev <= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MIN(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  return false;
}

// ******************************************************************
// *                                                                *
// *                  union_min_evplus_mxd  class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::union_min_evplus_mxd : public generic_binary_evplus_mxd {
    public:
        union_min_evplus_mxd(forest* arg1, forest* arg2, forest* res);

        virtual ct_entry_key* findResult(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle &c);
        virtual void saveResult(ct_entry_key* key, long aev, node_handle a,
            long bev, node_handle b, long cev, node_handle c);

        virtual bool checkTerminals(long aev, node_handle a,
            long bev, node_handle b, long& cev, node_handle& c);
};

MEDDLY::union_min_evplus_mxd::union_min_evplus_mxd(
  forest* arg1, forest* arg2, forest* res)
  : generic_binary_evplus_mxd(UNION_cache, arg1, arg2, res)
{
    operationCommutes();

    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::EVPLUS);
}

MEDDLY::ct_entry_key* MEDDLY::union_min_evplus_mxd::findResult(long aev, node_handle a,
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
  MEDDLY_DCASSERT(cev == 0);
  c = resF->linkNode(CTresult[0].readN());
  if (c != 0) {
    cev = MIN(aev, bev);
  }
  CT0->recycle(CTsrch);
  return 0;
}

void MEDDLY::union_min_evplus_mxd::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  MEDDLY_DCASSERT(c == 0 || cev == MIN(aev, bev));
  CTresult[0].reset();
  CTresult[0].writeL(0);   // why always 0?
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

bool MEDDLY::union_min_evplus_mxd::checkTerminals(long aev, node_handle a, long bev, node_handle b,
  long& cev, node_handle& c)
{
  if (a == 0) {
    if (b == 0) {
      cev = 0;
      b = 0;
      return true;
    }
    else if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (b == 0) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg1F->isTerminalNode(a) && aev <= bev) {
    if (arg1F == resF) {
      cev = aev;
      c = resF->linkNode(a);
      return true;
    }
    else {
      return false;
    }
  }
  if (arg2F->isTerminalNode(b) && bev <= aev) {
    if (arg2F == resF) {
      cev = bev;
      c = resF->linkNode(b);
      return true;
    }
    else {
      return false;
    }
  }
  if (a == b) {
    if (arg1F == arg2F && arg2F == resF) {
      cev = MIN(aev, bev);
      c = resF->linkNode(a);
      return true;
    }
    else {
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
MEDDLY::UNION(forest* a, forest* b, forest* c)
{
    if (!a || !b || !c) {
        return nullptr;
    }
    binary_operation* bop =  UNION_cache.find(a, b, c);
    if (bop) {
        return bop;
    }

    if (c->getEdgeLabeling() == edge_labeling::MULTI_TERMINAL) {
        if (c->isForRelations()) {
            return UNION_cache.add(new union_mxd(a, b, c));
        } else {
            return UNION_cache.add(new union_mdd(a, b, c));
        }
    }

    if (c->getEdgeLabeling() == edge_labeling::EVPLUS) {
        if (c->isForRelations()) {
            return UNION_cache.add(new union_min_evplus_mxd(a, b, c));
        } else {
            return UNION_cache.add(new union_min_evplus(a, b, c));
        }
    }

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::UNION_init()
{
    UNION_cache.reset("Union");
}

void MEDDLY::UNION_done()
{
    MEDDLY_DCASSERT(UNION_cache.isEmpty());
}

