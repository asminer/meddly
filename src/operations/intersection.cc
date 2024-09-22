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

#include "../ops_builtin.h" // for COPY
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class inter_mt;

    class inter_mdd;
    class inter_mxd;
    class inter_max_evplus;

    binary_list INTER_cache;
};

#define NEW_INTER

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                         inter_mt class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::inter_mt : public binary_operation {
    public:
        inter_mt(forest* arg1, forest* arg2, forest* res);
        virtual ~inter_mt();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    protected:
        void _compute(int L, unsigned in,
                node_handle A, node_handle B, node_handle &C);

    private:
        inline int topLevelOf(int L, int Alevel, int Blevel) const
        {
            if (! resF->isForRelations()) {
                //
                // MDD; skip as much as we can
                //
                return MDD_levels::topLevel(Alevel, Blevel);
            }

            if (identity_pattern) {
                const int t = MXD_levels::topUnprimed(Alevel, Blevel);
                return (t == MXD_levels::upLevel(L)) ? L : t;
            } else {
                return MXD_levels::topLevel(Alevel, Blevel);
            }
        }

        inline void chainToLevel(node_handle &C, int Clevel, int L, unsigned in)
        {
            //
            // Add nodes from Clevel to L
            //
            if (identity_pattern) {
#ifdef TRACE
                out << "I chain to " << C << ", levels " << Clevel
                    << " to " << L << "\n";
#endif
                C = resF->makeIdentitiesTo(C, Clevel, L, in);
            } else {
#ifdef TRACE
                out << "X chain to " << C << ", levels " << Clevel
                    << " to " << L << "\n";
#endif
                C = resF->makeRedundantsTo(C, Clevel, L);
            }
        }

    private:
        ct_entry_type* ct;
        ct_entry_type* ct_primed;
        unary_operation* copy_arg1res;
        unary_operation* copy_arg2res;
        bool identity_pattern;

#ifdef TRACE
        ostream_output out;
#endif
};


// ******************************************************************

MEDDLY::inter_mt::inter_mt(forest* arg1, forest* arg2, forest* res)
  : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    //
    // Quick copying, even across forests, for terminal cases :)
    //
    copy_arg1res = COPY(arg1, res);
    copy_arg2res = COPY(arg2, res);

    //
    // How to handle different reductions
    //
    //      quasi ^ quasi   :   no level skipping will occur
    //      quasi ^ fully   :
    //      quasi ^ ident   :
    //      fully ^ quasi   :
    //      ident ^ quasi   :
    //
    //      ident ^ ident   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q 0 0]   [ p^q  0   0  ]
    //                          [0 p 0] ^ [0 q 0] = [  0  p^q  0  ]
    //                          [0 0 p]   [0 0 q]   [  0   0  p^q ]
    //
    //      ident ^ fully   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q q q]   [ p^q  0   0  ]
    //                          [0 p 0] ^ [q q q] = [  0  p^q  0  ]
    //                          [0 0 p]   [q q q]   [  0   0  p^q ]
    //
    //      fully ^ ident   :   identity pattern: skip by unprimed
    //
    //                          [p p p]   [q 0 0]   [ p^q  0   0  ]
    //                          [p p p] ^ [0 q 0] = [  0  p^q  0  ]
    //                          [p p p]   [0 0 q]   [  0   0  p^q ]
    //
    //      fully ^ fully   :   fully pattern: skip by top level
    //
    //                          [p p p]   [q q q]   [ p^q p^q p^q ]
    //                          [p p p] ^ [q q q] = [ p^q p^q p^q ]
    //                          [p p p]   [q q q]   [ p^q p^q p^q ]
    //
    //
    identity_pattern = (arg1->isIdentityReduced() && !arg2->isQuasiReduced())
                        ||
                       (arg2->isIdentityReduced() && !arg1->isQuasiReduced());

    ct = new ct_entry_type("intersection");
    // CT key:      node from forest arg1, node from forest arg2
    // CT result:   node from forest res
    ct->setFixed(arg1, arg2);
    ct->setResult(res);
    ct->doneBuilding();

    //
    // If we're skipping by unprimed levels,
    // keep a second CT for the primed level computations
    // we are able to save. This is needed to differentiate
    // the unprimed level result and the primed level result.
    ct_primed = nullptr;
    if (identity_pattern) {
        ct_primed = new ct_entry_type("difference_pr");
        ct_primed->setFixed(arg1, arg2);
        ct_primed->setResult(res);
        ct_primed->doneBuilding();
    }
}

MEDDLY::inter_mt::~inter_mt()
{
    ct->markForDestroy();
    if (ct_primed) ct_primed->markForDestroy();
}

void MEDDLY::inter_mt::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
#ifdef TRACE
    out.indentation(0);
#endif
    _compute(L, in, ap, bp, cp);
    cv.set();
}

void MEDDLY::inter_mt::_compute(int L, unsigned in,
            node_handle A, node_handle B, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (A==0 || B==0) {
        C = 0;
        return;
    }

    if (arg1F->isTerminalNode(A)) {
        if (L==0 || arg1F->isFullyReduced()) {
            // TRUE and B = B
            edge_value dummy;
            dummy.set();
            MEDDLY_DCASSERT(copy_arg2res);
            copy_arg2res->compute(L, in, dummy, B, dummy, C);
            MEDDLY_DCASSERT(dummy.isVoid());
            return;
        }
    }

    if (arg2F->isTerminalNode(B)) {
        if (L==0 || arg2F->isFullyReduced()) {
            // A and TRUE = A
            edge_value dummy;
            dummy.set();
            MEDDLY_DCASSERT(copy_arg1res);
            copy_arg1res->compute(L, in, dummy, A, dummy, C);
            MEDDLY_DCASSERT(dummy.isVoid());
            return;
        }
    }

    if ((A == B) && (arg1F==arg2F)) {
        // A and A = A
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in, dummy, A, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
        return;
    }

    //
    // Reorder A and B if they commute (same forest)
    //
    if (arg1F == arg2F) {
        if (A > B) {
            SWAP(A, B);
        }
    }

    // **************************************************************
    //
    // Determine level information
    //
    // **************************************************************
    const int Alevel = arg1F->getNodeLevel(A);
    const int Blevel = arg2F->getNodeLevel(B);
    const int Clevel = topLevelOf(L, Alevel, Blevel);
    const int Cnextlevel = resF->isForRelations()
        ? MXD_levels::downLevel(Clevel)
        : MDD_levels::downLevel(Clevel)
    ;
    const bool useCT = !identity_pattern || Clevel>0;
    const bool useCTpr = !useCT && (L<0) &&
        (Alevel == L) && (Blevel == L) && (Clevel == L);
    MEDDLY_DCASSERT(!useCTpr || ct_primed);

#ifdef TRACE
    out << "inter_mt::_compute(" << L << ", " << in << ", " << A << ", "
        << B << ")\n";
    out << A << " level " << Alevel << "\n";
    out << B << " level " << Blevel << "\n";
    out << "result level " << Clevel << "\n";
#endif


    // **************************************************************
    //
    // Check the compute table or primed compute table
    //
    // **************************************************************
    ct_vector key(2);
    ct_vector res(1);
    key[0].setN(A);
    key[1].setN(B);


    if (useCT && ct->findCT(key, res)) {
        //
        // 'main' compute table hit
        //
        C = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "CT hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        if (L == resF->getNodeLevel(C)) {
            // Make sure we don't point to a singleton
            // from the same index.
            C = resF->redirectSingleton(in, C);
        } else {
            chainToLevel(C, Clevel, L, in);
        }
        return;
        //
        // done 'main' compute table hit
        //
    }
    if (useCTpr && ct_primed->findCT(key, res)) {
        //
        // 'primed' compute table hit
        //
        C = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "CT' hit ";
        key.show(out);
        out << " -> ";
        res.show(out);
        out << "\n";
#endif
        if (L == resF->getNodeLevel(C)) {
            // Make sure we don't point to a singleton
            // from the same index.
            C = resF->redirectSingleton(in, C);
        }
        return;
        //
        // done 'primed' compute table hit
        //
    }

    // **************************************************************
    //
    // Compute table 'miss'; do computation
    //
    // **************************************************************

    //
    // Set up unpacked nodes
    //

    unpacked_node* Au;
    if (Alevel != Clevel) {
        if (arg1F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            // ^ if we skip too much, in is irrelevant
            Au = unpacked_node::newIdentity(arg1F, Clevel, in, A, SPARSE_ONLY);
            MEDDLY_DCASSERT(Au->wasIdentity());
        } else {
            Au = unpacked_node::newRedundant(arg1F, Clevel, A, SPARSE_ONLY);
            MEDDLY_DCASSERT(!Au->wasIdentity());
        }
    } else {
        Au = arg1F->newUnpacked(A, SPARSE_ONLY);
        MEDDLY_DCASSERT(!Au->wasIdentity());
    }

    unpacked_node* Bu;
    if (Blevel != Clevel) {
        if (arg2F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            Bu = unpacked_node::newIdentity(arg2F, Clevel, in, B, SPARSE_ONLY);
            MEDDLY_DCASSERT(Bu->wasIdentity());
        } else {
            Bu = unpacked_node::newRedundant(arg2F, Clevel, B, SPARSE_ONLY);
            MEDDLY_DCASSERT(!Bu->wasIdentity());
        }
    } else {
        Bu = arg2F->newUnpacked(B, SPARSE_ONLY);
        MEDDLY_DCASSERT(!Bu->wasIdentity());
    }

    unpacked_node* Cu = unpacked_node::newSparse(resF, Clevel,
            MIN(Au->getSize(), Bu->getSize()));

#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
    out << "\nB: ";
    Bu->show(out, true);
    out.indent_more();
    out.put('\n');
#endif

    //
    // Recurse
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
        const unsigned i = Au->index(za);
        MEDDLY_DCASSERT(i == Bu->index(zb));

        node_handle cd;
        _compute(Cnextlevel, i, Au->down(za), Bu->down(zb), cd);
        if (cd) {
            Cu->setSparse(zc++, i, cd);
        }
        za++;
        zb++;
    }
    Cu->resize(zc);
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "inter_mt::_compute(" << L << ", " << in << ", " << A << ", "
        << B << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  B: ";
    Bu->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    edge_value dummy;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT, if we can
    //
    if (useCT) {
        if (Au->wasIdentity() || Bu->wasIdentity()) {
            ct->noaddCT(key);
        } else {
            res[0].setN(C);
            ct->addCT(key, res);
        }
    } else if (useCTpr) {
        MEDDLY_DCASSERT(!Au->wasIdentity());
        MEDDLY_DCASSERT(!Bu->wasIdentity());
        res[0].setN(C);
        ct_primed->addCT(key, res);
    } else {
        // recycle the key
        ct->noaddCT(key);
    }

    //
    // Cleanup
    //
    unpacked_node::Recycle(Bu);
    unpacked_node::Recycle(Au);

    //
    // Adjust result
    //
    if (L == resF->getNodeLevel(C)) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        C = resF->redirectSingleton(in, C);
    } else {
        //
        // Add nodes but only from Clevel;
        // if the actual level of C is below Clevel it means
        // nodes were eliminated in the result forest.
        //
        chainToLevel(C, Clevel, L, in);
    }
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
        return INTER_cache.add(new inter_mt(a, b, c));
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

