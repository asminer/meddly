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

#include "../ops_builtin.h" // for COPY
#include "../oper_binary.h"
#include "../oper_unary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class union_mdd;
    class union_mxd;

    class union_mt;

    class union_min_evplus;
    class union_min_evplus_mxd;

    binary_list UNION_cache;
};

// #define TRACE
#define USE_PRIMED_CACHE

#ifdef TRACE
#include "../operators.h"
#endif

#define NEW_UNION

// ******************************************************************
// *                                                                *
// *                         union_mt class                         *
// *                                                                *
// ******************************************************************

/**
    Union operation for multi-terminal forests.
*/

class MEDDLY::union_mt : public binary_operation {
    public:
        union_mt(forest* arg1, forest* arg2, forest* res);
        virtual ~union_mt();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp);

    private:
        void _compute(int L, unsigned in,
                node_handle A, node_handle B, node_handle &C);

    private:
        // inlined helpers

        inline int topLevelOf(int L, int alevel, int blevel) const
        {
            if (by_levels) return L;
            if (resF->isForRelations()) {
                if (both_identity) {
                    const int t = MXD_levels::topUnprimed(alevel, blevel);
                    return (t == MXD_levels::upLevel(L)) ? L : t;
                } else {
                    return MXD_levels::topLevel(alevel, blevel);
                }
            } else {
                return MDD_levels::topLevel(alevel, blevel);
            }
        }

        inline void chainToLevel(node_handle &C, int Clevel, int L, unsigned in)
        {
            //
            // Add nodes from Clevel to L
            //
            if (both_identity) {
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
#ifdef USE_PRIMED_CACHE
        ct_entry_type* ct_primed;
#endif
        unary_operation* copy_arg1res;
        unary_operation* copy_arg2res;

#ifdef TRACE
        ostream_output out;
        unsigned top_count;
#endif

        bool by_levels;
        bool forced_by_levels;
        bool both_fully;
        bool both_identity;
};

// ******************************************************************

MEDDLY::union_mt::union_mt(forest* arg1, forest* arg2, forest* res)
    : binary_operation(arg1, arg2, res)
#ifdef TRACE
      , out(std::cout), top_count(0)
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
    //      quasi v quasi   :   no level skipping will occur
    //      quasi v fully   :
    //      quasi v ident   :
    //      fully v quasi   :
    //      ident v quasi   :
    //
    //      ident v ident   :   identity pattern: skip by unprimed
    //
    //                          [p 0 0]   [q 0 0]   [ pvq  0   0  ]
    //                          [0 p 0] v [0 q 0] = [  0  pvq  0  ]
    //                          [0 0 p]   [0 0 q]   [  0   0  pvq ]
    //
    //      ident v fully   :   force by levels, and store levels in
    //      fully v ident   :   the compute table entry
    //
    //                          [p p p]   [q 0 0]   [ pvq  p   p  ]
    //                          [p p p] v [0 q 0] = [  p  pvq  p  ]
    //                          [p p p]   [0 0 q]   [  p   p  pvq ]
    //
    //      fully v fully   :   fully pattern: skip by top level
    //
    //                          [p p p]   [q q q]   [ pvq pvq pvq ]
    //                          [p p p] v [q q q] = [ pvq pvq pvq ]
    //                          [p p p]   [q q q]   [ pvq pvq pvq ]
    //
    //

    both_fully = arg1F->isFullyReduced() && arg2F->isFullyReduced();
    both_identity = arg1F->isIdentityReduced() && arg2F->isIdentityReduced();

    by_levels = !(both_fully || both_identity);
    forced_by_levels =
        (arg1F->isFullyReduced() && arg2F->isIdentityReduced())
        ||
        (arg2F->isFullyReduced() && arg1F->isIdentityReduced())
    ;

    ct = new ct_entry_type("union");
    if (forced_by_levels) {
        ct->setFixed('I', arg1, arg2);
    } else {
        ct->setFixed(arg1, arg2);
    }
    ct->setResult(res);
    ct->doneBuilding();

#ifdef USE_PRIMED_CACHE
    //
    // If we're skipping by unprimed levels,
    // keep a second CT for the primed level computations
    // we are able to save. This is needed to differentiate
    // the unprimed level result and the primed level result.
    ct_primed = nullptr;
    if (both_identity) {
        ct_primed = new ct_entry_type("union_pr");
        ct_primed->setFixed(arg1, arg2);
        ct_primed->setResult(res);
        ct_primed->doneBuilding();
    }
#endif // USE_PRIMED_CACHE
}

MEDDLY::union_mt::~union_mt()
{
    ct->markForDestroy();
#ifdef USE_PRIMED_CACHE
    if (ct_primed) ct_primed->markForDestroy();
#endif
}

void MEDDLY::union_mt::compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                const edge_value &bv, node_handle bp,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    MEDDLY_DCASSERT(bv.isVoid());
#ifdef TRACE
    out.indentation(0);
    ++top_count;
    out << "Union #" << top_count << " begin\n";
#endif
    _compute(L, in, ap, bp, cp);
#ifdef TRACE
    out << "Union #" << top_count << " end\n";
#endif
    cv.set();
}

void MEDDLY::union_mt::_compute(int L, unsigned in,
        node_handle A, node_handle B, node_handle &C)
{
    // **************************************************************
    //
    // Check terminal cases
    //
    // **************************************************************
    if (0==A && 0==B) {
        // This is correct, even if we're forced to go by levels.
        C = 0;
        return;
    }

    if (0==A) {
        //
        // Result is B
        //
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg2res);
        copy_arg2res->compute(L, in, dummy, B, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
        return;
    }

    if ( (0 == B) || ((A==B)&&(arg1F==arg2F)) )
    {
        //
        // Result is A
        //
        edge_value dummy;
        dummy.set();
        MEDDLY_DCASSERT(copy_arg1res);
        copy_arg1res->compute(L, in, dummy, A, dummy, C);
        MEDDLY_DCASSERT(dummy.isVoid());
        return;
    }

    //
    // Both terminal one; result is fully-fully one unless
    // both argument forests are identity reduced; then it is I.
    //
    if (arg1F->isTerminalNode(A) && arg2F->isTerminalNode(B)) {
        terminal tt(true, resF->getTerminalType());
        if (both_identity) {
            C = resF->makeIdentitiesTo(tt.getHandle(), 0, L, in);
        } else {
            C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
        }
        return;
    }

    //
    // One argument is a fully-fully ONE.
    // Result is fully-fully one.
    //
    if ( (arg1F->isTerminalNode(A) && arg1F->isFullyReduced())
        || (arg2F->isTerminalNode(B) && arg2F->isFullyReduced()) )
    {
        terminal tt(true, resF->getTerminalType());
        C = resF->makeRedundantsTo(tt.getHandle(), 0, L);
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
    const bool useCT = !both_identity || Clevel>0;
#ifdef USE_PRIMED_CACHE
    const bool useCTpr = !useCT && (L<0) &&
        (Alevel == L) && (Blevel == L) && (Clevel == L);
    MEDDLY_DCASSERT(!useCTpr || ct_primed);
#endif

#ifdef TRACE
    out << "union_mt::_compute(" << L << ", " << in << ", " << A << ", "
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
    ct_vector key(forced_by_levels ? 3 : 2);
    ct_vector res(1);
    if (forced_by_levels) {
        key[0].setI(L);
        key[1].setN(A);
        key[2].setN(B);
    } else {
        key[0].setN(A);
        key[1].setN(B);
    }

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
#ifdef USE_PRIMED_CACHE
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
#endif // USE_PRIMED_CACHE


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
            Au = unpacked_node::newIdentity(arg1F, Clevel, in, A, FULL_ONLY);
            MEDDLY_DCASSERT(Au->wasIdentity());
        } else {
            Au = unpacked_node::newRedundant(arg1F, Clevel, A, FULL_ONLY);
            MEDDLY_DCASSERT(!Au->wasIdentity());
        }
    } else {
        Au = arg1F->newUnpacked(A, FULL_ONLY);
        MEDDLY_DCASSERT(!Au->wasIdentity());
    }

    unpacked_node* Bu;
    if (Blevel != Clevel) {
        if (arg2F->isIdentityReduced() && Clevel<0) {
            MEDDLY_DCASSERT(Clevel == L);
            Bu = unpacked_node::newIdentity(arg2F, Clevel, in, B, FULL_ONLY);
            MEDDLY_DCASSERT(Bu->wasIdentity());
        } else {
            Bu = unpacked_node::newRedundant(arg2F, Clevel, B, FULL_ONLY);
            MEDDLY_DCASSERT(!Bu->wasIdentity());
        }
    } else {
        Bu = arg2F->newUnpacked(B, FULL_ONLY);
        MEDDLY_DCASSERT(!Bu->wasIdentity());
    }

    MEDDLY_DCASSERT(Au->getSize() == Bu->getSize());
    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

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
    for (unsigned i=0; i<Cu->getSize(); i++) {
        node_handle cd;
        _compute(Cnextlevel, i, Au->down(i), Bu->down(i), cd);
        Cu->setFull(i, cd);
    }

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "union_mt::_compute(" << L << ", " << in << ", " << A << ", "
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
#ifdef USE_PRIMED_CACHE
    } else if (useCTpr) {
        MEDDLY_DCASSERT(!Au->wasIdentity());
        MEDDLY_DCASSERT(!Bu->wasIdentity());
        res[0].setN(C);
        ct_primed->addCT(key, res);
#endif // USE_PRIMED_CACHE
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
        return UNION_cache.add(new union_mt(a, b, c));
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

