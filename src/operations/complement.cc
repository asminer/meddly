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
#include "../forests/mt.h"
#include "complement.h"

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"

#include "../ct_vector.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../operators.h"

namespace MEDDLY {
    class compl_mdd;
    class compl_mxd;

    unary_list COMPL_cache;
};

// #define DEBUG_MXD_COMPL

// #define OLD_OPER

// #define TRACE

// ******************************************************************
// *                                                                *
// *                        compl_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mdd : public unary_operation {
    public:
        compl_mdd(forest* arg, forest* res);
        virtual ~compl_mdd();

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, int L);

    private:
        ct_entry_type* ct;
};

// ******************************************************************

MEDDLY::compl_mdd::compl_mdd(forest* arg, forest* res)
    : unary_operation(COMPL_cache, arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("mdd_complement");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::compl_mdd::~compl_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::compl_mdd::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    cp = _compute(ap, L);
}

MEDDLY::node_handle MEDDLY::compl_mdd::_compute(node_handle A, int L)
{
    //
    // Terminal cases
    //
    if (argF->isTerminalNode(A)) {
        bool ta;
        argF->getValueFromHandle(A, ta);
        return resF->makeRedundantsTo(resF->handleForValue(!ta), 0, L);
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        return resF->makeRedundantsTo(resF->linkNode(res[0].getN()), Alevel, L);
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = argF->newUnpacked(A, FULL_ONLY);
    unpacked_node* Cu = unpacked_node::newFull(resF, Alevel, Au->getSize());

    //
    // Build result node
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i, _compute(Au->down(i), Cu->getLevel()-1));
    }

    //
    // Reduce / cleanup
    //
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

    return resF->makeRedundantsTo(C, Alevel, L);
}

// ******************************************************************
// *                                                                *
// *                        compl_mxd  class                        *
// *                                                                *
// ******************************************************************

#ifdef OLD_OPER

class MEDDLY::compl_mxd : public unary_operation {
    public:
        compl_mxd(forest* arg, forest* res);

        virtual void computeDDEdge(const dd_edge& a, dd_edge& b, bool userFlag);

        node_handle compute_r(int in, int k, node_handle a);
};

MEDDLY::compl_mxd::compl_mxd(forest* arg, forest* res)
 : unary_operation(COMPL_cache, 1, arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct_entry_type* et = new ct_entry_type(COMPL_cache.getName(), "IN:N");
    et->setForestForSlot(1, arg);
    et->setForestForSlot(3, res);
    registerEntryType(0, et);
    buildCTs();
}

void MEDDLY::compl_mxd::computeDDEdge(const dd_edge& a, dd_edge& b, bool userFlag)
{
  node_handle result = compute_r(-1, argF->getMaxLevelIndex(), a.getNode());
  b.set(result);
}

MEDDLY::node_handle MEDDLY::compl_mxd::compute_r(int in, int k, node_handle a)
{
    // Check terminals
    if ( (0==k) || (argF->isTerminalNode(a) && resF->isFullyReduced()) ) {
        bool ta;
        argF->getValueFromHandle(a, ta);
        return argF->handleForValue(!ta);
    }

  // Check compute table
  ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeI(k);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    node_handle ans = CTresult[0].readN();
#ifdef DEBUG_MXD_COMPL
    fprintf(stderr, "\tin CT:   compl_mxd(%d, %d) : %d\n", ht, a, ans);
#endif
    CT0->recycle(CTsrch);
    return resF->linkNode(ans);
  }

#ifdef DEBUG_MXD_COMPL
  fprintf(stderr, "\tstarting compl_mxd(%d, %d)\n", ht, a);
#endif

  // Initialize unpacked node
  const unsigned size = unsigned(resF->getLevelSize(k));
  const int aLevel = argF->getNodeLevel(a);
  MEDDLY_DCASSERT(!isLevelAbove(aLevel, k));
  unpacked_node* A = unpacked_node::New(argF);
  bool canSave = true;
  if (aLevel == k) {
    argF->unpackNode(A, a, FULL_ONLY);
  } else if (k>0 || argF->isFullyReduced()) {
    A->initRedundant(argF, k, a, FULL_ONLY);
  } else {
    MEDDLY_DCASSERT(in>=0);
    A->initIdentity(argF, k, unsigned(in), a, FULL_ONLY);
    canSave = false;
  }
  unpacked_node* C = unpacked_node::newFull(resF, k, size);

  // recurse
  int nextLevel = argF->downLevel(k);
  unsigned nnz = 0;
  bool addRedundentNode=(resF->isQuasiReduced() && (k>0 || k<-1));

  // recurse
  for (unsigned i=0; i<size; i++) {
      node_handle cdi = compute_r(int(i), nextLevel, A->down(i));

      if (cdi != resF->getTransparentNode()) nnz++;

    if(addRedundentNode && resF->isTerminalNode(cdi) && cdi!=resF->getTransparentNode()){
      cdi=((mt_forest*)resF)->makeNodeAtLevel(nextLevel, cdi);
    }

    C->setFull(i, cdi);
  }

  // reduce, save in CT
  unpacked_node::Recycle(A);
  edge_value ev;
  node_handle result;
  resF->createReducedNode(C, ev, result, in);
  MEDDLY_DCASSERT(ev.isVoid());
  if (k<0 && 1==nnz) canSave = false;
  if (canSave) {
    CTresult[0].reset();
    CTresult[0].writeN(result);
    CT0->addEntry(CTsrch, CTresult[0]);
  } else {
    CT0->recycle(CTsrch);
  }
  return result;
}

#else

class MEDDLY::compl_mxd : public unary_operation {
    public:
        compl_mxd(forest* arg, forest* res);
        virtual ~compl_mxd();

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, int L);

        node_handle _compute_primed(int in, node_handle A, int L);

        /*
            Build a complemented identity pattern:
                [ p 1 1 ... 1 ]
                [ 1 p 1 ... 1 ]
                [ 1 1 p ... 1 ]
                [ :         : ]
                [ 1 1 1 ... p ]


                @param  p   Diagonal pointer
                @param  K   Start the pattern above this level
                @param  L   Stop the pattern at this level
        */
        node_handle _identity_complement(node_handle p, int K, int L);

    protected:
        inline unpacked_node* patternNode(forest* F, int L, int in,
                node_handle A, node_storage_flags fs)
        {
            MEDDLY_DCASSERT(L<0);
            if (0==A || F->isFullyReduced()) {
                return unpacked_node::newRedundant(F, L, A, fs);
            }
            if (F->isIdentityReduced()) {
                return unpacked_node::newIdentity(F, L, in, A, fs);
            }
            std::cout << "L: " << L << "\n";
            std::cout << "A: " << A << "\n";
            std::cout << "A.level: " << F->getNodeLevel(A) << "\n";
            MEDDLY_DCASSERT(false);
        }

    protected:
        inline node_handle identity_complement(node_handle p, int K, int L)
        {
            MEDDLY_DCASSERT(L>=0);
            MEDDLY_DCASSERT(K>=0);
            if (0==L) return p;
            if (K==L) return p;
            MEDDLY_DCASSERT(L>=K);
            return _identity_complement(p, K, L);
        }

    private:
        ct_entry_type* ct;

        ostream_output out;
};

// ******************************************************************

MEDDLY::compl_mxd::compl_mxd(forest* arg, forest* res)
    : unary_operation(COMPL_cache, arg, res), out(std::cout)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("mxd_complement");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::compl_mxd::~compl_mxd()
{
    ct->markForDestroy();
}

void MEDDLY::compl_mxd::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    out.indentation(0);
    cp = _compute(ap, L);
}

MEDDLY::node_handle MEDDLY::compl_mxd::_compute(node_handle A, int L)
{
    //
    // Terminal cases
    //
    if (0==A) {
        // Terminal zero;
        // that's a matrix of all zeroes.
        return resF->makeRedundantsTo(resF->handleForValue(true), 0, L);
    }

    if (A < 0) {
        // Terminal one; that's either a matrix of all ones
        // or an identity matrix.

        if (0==L || argF->isFullyReduced()) {
            return resF->makeRedundantsTo(resF->handleForValue(false), 0, L);
        }

        MEDDLY_DCASSERT(argF->isIdentityReduced());
        return identity_complement(0, 0, L);
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
    const int Clevel = ABS(Alevel);

#ifdef TRACE
    out << "compl_mxd::_compute(" << A << ", " << L << ")\n";
    out << A << " level " << Alevel << "\n";
    out << "result level " << Clevel << " before chain\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        node_handle C = resF->linkNode(res[0].getN());
        if (argF->isIdentityReduced()) {
            C = identity_complement(C, Clevel, L);
        } else {
            C = resF->makeRedundantsTo(C, Clevel, L);
        }
#ifdef TRACE
        out << "CT hit " << res[0].getN() << "\n";
        out << "after chain " << C << "\n";
#endif
        return C;
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   unpacked_node::newRedundant(argF, Clevel, A, FULL_ONLY)
        :   argF->newUnpacked(A, FULL_ONLY);
    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i,
            _compute_primed(int(i), Au->down(i), forest::downLevel(Clevel))
        );
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "compl_mxd::_compute(" << A << ", " << L
              << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);

    if (argF->isIdentityReduced()) {
        C = identity_complement(C, Clevel, L);
    } else {
        C = resF->makeRedundantsTo(C, Clevel, L);
    }

#ifdef TRACE
    out << "chain to level " << L << " = " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif
    return C;
}

MEDDLY::node_handle MEDDLY::compl_mxd::_compute_primed(int in,
            node_handle A, const int Clevel)
{
    MEDDLY_DCASSERT(Clevel<0);

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
#ifdef TRACE
    out << "compl_mxd::_compute_primed(" << in << ", " << A << ", " << Clevel << ")\n";
    out << A << " level " << Alevel << "\n";
#endif

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   patternNode(argF, Clevel, in, A, FULL_ONLY)
        :   argF->newUnpacked(A, FULL_ONLY);

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i, _compute(Au->down(i), forest::downLevel(Clevel)));
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "compl_mxd::_compute_primed(" << in << ", " << A << ", " << Clevel << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C, in);
    MEDDLY_DCASSERT(dummy.isVoid());

#ifdef TRACE
    out << "compl_mxd::_compute_primed(" << in << ", " << A << ", "
              << Clevel << ") = " << C << "\n  ";
    resF->showNode(out, C, SHOW_DETAILS);
    out.put('\n');
#endif
    return C;
}

MEDDLY::node_handle MEDDLY::compl_mxd::_identity_complement(node_handle p,
    int K, int L)
{
    MEDDLY_DCASSERT(L>0);
    MEDDLY_DCASSERT(K>=0);
    unpacked_node* Uun;
    unpacked_node* Upr;
    edge_value ev;

    terminal ONE(true);
    node_handle chain_to_one = resF->makeRedundantsTo(ONE.getHandle(), 0, K);

    for (K++; K<=L; K++) {
        Uun = unpacked_node::newFull(resF, K, resF->getLevelSize(K));

        for (unsigned i=0; i<Uun->getSize(); i++) {
            Upr = unpacked_node::newFull(resF, -K, Uun->getSize());
            for (unsigned j=0; j<Uun->getSize(); j++) {
                Upr->setFull(j, resF->linkNode(chain_to_one));
            }
            Upr->setFull(i, (i ? resF->linkNode(p) : p));
            node_handle h;
            resF->createReducedNode(Upr, ev, h);
            MEDDLY_DCASSERT(ev.isVoid());
            Uun->setFull(i, h);
        }

        resF->createReducedNode(Uun, ev, p);
        MEDDLY_DCASSERT(ev.isVoid());

        chain_to_one = resF->makeRedundantsTo(chain_to_one, K-1, K);
    } // for k
    resF->unlinkNode(chain_to_one);
    return p;
}

#endif

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_operation* MEDDLY::COMPLEMENT(forest* arg, forest* res)
{
    if (!arg || !res) return nullptr;

    unary_operation* uop = COMPL_cache.find(arg, res);
    if (uop) {
        return uop;
    }

    if (arg->isForRelations()) {
        return COMPL_cache.add(new compl_mxd(arg, res));
    } else {
        return COMPL_cache.add(new compl_mdd(arg, res));
    }
}

void MEDDLY::COMPLEMENT_init()
{
    COMPL_cache.reset("Complement");
}

void MEDDLY::COMPLEMENT_done()
{
    MEDDLY_DCASSERT( COMPL_cache.isEmpty() );
}

