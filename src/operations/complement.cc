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

namespace MEDDLY {
    class compl_mdd;
    class compl_mxd;

    unary_list COMPL_cache;
};

// #define OLD_COMP

// #define DEBUG_MXD_COMPL

// ******************************************************************
// *                                                                *
// *                        compl_mdd  class                        *
// *                                                                *
// ******************************************************************

#ifdef OLD_COMP

class MEDDLY::compl_mdd : public unary_operation {
    public:
        compl_mdd(forest* arg, forest* res);

        virtual void computeDDEdge(const dd_edge& a, dd_edge& b, bool userFlag);

        node_handle compute_r(node_handle a);

        inline ct_entry_key*
        findResult(node_handle a, node_handle &b)
        {
            ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
            MEDDLY_DCASSERT(CTsrch);
            CTsrch->writeN(a);
            CT0->find(CTsrch, CTresult[0]);
            if (!CTresult[0]) return CTsrch;
            b = resF->linkNode(CTresult[0].readN());
            CT0->recycle(CTsrch);
            return 0;
        }
        inline node_handle saveResult(ct_entry_key* Key,
            node_handle a, node_handle b)
        {
            CTresult[0].reset();
            CTresult[0].writeN(b);
            CT0->addEntry(Key, CTresult[0]);
            return b;
        }
};

MEDDLY::compl_mdd::compl_mdd(forest* arg, forest* res)
    : unary_operation(COMPL_cache, 1, arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct_entry_type* et = new ct_entry_type(COMPL_cache.getName(), "N:N");
    et->setForestForSlot(0, arg);
    et->setForestForSlot(2, res);
    registerEntryType(0, et);
    buildCTs();
}

void MEDDLY::compl_mdd::computeDDEdge(const dd_edge& a, dd_edge& b, bool userFlag)
{
  node_handle result = compute_r(a.getNode());
  const int num_levels = resF->getMaxLevelIndex();
  if (userFlag && result != resF->getTransparentNode() && resF->isQuasiReduced() &&
      resF->getNodeLevel(result) < num_levels) {
    node_handle temp = ((mt_forest*)resF)->makeNodeAtLevel(num_levels, result);
    resF->unlinkNode(result);
    result = temp;
  }
  b.set(result);
}

MEDDLY::node_handle MEDDLY::compl_mdd::compute_r(node_handle a)
{
    // Check terminals
    if (argF->isTerminalNode(a)) {
        bool ta;
        argF->getValueFromHandle(a, ta);
        return argF->handleForValue(!ta);
    }

  // Check compute table
  node_handle b;
  ct_entry_key* Key = findResult(a, b);
  if (0==Key) return b;

  const int level = argF->getNodeLevel(a);
  const unsigned size = unsigned(resF->getLevelSize(level));
  bool addRedundentNode=(resF->isQuasiReduced() && level>1);

  // Initialize unpacked nodes
  unpacked_node* A = argF->newUnpacked(a, FULL_ONLY);
  unpacked_node* C = unpacked_node::newFull(resF, level, size);

  // recurse
  for (unsigned i=0; i<size; i++) {

    node_handle cdi = compute_r(A->down(i));

    if(addRedundentNode && resF->isTerminalNode(cdi) && cdi!=resF->getTransparentNode()){
    	cdi =((mt_forest*)resF)->makeNodeAtLevel(level-1, cdi);
    }

    C->setFull(i, cdi);

  }

  // cleanup and Reduce
  unpacked_node::Recycle(A);
  edge_value ev;
  resF->createReducedNode(C, ev, b);
  MEDDLY_DCASSERT(ev.isVoid());

  // Add to compute table
  return saveResult(Key, a, b);
}

#else

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

    ct = new ct_entry_type("union");
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
        return resF->makeRedundantsTo(resF->handleForValue(!ta), L);
    }

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        return resF->makeRedundantsTo(resF->linkNode(res[0].getN()), L);
    }

    //
    // Do computation
    //

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);

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

    return resF->makeRedundantsTo(C, L);
}

#endif

// ******************************************************************
// *                                                                *
// *                        compl_mxd  class                        *
// *                                                                *
// ******************************************************************

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

