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
#include "equant.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"
#include "../oper_binary.h"
#include "../ops_builtin.h"
namespace MEDDLY {
  class equant_opname;

  class equant_EV2EV;
};

// ******************************************************************
// *                                                                *
// *                       cycle_EV2MT  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::equant_EV2EV : public unary_operation {
  public:
    equant_EV2EV(unary_opname* oc, expert_forest* arg, expert_forest* res);

    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);

  protected:
    virtual void compute_r(long aev, node_handle a, int k, long& bev, node_handle& b);
    virtual node_handle compute_r( int k, node_handle a);
    binary_operation* mddUnion;

    inline ct_entry_key*
    findResult(long aev, node_handle a, long& bev, node_handle &b)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      bev = CTresult[0].readL();
      b = resF->linkNode(CTresult[0].readN());
      if (b != 0) {
        bev += aev;
      }
      CT0->recycle(CTsrch);
      return 0;
    }
    inline ct_entry_key*
    findResult(node_handle a, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      // CTsrch->writeN(b);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(ct_entry_key* Key,
      long aev, node_handle a, long bev, node_handle b)
    {
      CTresult[0].reset();
      CTresult[0].writeL(b == 0 ? 0L : bev - aev);
      CTresult[0].writeN(b);
      CT0->addEntry(Key, CTresult[0]);
      return b;
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

MEDDLY::equant_EV2EV::equant_EV2EV(unary_opname* oc, expert_forest* arg, expert_forest* res)
  : unary_operation(oc, 1, arg, res)
{
  MEDDLY_DCASSERT( argF->isForRelations());
  MEDDLY_DCASSERT(!resF->isForRelations());

  ct_entry_type* et = new ct_entry_type(oc->getName(), "N:N");
  et->setForestForSlot(0, arg);
  et->setForestForSlot(2, res);
  registerEntryType(0, et);
  buildCTs();
  dd_edge er(resF);
  mddUnion=getOperation(UNION, er,er,er);
  MEDDLY_DCASSERT(mddUnion);
}

void MEDDLY::equant_EV2EV::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{

  node_handle b = 0;
  b=compute_r(argF->getNumVariables(),arg.getNode());
  res.set(b);
}
MEDDLY::node_handle MEDDLY::equant_EV2EV::compute_r(int k, node_handle mxd){
    ostream_output out(std::cout);

    if (mxd == 0 ) return 0;
  //   // check the cache
  node_handle result = 0;

    ct_entry_key* key = findResult(mxd,result);
    if (key == 0) {
      return result;
    }
  //
    const int mxdLevel = argF->getNodeLevel(mxd);
    if(mxdLevel==0){
        if (mxd == 0 ) return 0;
        return resF->linkNode(mxd);
    }
    const int rLevel = ABS(mxdLevel);
    const unsigned rSize = unsigned(resF->getLevelSize(rLevel));
    unpacked_node* C = unpacked_node::newFull(resF, rLevel, rSize);
      // clear out result (important!)
      for (unsigned i=0; i<rSize; i++) C->d_ref(i) = 0;

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node *Ru = unpacked_node::New();
      unpacked_node *Rp = unpacked_node::New();
      if (mxdLevel < 0) {
        Ru->initRedundant(argF, rLevel, mxd, false);
      } else {
        argF->unpackNode(Ru, mxd, SPARSE_ONLY);
      }
      dd_edge newstatesE(resF), cdi(resF);

      // loop over mxd "rows"
      for (unsigned iz=0; iz<Ru->getNNZs(); iz++) {
        unsigned i = Ru->i(iz);
        if (isLevelAbove(-rLevel, argF->getNodeLevel(Ru->d(iz)))) {
          Rp->initIdentity(argF, rLevel, i, Ru->d(iz), false);
        } else {
          argF->unpackNode(Rp, Ru->d(iz), SPARSE_ONLY);
        }
        // loop over mxd "columns"
        for (unsigned jz=0; jz<Rp->getNNZs(); jz++) {
          unsigned j = Rp->i(jz);
          unsigned ix=i;
          // if (0==Rp->d(jz))   continue;
          // ok, there is an i->j "edge".
          // determine new states to be added (recursively)
          // and add them;
          node_handle newstates = compute_r(k-1, Rp->d(jz));
          if (0==newstates) continue;
          if (0==C->d(ix)) {
            C->d_ref(ix) = newstates;
            continue;
          }
           // there's new states and existing states; union them.
          newstatesE.set(newstates);
          cdi.set(C->d(ix));
           mddUnion->computeTemp(newstatesE, cdi, cdi);
          C->set_d(ix, cdi);
        } // for j
      } // for i
      unpacked_node::recycle(Rp);
      unpacked_node::recycle(Ru);
    result = resF->createReducedNode(rLevel, C);
    return saveResult(key, mxd, result);
}

void MEDDLY::equant_EV2EV::compute_r(long aev, node_handle a, int k, long& bev, node_handle& b)
{
  // if ((!resF->isQuasiReduced() || k == 0) && argF->isTerminalNode(a)) {
  //   if (a == 0) {
  //     bev = 0;
  //     b = 0;
  //   }
  //   else {
  //     bev = aev;
  //     b = a;
  //   }
  //   return;
  // }
  //
  // if (resF->isQuasiReduced() && ABS(argF->getNodeLevel(a)) < k) {
  //   int size = resF->getLevelSize(k);
  //   unpacked_node* T = unpacked_node::newFull(resF, k, size);
  //   long tev = Inf<long>();
  //   node_handle t = 0;
  //   compute_r(aev, a, k - 1, tev, t);
  //   for (int i = 0; i < size; i++) {
  //     T->setEdge(i, tev);
  //     T->d_ref(i) = resF->linkNode(t);
  //   }
  //   resF->unlinkNode(t);
  //   resF->createReducedNode(-1, T, bev, b);
  //   return;
  // }
  //
  // // check the cache
  // ct_entry_key* key = findResult(aev, a, bev, b);
  // if (key == 0) {
  //   return;
  // }
  //
  // const int aLevel = argF->getNodeLevel(a);
  // const int level = ABS(aLevel);
  // const int size = resF->getLevelSize(level);
  //
  // unpacked_node* A = aLevel < 0
  //   ? unpacked_node::newRedundant(argF, level, 0L, a, true)
  //   : argF->newUnpacked(a, FULL_ONLY);
  // unpacked_node* T = unpacked_node::newFull(resF, level, size);
  // for (int i = 0; i < size; i++) {
  //   unpacked_node* B = isLevelAbove(-level, argF->getNodeLevel(A->d(i)))
  //     ? (argF->isIdentityReduced()
  //       ? unpacked_node::newIdentity(argF, -level, i, 0L, A->d(i), true)
  //       : unpacked_node::newRedundant(argF, -level, 0L, A->d(i), true))
  //     : argF->newUnpacked(A->d(i), FULL_ONLY);
  //
  //   long tev = Inf<long>();
  //   node_handle t = 0;
  //   compute_r(aev + A->ei(i) + B->ei(i), B->d(i), level - 1, tev, t);
  //   T->setEdge(i, tev);
  //   T->d_ref(i) = t;
  //
  //   unpacked_node::recycle(B);
  // }
  //
  // unpacked_node::recycle(A);
  //
  // resF->createReducedNode(-1, T, bev, b);
  // saveResult(key, aev, a, bev, b);
}

// ******************************************************************
// *                                                                *
// *                         equant_opname                           *
// *                                                                *
// ******************************************************************

class MEDDLY::equant_opname : public unary_opname {
public:
  equant_opname();
  virtual unary_operation* buildOperation(expert_forest* arg, expert_forest* res);
};

MEDDLY::equant_opname::equant_opname()
  : unary_opname("existential quantifier")
{
}

MEDDLY::unary_operation* MEDDLY::equant_opname::buildOperation(
  expert_forest* arg, expert_forest* res)
{
  unary_operation* op = 0;
  if ( arg->isForRelations() && !res->isForRelations()) {
    op = new equant_EV2EV(this, arg, res);
  }
  else {
    throw error(error::NOT_IMPLEMENTED);
  }
  return op;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeEquant()
{
  return new equant_opname;
}
