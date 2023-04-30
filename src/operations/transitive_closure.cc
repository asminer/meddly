
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

#include <deque>
#include <vector>

#include "defines.h"
#include "transitive_closure.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"
#include "../opname_satur.h"

// ******************************************************************
// *                                                                *
// *                       common_constraint                        *
// *                                                                *
// ******************************************************************

MEDDLY::common_transitive_closure::common_transitive_closure(const constrained_opname* code,
  unsigned slots,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : specialized_operation(code, slots)
{
  MEDDLY_DCASSERT(cons->isEVPlus() && !cons->isForRelations());
  MEDDLY_DCASSERT(tc->isEVPlus() && tc->isForRelations());
  MEDDLY_DCASSERT(trans->isMultiTerminal() && trans->isForRelations());
  MEDDLY_DCASSERT(res->isEVPlus() && res->isForRelations());

  consF = cons;
  tcF = tc;
  transF = trans;
  resF = res;

  registerInForest(consF);
  registerInForest(tcF);
  registerInForest(transF);
  registerInForest(resF);
}

MEDDLY::common_transitive_closure::~common_transitive_closure()
{
  unregisterInForest(consF);
  unregisterInForest(tcF);
  unregisterInForest(transF);
  unregisterInForest(resF);
}

bool MEDDLY::common_transitive_closure::checkForestCompatibility() const
{
  auto o1 = consF->variableOrder();
  auto o2 = tcF->variableOrder();
  auto o3 = transF->variableOrder();
  auto o4 = resF->variableOrder();
  return o1->is_compatible_with(*o2)
    && o1->is_compatible_with(*o3)
    && o1->is_compatible_with(*o4);
}

// ******************************************************************
// *                                                                *
// *                transitive_closure_forwd_bfs                    *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_forwd_bfs::transitive_closure_forwd_bfs(const constrained_opname* code,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : common_transitive_closure(code, 0, cons, tc, trans, res)
{
  if (resF->getRangeType() == range_type::INTEGER && resF->isForRelations()) {
    plusOp = getOperation(POST_PLUS, resF, consF, resF);
    minOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(TC_POST_IMAGE, tcF, transF, resF);
}

void MEDDLY::transitive_closure_forwd_bfs::compute(const dd_edge &a, const dd_edge &b, const dd_edge &r, dd_edge &res)
{
  MEDDLY_DCASSERT(consF == a.getForest());
  MEDDLY_DCASSERT(tcF == b.getForest());
  MEDDLY_DCASSERT(transF == r.getForest());
  MEDDLY_DCASSERT(resF == res.getForest());

  /*
  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  iterate(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, cnode);

  res.set(cnode, cev);
  */
  iterate(a, b, r, res);
}

/*
void MEDDLY::transitive_closure_forwd_bfs::iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  cev = bev;
  MEDDLY_DCASSERT(tcF == resF);
  c = resF->linkNode(b);

  node_handle prev = 0;
  while (prev != c) {
    resF->unlinkNode(prev);
    prev = c;
    long tev = Inf<long>();
    node_handle t = 0;
    imageOp->computeTemp(cev, c, r, tev, t);
    // tev was increased by one step
    plusOp->computeTemp(tev - 1, t, aev, a, tev, t);
    minOp->computeTemp(cev, c, tev, t, cev, c);
    resF->unlinkNode(t);
  }
  resF->unlinkNode(prev);
}
*/

void MEDDLY::transitive_closure_forwd_bfs::iterate(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& c)
{
  c = b;
  dd_edge prev(resF), t(resF);
  while (prev != c) {
    prev = c;
    imageOp->computeTemp(c, r, t);

    // t was increased by one step; roll it back one
    long tev;
    t.getEdgeValue(tev);
    tev--;
    t.setEdgeValue(tev);

    plusOp->computeTemp(t, a, t);
    minOp->computeTemp(c, t, c);
  }
}

// ******************************************************************
// *                                                                *
// *                     constraint_dfs_opname                      *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_dfs_opname::transitive_closure_dfs_opname()
 : constrained_opname("Transitive Closure")
{
}

MEDDLY::specialized_operation* MEDDLY::transitive_closure_dfs_opname::buildOperation(arguments* a) const
{
  constrained_opname::constrained_args* args = dynamic_cast<constrained_opname::constrained_args*>(a);
  return new transitive_closure_forwd_dfs(this,
    static_cast<expert_forest*>(args->consForest),
    static_cast<expert_forest*>(args->inForest),
    static_cast<expert_forest*>(args->relForest),
    static_cast<expert_forest*>(args->outForest));
}

// ******************************************************************
// *                                                                *
// *                  transitive_closure_dfs                        *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_dfs::transitive_closure_dfs(const constrained_opname* code,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : common_transitive_closure(code, 1, cons, tc, trans, res)
{
  mxdIntersectionOp = getOperation(INTERSECTION, transF, transF, transF);
  mxdDifferenceOp = getOperation(DIFFERENCE, transF, transF, transF);
  minOp = getOperation(UNION, resF, resF, resF);

  splits = nullptr;

  ct_entry_type* et = new ct_entry_type(code->getName(), "LNNN:LN");
  et->setForestForSlot(1, cons);
  et->setForestForSlot(2, tc);
  et->setForestForSlot(3, trans);
  et->setForestForSlot(6, res);
  registerEntryType(0, et);
  buildCTs();
}

bool MEDDLY::transitive_closure_dfs::checkTerminals(int aev, node_handle a, int bev, node_handle b, node_handle c,
  long& dev, node_handle& d)
{
  if (a == -1 && b == -1 && c == -1) {
    d = -1;
    dev = bev;
    return true;
  }
  if (a == 0 || b == 0 || c == 0) {
    d = 0;
    dev = 0;
    return true;
  }
  return false;
}

MEDDLY::ct_entry_key* MEDDLY::transitive_closure_dfs::findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle c, long& dev, node_handle &d)
{
  ct_entry_key* key = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(key);
  key->writeL(aev);
  key->writeN(a);
  key->writeN(b);
  key->writeN(c);

  CT0->find(key, CTresult[0]);
  if (!CTresult[0]) {
    return key;
  }

  dev = CTresult[0].readL();
  d = resF->linkNode(CTresult[0].readN());

  if (d != 0) {
    dev += bev;
  }
  else {
    MEDDLY_DCASSERT(dev == 0);
  }
  MEDDLY_DCASSERT(dev >= 0);

  CT0->recycle(key);
  return 0;
}

void MEDDLY::transitive_closure_dfs::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, node_handle c, long dev, node_handle d)
{
  CTresult[0].reset();
  if (d == 0) {
    // Write long
    CTresult[0].writeL(0);
  }
  else {
    MEDDLY_DCASSERT(dev - bev >= 0);
    CTresult[0].writeL(dev - bev);
  }
  CTresult[0].writeN(d);
  CT0->addEntry(key, CTresult[0]);
}


// Partition the nsf based on "top level"
void MEDDLY::transitive_closure_dfs::splitMxd(const dd_edge& mxd)
{
  MEDDLY_DCASSERT(transF);
  MEDDLY_DCASSERT(splits == nullptr);

  splits = new dd_edge[transF->getNumVariables() + 1];

  dd_edge root(mxd), maxDiag, Rpdi;
  maxDiag.setForest(transF);
  Rpdi.setForest(transF);

  // Build from top down
  for (int level = transF->getNumVariables(); level > 0; level--) {
    splits[level].setForest(transF);
    if (root.getNode() == 0) {
      // common and easy special case
      continue;
    }

    int mxdLevel = root.getLevel(); // transF->getNodeLevel(mxd);
    MEDDLY_DCASSERT(ABS(mxdLevel) <= level);

    // Initialize readers
    unpacked_node* Ru = isLevelAbove(level, mxdLevel)
      ? unpacked_node::newRedundant(transF, level, root.getNode(), true)
      : transF->newUnpacked(root.getNode(), FULL_ONLY);

    bool first = true;

    // Read "rows"
    for (int i = 0; i < Ru->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = transF->getNodeLevel(Ru->d(i));
      unpacked_node* Rp = isLevelAbove(-level, mxdPLevel)
        ? unpacked_node::newIdentity(transF, -level, i, Ru->d(i), true)
        : transF->newUnpacked(Ru->d(i), FULL_ONLY);

      // Intersect along the diagonal
      if (first) {
        maxDiag.set( transF->linkNode(Rp->d(i)) );
        first = false;
      } else {
        Rpdi.set( transF->linkNode(Rp->d(i)) );
        mxdIntersectionOp->computeTemp(maxDiag, Rpdi, maxDiag);
      }

      // cleanup
      unpacked_node::recycle(Rp);
    } // for i

    // maxDiag is what we can split from here
    mxdDifferenceOp->computeTemp(root, maxDiag, splits[level]);
    root = maxDiag;

    // Cleanup
    unpacked_node::recycle(Ru);
  } // for level

#ifdef DEBUG_SPLIT
  printf("After splitting monolithic event in msat\n");
  printf("splits array: [");
  for (int k = 0; k <= transF->getNumVariables(); k++) {
    if (k) printf(", ");
    printf("%d", splits[k]);
  }
  printf("]\n");
#endif
}

void MEDDLY::transitive_closure_dfs::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(consF == a.getForest());
  MEDDLY_DCASSERT(tcF == b.getForest());
  MEDDLY_DCASSERT(transF == r.getForest());
  MEDDLY_DCASSERT(resF == res.getForest());

  // Partition NSF by levels
  splitMxd(r);

  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle c = 0;
  _compute(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, c);

  res.set(c, cev);
}

void MEDDLY::transitive_closure_dfs::_compute(int aev, node_handle a, int bev, node_handle b, node_handle r,
  long& cev, node_handle& c)
{
  // Execute saturation operation
  transitive_closure_evplus* tcSatOp = new transitive_closure_evplus(this, consF, tcF, resF);
  tcSatOp->saturate(aev, a, bev, b, cev, c);

  // Cleanup
  //tcSatOp->removeAllComputeTableEntries();
  // delete tcSatOp;
  //minOp->removeAllComputeTableEntries();
  //removeAllComputeTableEntries();

  // for (int i = transF->getNumVariables(); i > 0; i--) {
    // transF->unlinkNode(splits[i]);
  // }
  delete[] splits;
  splits = nullptr;
}

// ******************************************************************
// *                                                                *
// *                transitive_closure_forwd_dfs                    *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_forwd_dfs::transitive_closure_forwd_dfs(const constrained_opname* code,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : transitive_closure_dfs(code, cons, tc, trans, res)
{
}

void MEDDLY::transitive_closure_forwd_dfs::saturateHelper(long aev, node_handle a, int in, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);
  // nb is at a primed level
  MEDDLY_DCASSERT(nb.getLevel() < 0);

  const dd_edge& mxd = splits[-nb.getLevel()];
  if (mxd.getNode() == 0) {
    return;
  }

  // const int mxdLevel = transF->getNodeLevel(mxd);
  const int mxdLevel = mxd.getLevel();
  MEDDLY_DCASSERT(ABS(mxdLevel) == -nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(transF, -nb.getLevel(), mxd.getNode(), true)
    : transF->newUnpacked(mxd.getNode(), FULL_ONLY);
  unpacked_node* Rp = unpacked_node::New();

  unpacked_node* A = isLevelAbove(-nb.getLevel(), consF->getNodeLevel(a))
    ? unpacked_node::newRedundant(consF, -nb.getLevel(), 0L, a, true)
    : consF->newUnpacked(a, FULL_ONLY);

  dd_edge nbdj(resF), newst(resF);

  // indices to explore
  int size = nb.getSize();
  std::deque<int> queue;
  std::vector<bool> waiting(size, false);
  for (int ip = 0; ip < size; ip++) {
    if (nb.d(ip) != 0 && Ru->d(ip) != 0) {
      queue.push_back(ip);
      waiting[ip] = true;
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int ip = queue.front();
    queue.pop_front();
    waiting[ip] = false;

    MEDDLY_DCASSERT(nb.d(ip) != 0);
    MEDDLY_DCASSERT(Ru->d(ip) != 0);

    const int dlevel = transF->getNodeLevel(Ru->d(ip));
    if (dlevel == -Ru->getLevel()) {
      transF->unpackNode(Rp, Ru->d(ip), SPARSE_ONLY);
    }
    else {
      Rp->initIdentity(transF, -Ru->getLevel(), ip, Ru->d(ip), false);
    }

    for (int jpz = 0; jpz < Rp->getNNZs(); jpz++) {
      const int jp = Rp->i(jpz);
      if (A->d(jp) == 0) {
        MEDDLY_DCASSERT(A->ei(jp) == 0);
        continue;
      }

      long recev = 0;
      node_handle rec = 0;
      recFire(aev + A->ei(jp), A->d(jp), nb.ei(ip), nb.d(ip), Rp->d(jpz), recev, rec);
      MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

      if (rec == 0) {
        MEDDLY_DCASSERT(recev == 0);
        continue;
      }

      if (rec == nb.d(jp)) {
        // Compute the minimum
        if (recev < nb.ei(jp)) {
          nb.setEdge(jp, recev);
        }
        resF->unlinkNode(rec);
        continue;
      }

      bool updated = true;

      if (nb.d(jp) == 0) {
        MEDDLY_DCASSERT(nb.ei(jp) == 0);
        nb.setEdge(jp, recev);
        nb.d_ref(jp) = rec;
      }
      else {
        nbdj.set(nb.d(jp), nb.ei(jp));
        newst.set(rec, recev);
        minOp->computeTemp(nbdj, newst, nbdj);
        updated = (nbdj.getNode() != nb.d(jp));
        nb.set_de(jp, nbdj);

        // OLD
        /*
        long accev = Inf<long>();
        node_handle acc = 0;
        minOp->computeTemp(nb.ei(jp), nb.d(jp), recev, rec, accev, acc);

        MEDDLY_DCASSERT(acc != 0);
        resF->unlinkNode(rec);
        if (acc != nb.d(jp)) {
          resF->unlinkNode(nb.d(jp));
          nb.setEdge(jp, accev);
          nb.d_ref(jp) = acc;
        }
        else {
          MEDDLY_DCASSERT(accev == nb.ei(jp));
          resF->unlinkNode(acc);
          updated = false;
        }
        */

      }

      if (updated) {
        if (jp == ip) {
          // Restart inner for-loop.
          jpz = -1;
        }
        else {
          if (!waiting[jp] && Ru->d(jp) != 0) {
            MEDDLY_DCASSERT(A->d(jp) != 0);
            queue.push_back(jp);
            waiting[jp] = true;
          }
        }
      }
    }
  }

  unpacked_node::recycle(Ru);
  unpacked_node::recycle(Rp);
  unpacked_node::recycle(A);
}

void MEDDLY::transitive_closure_forwd_dfs::recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  // termination conditions
  if (a == 0 || b == 0 || r == 0) {
    cev = 0;
    c = 0;
    return;
  }
  if (a == -1 && r == -1 && tcF == resF) {
    cev = aev + bev;
    c = resF->linkNode(b);
    return;
  }

  // check the cache
  ct_entry_key* key = findResult(aev, a, bev, b, r, cev, c);
  if (key == 0) {
    MEDDLY_DCASSERT(cev >= 0);
    return;
  }

  // check if mxd and evmdd are at the same level
  const int aLevel = consF->getNodeLevel(a);
  MEDDLY_DCASSERT(aLevel >= 0);
  const int bLevel = tcF->getNodeLevel(b);
  const int rLevel = transF->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), ABS(bLevel)), aLevel);
  const int size = resF->getLevelSize(level);

  dd_edge Tdj(resF), newst(resF);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, 0L, a, true)
    : consF->newUnpacked(a, FULL_ONLY);
  // Initialize evmxd reader
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(tcF, level, 0L, b, true)
    : tcF->newUnpacked(b, FULL_ONLY);
  unpacked_node* D = unpacked_node::New();

  unpacked_node* T = unpacked_node::newFull(resF, level, size);

  for (int i = 0; i < size; i++) {
    if (B->d(i) == 0) {
      T->setEdge(i, 0L);
      T->d_ref(i) = 0;
      continue;
    }

    if (isLevelAbove(-level, tcF->getNodeLevel(B->d(i)))) {
      D->initIdentity(tcF, -level, i, 0L, B->d(i), true);
    }
    else {
      tcF->unpackNode(D, B->d(i), FULL_ONLY);
    }

    unpacked_node* Tp = unpacked_node::newFull(resF, -level, size);

    if (ABS(rLevel) < level) {
      // Assume identity reduction
      MEDDLY_DCASSERT(transF->isIdentityReduced());
      for (int ip = 0; ip < size; ip++) {
        long tev = 0;
        node_handle t = 0;
        recFire(aev + A->ei(ip), A->d(ip), bev + B->ei(i) + D->ei(ip), D->d(ip), r, tev, t);
        Tp->setEdge(ip, tev);
        Tp->d_ref(ip) = t;
      }
    }
    else {
      MEDDLY_DCASSERT(ABS(rLevel) == level);

      //
      // Need to process this level in the MXD.

      // clear out result (important!)
      for (int ip = 0; ip < size; ip++) {
        Tp->setEdge(ip, 0L);
        Tp->d_ref(ip) = 0;
      }

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node* Ru = (rLevel < 0)
        ? unpacked_node::newRedundant(transF, level, r, false)
        : transF->newUnpacked(r, SPARSE_ONLY);
      unpacked_node* Rp = unpacked_node::New();

      // loop over mxd "rows"
      for (int ipz = 0; ipz < Ru->getNNZs(); ipz++) {
        const int ip = Ru->i(ipz);
        if (D->d(ip) == 0) {
          continue;
        }

        if (isLevelAbove(-level, transF->getNodeLevel(Ru->d(ipz)))) {
          Rp->initIdentity(transF, -level, ip, Ru->d(ipz), false);
        }
        else {
          transF->unpackNode(Rp, Ru->d(ipz), SPARSE_ONLY);
        }

        // loop over mxd "columns"
        for (int jpz = 0; jpz < Rp->getNNZs(); jpz++) {
          const int jp = Rp->i(jpz);
          if (A->d(jp) == 0) {
            MEDDLY_DCASSERT(A->ei(jp) == 0);
            continue;
          }

          // ok, there is an ip->jp "edge".
          // determine new states to be added (recursively)
          // and add them
          long nev = 0;
          node_handle n = 0;
          recFire(aev + A->ei(jp), A->d(jp), bev + B->ei(i) + D->ei(ip), D->d(ip), Rp->d(jpz), nev, n);

          if (n == 0) {
            MEDDLY_DCASSERT(nev == 0);
            continue;
          }

          if (Tp->d(jp) == 0) {
            MEDDLY_DCASSERT(Tp->ei(jp) == 0);
            Tp->setEdge(jp, nev);
            Tp->d_ref(jp) = n;
            continue;
          }

          // there's new states and existing states; union them.
          /*
          const node_handle oldjp = Tp->d(jp);
          long newev = Inf<long>();
          node_handle newstates = 0;
          minOp->computeTemp(nev, n, Tp->ei(jp), oldjp, newev, newstates);
          Tp->setEdge(jp, newev);
          Tp->d_ref(jp) = newstates;

          resF->unlinkNode(oldjp);
          resF->unlinkNode(n);
          */

          // there's new states and existing states; union them.
          Tdj.set(Tp->d(jp), Tp->ei(jp));
          newst.set(n, nev);
          minOp->computeTemp(newst, Tdj, Tdj);
          Tp->set_de(jp, Tdj);
        } // for j
      } // for i

      unpacked_node::recycle(Ru);
      unpacked_node::recycle(Rp);
    }

    saturateHelper(aev, a, i, *Tp);

    long tpev = Inf<long>();
    node_handle tp = 0;
    resF->createReducedNode(i, Tp, tpev, tp);
    T->setEdge(i, tpev);
    T->d_ref(i) = tp;
  }

  // cleanup mdd reader
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);
  unpacked_node::recycle(D);

  resF->createReducedNode(-1, T, cev, c);
  MEDDLY_DCASSERT(cev >= 0);

  saveResult(key, aev, a, bev, b, r, cev, c);
}

// ******************************************************************
// *                                                                *
// *                  transitive_closure_evplus                     *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_evplus::transitive_closure_evplus(transitive_closure_dfs* p,
  expert_forest* cons, expert_forest* tc, expert_forest* res)
  : specialized_operation(nullptr, 1)
{
  MEDDLY_DCASSERT(cons->isEVPlus() && !cons->isForRelations());
  MEDDLY_DCASSERT(tc->isEVPlus() && tc->isForRelations());
  MEDDLY_DCASSERT(res->isEVPlus() && res->isForRelations());

  parent = p;
  consF = cons;
  tcF = tc;
  resF = res;

  registerInForest(consF);
  registerInForest(tcF);
  registerInForest(resF);

  ct_entry_type* et;

  if (tcF->isFullyReduced() || tcF->isIdentityReduced()) {
    // CT entry includes level info
    et = new ct_entry_type("transitive_closure_evplus", "LNNI:LN");
    et->setForestForSlot(1, cons);
    et->setForestForSlot(2, tc);
    et->setForestForSlot(6, res);
  } else {
    et = new ct_entry_type("transitive_closure_evplus", "LNN:LN");
    et->setForestForSlot(1, cons);
    et->setForestForSlot(2, tc);
    et->setForestForSlot(5, res);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::transitive_closure_evplus::~transitive_closure_evplus()
{
  unregisterInForest(consF);
  unregisterInForest(tcF);
  unregisterInForest(resF);
}

bool MEDDLY::transitive_closure_evplus::checkForestCompatibility() const
{
  auto o1 = consF->variableOrder();
  auto o2 = tcF->variableOrder();
  auto o3 = resF->variableOrder();
  return o1->is_compatible_with(*o2)
    && o1->is_compatible_with(*o3);
}

bool MEDDLY::transitive_closure_evplus::checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c)
{
  if (a == -1 && b == -1) {
    c = -1;
    cev = bev;
    return true;
  }
  if (a == 0 || b == 0) {
    c = 0;
    cev = 0;
    return true;
  }
  return false;
}

MEDDLY::ct_entry_key* MEDDLY::transitive_closure_evplus::findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c)
{
  ct_entry_key* key = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(key);
  key->writeL(aev);
  key->writeN(a);
  key->writeN(b);
  if(tcF->isFullyReduced() || tcF->isIdentityReduced()) {
    // Level is part of key for fully-reduced or identity-reduced forest
    key->writeI(level);
  }

  CT0->find(key, CTresult[0]);
  if (!CTresult[0]) return key;

  cev = CTresult[0].readL();
  c = resF->linkNode(CTresult[0].readN());

  if (c != 0) {
    cev += bev;
  }
  else {
    MEDDLY_DCASSERT(cev == 0);
  }
  CT0->recycle(key);
  return 0;
}

void MEDDLY::transitive_closure_evplus::saveResult(ct_entry_key* key,
  long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c)
{
  CTresult[0].reset();
  if (c == 0) {
    // Write long
    CTresult[0].writeL(0);
  }
  else {
    MEDDLY_DCASSERT(cev - bev >= 0);
    CTresult[0].writeL(cev - bev);
  }
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

void MEDDLY::transitive_closure_evplus::saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c)
{
  saturate(aev, a, bev, b, tcF->getNumVariables(), cev, c);
}

void MEDDLY::transitive_closure_evplus::saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c)
{
  MEDDLY_DCASSERT(level >= 0);

  if (checkTerminals(aev, a, bev, b, cev, c)) {
    return;
  }

  ct_entry_key* key = findResult(aev, a, bev, b, level, cev, c);
  if (key == 0) {
    return;
  }

  const int sz = tcF->getLevelSize(level);
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = tcF->getNodeLevel(b);

  MEDDLY_DCASSERT(aLevel >= 0);

  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, 0L, a, true)
    : consF->newUnpacked(a, FULL_ONLY);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(tcF, level, 0L, b, true)
    : tcF->newUnpacked(b, FULL_ONLY);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
    if (B->d(i) == 0) {
      MEDDLY_DCASSERT(resF == tcF);
      T->setEdge(i, 0L);
      T->d_ref(i) = 0;
    }
    else {
      unpacked_node* D = isLevelAbove(-level, tcF->getNodeLevel(B->d(i)))
        ? unpacked_node::newIdentity(tcF, -level, i, 0L, B->d(i), true)
        : tcF->newUnpacked(B->d(i), FULL_ONLY);
      unpacked_node* Tp = unpacked_node::newFull(resF, -level, sz);
      for (int j = 0; j < sz; j++) {
        if (A->d(j) == 0 || D->d(j) == 0) {
          Tp->setEdge(j, 0L);
          Tp->d_ref(j) = 0;
        }
        else {
          long tpev = Inf<long>();
          node_handle tp = 0;
          saturate(aev + A->ei(j), A->d(j), bev + B->ei(i) + D->ei(j), D->d(j), level - 1, tpev, tp);
          Tp->setEdge(j, tpev);
          Tp->d_ref(j) = tp;
        }
      }
      unpacked_node::recycle(D);

      parent->saturateHelper(aev, a, i, *Tp);

      long tev = Inf<long>();
      node_handle t = 0;
      resF->createReducedNode(i, Tp, tev, t);
      T->setEdge(i, tev);
      T->d_ref(i) = t;
    }
  }

  // Cleanup
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);

  resF->createReducedNode(-1, T, cev, c);

  // save in compute table
  saveResult(key, aev, a, bev, b, level, cev, c);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_opname* MEDDLY::initTransitiveClosureDFS()
{
  return new transitive_closure_dfs_opname();
}
