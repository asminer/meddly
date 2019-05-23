
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
#include "constrained.h"

// ******************************************************************
// *                                                                *
// *                       common_constraint                        *
// *                                                                *
// ******************************************************************

MEDDLY::common_constrained::common_constrained(const constrained_opname* code,
  unsigned slots,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : specialized_operation(code, slots)
{
  consF = cons;
  argF = arg;
  transF = trans;
  resF = res;

  registerInForest(consF);
  registerInForest(argF);
  registerInForest(transF);
  registerInForest(resF);
}

MEDDLY::common_constrained::~common_constrained()
{
  unregisterInForest(consF);
  unregisterInForest(argF);
  unregisterInForest(transF);
  unregisterInForest(resF);
}

bool MEDDLY::common_constrained::checkForestCompatibility() const
{
  auto o1 = consF->variableOrder();
  auto o2 = argF->variableOrder();
  auto o3 = transF->variableOrder();
  auto o4 = resF->variableOrder();
  return o1->is_compatible_with(*o2)
    && o1->is_compatible_with(*o3)
    && o1->is_compatible_with(*o4);
}

// ******************************************************************
// *                                                                *
// *                     constraint_bfs_opname                      *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_bfs_opname::constrained_bfs_opname(bool fwd)
 : constrained_opname(fwd ? "PostImage" : "PreImage")
{
  forward = fwd;
}

MEDDLY::specialized_operation* MEDDLY::constrained_bfs_opname::buildOperation(arguments* a) const
{
  constrained_opname::constrained_args* args = dynamic_cast<constrained_opname::constrained_args*>(a);
  specialized_operation* op = 0;
  if (forward) {
    throw error(error::NOT_IMPLEMENTED);
  }
  else {
    op = new constrained_bckwd_bfs_evplus(this,
      static_cast<expert_forest*>(args->consForest),
      static_cast<expert_forest*>(args->inForest),
      static_cast<expert_forest*>(args->relForest),
      static_cast<expert_forest*>(args->outForest));
  }
  return op;
}

// ******************************************************************
// *                                                                *
// *                constraint_bckwd_bfs_evplus                     *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_bckwd_bfs_evplus::constrained_bckwd_bfs_evplus(const constrained_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : common_constrained(code, 0, cons, arg, trans, res)
{
  if (resF->getRangeType() == forest::INTEGER) {
    plusOp = getOperation(PLUS, resF, consF, resF);
    minOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(PRE_IMAGE, argF, transF, resF);
}

void MEDDLY::constrained_bckwd_bfs_evplus::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(res.getForest() == resF);

  if (resF->getRangeType() == forest::INTEGER) {
    plusOp = getOperation(PLUS, resF, consF, resF);
    minOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(PRE_IMAGE, argF, transF, resF);

  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  iterate(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, cnode);

  res.set(cnode, cev);
}

void MEDDLY::constrained_bckwd_bfs_evplus::iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  cev = bev;
  MEDDLY_DCASSERT(argF == resF);
  c = resF->linkNode(b);

  node_handle prev = 0;
  while (prev != c) {
    resF->unlinkNode(prev);
    prev = c;
    long tev = Inf<long>();
    node_handle t = 0;
    imageOp->compute(cev, c, r, tev, t);
    node_handle oldt = t;
    plusOp->compute(tev, oldt, aev, a, tev, t);
    minOp->compute(cev, c, tev, t, cev, c);
    resF->unlinkNode(oldt);
    resF->unlinkNode(t);
  }
  resF->unlinkNode(prev);
}

// ******************************************************************
// *                                                                *
// *                     constraint_dfs_opname                      *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_dfs_opname::constrained_dfs_opname(bool fwd)
 : constrained_opname(fwd ? "PostImage" : "PreImage")
{
  forward = fwd;
}

MEDDLY::specialized_operation* MEDDLY::constrained_dfs_opname::buildOperation(arguments* a) const
{
  constrained_opname::constrained_args* args = dynamic_cast<constrained_opname::constrained_args*>(a);
  specialized_operation* op = 0;
  if (forward) {
    if (args->consForest->isMultiTerminal() && args->inForest->isMultiTerminal()
      && args->relForest->isMultiTerminal() && args->outForest->isMultiTerminal()) {
      op = new constrained_forwd_dfs_mt(this,
        static_cast<expert_forest*>(args->consForest),
        static_cast<expert_forest*>(args->inForest),
        static_cast<expert_forest*>(args->relForest),
        static_cast<expert_forest*>(args->outForest));
    }
    else {
      throw error(error::NOT_IMPLEMENTED);
    }
  }
  else {
    if (args->consForest->isMultiTerminal() && args->inForest->isMultiTerminal()
      && args->relForest->isMultiTerminal() && args->outForest->isMultiTerminal()) {
      op = new constrained_bckwd_dfs_mt(this,
        static_cast<expert_forest*>(args->consForest),
        static_cast<expert_forest*>(args->inForest),
        static_cast<expert_forest*>(args->relForest),
        static_cast<expert_forest*>(args->outForest));
    }
    else if (args->consForest->isEVPlus() && args->inForest->isEVPlus()
      && args->relForest->isMultiTerminal() && args->outForest->isEVPlus()) {
      op = new constrained_bckwd_dfs_evplus(this,
        static_cast<expert_forest*>(args->consForest),
        static_cast<expert_forest*>(args->inForest),
        static_cast<expert_forest*>(args->relForest),
        static_cast<expert_forest*>(args->outForest));
    }
    else {
      throw error(error::NOT_IMPLEMENTED);
    }
  }
  return op;
}


// ******************************************************************
// *                                                                *
// *                  constraint_bckwd_dfs_mt                       *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_dfs_mt::constrained_dfs_mt(const constrained_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : common_constrained(code, 1, cons, arg, trans, res)
{
  MEDDLY_DCASSERT(cons->isMultiTerminal() && !cons->isForRelations());
  MEDDLY_DCASSERT(arg->isMultiTerminal() && !arg->isForRelations());
  MEDDLY_DCASSERT(trans->isMultiTerminal() && trans->isForRelations());
  MEDDLY_DCASSERT(res->isMultiTerminal() && !res->isForRelations());

  mxdIntersectionOp = getOperation(INTERSECTION, transF, transF, transF);
  mxdDifferenceOp = getOperation(DIFFERENCE, transF, transF, transF);
  unionOp = getOperation(UNION, resF, resF, resF);

  splits = nullptr;

  compute_table::entry_type* et = new compute_table::entry_type(code->getName(), "NNN:N");
  et->setForestForSlot(0, cons);
  et->setForestForSlot(1, arg);
  et->setForestForSlot(2, trans);
  et->setForestForSlot(4, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::compute_table::entry_key* MEDDLY::constrained_dfs_mt::findResult(
    node_handle a, node_handle b, node_handle r, node_handle &c)
{
  compute_table::entry_key* key = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(key);
  key->writeN(a);
  key->writeN(b);
  key->writeN(r);

  CT0->find(key, CTresult[0]);
  if (!CTresult[0]) {
    return key;
  }

  c = resF->linkNode(CTresult[0].readN());

  CT0->recycle(key);
  return 0;
}

void MEDDLY::constrained_dfs_mt::saveResult(compute_table::entry_key* key,
  node_handle a, node_handle b, node_handle r, node_handle c)
{
  CTresult[0].reset();
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}


// Partition the nsf based on "top level"
void MEDDLY::constrained_dfs_mt::splitMxd(node_handle mxd)
{
  MEDDLY_DCASSERT(transF);
  MEDDLY_DCASSERT(splits == nullptr);

  splits = new node_handle[transF->getNumVariables() + 1];
  splits[0] = 0;

  // we'll be unlinking later, so...
  transF->linkNode(mxd);

  // Build from top down
  for (int level = transF->getNumVariables(); level > 0; level--) {
    if (mxd == 0) {
      // common and easy special case
      splits[level] = 0;
      continue;
    }

    int mxdLevel = transF->getNodeLevel(mxd);
    MEDDLY_DCASSERT(ABS(mxdLevel) <= level);

    // Initialize readers
    unpacked_node* Ru = isLevelAbove(level, mxdLevel)
      ? unpacked_node::newRedundant(transF, level, mxd, true)
      : unpacked_node::newFromNode(transF, mxd, true);

    bool first = true;
    node_handle maxDiag;

    // Read "rows"
    for (int i = 0; i < Ru->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = transF->getNodeLevel(Ru->d(i));
      unpacked_node* Rp = isLevelAbove(-level, mxdPLevel)
        ? unpacked_node::newIdentity(transF, -level, i, Ru->d(i), true)
        : unpacked_node::newFromNode(transF, Ru->d(i), true);

      // Intersect along the diagonal
      if (first) {
        maxDiag = transF->linkNode(Rp->d(i));
        first = false;
      } else {
        node_handle nmd = mxdIntersectionOp->compute(maxDiag, Rp->d(i));
        transF->unlinkNode(maxDiag);
        maxDiag = nmd;
      }

      // cleanup
      unpacked_node::recycle(Rp);
    } // for i

    // maxDiag is what we can split from here
    splits[level] = mxdDifferenceOp->compute(mxd, maxDiag);
    transF->unlinkNode(mxd);
    mxd = maxDiag;

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

void MEDDLY::constrained_dfs_mt::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(res.getForest() == resF);

  node_handle c = 0;
  if (a.getNode() == 0) {
    c = smart_cast<expert_forest*>(resF)->linkNode(b.getNode());
  } else {
    _compute(a.getNode(), b.getNode(), r.getNode(), c);
  }
  res.set(c);
}

void MEDDLY::constrained_dfs_mt::_compute(node_handle a, node_handle b, node_handle r,
  node_handle& c)
{
  // Partition NSF by levels
  splitMxd(r);

  // Execute saturation operation
  constrained_saturation_mt* bckwdSatOp = new constrained_saturation_mt(this, consF, argF, resF);
  bckwdSatOp->saturate(a, b, c);

  // Cleanup
//  bckwdSatOp->removeAllComputeTableEntries();
  //delete bckwdSatOp;
//  removeAllComputeTableEntries();
  for (int i = transF->getNumVariables(); i > 0; i--) {
    transF->unlinkNode(splits[i]);
  }
  delete[] splits;
  splits = nullptr;
}

// ******************************************************************
// *                                                                *
// *                   constrained_forwd_dfs_mt                     *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_forwd_dfs_mt::constrained_forwd_dfs_mt(const constrained_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : constrained_dfs_mt(code, cons, arg, trans, res)
{
}

void MEDDLY::constrained_forwd_dfs_mt::saturateHelper(node_handle a, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);

  node_handle mxd = splits[nb.getLevel()];
  if (mxd == 0) {
    return;
  }

  const int mxdLevel = transF->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(transF, nb.getLevel(), mxd, true)
    : unpacked_node::newFromNode(transF, mxd, true);

  unpacked_node* A = isLevelAbove(nb.getLevel(), consF->getNodeLevel(a))
    ? unpacked_node::newRedundant(consF, nb.getLevel(), a, true)
    : unpacked_node::newFromNode(consF, a, true);

  // indices to explore
  std::deque<int> queue;
  std::vector<bool> waiting(nb.getSize(), false);
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.d(i) != 0 && Ru->d(i) != 0) {
      queue.push_back(i);
      waiting[i] = true;
    }
  }

  unpacked_node** Rps = new unpacked_node*[Ru->getSize()];
  for (int i = 0; i < Ru->getSize(); i++) {
    if (Ru->d(i) == 0) {
      Rps[i] = nullptr;
    }
    else {
      Rps[i] = (transF->getNodeLevel(Ru->d(i)) == -nb.getLevel())
        ? unpacked_node::newFromNode(transF, Ru->d(i), false)
        : unpacked_node::newIdentity(transF, -nb.getLevel(), i, Ru->d(i), false);
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int i = queue.front();
    queue.pop_front();
    waiting[i] = false;

    MEDDLY_DCASSERT(nb.d(i) != 0 && Rps[i] != nullptr);

    for (int jz = 0; jz < Rps[i]->getNNZs(); jz++) {
      const int j = Rps[i]->i(jz);
      if (A->d(j) == 0) {
        continue;
      }

      if (Rps[i]->d(jz) != 0) {
        node_handle rec = 0;
        recFire(A->d(j), nb.d(i), Rps[i]->d(jz), rec);
        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == 0) {
          continue;
        }

        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == nb.d(j)) {
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (nb.d(j) == 0) {
          nb.d_ref(j) = rec;
        }
        else {
          node_handle acc = unionOp->compute(nb.d(j), rec);
          resF->unlinkNode(rec);
          if (acc != nb.d(j)) {
            resF->unlinkNode(nb.d(j));
            nb.d_ref(j) = acc;
          }
          else {
            resF->unlinkNode(acc);
            updated = false;
          }
        }

        if (updated) {
          if (j == i) {
            // Restart inner for-loop.
            jz = -1;
          }
          else {
            if (!waiting[j] && Rps[j] != nullptr) {
              queue.push_back(j);
              waiting[j] = true;
            }
          }
        }
      }
    }
  }

  unpacked_node::recycle(Ru);
  for (int i = 0; i < Ru->getSize(); i++) {
    if (Rps[i] != nullptr) {
      unpacked_node::recycle(Rps[i]);
    }
  }
  delete[] Rps;
  unpacked_node::recycle(A);
}

void MEDDLY::constrained_forwd_dfs_mt::recFire(node_handle a, node_handle b, node_handle r, node_handle& c)
{
  // termination conditions
  if (a == 0 || b == 0 || r == 0) {
    c = 0;
    return;
  }
  if ((a == -1 || a == b) && r == -1) {
    c = resF->linkNode(b);
    return;
  }

  // check the cache
  compute_table::entry_key* key = findResult(a, b, r, c);
  if (key == 0) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);
  const int rLevel = transF->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), bLevel), aLevel);
  const int size = resF->getLevelSize(level);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(argF, level, b, true)
    : unpacked_node::newFromNode(argF, b, true);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  if (ABS(rLevel) < level) {
    // Assume identity reduction
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i = 0; i < size; i++) {
      node_handle t = 0;
      recFire(A->d(i), B->d(i), r, t);
      T->d_ref(i) = t;
    }
  }
  else {
    //
    // Need to process this level in the MXD.

    // clear out result (important!)
    for (int i = 0; i < size; i++) {
      T->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node* Ru = (rLevel < 0)
      ? unpacked_node::newRedundant(transF, -rLevel, r, false)
      : unpacked_node::newFromNode(transF, r, false);

    // loop over mxd "rows"
    for (int iz = 0; iz < Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);

      unpacked_node* Rp = isLevelAbove(-level, transF->getNodeLevel(Ru->d(iz)))
        ? unpacked_node::newIdentity(transF, -level, i, Ru->d(iz), false)
        : unpacked_node::newFromNode(transF, Ru->d(iz), false);

      // loop over mxd "columns"
      for (int jz = 0; jz < Rp->getNNZs(); jz++) {
        const int j = Rp->i(jz);
        if (A->d(j) == 0) {
          continue;
        }

        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle n = 0;
        recFire(A->d(j), B->d(i), Rp->d(jz), n);

        if (n == 0) {
          continue;
        }

        if (T->d(j) == 0) {
          T->d_ref(j) = n;
          continue;
        }

        // there's new states and existing states; union them.
        const node_handle oldj = T->d(j);
        node_handle newstates = unionOp->compute(n, oldj);
        T->d_ref(j) = newstates;

        resF->unlinkNode(oldj);
        resF->unlinkNode(n);
      } // for j

      unpacked_node::recycle(Rp);
    } // for i

    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);

  saturateHelper(a, *T);
  c = resF->createReducedNode(-1, T);

  saveResult(key, a, b, r, c);
}

// ******************************************************************
// *                                                                *
// *                   constrained_bckwd_dfs_mt                     *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_bckwd_dfs_mt::constrained_bckwd_dfs_mt(const constrained_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : constrained_dfs_mt(code, cons, arg, trans, res)
{
}

void MEDDLY::constrained_bckwd_dfs_mt::saturateHelper(node_handle a, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);

  node_handle mxd = splits[nb.getLevel()];
  if (mxd == 0) {
    return;
  }

  const int mxdLevel = transF->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(transF, nb.getLevel(), mxd, false)
    : unpacked_node::newFromNode(transF, mxd, false);

  unpacked_node* A = isLevelAbove(nb.getLevel(), consF->getNodeLevel(a))
    ? unpacked_node::newRedundant(consF, nb.getLevel(), a, true)
    : unpacked_node::newFromNode(consF, a, true);

  // indices to explore
  std::deque<int> queue;
  std::vector<bool> waiting(nb.getSize(), false);
  for (int j = 0; j < nb.getSize(); j++) {
    if (nb.d(j) != 0) {
      queue.push_back(j);
      waiting[j] = true;
    }
  }

  unpacked_node** Rps = new unpacked_node*[Ru->getNNZs()];
  for (int iz = 0; iz < Ru->getNNZs(); iz++) {
    const int i = Ru->i(iz);
    if (A->d(i) == 0) {
      Rps[iz] = nullptr;
    }
    else {
      Rps[iz] = (transF->getNodeLevel(Ru->d(iz)) == -nb.getLevel())
        ? unpacked_node::newFromNode(transF, Ru->d(iz), true)
        : unpacked_node::newIdentity(transF, -nb.getLevel(), i, Ru->d(iz), true);
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int j = queue.front();
    queue.pop_front();
    waiting[j] = false;

    MEDDLY_DCASSERT(nb.d(j) != 0);

    for (int iz = 0; iz < Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      if (A->d(i) == 0) {
        continue;
      }

      MEDDLY_DCASSERT(Rps[iz] != nullptr);

      if (Rps[iz]->d(j) != 0) {
        node_handle rec = 0;
        recFire(A->d(i), nb.d(j), Rps[iz]->d(j), rec);
        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == 0) {
          continue;
        }

        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == nb.d(i)) {
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (nb.d(i) == 0) {
          nb.d_ref(i) = rec;
        }
        else {
          node_handle acc = unionOp->compute(nb.d(i), rec);
          resF->unlinkNode(rec);
          if (acc != nb.d(i)) {
            resF->unlinkNode(nb.d(i));
            nb.d_ref(i) = acc;
          }
          else {
            resF->unlinkNode(acc);
            updated = false;
          }
        }

        if (updated) {
          if (j == i) {
            // Restart inner for-loop.
            iz = -1;
          }
          else {
            if (!waiting[i]) {
              queue.push_back(i);
              waiting[i] = true;
            }
          }
        }
      }
    }
  }

  unpacked_node::recycle(Ru);
  for (int iz = 0; iz < Ru->getNNZs(); iz++) {
    unpacked_node::recycle(Rps[iz]);
  }
  delete[] Rps;
  unpacked_node::recycle(A);
}

void MEDDLY::constrained_bckwd_dfs_mt::recFire(node_handle a, node_handle b, node_handle r, node_handle& c)
{
  // termination conditions
  if (a == 0 || b == 0 || r == 0) {
    c = 0;
    return;
  }
  if ((a == -1 || a == b) && r == -1) {
    c = resF->linkNode(b);
    return;
  }

  // check the cache
  compute_table::entry_key* key = findResult(a, b, r, c);
  if (key == 0) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);
  const int rLevel = transF->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), bLevel), aLevel);
  const int size = resF->getLevelSize(level);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(argF, level, b, true)
    : unpacked_node::newFromNode(argF, b, true);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  if (ABS(rLevel) < level) {
    // Assume identity reduction
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i = 0; i < size; i++) {
      node_handle t = 0;
      recFire(A->d(i), B->d(i), r, t);
      T->d_ref(i) = t;
    }
  }
  else {
    //
    // Need to process this level in the MXD.

    // clear out result (important!)
    for (int i = 0; i < size; i++) {
      T->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node* Ru = (rLevel < 0)
      ? unpacked_node::newRedundant(transF, -rLevel, r, false)
      : unpacked_node::newFromNode(transF, r, false);

    // loop over mxd "rows"
    for (int iz = 0; iz < Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      if (A->d(i) == 0) {
        continue;
      }

      unpacked_node* Rp = isLevelAbove(-level, transF->getNodeLevel(Ru->d(iz)))
        ? unpacked_node::newIdentity(transF, -level, i, Ru->d(iz), false)
        : unpacked_node::newFromNode(transF, Ru->d(iz), false);

      // loop over mxd "columns"
      for (int jz = 0; jz < Rp->getNNZs(); jz++) {
        const int j = Rp->i(jz);

        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle n = 0;
        recFire(A->d(i), B->d(j), Rp->d(jz), n);

        if (n == 0) {
          continue;
        }

        if (T->d(i) == 0) {
          T->d_ref(i) = n;
          continue;
        }

        // there's new states and existing states; union them.
        const node_handle oldi = T->d(i);
        node_handle newstates = unionOp->compute(n, oldi);
        T->d_ref(i) = newstates;

        resF->unlinkNode(oldi);
        resF->unlinkNode(n);
      } // for j

      unpacked_node::recycle(Rp);
    } // for i

    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);

  saturateHelper(a, *T);
  c = resF->createReducedNode(-1, T);

  saveResult(key, a, b, r, c);
}

// ******************************************************************
// *                                                                *
// *                  constrained_saturation_mt                     *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_saturation_mt::constrained_saturation_mt(constrained_dfs_mt* p,
  expert_forest* cons, expert_forest* arg, expert_forest* res)
  : specialized_operation(nullptr, 1)
{
  MEDDLY_DCASSERT(cons->isMultiTerminal() && !cons->isForRelations());
  MEDDLY_DCASSERT(arg->isMultiTerminal() && !arg->isForRelations());
  MEDDLY_DCASSERT(res->isMultiTerminal() && !res->isForRelations());

  parent = p;
  consF = cons;
  argF = arg;
  resF = res;

  registerInForest(consF);
  registerInForest(argF);
  registerInForest(resF);

  compute_table::entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new compute_table::entry_type("constrained_saturation_mt", "NNI:N");
    et->setForestForSlot(0, cons);
    et->setForestForSlot(1, arg);
    et->setForestForSlot(4, res);
  } else {
    et = new compute_table::entry_type("constrained_saturation_mt", "NN:N");
    et->setForestForSlot(0, cons);
    et->setForestForSlot(1, arg);
    et->setForestForSlot(3, res);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::constrained_saturation_mt::~constrained_saturation_mt()
{
  unregisterInForest(consF);
  unregisterInForest(argF);
  unregisterInForest(resF);
}

bool MEDDLY::constrained_saturation_mt::checkForestCompatibility() const
{
  auto o1 = consF->variableOrder();
  auto o2 = argF->variableOrder();
  auto o3 = resF->variableOrder();
  return o1->is_compatible_with(*o2)
    && o1->is_compatible_with(*o3);
}

bool MEDDLY::constrained_saturation_mt::checkTerminals(node_handle a, node_handle b, node_handle& c)
{
  if (argF->isTerminalNode(b)) {
    c = b;
    return true;
  }
  return false;
}

MEDDLY::compute_table::entry_key* MEDDLY::constrained_saturation_mt::findResult(
    node_handle a, node_handle b, int level, node_handle &c)
{
  compute_table::entry_key* key = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(key);
  key->writeN(a);
  key->writeN(b);
  if(argF->isFullyReduced()) {
    // Level is part of key for fully-reduced forest
    key->writeI(level);
  }

  CT0->find(key, CTresult[0]);
  if (!CTresult[0]) return key;
  c = resF->linkNode(CTresult[0].readN());

  CT0->recycle(key);
  return 0;
}

void MEDDLY::constrained_saturation_mt::saveResult(compute_table::entry_key* key,
  node_handle a, node_handle b, int level, node_handle c)
{
  CTresult[0].reset();
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

void MEDDLY::constrained_saturation_mt::saturate(node_handle a, node_handle b, node_handle& c)
{
  saturate(a, b, argF->getNumVariables(), c);
}

void MEDDLY::constrained_saturation_mt::saturate(node_handle a, node_handle b, int level, node_handle& c)
{
  if (checkTerminals(a, b, c)) {
    return;
  }

  compute_table::entry_key* key = findResult(a, b, level, c);
  if (key == 0) {
    return;
  }

  const int sz = argF->getLevelSize(level);
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);

  unpacked_node* A = (aLevel < level)
    ? unpacked_node::newRedundant(consF, level, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = (bLevel < level)
    ? unpacked_node::newRedundant(argF, level, b, true)
    : unpacked_node::newFromNode(argF, b, true);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
    if (A->d(i) == 0) {
      MEDDLY_DCASSERT(resF == argF);
      T->d_ref(i) = resF->linkNode(B->d(i));
    }
    else {
      node_handle t = 0;
      saturate(A->d(i), B->d(i), level - 1, t);
      T->d_ref(i) = t;
    }
  }

  // Cleanup
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);

  parent->saturateHelper(a, *T);
  c = resF->createReducedNode(-1, T);

  // save in compute table
  saveResult(key, a, b, level, c);
}

// ******************************************************************
// *                                                                *
// *                constraint_bckwd_dfs_evplus                     *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_bckwd_dfs_evplus::constrained_bckwd_dfs_evplus(const constrained_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : common_constrained(code, 1, cons, arg, trans, res)
{
  MEDDLY_DCASSERT(cons->isEVPlus() && !cons->isForRelations());
  MEDDLY_DCASSERT(arg->isEVPlus() && !arg->isForRelations());
  MEDDLY_DCASSERT(trans->isMultiTerminal() && trans->isForRelations());
  MEDDLY_DCASSERT(res->isEVPlus() && !res->isForRelations());

  mxdIntersectionOp = getOperation(INTERSECTION, transF, transF, transF);
  mxdDifferenceOp = getOperation(DIFFERENCE, transF, transF, transF);
  minOp = getOperation(UNION, resF, resF, resF);

  splits = nullptr;

  compute_table::entry_type* et = new compute_table::entry_type(code->getName(), "LNNN:LN");
  et->setForestForSlot(1, cons);
  et->setForestForSlot(2, arg);
  et->setForestForSlot(3, trans);
  et->setForestForSlot(6, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::compute_table::entry_key* MEDDLY::constrained_bckwd_dfs_evplus::findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle r, long& cev, node_handle &c)
{
  compute_table::entry_key* key = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(key);
  key->writeL(aev);
  key->writeN(a);
  key->writeN(b);
  key->writeN(r);

  CT0->find(key, CTresult[0]);
  if (!CTresult[0]) {
    return key;
  }

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

void MEDDLY::constrained_bckwd_dfs_evplus::saveResult(compute_table::entry_key* key,
  long aev, node_handle a, long bev, node_handle b, node_handle r, long cev, node_handle c)
{
  CTresult[0].reset();
  if (c == 0) {
    // Write long
    CTresult[0].writeL(0);
  }
  else {
    CTresult[0].writeL(cev - bev);
  }
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}


// Partition the nsf based on "top level"
void MEDDLY::constrained_bckwd_dfs_evplus::splitMxd(node_handle mxd)
{
  MEDDLY_DCASSERT(transF);
  MEDDLY_DCASSERT(splits == nullptr);

  splits = new node_handle[transF->getNumVariables() + 1];
  splits[0] = 0;

  // we'll be unlinking later, so...
  transF->linkNode(mxd);

  // Build from top down
  for (int level = transF->getNumVariables(); level > 0; level--) {
    if (mxd == 0) {
      // common and easy special case
      splits[level] = 0;
      continue;
    }

    int mxdLevel = transF->getNodeLevel(mxd);
    MEDDLY_DCASSERT(ABS(mxdLevel) <= level);

    // Initialize readers
    unpacked_node* Ru = isLevelAbove(level, mxdLevel)
      ? unpacked_node::newRedundant(transF, level, mxd, true)
      : unpacked_node::newFromNode(transF, mxd, true);

    bool first = true;
    node_handle maxDiag;

    // Read "rows"
    for (int i = 0; i < Ru->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = transF->getNodeLevel(Ru->d(i));
      unpacked_node* Rp = isLevelAbove(-level, mxdPLevel)
        ? unpacked_node::newIdentity(transF, -level, i, Ru->d(i), true)
        : unpacked_node::newFromNode(transF, Ru->d(i), true);

      // Intersect along the diagonal
      if (first) {
        maxDiag = transF->linkNode(Rp->d(i));
        first = false;
      } else {
        node_handle nmd = mxdIntersectionOp->compute(maxDiag, Rp->d(i));
        transF->unlinkNode(maxDiag);
        maxDiag = nmd;
      }

      // cleanup
      unpacked_node::recycle(Rp);
    } // for i

    // maxDiag is what we can split from here
    splits[level] = mxdDifferenceOp->compute(mxd, maxDiag);
    transF->unlinkNode(mxd);
    mxd = maxDiag;

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

void MEDDLY::constrained_bckwd_dfs_evplus::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(res.getForest() == resF);

  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle c = 0;
  _compute(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, c);

  res.set(c, cev);
}

void MEDDLY::constrained_bckwd_dfs_evplus::_compute(int aev, node_handle a, int bev, node_handle b, node_handle r,
  long& cev, node_handle& c)
{
  // Partition NSF by levels
  splitMxd(r);

  // Execute saturation operation
  constrained_saturation_evplus* bckwdSatOp = new constrained_saturation_evplus(this, consF, argF, resF);
  bckwdSatOp->saturate(aev, a, bev, b, cev, c);

  // Cleanup
//  bckwdSatOp->removeAllComputeTableEntries();
  //delete bckwdSatOp;
//  removeAllComputeTableEntries();
  for (int i = transF->getNumVariables(); i > 0; i--) {
    transF->unlinkNode(splits[i]);
  }
  delete[] splits;
  splits = nullptr;
}

void MEDDLY::constrained_bckwd_dfs_evplus::saturateHelper(long aev, node_handle a, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);

  node_handle mxd = splits[nb.getLevel()];
  if (mxd == 0) {
    return;
  }

  const int mxdLevel = transF->getNodeLevel(mxd);
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(transF, nb.getLevel(), mxd, false)
    : unpacked_node::newFromNode(transF, mxd, false);

  unpacked_node* A = isLevelAbove(nb.getLevel(), consF->getNodeLevel(a))
    ? unpacked_node::newRedundant(consF, nb.getLevel(), 0L, a, true)
    : unpacked_node::newFromNode(consF, a, true);

  // indices to explore
  std::deque<int> queue;
  std::vector<bool> waiting(nb.getSize(), false);
  for (int j = 0; j < nb.getSize(); j++) {
    if (nb.d(j) != 0) {
      queue.push_back(j);
      waiting[j] = true;
    }
  }

  unpacked_node** Rps = new unpacked_node*[Ru->getNNZs()];
  for (int iz = 0; iz < Ru->getNNZs(); iz++) {
    const int i = Ru->i(iz);
    if (A->d(i) == 0) {
      Rps[iz] = nullptr;
    }
    else {
      Rps[iz] = (transF->getNodeLevel(Ru->d(iz)) == -nb.getLevel())
        ? unpacked_node::newFromNode(transF, Ru->d(iz), true)
        : unpacked_node::newIdentity(transF, -nb.getLevel(), i, Ru->d(iz), true);
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int j = queue.front();
    queue.pop_front();
    waiting[j] = false;

    MEDDLY_DCASSERT(nb.d(j) != 0);

    for (int iz = 0; iz < Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      if (A->d(i) == 0) {
        MEDDLY_DCASSERT(A->ei(i) == 0);
        continue;
      }

      MEDDLY_DCASSERT(Rps[iz] != nullptr);

      if (Rps[iz]->d(j) != 0) {
        long recev = 0;
        node_handle rec = 0;
        recFire(aev + A->ei(i), A->d(i), nb.ei(j), nb.d(j), Rps[iz]->d(j), recev, rec);
        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == 0) {
          MEDDLY_DCASSERT(recev == 0);
          continue;
        }

        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == nb.d(i)) {
          // Compute the minimum
          if (recev < nb.ei(i)) {
            nb.setEdge(i, recev);
          }
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (nb.d(i) == 0) {
          MEDDLY_DCASSERT(nb.ei(i) == 0);
          nb.setEdge(i, recev);
          nb.d_ref(i) = rec;
        }
        else {
          long accev = Inf<long>();
          node_handle acc = 0;
          minOp->compute(nb.ei(i), nb.d(i), recev, rec, accev, acc);
          resF->unlinkNode(rec);
          if (acc != nb.d(i)) {
            resF->unlinkNode(nb.d(i));
            nb.setEdge(i, accev);
            nb.d_ref(i) = acc;
          }
          else {
            MEDDLY_DCASSERT(accev == nb.ei(i));
            resF->unlinkNode(acc);
            updated = false;
          }
        }

        if (updated) {
          if (j == i) {
            // Restart inner for-loop.
            iz = -1;
          }
          else {
            if (!waiting[i]) {
              queue.push_back(i);
              waiting[i] = true;
            }
          }
        }
      }
    }
  }

  unpacked_node::recycle(Ru);
  for (int iz = 0; iz < Ru->getNNZs(); iz++) {
    unpacked_node::recycle(Rps[iz]);
  }
  delete[] Rps;
  unpacked_node::recycle(A);
}

void MEDDLY::constrained_bckwd_dfs_evplus::recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  // termination conditions
  if (a == 0 || b == 0 || r == 0) {
    cev = 0;
    c = 0;
    return;
  }
  if ((a == -1 || a == b) && r == -1) {
    cev = aev + bev;
    c = resF->linkNode(b);
    return;
  }

  // check the cache
  compute_table::entry_key* key = findResult(aev, a, bev, b, r, cev, c);
  if (key == 0) {
    MEDDLY_DCASSERT(cev >= 0);
    return;
  }

  // check if mxd and evmdd are at the same level
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);
  const int rLevel = transF->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), bLevel), aLevel);
  const int size = resF->getLevelSize(level);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, 0L, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(argF, level, 0L, b, true)
    : unpacked_node::newFromNode(argF, b, true);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  if (ABS(rLevel) < level) {
    // Assume identity reduction
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i = 0; i < size; i++) {
      long tev = 0;
      node_handle t = 0;
      recFire(aev + A->ei(i), A->d(i), bev + B->ei(i), B->d(i), r, tev, t);
      T->setEdge(i, tev);
      T->d_ref(i) = t;
    }
  }
  else {
    //
    // Need to process this level in the MXD.

    // clear out result (important!)
    for (int i = 0; i < size; i++) {
      T->setEdge(i, 0L);
      T->d_ref(i) = 0;
    }

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node* Ru = (rLevel < 0)
      ? unpacked_node::newRedundant(transF, -rLevel, r, false)
      : unpacked_node::newFromNode(transF, r, false);

    // loop over mxd "rows"
    for (int iz = 0; iz < Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      if (A->d(i) == 0) {
        MEDDLY_DCASSERT(A->ei(i) == 0);
        continue;
      }

      unpacked_node* Rp = isLevelAbove(-level, transF->getNodeLevel(Ru->d(iz)))
        ? unpacked_node::newIdentity(transF, -level, i, Ru->d(iz), false)
        : unpacked_node::newFromNode(transF, Ru->d(iz), false);

      // loop over mxd "columns"
      for (int jz = 0; jz < Rp->getNNZs(); jz++) {
        int j = Rp->i(jz);

        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = 0;
        node_handle n = 0;
        recFire(aev + A->ei(i), A->d(i), bev + B->ei(j), B->d(j), Rp->d(jz), nev, n);

        if (n == 0) {
          MEDDLY_DCASSERT(nev == 0);
          continue;
        }

        if (T->d(i) == 0) {
          MEDDLY_DCASSERT(T->ei(i) == 0);
          T->setEdge(i, nev);
          T->d_ref(i) = n;
          continue;
        }

        // there's new states and existing states; union them.
        const node_handle oldi = T->d(i);
        long newev = Inf<long>();
        node_handle newstates = 0;
        minOp->compute(nev, n, T->ei(i), oldi, newev, newstates);
        T->setEdge(i, newev);
        T->d_ref(i) = newstates;

        resF->unlinkNode(oldi);
        resF->unlinkNode(n);
      } // for j

      unpacked_node::recycle(Rp);
    } // for i

    unpacked_node::recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);

  saturateHelper(aev, a, *T);
  resF->createReducedNode(-1, T, cev, c);
  MEDDLY_DCASSERT(cev >= 0);

  saveResult(key, aev, a, bev, b, r, cev, c);
}

// ******************************************************************
// *                                                                *
// *                  constraint_bckwd_sat_evplus                   *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_saturation_evplus::constrained_saturation_evplus(constrained_bckwd_dfs_evplus* p,
  expert_forest* cons, expert_forest* arg, expert_forest* res)
  : specialized_operation(nullptr, 1)
{
  MEDDLY_DCASSERT(cons->isEVPlus());
  MEDDLY_DCASSERT(arg->isEVPlus());
  MEDDLY_DCASSERT(res->isEVPlus());

  parent = p;
  consF = cons;
  argF = arg;
  resF = res;

  registerInForest(consF);
  registerInForest(argF);
  registerInForest(resF);

  compute_table::entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new compute_table::entry_type("constrained_saturation_evplus", "LNNI:LN");
    et->setForestForSlot(1, cons);
    et->setForestForSlot(2, arg);
    et->setForestForSlot(6, res);
  } else {
    et = new compute_table::entry_type("constrained_saturation_evplus", "LNN:LN");
    et->setForestForSlot(1, cons);
    et->setForestForSlot(2, arg);
    et->setForestForSlot(5, res);
  }
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::constrained_saturation_evplus::~constrained_saturation_evplus()
{
  unregisterInForest(consF);
  unregisterInForest(argF);
  unregisterInForest(resF);
}

bool MEDDLY::constrained_saturation_evplus::checkForestCompatibility() const
{
  auto o1 = consF->variableOrder();
  auto o2 = argF->variableOrder();
  auto o3 = resF->variableOrder();
  return o1->is_compatible_with(*o2)
    && o1->is_compatible_with(*o3);
}

bool MEDDLY::constrained_saturation_evplus::checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c)
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

MEDDLY::compute_table::entry_key* MEDDLY::constrained_saturation_evplus::findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c)
{
  compute_table::entry_key* key = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(key);
  key->writeL(aev);
  key->writeN(a);
//  key->write(bev);
  key->writeN(b);
  if(argF->isFullyReduced()) {
    // Level is part of key for fully-reduced forest
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

void MEDDLY::constrained_saturation_evplus::saveResult(compute_table::entry_key* key,
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

void MEDDLY::constrained_saturation_evplus::saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c)
{
  saturate(aev, a, bev, b, argF->getNumVariables(), cev, c);
}

void MEDDLY::constrained_saturation_evplus::saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c)) {
    return;
  }

  compute_table::entry_key* key = findResult(aev, a, bev, b, level, cev, c);
  if (key == 0) {
    return;
  }

  const int sz = argF->getLevelSize(level);
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);

  unpacked_node* A = (aLevel < level)
    ? unpacked_node::newRedundant(consF, level, 0L, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = (bLevel < level)
    ? unpacked_node::newRedundant(argF, level, 0L, b, true)
    : unpacked_node::newFromNode(argF, b, true);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
    if (A->d(i) == 0) {
      MEDDLY_DCASSERT(resF == argF);
      T->setEdge(i, B->ei(i));
      T->d_ref(i) = resF->linkNode(B->d(i));
    }
    else {
      long tev = Inf<long>();
      node_handle t = 0;
      saturate(aev + A->ei(i), A->d(i), B->ei(i), B->d(i), level - 1, tev, t);
      T->setEdge(i, tev);
      T->d_ref(i) = t;
    }
  }

  // Cleanup
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);

  parent->saturateHelper(aev, a, *T);
  resF->createReducedNode(-1, T, cev, c);
  cev += bev;

  // save in compute table
  saveResult(key, aev, a, bev, b, level, cev, c);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_opname* MEDDLY::initConstrainedBFSBackward()
{
  return new constrained_bfs_opname(false);
}

MEDDLY::constrained_opname* MEDDLY::initConstrainedDFSForward()
{
  return new constrained_dfs_opname(true);
}

MEDDLY::constrained_opname* MEDDLY::initConstrainedDFSBackward()
{
  return new constrained_dfs_opname(false);
}
