
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
#include "minimum_witness.h"

// ******************************************************************
// *                                                                *
// *                       common_constraint                        *
// *                                                                *
// ******************************************************************

MEDDLY::common_constraint::common_constraint(const minimum_witness_opname* code,
  int kl, int al,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : specialized_operation(code, kl, al)
{
  MEDDLY_DCASSERT(cons->isEVPlus());
  MEDDLY_DCASSERT(arg->isEVPlus());
  MEDDLY_DCASSERT(trans->isMultiTerminal() && trans->isForRelations());
  MEDDLY_DCASSERT(res->isEVPlus());

  consF = cons;
  argF = arg;
  transF = trans;
  resF = res;

  registerInForest(consF);
  registerInForest(argF);
  registerInForest(transF);
  registerInForest(resF);

  setAnswerForest(resF);
}

MEDDLY::common_constraint::~common_constraint()
{
  unregisterInForest(consF);
  unregisterInForest(argF);
  unregisterInForest(transF);
  unregisterInForest(resF);
}

bool MEDDLY::common_constraint::checkForestCompatibility() const
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

MEDDLY::constraint_bfs_opname::constraint_bfs_opname(bool fwd)
 : minimum_witness_opname(fwd ? "PostImage" : "PreImage")
{
  forward = fwd;
}

MEDDLY::specialized_operation* MEDDLY::constraint_bfs_opname::buildOperation(arguments* a) const
{
  minimum_witness_opname::minimum_witness_args* args = dynamic_cast<minimum_witness_opname::minimum_witness_args*>(a);
  specialized_operation* op = 0;
  if (forward) {
    throw error(error::NOT_IMPLEMENTED);
  }
  else {
    op = new constraint_bckwd_bfs(this,
      static_cast<expert_forest*>(args->consForest),
      static_cast<expert_forest*>(args->inForest),
      static_cast<expert_forest*>(args->relForest),
      static_cast<expert_forest*>(args->outForest));
  }
  return op;
}

// ******************************************************************
// *                                                                *
// *                    constraint_bckwd_bfs                        *
// *                                                                *
// ******************************************************************

MEDDLY::constraint_bckwd_bfs::constraint_bckwd_bfs(const minimum_witness_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : common_constraint(code, 0, 0, cons, arg, trans, res)
{
  if (resF->getRangeType() == forest::INTEGER) {
    plusOp = getOperation(PLUS, resF, consF, resF);
    minOp = getOperation(UNION, resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION);
  }
  imageOp = getOperation(PRE_IMAGE, argF, transF, resF);
}

void MEDDLY::constraint_bckwd_bfs::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
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

void MEDDLY::constraint_bckwd_bfs::iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
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
    // tev was increased by one step
    node_handle oldt = t;
    plusOp->compute(tev - 1, oldt, aev, a, tev, t);
    minOp->compute(cev, c, tev, t, cev, c);
    resF->unlinkNode(oldt);
    resF->unlinkNode(t);
  }
  resF->unlinkNode(prev);
}

bool MEDDLY::constraint_bckwd_bfs::isStaleEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::constraint_bckwd_bfs::discardEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::constraint_bckwd_bfs::showEntry(output &strm, const node_handle* entryData) const
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

// ******************************************************************
// *                                                                *
// *                     constraint_dfs_opname                      *
// *                                                                *
// ******************************************************************

MEDDLY::constraint_dfs_opname::constraint_dfs_opname(bool fwd)
 : minimum_witness_opname(fwd ? "PostImage" : "PreImage")
{
  forward = fwd;
}

MEDDLY::specialized_operation* MEDDLY::constraint_dfs_opname::buildOperation(arguments* a) const
{
  minimum_witness_opname::minimum_witness_args* args = dynamic_cast<minimum_witness_opname::minimum_witness_args*>(a);
  specialized_operation* op = 0;
  if (forward) {
    throw error(error::NOT_IMPLEMENTED);
  }
  else {
    op = new constraint_bckwd_dfs(this,
      static_cast<expert_forest*>(args->consForest),
      static_cast<expert_forest*>(args->inForest),
      static_cast<expert_forest*>(args->relForest),
      static_cast<expert_forest*>(args->outForest));
  }
  return op;
}

// ******************************************************************
// *                                                                *
// *                    constraint_bckwd_dfs                        *
// *                                                                *
// ******************************************************************

const int MEDDLY::constraint_bckwd_dfs::NODE_INDICES_IN_KEY[4] = {
  sizeof(long) / sizeof(node_handle),
  (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
  (sizeof(long) + 2 * sizeof(node_handle)) / sizeof(node_handle),
  (3 * sizeof(node_handle) + 2 * sizeof(long)) / sizeof(node_handle)
};

MEDDLY::constraint_bckwd_dfs::constraint_bckwd_dfs(const minimum_witness_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : common_constraint(code,
      (sizeof(long) + 3 * sizeof(node_handle)) / sizeof(node_handle),
      (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
      cons, arg, trans, res)
{
  mxdIntersectionOp = getOperation(INTERSECTION, transF, transF, transF);
  mxdDifferenceOp = getOperation(DIFFERENCE, transF, transF, transF);
  plusOp = getOperation(PLUS, resF, consF, resF);
  minOp = getOperation(UNION, resF, resF, resF);

  splits = nullptr;
}

MEDDLY::compute_table::search_key* MEDDLY::constraint_bckwd_dfs::findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle c, long& dev, node_handle &d)
{
  compute_table::search_key* key = useCTkey();
  MEDDLY_DCASSERT(key);
  key->reset();
  key->write(aev);
  key->writeNH(a);
  key->writeNH(b);
  key->writeNH(c);

  compute_table::search_result& cacheFind = CT->find(key);
  if (!cacheFind) {
    return key;
  }

  cacheFind.read(dev);
  d = resF->linkNode(cacheFind.readNH());
  if (d != 0) {
    dev += bev;
  }
  else {
    MEDDLY_DCASSERT(dev == 0);
  }

  doneCTkey(key);
  return 0;
}

void MEDDLY::constraint_bckwd_dfs::saveResult(compute_table::search_key* key,
  long aev, node_handle a, long bev, node_handle b, node_handle c, long dev, node_handle d)
{
  consF->cacheNode(a);
  argF->cacheNode(b);
  transF->cacheNode(c);
  compute_table::entry_builder& entry = CT->startNewEntry(key);
  if (d == 0) {
    // Write long
    entry.writeResult(0L);
  }
  else {
    entry.writeResult(dev - bev);
  }
  entry.writeResultNH(resF->cacheNode(d));
  CT->addEntry();
}

bool MEDDLY::constraint_bckwd_dfs::isStaleEntry(const node_handle* data)
{
  return consF->isStale(data[NODE_INDICES_IN_KEY[0]])
    || argF->isStale(data[NODE_INDICES_IN_KEY[1]])
    || transF->isStale(data[NODE_INDICES_IN_KEY[2]])
    || resF->isStale(data[NODE_INDICES_IN_KEY[3]]);
}

void MEDDLY::constraint_bckwd_dfs::discardEntry(const node_handle* data)
{
  consF->uncacheNode(data[NODE_INDICES_IN_KEY[0]]);
  argF->uncacheNode(data[NODE_INDICES_IN_KEY[1]]);
  transF->uncacheNode(data[NODE_INDICES_IN_KEY[2]]);
  resF->uncacheNode(data[NODE_INDICES_IN_KEY[3]]);
}

void MEDDLY::constraint_bckwd_dfs::showEntry(output &strm, const node_handle* data) const
{
  strm << "[" << getName()
    << "(" << long(data[NODE_INDICES_IN_KEY[0]])
    << ", " << long(data[NODE_INDICES_IN_KEY[1]])
    << ", " << long(data[NODE_INDICES_IN_KEY[2]])
    << "): " << long(data[NODE_INDICES_IN_KEY[3]])
    << "]";
}

// Partition the nsf based on "top level"
void MEDDLY::constraint_bckwd_dfs::splitMxd(node_handle mxd)
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

void MEDDLY::constraint_bckwd_dfs::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(res.getForest() == resF);

  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle c = 0;
  compute(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, c);

  res.set(c, cev);
}

void MEDDLY::constraint_bckwd_dfs::compute(int aev, node_handle a, int bev, node_handle b, node_handle r,
  long& cev, node_handle& c)
{
  // Partition NSF by levels
  splitMxd(r);

  // Execute saturation operation
  constraint_saturation* bckwdSatOp = new constraint_saturation(this, consF, argF, resF);
  bckwdSatOp->saturate(aev, a, bev, b, cev, c);

  // Cleanup
  bckwdSatOp->removeAllComputeTableEntries();
  //delete bckwdSatOp;
  removeAllComputeTableEntries();
  for (int i = transF->getNumVariables(); i > 0; i--) {
    transF->unlinkNode(splits[i]);
  }
  delete[] splits;
  splits = nullptr;
}

void MEDDLY::constraint_bckwd_dfs::saturateHelper(long aev, node_handle a, unpacked_node& nb)
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

        // Increase the distance
        node_handle oldrec = rec;
        plusOp->compute(recev, oldrec, aev + A->ei(i), A->d(i), recev, rec);
        resF->unlinkNode(oldrec);
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

void MEDDLY::constraint_bckwd_dfs::recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  // termination conditions
  if (a == 0 || b == 0 || r == 0) {
    cev = 0;
    c = 0;
    return;
  }
  if ((a == -1 || a == b) && r == -1) {
    cev = bev;
    c = resF->linkNode(b);
    return;
  }

  // check the cache
  compute_table::search_key* key = findResult(aev, a, bev, b, r, cev, c);
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
// *                     constraint_bckwd_sat                       *
// *                                                                *
// ******************************************************************

MEDDLY::constraint_saturation::constraint_saturation(constraint_bckwd_dfs* p,
  expert_forest* cons, expert_forest* arg, expert_forest* res)
  : specialized_operation(nullptr,
      (arg->isFullyReduced()
          ? (sizeof(long) + 2 * sizeof(node_handle) + sizeof(int)) / sizeof(node_handle)
          : (sizeof(long) + 2 * sizeof(node_handle)) / sizeof(node_handle)),
      (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle))
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

  setAnswerForest(resF);

  if (argF->isFullyReduced()) {
    NODE_INDICES_IN_KEY[0] = sizeof(long) / sizeof(node_handle);
    NODE_INDICES_IN_KEY[1] = (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle);
    // Store level in key for fully-reduced forest
    NODE_INDICES_IN_KEY[2] = (2 * sizeof(long) + 2 * sizeof(node_handle) + sizeof(int)) / sizeof(node_handle);
  }
  else {
    NODE_INDICES_IN_KEY[0] = sizeof(long) / sizeof(node_handle);
    NODE_INDICES_IN_KEY[1] = (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle);
    NODE_INDICES_IN_KEY[2] = 2 * (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle);
  }
}

MEDDLY::constraint_saturation::~constraint_saturation()
{
  unregisterInForest(consF);
  unregisterInForest(argF);
  unregisterInForest(resF);
}

bool MEDDLY::constraint_saturation::checkForestCompatibility() const
{
  auto o1 = consF->variableOrder();
  auto o2 = argF->variableOrder();
  auto o3 = resF->variableOrder();
  return o1->is_compatible_with(*o2)
    && o1->is_compatible_with(*o3);
}

bool MEDDLY::constraint_saturation::checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c)
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

MEDDLY::compute_table::search_key* MEDDLY::constraint_saturation::findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c)
{
  compute_table::search_key* key = useCTkey();
  MEDDLY_DCASSERT(key);
  key->reset();
  key->write(aev);
  key->writeNH(a);
//  key->write(bev);
  key->writeNH(b);
  if(argF->isFullyReduced()) {
    // Level is part of key for fully-reduced forest
    key->write(level);
  }

  compute_table::search_result &cacheFind = CT->find(key);
  if (!cacheFind) return key;
  cacheFind.read(cev);
  c = resF->linkNode(cacheFind.readNH());
  if (c != 0) {
    cev += bev;
  }
  else {
    MEDDLY_DCASSERT(cev == 0);
  }

  doneCTkey(key);
  return 0;
}

void MEDDLY::constraint_saturation::saveResult(compute_table::search_key* key,
  long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c)
{
  consF->cacheNode(a);
  argF->cacheNode(b);
  compute_table::entry_builder &entry = CT->startNewEntry(key);
  if (c == 0) {
    // Write long
    entry.writeResult(0L);
  }
  else {
    MEDDLY_DCASSERT(cev - bev >= 0);
    entry.writeResult(cev - bev);
  }
  entry.writeResultNH(resF->cacheNode(c));
  CT->addEntry();
}

bool MEDDLY::constraint_saturation::isStaleEntry(const node_handle* data)
{
  return consF->isStale(data[NODE_INDICES_IN_KEY[0]])
    || argF->isStale(data[NODE_INDICES_IN_KEY[1]])
    || resF->isStale(data[NODE_INDICES_IN_KEY[2]]);
}

void MEDDLY::constraint_saturation::discardEntry(const node_handle* data)
{
  consF->uncacheNode(data[NODE_INDICES_IN_KEY[0]]);
  argF->uncacheNode(data[NODE_INDICES_IN_KEY[1]]);
  resF->uncacheNode(data[NODE_INDICES_IN_KEY[2]]);
}

void MEDDLY::constraint_saturation::showEntry(output &strm, const node_handle* data) const
{
  strm << "[" << getName()
    << "(" << long(data[NODE_INDICES_IN_KEY[0]])
    << ", " << long(data[NODE_INDICES_IN_KEY[1]])
    << "): " << long(data[NODE_INDICES_IN_KEY[2]])
    << "]";
}

void MEDDLY::constraint_saturation::saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c)
{
  saturate(aev, a, bev, b, argF->getNumVariables(), cev, c);
}

void MEDDLY::constraint_saturation::saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c)) {
    return;
  }

  compute_table::search_key* key = findResult(aev, a, bev, b, level, cev, c);
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

MEDDLY::minimum_witness_opname* MEDDLY::initConstraintBFSBackward()
{
  return new constraint_bfs_opname(false);
}

MEDDLY::minimum_witness_opname* MEDDLY::initConstraintDFSBackward()
{
  return new constraint_dfs_opname(false);
}
