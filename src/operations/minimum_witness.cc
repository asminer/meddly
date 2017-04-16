// $Id$

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
// *                    constraint_image_opname                     *
// *                                                                *
// ******************************************************************

MEDDLY::constraint_dfs_opname::constraint_dfs_opname(bool fwd)
 : minimum_witness_opname(fwd ? "PostImage" : "PreImage")
{
  forward = fwd;
}

MEDDLY::specialized_operation* MEDDLY::constraint_dfs_opname::buildOperation(
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res) const
{
  specialized_operation* op = 0;
  if (forward) {
    throw error(error::NOT_IMPLEMENTED);
  }
  else {
    op = new constraint_bckwd_dfs(this, cons, arg, trans, res);
  }
  return op;
}

// ******************************************************************
// *                                                                *
// *                     constraint_preimage                        *
// *                                                                *
// ******************************************************************

const int MEDDLY::constraint_bckwd_dfs::NODE_INDICES_IN_KEY[4] = {
  0,
  sizeof(node_handle) / sizeof(node_handle),
  2 * sizeof(node_handle) / sizeof(node_handle),
  (3 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle)
};

MEDDLY::constraint_bckwd_dfs::constraint_bckwd_dfs(const minimum_witness_opname* code,
  expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res)
  : specialized_operation(code,
      3 * (sizeof(node_handle)) / sizeof(node_handle),
      (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle))
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

  mxdIntersectionOp = getOperation(INTERSECTION, transF, transF, transF);
  mxdDifferenceOp = getOperation(DIFFERENCE, transF, transF, transF);
  plusOp = getOperation(PLUS, resF, resF, resF);
  minOp = getOperation(UNION, resF, resF, resF);

  splits = nullptr;
}

MEDDLY::constraint_bckwd_dfs::~constraint_bckwd_dfs()
{
  unregisterInForest(consF);
  unregisterInForest(argF);
  unregisterInForest(transF);
  unregisterInForest(resF);

  //delete mxdDifferenceOp;
  //delete mxdIntersectionOp;
  //delete plusOp;
  //delete minOp;
}

bool MEDDLY::constraint_bckwd_dfs::checkForestCompatibility() const
{
  auto o1 = consF->variableOrder();
  auto o2 = argF->variableOrder();
  auto o3 = transF->variableOrder();
  auto o4 = resF->variableOrder();
  return o1->is_compatible_with(*o2)
    && o1->is_compatible_with(*o3)
    && o1->is_compatible_with(*o4);
}

bool MEDDLY::constraint_bckwd_dfs::checkTerminals(int aev, node_handle a, int bev, node_handle b, node_handle c,
  long& dev, node_handle& d)
{
  if (a == -1 && b == -1 && c == -1) {
    d = -1;
    dev = bev;
    return true;
  }
  if (aev == Inf<long>() || bev == Inf<long>() || c == 0) {
//  if (a == 0 || b == 0 || c == 0) {
    d = 0;
    // XXX: 0 or infinity???
    dev = Inf<long>();
    return true;
  }
  return false;
}

MEDDLY::compute_table::search_key* MEDDLY::constraint_bckwd_dfs::findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle c, long& dev, node_handle &d)
{
  compute_table::search_key* key = useCTkey();
  MEDDLY_DCASSERT(key);
  key->reset();
  key->writeNH(a);
  key->writeNH(b);
  key->writeNH(c);

  compute_table::search_result& cacheFind = CT->find(key);
  if (!cacheFind) {
    return key;
  }
  cacheFind.read(dev);
  dev += aev + bev;
  d = resF->linkNode(cacheFind.readNH());
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
  entry.writeResult(dev - aev - bev);
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
  for (int i = transF->getNumVariables(); i > 0; i--) {
    transF->unlinkNode(splits[i]);
  }
  delete[] splits;
  splits = nullptr;
}

void MEDDLY::constraint_bckwd_dfs::saturateHelper(long aev, node_handle a, unpacked_node& nb)
{
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
    ? unpacked_node::newRedundant(consF, nb.getLevel(), 0, a, true)
    : unpacked_node::newFromNode(consF, a, true);

  // indices to explore
  std::deque<int> queue;
  std::vector<bool> waiting(nb.getSize(), false);
  for (int j = 0; j < nb.getSize(); j++) {
    if (nb.ei(j) != Inf<long>() && A->ei(j) != Inf<long>()) {
      queue.push_back(j);
      waiting[j] = true;
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int j = queue.front();
    queue.pop_front();
    waiting[j] = false;

    MEDDLY_DCASSERT(nb.d(j) && A->ei(j) != Inf<long>());

    for (int iz = 0; iz < Ru->getNNZs(); iz++) {
      const int i = Ru->i(iz);
      const int dlevel = transF->getNodeLevel(Ru->d(iz));
      unpacked_node* Rp = (dlevel == -nb.getLevel())
        ? unpacked_node::newFromNode(transF, Ru->d(iz), true)
        : unpacked_node::newIdentity(transF, -nb.getLevel(), i, Ru->d(iz), true);
      if (Rp->d(j) != 0) {
        long recev = Inf<long>();
        node_handle rec = 0;
        recFire(aev + A->ei(j), A->d(j), nb.ei(j), nb.d(j), Rp->d(j), recev, rec);

//        if (recev == Inf<long>()) {
        if (rec == 0) {
          continue;
        }

        MEDDLY_DCASSERT(recev != Inf<long>());
        // Increase the distance
        plusOp->compute(recev, rec, aev, a, recev, rec);

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
            if (!waiting[i] && A->ei(i) != Inf<long>()) {
              queue.push_back(i);
              waiting[i] = true;
            }
          }
        }
      }

      unpacked_node::recycle(Rp);
    }
  }

  unpacked_node::recycle(Ru);
  unpacked_node::recycle(A);
}

void MEDDLY::constraint_bckwd_dfs::recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  // termination conditions
  if (aev == Inf<long>() || r == 0) {
    cev = Inf<long>();
    c = 0;
    return;
  }
  if (a == -1 && r == -1) {
    // XXX: a == b ???
    cev = bev;
    c = b;
    return;
  }

  // check the cache
  compute_table::search_key* key = findResult(aev, a, bev, b, r, cev, c);
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
    ? unpacked_node::newRedundant(consF, level, 0, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(argF, level, 0, b, true)
    : unpacked_node::newFromNode(argF, b, true);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  if (ABS(rLevel) < level) {
    // Assume identity reduction
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i = 0; i < size; i++) {
      long tev = Inf<long>();
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
      T->setEdge(i, Inf<long>());
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
        int j = Rp->i(jz);

        if (A->ei(j) == Inf<long>() || Rp->d(jz) == 0) {
          continue;
        }

        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = Inf<long>();
        node_handle n = 0;
        recFire(aev + A->ei(i), A->d(i), bev + B->ei(j), B->d(j), Rp->d(jz), nev, n);

        if (nev == Inf<long>()) {
          continue;
        }
        if (T->ei(i) == Inf<long>()) {
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

  saveResult(key, aev, a, bev, b, r, cev, c);
}

// ******************************************************************
// *                                                                *
// *                     constraint_bckwd_sat                       *
// *                                                                *
// ******************************************************************

const int MEDDLY::constraint_saturation::NODE_INDICES_IN_KEY[3] = {
  sizeof(long) / sizeof(node_handle),
  (2 * sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
  (3 * sizeof(long) + 2 * sizeof(node_handle)) / sizeof(node_handle)
};

MEDDLY::constraint_saturation::constraint_saturation(constraint_bckwd_dfs* p,
  expert_forest* cons, expert_forest* arg, expert_forest* res)
  : specialized_operation(nullptr,
      (arg->isFullyReduced()
          ? (2 * (sizeof(long) + sizeof(node_handle)) + sizeof(int)) / sizeof(node_handle)
          : 2 * (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle)),
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

  minOp = getOperation(UNION, resF, resF, resF);
  plusOp = getOperation(PLUS, resF, resF, resF);
}

MEDDLY::constraint_saturation::~constraint_saturation()
{
  unregisterInForest(consF);
  unregisterInForest(argF);
  unregisterInForest(resF);

//  delete minOp;
//  delete plusOp;
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
    // XXX: 0 or infinity???
    cev = Inf<long>();
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
  key->write(bev);
  key->writeNH(b);
  if(argF->isFullyReduced()) {
    // Level is part of key for fully-reduced forest
    key->write(level);
  }

  compute_table::search_result &cacheFind = CT->find(key);
  if (!cacheFind) return key;
  cacheFind.read(cev);
  c = resF->linkNode(cacheFind.readNH());
  doneCTkey(key);
  return 0;
}

void MEDDLY::constraint_saturation::saveResult(compute_table::search_key* key,
  long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
{
  consF->cacheNode(a);
  argF->cacheNode(b);
  compute_table::entry_builder &entry = CT->startNewEntry(key);
  entry.writeResult(cev);
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
    ? unpacked_node::newRedundant(consF, level, 0, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = (bLevel < level)
    ? unpacked_node::newRedundant(argF, level, 0, b, true)
    : unpacked_node::newFromNode(argF, b, true);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
    long tev = Inf<long>();
    node_handle t = 0;
    saturate(A->ei(i), A->d(i), B->ei(i), B->d(i), level - 1, tev, t);
    T->setEdge(i, tev);
    T->d_ref(i) = t;
  }

  // Cleanup
  unpacked_node::recycle(A);
  unpacked_node::recycle(B);

  parent->saturateHelper(aev, a, *T);
  resF->createReducedNode(-1, T, cev, c);
  cev += bev;

  // save in compute table
  saveResult(key, aev, a, bev, b, cev, c);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::minimum_witness_opname* MEDDLY::initConstraintDFSBackward()
{
  return new constraint_dfs_opname(false);
}
