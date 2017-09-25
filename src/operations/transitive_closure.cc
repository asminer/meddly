
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

// XXX: DEBUG
void MEDDLY::printAssignmentsDEBUG(const dd_edge& res)
{
  int numVariables = ((expert_forest*)res.getForest())->getNumVariables();
  if (res.getForest()->isForRelations()) {
    unsigned counter = 0;
    for (enumerator iter(res); iter; ++iter, ++counter)
    {
      const int* element = iter.getAssignments();
      const int* curr = element + numVariables;
      const int* end = element;
      printf("%d: [%d", counter, *curr--);
      while (curr != end) {
        printf(" %d", *curr--);
      }
      printf("] -> [");
      element = iter.getPrimedAssignments();
      curr = element + numVariables;
      end = element;
      while (curr != end) {
        printf(" %d", *curr--);
      }
      long ev = Inf<long>();
      iter.getValue(ev);
      printf("]: %ld\n", ev);
    }
  }
  else {
    unsigned counter = 0;
    for (enumerator iter(res); iter; ++iter, ++counter)
    {
      const int* element = iter.getAssignments();
      const int* curr = element + numVariables;
      const int* end = element;
      printf("%d: [%d", counter, *curr--);
      while (curr != end) {
        printf(" %d", *curr--);
      }
      long ev = Inf<long>();
      iter.getValue(ev);
      printf("]: %ld\n", ev);
    }
  }
}






// ******************************************************************
// *                                                                *
// *                       common_constraint                        *
// *                                                                *
// ******************************************************************

MEDDLY::common_transitive_closure::common_transitive_closure(const minimum_witness_opname* code,
  int kl, int al,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : specialized_operation(code, kl, al)
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

  setAnswerForest(resF);
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

MEDDLY::transitive_closure_forwd_bfs::transitive_closure_forwd_bfs(const minimum_witness_opname* code,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : common_transitive_closure(code, 0, 0, cons, tc, trans, res)
{
  if (resF->getRangeType() == forest::INTEGER && resF->isForRelations()) {
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

  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  iterate(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, cnode);

  res.set(cnode, cev);
}

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
    imageOp->compute(cev, c, r, tev, t);
    // tev was increased by one step
    plusOp->compute(tev - 1, t, aev, a, tev, t);
    minOp->compute(cev, c, tev, t, cev, c);
    resF->unlinkNode(t);
  }
  resF->unlinkNode(prev);
}

bool MEDDLY::transitive_closure_forwd_bfs::isStaleEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::transitive_closure_forwd_bfs::discardEntry(const node_handle* entryData)
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

void MEDDLY::transitive_closure_forwd_bfs::showEntry(output &strm, const node_handle* entryData) const
{
  throw error(error::MISCELLANEOUS);
  // this operation won't add any CT entries.
}

// ******************************************************************
// *                                                                *
// *                     constraint_dfs_opname                      *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_dfs_opname::transitive_closure_dfs_opname()
 : minimum_witness_opname("Transitive Closure")
{
}

MEDDLY::specialized_operation* MEDDLY::transitive_closure_dfs_opname::buildOperation(arguments* a) const
{
  minimum_witness_opname::minimum_witness_args* args = dynamic_cast<minimum_witness_opname::minimum_witness_args*>(a);
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

const int MEDDLY::transitive_closure_dfs::NODE_INDICES_IN_KEY[4] = {
  sizeof(long) / sizeof(node_handle),
  (sizeof(node_handle) + sizeof(long)) / sizeof(node_handle),
  (2 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle),
  (3 * sizeof(node_handle) + 2 * sizeof(long)) / sizeof(node_handle)
};

MEDDLY::transitive_closure_dfs::transitive_closure_dfs(const minimum_witness_opname* code,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : common_transitive_closure(code,
      (3 * sizeof(node_handle) + sizeof(long)) / sizeof(node_handle),
      (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle),
      cons, tc, trans, res)
{
  mxdIntersectionOp = getOperation(INTERSECTION, transF, transF, transF);
  mxdDifferenceOp = getOperation(DIFFERENCE, transF, transF, transF);
  minOp = getOperation(UNION, resF, resF, resF);

  splits = nullptr;
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

MEDDLY::compute_table::search_key* MEDDLY::transitive_closure_dfs::findResult(long aev, node_handle a,
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
  MEDDLY_DCASSERT(dev >= 0);

  doneCTkey(key);
  return 0;
}

void MEDDLY::transitive_closure_dfs::saveResult(compute_table::search_key* key,
  long aev, node_handle a, long bev, node_handle b, node_handle c, long dev, node_handle d)
{
  consF->cacheNode(a);
  tcF->cacheNode(b);
  transF->cacheNode(c);
  compute_table::entry_builder& entry = CT->startNewEntry(key);
  if (d == 0) {
    // Write long
    entry.writeResult(0L);
  }
  else {
    MEDDLY_DCASSERT(dev - bev >= 0);
    entry.writeResult(dev - bev);
  }
  entry.writeResultNH(resF->cacheNode(d));
  CT->addEntry();
}

bool MEDDLY::transitive_closure_dfs::isStaleEntry(const node_handle* data)
{
  return consF->isStale(data[NODE_INDICES_IN_KEY[0]])
    || tcF->isStale(data[NODE_INDICES_IN_KEY[1]])
    || transF->isStale(data[NODE_INDICES_IN_KEY[2]])
    || resF->isStale(data[NODE_INDICES_IN_KEY[3]]);
}

void MEDDLY::transitive_closure_dfs::discardEntry(const node_handle* data)
{
  consF->uncacheNode(data[NODE_INDICES_IN_KEY[0]]);
  tcF->uncacheNode(data[NODE_INDICES_IN_KEY[1]]);
  transF->uncacheNode(data[NODE_INDICES_IN_KEY[2]]);
  resF->uncacheNode(data[NODE_INDICES_IN_KEY[3]]);
}

void MEDDLY::transitive_closure_dfs::showEntry(output &strm, const node_handle* data) const
{
  strm << "[" << getName()
    << "(" << long(data[NODE_INDICES_IN_KEY[0]])
    << ", " << long(data[NODE_INDICES_IN_KEY[1]])
    << ", " << long(data[NODE_INDICES_IN_KEY[2]])
    << "): " << long(data[NODE_INDICES_IN_KEY[3]])
    << "]";
}

// Partition the nsf based on "top level"
void MEDDLY::transitive_closure_dfs::splitMxd(node_handle mxd)
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

void MEDDLY::transitive_closure_dfs::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(consF == a.getForest());
  MEDDLY_DCASSERT(tcF == b.getForest());
  MEDDLY_DCASSERT(transF == r.getForest());
  MEDDLY_DCASSERT(resF == res.getForest());

  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle c = 0;
  compute(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, c);

  res.set(c, cev);
}

void MEDDLY::transitive_closure_dfs::compute(int aev, node_handle a, int bev, node_handle b, node_handle r,
  long& cev, node_handle& c)
{
  // Partition NSF by levels
  splitMxd(r);

  // Execute saturation operation
  transitive_closure_evplus* tcSatOp = new transitive_closure_evplus(this, consF, tcF, resF);
  tcSatOp->saturate(aev, a, bev, b, cev, c);

  // Cleanup
  tcSatOp->removeAllComputeTableEntries();
  for (int i = transF->getNumVariables(); i > 0; i--) {
    transF->unlinkNode(splits[i]);
  }
  delete[] splits;
  splits = nullptr;
}

// ******************************************************************
// *                                                                *
// *                transitive_closure_forwd_dfs                    *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_forwd_dfs::transitive_closure_forwd_dfs(const minimum_witness_opname* code,
  expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res)
  : transitive_closure_dfs(code, cons, tc, trans, res)
{
}

void MEDDLY::transitive_closure_forwd_dfs::saturateHelper(long aev, node_handle a, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);
  MEDDLY_DCASSERT(nb.getLevel() > 0);

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
    ? unpacked_node::newRedundant(consF, nb.getLevel(), 0L, a, true)
    : unpacked_node::newFromNode(consF, a, true);

  int size = nb.getSize();
  for (int i = 0; i < size; i++) {
    if (nb.d(i) == 0) {
      continue;
    }

    unpacked_node* D = isLevelAbove(Ru->getLevel(), ABS(resF->getNodeLevel(nb.d(i))))
      ? unpacked_node::newIdentity(resF, -nb.getLevel(), i, 0L, nb.d(i), true)
      : unpacked_node::newFromNode(resF, nb.d(i), true);

    // indices to explore
    std::deque<int> queue;
    std::vector<bool> waiting(size, false);
    for (int ip = 0; ip < size; ip++) {
      // Must link the node because D will be updated and reduced to a new node
      resF->linkNode(D->d(ip));

      if (D->d(ip) != 0) {
        D->setEdge(ip, nb.ei(i) + D->ei(ip));

        if (Ru->d(ip) != 0) {
          queue.push_back(ip);
          waiting[ip] = true;
        }
      }
    }

    // explore indexes
    while (!queue.empty()) {
      const int ip = queue.front();
      queue.pop_front();
      waiting[ip] = false;

      MEDDLY_DCASSERT(D->d(ip) != 0);
      MEDDLY_DCASSERT(Ru->d(ip) != 0);

      const int dlevel = transF->getNodeLevel(Ru->d(ip));
      unpacked_node* Rp = (dlevel == -Ru->getLevel())
        ? unpacked_node::newFromNode(transF, Ru->d(ip), false)
        : unpacked_node::newIdentity(transF, -Ru->getLevel(), ip, Ru->d(ip), false);

      for (int jpz = 0; jpz < Rp->getNNZs(); jpz++) {
        const int jp = Rp->i(jpz);
        if (A->d(jp) == 0) {
          MEDDLY_DCASSERT(A->ei(jp) == 0);
          continue;
        }

        long recev = 0;
        node_handle rec = 0;
        recFire(aev + A->ei(jp), A->d(jp), D->ei(ip), D->d(ip), Rp->d(jpz), recev, rec);
        MEDDLY_DCASSERT(isLevelAbove(D->getLevel(), resF->getNodeLevel(rec)));

        if (rec == 0) {
          MEDDLY_DCASSERT(recev == 0);
          continue;
        }

        MEDDLY_DCASSERT(isLevelAbove(D->getLevel(), resF->getNodeLevel(rec)));

        if (rec == D->d(jp)) {
          // Compute the minimum
          if (recev < D->ei(jp)) {
            D->setEdge(jp, recev);
          }
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (D->d(jp) == 0) {
          MEDDLY_DCASSERT(D->ei(jp) == 0);
          D->setEdge(jp, recev);
          D->d_ref(jp) = rec;
        }
        else {
          long accev = Inf<long>();
          node_handle acc = 0;
          minOp->compute(D->ei(jp), D->d(jp), recev, rec, accev, acc);

          MEDDLY_DCASSERT(acc != 0);
          resF->unlinkNode(rec);
          if (acc != D->d(jp)) {
            resF->unlinkNode(D->d(jp));
            D->setEdge(jp, accev);
            D->d_ref(jp) = acc;
          }
          else {
            MEDDLY_DCASSERT(accev == D->ei(jp));
            resF->unlinkNode(acc);
            updated = false;
          }
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

      unpacked_node::recycle(Rp);
    }

    long tpev = Inf<long>();
    node_handle tp = 0;
    resF->createReducedNode(i, D, tpev, tp);

    resF->unlinkNode(nb.d(i));
    nb.setEdge(i, tpev);
    nb.d_ref(i) = tp;
  }

  unpacked_node::recycle(Ru);
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
  if (a == -1 && r == -1) {
    cev = aev + bev;
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
  MEDDLY_DCASSERT(aLevel >= 0);
  const int bLevel = tcF->getNodeLevel(b);
  const int rLevel = transF->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), ABS(bLevel)), aLevel);
  const int size = resF->getLevelSize(level);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, 0L, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  // Initialize evmxd reader
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(resF, level, 0L, b, true)
    : unpacked_node::newFromNode(resF, b, true);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);

  for (int i = 0; i < size; i++) {
    unpacked_node* D = isLevelAbove(-level, resF->getNodeLevel(B->d(i)))
      ? unpacked_node::newIdentity(resF, -level, i, 0L, B->d(i), true)
      : unpacked_node::newFromNode(resF, B->d(i), true);

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
        : unpacked_node::newFromNode(transF, r, false);

      // loop over mxd "rows"
      for (int ipz = 0; ipz < Ru->getNNZs(); ipz++) {
        const int ip = Ru->i(ipz);

        unpacked_node* Rp = isLevelAbove(-level, transF->getNodeLevel(Ru->d(ipz)))
          ? unpacked_node::newIdentity(transF, -level, ip, Ru->d(ipz), false)
          : unpacked_node::newFromNode(transF, Ru->d(ipz), false);

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
          const node_handle oldjp = Tp->d(jp);
          long newev = Inf<long>();
          node_handle newstates = 0;
          minOp->compute(nev, n, Tp->ei(jp), oldjp, newev, newstates);
          Tp->setEdge(jp, newev);
          Tp->d_ref(jp) = newstates;

          resF->unlinkNode(oldjp);
          resF->unlinkNode(n);
        } // for j

        unpacked_node::recycle(Rp);
      } // for i

      unpacked_node::recycle(Ru);
    }

    long tpev = Inf<long>();
    node_handle tp = 0;
    resF->createReducedNode(i, Tp, tpev, tp);
    T->setEdge(i, tpev);
    T->d_ref(i) = tp;
  }

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
// *                  transitive_closure_evplus                     *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_evplus::transitive_closure_evplus(transitive_closure_dfs* p,
  expert_forest* cons, expert_forest* tc, expert_forest* res)
  : specialized_operation(nullptr,
      ((tc->isFullyReduced() || tc->isIdentityReduced())
          ? (sizeof(long) + 2 * sizeof(node_handle) + sizeof(int)) / sizeof(node_handle)
          : (sizeof(long) + 2 * sizeof(node_handle)) / sizeof(node_handle)),
      (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle))
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

  setAnswerForest(resF);

  if (tcF->isFullyReduced() || tcF->isIdentityReduced()) {
    NODE_INDICES_IN_KEY[0] = sizeof(long) / sizeof(node_handle);
    NODE_INDICES_IN_KEY[1] = (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle);
    // Store level in key for fully-reduced or identity-reduced forest
    NODE_INDICES_IN_KEY[2] = (2 * (sizeof(long) + sizeof(node_handle)) + sizeof(int)) / sizeof(node_handle);
  }
  else {
    NODE_INDICES_IN_KEY[0] = sizeof(long) / sizeof(node_handle);
    NODE_INDICES_IN_KEY[1] = (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle);
    NODE_INDICES_IN_KEY[2] = 2 * (sizeof(long) + sizeof(node_handle)) / sizeof(node_handle);
  }
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

MEDDLY::compute_table::search_key* MEDDLY::transitive_closure_evplus::findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c)
{
  compute_table::search_key* key = useCTkey();
  MEDDLY_DCASSERT(key);
  key->reset();
  key->write(aev);
  key->writeNH(a);
  key->writeNH(b);
  if(tcF->isFullyReduced() || tcF->isIdentityReduced()) {
    // Level is part of key for fully-reduced or identity-reduced forest
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

void MEDDLY::transitive_closure_evplus::saveResult(compute_table::search_key* key,
  long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c)
{
  consF->cacheNode(a);
  tcF->cacheNode(b);
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

bool MEDDLY::transitive_closure_evplus::isStaleEntry(const node_handle* data)
{
  return consF->isStale(data[NODE_INDICES_IN_KEY[0]])
    || tcF->isStale(data[NODE_INDICES_IN_KEY[1]])
    || resF->isStale(data[NODE_INDICES_IN_KEY[2]]);
}

void MEDDLY::transitive_closure_evplus::discardEntry(const node_handle* data)
{
  consF->uncacheNode(data[NODE_INDICES_IN_KEY[0]]);
  tcF->uncacheNode(data[NODE_INDICES_IN_KEY[1]]);
  resF->uncacheNode(data[NODE_INDICES_IN_KEY[2]]);
}

void MEDDLY::transitive_closure_evplus::showEntry(output &strm, const node_handle* data) const
{
  strm << "[" << getName()
    << "(" << long(data[NODE_INDICES_IN_KEY[0]])
    << ", " << long(data[NODE_INDICES_IN_KEY[1]])
    << "): " << long(data[NODE_INDICES_IN_KEY[2]])
    << "]";
}

void MEDDLY::transitive_closure_evplus::saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c)
{
  saturate(aev, a, bev, b, tcF->getNumVariables(), cev, c);
}

void MEDDLY::transitive_closure_evplus::saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c)
{
//  if (level == 0 && a == -1 && b == -1) {
//    cev = bev;
//    c = -1;
//    return;
//  }

  if (checkTerminals(aev, a, bev, b, cev, c)) {
    return;
  }

  compute_table::search_key* key = findResult(aev, a, bev, b, level, cev, c);
  if (key == 0) {
    return;
  }

  const int sz = tcF->getLevelSize(level);
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = tcF->getNodeLevel(b);

  MEDDLY_DCASSERT(aLevel >= 0);

  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, 0L, a, true)
    : unpacked_node::newFromNode(consF, a, true);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(tcF, level, 0L, b, true)
    : unpacked_node::newFromNode(tcF, b, true);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
    if (B->d(i) == 0) {
      MEDDLY_DCASSERT(resF == tcF);
      T->setEdge(i, 0L);
      T->d_ref(i) = 0;
    }
    else {
      int dLevel = tcF->getNodeLevel(B->d(i));
      unpacked_node* D = isLevelAbove(level, ABS(dLevel))
        ? unpacked_node::newIdentity(tcF, -level, i, 0L, B->d(i), true)
        : unpacked_node::newFromNode(tcF, B->d(i), true);
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

  parent->saturateHelper(aev, a, *T);
  resF->createReducedNode(-1, T, cev, c);

  // save in compute table
  saveResult(key, aev, a, bev, b, level, cev, c);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::minimum_witness_opname* MEDDLY::initTransitiveClosureDFS()
{
  return new transitive_closure_dfs_opname();
}
