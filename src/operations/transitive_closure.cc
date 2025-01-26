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

#include "../defines.h"
#include "transitive_closure.h"

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"
#include "../oper_ternary.h"
#include "../ops_builtin.h"

namespace MEDDLY {

    class transitive_closure_forwd_bfs;

    class transitive_closure_dfs;
    class transitive_closure_forwd_dfs;

    class transitive_closure_evplus;

    ternary_list TRANSITIVE_CLOSURE_DFS_cache;
}

// ***********************************************************************
// ***********************************************************************
// ***                                                                 ***
// ***                             Classes                             ***
// ***                                                                 ***
// ***********************************************************************
// ***********************************************************************


// ******************************************************************
// *                                                                *
// *               transitive_closure_forwd_bfs class               *
// *                                                                *
// ******************************************************************

class MEDDLY::transitive_closure_forwd_bfs: public ternary_operation
{
protected:
  binary_operation* imageOp;
  binary_operation* plusOp;
  binary_operation* minOp;

  // void iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);
  void iterate(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& c);

public:
  transitive_closure_forwd_bfs(ternary_list &c, forest* cons, forest* tc,
          forest* trans, forest* res);

  virtual void compute(const dd_edge &a, const dd_edge &b, const dd_edge &r, dd_edge &res);
};


// ******************************************************************
// *                                                                *
// *                  transitive_closure_dfs class                  *
// *                                                                *
// ******************************************************************

class MEDDLY::transitive_closure_dfs: public ternary_operation
{
protected:
  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* minOp;

  dd_edge* splits;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, node_handle c, long& dev, node_handle& d);

  ct_entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle c, long& dev, node_handle& d);
  void saveResult(ct_entry_key* key,
    long aev, node_handle a, long bev, node_handle b, node_handle c, long dev, node_handle d);

  void splitMxd(const dd_edge& mxd);
  virtual void recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c) = 0;

public:
  transitive_closure_dfs(ternary_list &c, forest* cons, forest* tc,
          forest* trans, forest* res);

  virtual void computeDDEdge(const dd_edge &a, const dd_edge &b, const dd_edge &r, dd_edge &res, bool uf);
  void _compute(int aev, node_handle a, int bev, node_handle b, node_handle r, long& cev, node_handle& c);

  virtual void saturateHelper(long aev, node_handle a, int in, unpacked_node& nb) = 0;
};

// ******************************************************************
// *                                                                *
// *               transitive_closure_forwd_dfs class               *
// *                                                                *
// ******************************************************************


class MEDDLY::transitive_closure_forwd_dfs: public transitive_closure_dfs
{
protected:
  virtual void recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);

public:
  transitive_closure_forwd_dfs(ternary_list &c,
    forest* cons, forest* tc, forest* trans, forest* res);

  virtual void saturateHelper(long aev, node_handle a, int in, unpacked_node& nb);
};


// ******************************************************************
// *                                                                *
// *                transitive_closure_evplus  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::transitive_closure_evplus: public operation
{
protected:
  transitive_closure_dfs* parent;

  forest* tcF;
  forest* consF;
  forest* resF;

  virtual ~transitive_closure_evplus();

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  ct_entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c);
  void saveResult(ct_entry_key* Key,
    long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c);

public:
  transitive_closure_evplus(transitive_closure_dfs* p,
    forest* cons, forest* tc, forest* res);

  // high-level front-end
  void saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  void saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c);
};


// ***********************************************************************
// ***********************************************************************
// ***                                                                 ***
// ***                             Methods                             ***
// ***                                                                 ***
// ***********************************************************************
// ***********************************************************************


// ******************************************************************
// *                                                                *
// *              transitive_closure_forwd_bfs methods              *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_forwd_bfs
::transitive_closure_forwd_bfs(ternary_list &c, forest* cons, forest* tc,
    forest* trans, forest* res)
    : ternary_operation(c, 0, cons, tc, trans, res)
{
  if (resF->getRangeType() == range_type::INTEGER && resF->isForRelations()) {
    plusOp = POST_PLUS(resF, arg1F, resF);
    minOp = UNION(resF, resF, resF);
  } else {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  imageOp = TC_POST_IMAGE(arg2F, arg3F, resF);
}

void MEDDLY::transitive_closure_forwd_bfs::compute(const dd_edge &a, const dd_edge &b, const dd_edge &r, dd_edge &res)
{
  MEDDLY_DCASSERT(arg1F == a.getForest());
  MEDDLY_DCASSERT(arg2F == b.getForest());
  MEDDLY_DCASSERT(arg3F == r.getForest());
  MEDDLY_DCASSERT(resF == res.getForest());

  /*
  long aev = Inf<long>();
  a.getEdgeValue(aev);
  long bev = Inf<long>();
  b.getEdgeValue(bev);
  long cev = Inf<long>();
  node_handle cnode = 0;
  iterate(aev, a.getNode(), bev, b.getNode(), r.getNode(), cev, cnode);

  res.set(cev, cnode);
  */
  iterate(a, b, r, res);
}

/*
void MEDDLY::transitive_closure_forwd_bfs::iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  cev = bev;
  MEDDLY_DCASSERT(arg2F == resF);
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
// *                 transitive_closure_dfs methods                 *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_dfs::transitive_closure_dfs(ternary_list &c,
  forest* cons, forest* tc, forest* trans, forest* res)
  : ternary_operation(c, 1, cons, tc, trans, res)
{
  mxdIntersectionOp = INTERSECTION(arg3F, arg3F, arg3F);
  mxdDifferenceOp = DIFFERENCE(arg3F, arg3F, arg3F);
  minOp = UNION(resF, resF, resF);

  splits = nullptr;

  ct_entry_type* et = new ct_entry_type(c.getName(), "LNNN:LN");
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
  MEDDLY_DCASSERT(arg3F);
  MEDDLY_DCASSERT(splits == nullptr);

  splits = new dd_edge[arg3F->getNumVariables() + 1];

  dd_edge root(mxd), maxDiag(arg3F), Rpdi(arg3F);

  // Build from top down
  for (int level = arg3F->getMaxLevelIndex(); level > 0; level--) {
    splits[level].attach(arg3F);
    if (root.getNode() == 0) {
      // common and easy special case
      continue;
    }

    int mxdLevel = root.getLevel(); // arg3F->getNodeLevel(mxd);
    MEDDLY_DCASSERT(ABS(mxdLevel) <= level);

    // Initialize readers
    unpacked_node* Ru = isLevelAbove(level, mxdLevel)
      ? unpacked_node::newRedundant(arg3F, level, root.getNode(), FULL_ONLY)
      : arg3F->newUnpacked(root.getNode(), FULL_ONLY);

    bool first = true;

    // Read "rows"
    for (int i = 0; i < Ru->getSize(); i++) {
      // Initialize column reader
      int mxdPLevel = arg3F->getNodeLevel(Ru->down(i));
      unpacked_node* Rp = isLevelAbove(-level, mxdPLevel)
        ? unpacked_node::newIdentity(arg3F, -level, i, Ru->down(i), FULL_ONLY)
        : arg3F->newUnpacked(Ru->down(i), FULL_ONLY);

      // Intersect along the diagonal
      if (first) {
        maxDiag.set( arg3F->linkNode(Rp->down(i)) );
        first = false;
      } else {
        Rpdi.set( arg3F->linkNode(Rp->down(i)) );
        mxdIntersectionOp->computeTemp(maxDiag, Rpdi, maxDiag);
      }

      // cleanup
      unpacked_node::Recycle(Rp);
    } // for i

    // maxDiag is what we can split from here
    mxdDifferenceOp->computeTemp(root, maxDiag, splits[level]);
    root = maxDiag;

    // Cleanup
    unpacked_node::Recycle(Ru);
  } // for level

#ifdef DEBUG_SPLIT
  printf("After splitting monolithic event in msat\n");
  printf("splits array: [");
  for (unsigned k = 0; k <= arg3F->getNumVariables(); k++) {
    if (k) printf(", ");
    printf("%d", splits[k]);
  }
  printf("]\n");
#endif
}

void MEDDLY::transitive_closure_dfs::computeDDEdge(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res, bool uf)
{
  MEDDLY_DCASSERT(arg1F == a.getForest());
  MEDDLY_DCASSERT(arg2F == b.getForest());
  MEDDLY_DCASSERT(arg3F == r.getForest());
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

  res.set(cev, c);
}

void MEDDLY::transitive_closure_dfs::_compute(int aev, node_handle a, int bev, node_handle b, node_handle r,
  long& cev, node_handle& c)
{
  // Execute saturation operation
  transitive_closure_evplus* tcSatOp = new transitive_closure_evplus(this, arg1F, arg2F, resF);
  tcSatOp->saturate(aev, a, bev, b, cev, c);

  // Cleanup
  //tcSatOp->removeAllComputeTableEntries();
  // delete tcSatOp;
  //minOp->removeAllComputeTableEntries();
  //removeAllComputeTableEntries();

  // for (int i = arg3F->getNumVariables(); i > 0; i--) {
    // arg3F->unlinkNode(splits[i]);
  // }
  delete[] splits;
  splits = nullptr;
}

// ******************************************************************
// *                                                                *
// *              transitive_closure_forwd_dfs methods              *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_forwd_dfs::transitive_closure_forwd_dfs(
    ternary_list &c, forest* cons, forest* tc, forest* trans, forest* res)
  : transitive_closure_dfs(c, cons, tc, trans, res)
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

  // const int mxdLevel = arg3F->getNodeLevel(mxd);
  const int mxdLevel = mxd.getLevel();
  MEDDLY_DCASSERT(ABS(mxdLevel) == -nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(arg3F, -nb.getLevel(), mxd.getNode(), FULL_ONLY)
    : arg3F->newUnpacked(mxd.getNode(), FULL_ONLY);
  unpacked_node* Rp = unpacked_node::New(arg3F, SPARSE_ONLY);

  unpacked_node* A = isLevelAbove(-nb.getLevel(), arg1F->getNodeLevel(a))
    ? unpacked_node::newRedundant(arg1F, -nb.getLevel(), 0L, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY);

  dd_edge nbdj(resF), newst(resF);

  // indices to explore
  int size = nb.getSize();
  std::deque<int> queue;
  std::vector<bool> waiting(size, false);
  for (int ip = 0; ip < size; ip++) {
    if (nb.down(ip) != 0 && Ru->down(ip) != 0) {
      queue.push_back(ip);
      waiting[ip] = true;
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int ip = queue.front();
    queue.pop_front();
    waiting[ip] = false;

    MEDDLY_DCASSERT(nb.down(ip) != 0);
    MEDDLY_DCASSERT(Ru->down(ip) != 0);

    const int dlevel = arg3F->getNodeLevel(Ru->down(ip));
    if (dlevel == -Ru->getLevel()) {
      arg3F->unpackNode(Rp, Ru->down(ip), SPARSE_ONLY);
    }
    else {
      Rp->initIdentity(arg3F, -Ru->getLevel(), ip, Ru->down(ip), SPARSE_ONLY);
    }

    for (int jpz = 0; jpz < Rp->getSize(); jpz++) {
      const int jp = Rp->index(jpz);
      if (A->down(jp) == 0) {
        MEDDLY_DCASSERT(A->edge_long(jp) == 0);
        continue;
      }

      long recev = 0;
      node_handle rec = 0;
      recFire(aev + A->edge_long(jp), A->down(jp), nb.edge_long(ip), nb.down(ip), Rp->down(jpz), recev, rec);
      MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

      if (rec == 0) {
        MEDDLY_DCASSERT(recev == 0);
        continue;
      }

      if (rec == nb.down(jp)) {
        // Compute the minimum
        if (recev < nb.edge_long(jp)) {
          nb.setEdgeval(jp, recev);
        }
        resF->unlinkNode(rec);
        continue;
      }

      bool updated = true;

      if (nb.down(jp) == 0) {
        MEDDLY_DCASSERT(nb.edge_long(jp) == 0);
        nb.setFull(jp, edge_value(recev), rec);
        // nb.setEdge(jp, recev);
        // nb.d_ref(jp) = rec;
      }
      else {
        nbdj.set(nb.edge_long(jp), nb.down(jp));
        newst.set(recev, rec);
        minOp->computeTemp(nbdj, newst, nbdj);
        updated = (nbdj.getNode() != nb.down(jp));
        nb.setFull(jp, nbdj);

        // OLD
        /*
        long accev = Inf<long>();
        node_handle acc = 0;
        minOp->computeTemp(nb.edge_long(jp), nb.down(jp), recev, rec, accev, acc);

        MEDDLY_DCASSERT(acc != 0);
        resF->unlinkNode(rec);
        if (acc != nb.down(jp)) {
          resF->unlinkNode(nb.down(jp));
          nb.setEdge(jp, accev);
          nb.d_ref(jp) = acc;
        }
        else {
          MEDDLY_DCASSERT(accev == nb.edge_long(jp));
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
          if (!waiting[jp] && Ru->down(jp) != 0) {
            MEDDLY_DCASSERT(A->down(jp) != 0);
            queue.push_back(jp);
            waiting[jp] = true;
          }
        }
      }
    }
  }

  unpacked_node::Recycle(Ru);
  unpacked_node::Recycle(Rp);
  unpacked_node::Recycle(A);
}

void MEDDLY::transitive_closure_forwd_dfs::recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
{
  // termination conditions
  if (a == 0 || b == 0 || r == 0) {
    cev = 0;
    c = 0;
    return;
  }
  if (a == -1 && r == -1 && arg2F == resF) {
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
  const int aLevel = arg1F->getNodeLevel(a);
  MEDDLY_DCASSERT(aLevel >= 0);
  const int bLevel = arg2F->getNodeLevel(b);
  const int rLevel = arg3F->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), ABS(bLevel)), aLevel);
  const int size = resF->getLevelSize(level);

  dd_edge Tdj(resF), newst(resF);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(arg1F, level, 0L, a, FULL_ONLY)
    : arg1F->newUnpacked(a, FULL_ONLY);
  // Initialize evmxd reader
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(arg2F, level, 0L, b, FULL_ONLY)
    : arg2F->newUnpacked(b, FULL_ONLY);
  unpacked_node* D = unpacked_node::New(arg2F, FULL_ONLY);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);

  for (int i = 0; i < size; i++) {
    if (B->down(i) == 0) {
      // T->setEdge(i, 0L);
      // T->d_ref(i) = 0;
      // T is already zeroed
      continue;
    }

    if (isLevelAbove(-level, arg2F->getNodeLevel(B->down(i)))) {
      D->initIdentity(arg2F, -level, i, 0L, B->down(i), FULL_ONLY);
    }
    else {
      arg2F->unpackNode(D, B->down(i), FULL_ONLY);
    }

    unpacked_node* Tp = unpacked_node::newFull(resF, -level, size);

    if (ABS(rLevel) < level) {
      // Assume identity reduction
      MEDDLY_DCASSERT(arg3F->isIdentityReduced());
      for (int ip = 0; ip < size; ip++) {
        long tev = 0;
        node_handle t = 0;
        recFire(aev + A->edge_long(ip), A->down(ip), bev + B->edge_long(i) + D->edge_long(ip), D->down(ip), r, tev, t);
        Tp->setFull(ip, edge_value(tev), t);
        // Tp->setEdge(ip, tev);
        // Tp->d_ref(ip) = t;
      }
    }
    else {
      MEDDLY_DCASSERT(ABS(rLevel) == level);

      //
      // Need to process this level in the MXD.

      // clear out result (important!)
      // Tp is already zeroed
      /*
      for (int ip = 0; ip < size; ip++) {
        Tp->setEdge(ip, 0L);
        Tp->d_ref(ip) = 0;
      }
      */

      // Initialize mxd readers, note we might skip the unprimed level
      unpacked_node* Ru = (rLevel < 0)
        ? unpacked_node::newRedundant(arg3F, level, r, SPARSE_ONLY)
        : arg3F->newUnpacked(r, SPARSE_ONLY);
      unpacked_node* Rp = unpacked_node::New(arg3F, SPARSE_ONLY);

      // loop over mxd "rows"
      for (int ipz = 0; ipz < Ru->getSize(); ipz++) {
        const int ip = Ru->index(ipz);
        if (D->down(ip) == 0) {
          continue;
        }

        if (isLevelAbove(-level, arg3F->getNodeLevel(Ru->down(ipz)))) {
          Rp->initIdentity(arg3F, -level, ip, Ru->down(ipz), SPARSE_ONLY);
        }
        else {
          arg3F->unpackNode(Rp, Ru->down(ipz), SPARSE_ONLY);
        }

        // loop over mxd "columns"
        for (int jpz = 0; jpz < Rp->getSize(); jpz++) {
          const int jp = Rp->index(jpz);
          if (A->down(jp) == 0) {
            MEDDLY_DCASSERT(A->edge_long(jp) == 0);
            continue;
          }

          // ok, there is an ip->jp "edge".
          // determine new states to be added (recursively)
          // and add them
          long nev = 0;
          node_handle n = 0;
          recFire(aev + A->edge_long(jp), A->down(jp), bev + B->edge_long(i) + D->edge_long(ip), D->down(ip), Rp->down(jpz), nev, n);

          if (n == 0) {
            MEDDLY_DCASSERT(nev == 0);
            continue;
          }

          if (Tp->down(jp) == 0) {
            MEDDLY_DCASSERT(Tp->edge_long(jp) == 0);
            Tp->setFull(jp, edge_value(nev), n);
            // Tp->setEdge(jp, nev);
            // Tp->d_ref(jp) = n;
            continue;
          }

          // there's new states and existing states; union them.
          /*
          const node_handle oldjp = Tp->down(jp);
          long newev = Inf<long>();
          node_handle newstates = 0;
          minOp->computeTemp(nev, n, Tp->edge_long(jp), oldjp, newev, newstates);
          Tp->setEdge(jp, newev);
          Tp->d_ref(jp) = newstates;

          resF->unlinkNode(oldjp);
          resF->unlinkNode(n);
          */

          // there's new states and existing states; union them.
          Tdj.set(Tp->edge_long(jp), Tp->down(jp));
          newst.set(nev, n);
          minOp->computeTemp(newst, Tdj, Tdj);
          Tp->setFull(jp, Tdj);
        } // for j
      } // for i

      unpacked_node::Recycle(Ru);
      unpacked_node::Recycle(Rp);
    }

    saturateHelper(aev, a, i, *Tp);

    edge_value tpev;
    node_handle tp;
    resF->createReducedNode(Tp, tpev, tp, i);
    T->setFull(i, tpev, tp);
  }

  // cleanup mdd reader
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(B);
  unpacked_node::Recycle(D);

  edge_value ev;
  resF->createReducedNode(T, ev, c);
  cev = ev.getLong();
  MEDDLY_DCASSERT(cev >= 0);

  saveResult(key, aev, a, bev, b, r, cev, c);
}

// ******************************************************************
// *                                                                *
// *               transitive_closure_evplus  methods               *
// *                                                                *
// ******************************************************************

MEDDLY::transitive_closure_evplus::transitive_closure_evplus(transitive_closure_dfs* p,
  forest* cons, forest* tc, forest* res)
  : operation("transitive_closure_evplus", 1)
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
  saturate(aev, a, bev, b, tcF->getMaxLevelIndex(), cev, c);
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
    ? unpacked_node::newRedundant(consF, level, 0L, a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(tcF, level, 0L, b, FULL_ONLY)
    : tcF->newUnpacked(b, FULL_ONLY);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
      /* T is already zeroed
    if (B->down(i) == 0) {
      MEDDLY_DCASSERT(resF == tcF);
      T->setEdge(i, 0L);
      T->d_ref(i) = 0;
    }
    */
    if (B->down(i)) {
      unpacked_node* D = isLevelAbove(-level, tcF->getNodeLevel(B->down(i)))
        ? unpacked_node::newIdentity(tcF, -level, i, 0L, B->down(i), FULL_ONLY)
        : tcF->newUnpacked(B->down(i), FULL_ONLY);
      unpacked_node* Tp = unpacked_node::newFull(resF, -level, sz);
      for (int j = 0; j < sz; j++) {
          /*
        if (A->down(j) == 0 || D->down(j) == 0) {
          Tp->setEdge(j, 0L);
          Tp->d_ref(j) = 0;
        }
        */
        if (A->down(j) && D->down(j)) {
          long tpev = Inf<long>();
          node_handle tp = 0;
          saturate(aev + A->edge_long(j), A->down(j), bev + B->edge_long(i) + D->edge_long(j), D->down(j), level - 1, tpev, tp);
          Tp->setFull(j, edge_value(tpev), tp);
          // Tp->setEdge(j, tpev);
          // Tp->d_ref(j) = tp;
        }
      }
      unpacked_node::Recycle(D);

      parent->saturateHelper(aev, a, i, *Tp);

      edge_value tev;
      node_handle t;
      resF->createReducedNode(Tp, tev, t, i);
      T->setFull(i, tev, t);
    }
  }

  // Cleanup
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(B);

  edge_value ev;
  resF->createReducedNode(T, ev, c);
  cev = ev.getLong();

  // save in compute table
  saveResult(key, aev, a, bev, b, level, cev, c);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::ternary_operation*
MEDDLY::TRANSITIVE_CLOSURE_DFS(forest* consF, forest* inF,
        forest* relF, forest* outF)
{
    return new transitive_closure_forwd_dfs(TRANSITIVE_CLOSURE_DFS_cache,
        consF, inF, relF, outF);
}

void MEDDLY::TRANSITIVE_CLOSURE_DFS_init()
{
    TRANSITIVE_CLOSURE_DFS_cache.reset("trans_closure_dfs");
}

void MEDDLY::TRANSITIVE_CLOSURE_DFS_done()
{
    MEDDLY_DCASSERT( TRANSITIVE_CLOSURE_DFS_cache.isEmpty() );
}

