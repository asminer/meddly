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
#include "constrained.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"
#include "../oper_ternary.h"
#include "../ops_builtin.h"


namespace MEDDLY {
    class constrained_bckwd_bfs_evplus;

    class constrained_dfs_mt;
    class constrained_forwd_dfs_mt;
    class constrained_bckwd_dfs_mt;
    class constrained_bckwd_dfs_evplus;

    class constrained_saturation_mt;
    class constrained_saturation_evplus;



};

// ***********************************************************************
// ***********************************************************************
// ***                                                                 ***
// ***                             Classes                             ***
// ***                                                                 ***
// ***********************************************************************
// ***********************************************************************

// ******************************************************************
// *                                                                *
// *               constrained_bckwd_bfs_evplus class               *
// *                                                                *
// ******************************************************************

class MEDDLY::constrained_bckwd_bfs_evplus: public ternary_operation
{
    protected:
        binary_operation* imageOp;
        binary_operation* plusOp;
        binary_operation* minOp;

        void iterate(const dd_edge& a, const dd_edge& b,
                const dd_edge& r, dd_edge& c);

    public:
        constrained_bckwd_bfs_evplus(ternary_list &c, forest* cons,
                forest* arg, forest* trans, forest* res);

        virtual void compute(const dd_edge& a, const dd_edge& b,
                const dd_edge& r, dd_edge& res);
};

// ******************************************************************
// *                                                                *
// *                    constrained_dfs_mt class                    *
// *                                                                *
// ******************************************************************

class MEDDLY::constrained_dfs_mt: public ternary_operation
{
protected:
  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* unionOp;

  dd_edge* splits;

  ct_entry_key* findResult(node_handle a, node_handle b, node_handle r, node_handle& c);
  void saveResult(ct_entry_key* key,
    node_handle a, node_handle b, node_handle r, node_handle c);

  void splitMxd(const dd_edge& mxd);

public:
  constrained_dfs_mt(ternary_list &c, forest* cons, forest* arg,
          forest* trans, forest* res);

  virtual void compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res);
  void _compute(node_handle a, node_handle b, node_handle r, node_handle& c);

  virtual void saturateHelper(node_handle a, unpacked_node& nb) = 0;
};

// ******************************************************************
// *                                                                *
// *                 constrained_forwd_dfs_mt class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::constrained_forwd_dfs_mt: public constrained_dfs_mt
{
protected:
  void recFire(node_handle a, node_handle b, node_handle r, node_handle& c);

public:
  constrained_forwd_dfs_mt(ternary_list &c,
    forest* cons, forest* arg, forest* trans, forest* res);

  virtual void saturateHelper(node_handle a, unpacked_node& nb) override;
};

// ******************************************************************
// *                                                                *
// *                 constrained_bckwd_dfs_mt class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::constrained_bckwd_dfs_mt: public constrained_dfs_mt
{
protected:
  void recFire(node_handle a, node_handle b, node_handle r, node_handle& c);

public:
  constrained_bckwd_dfs_mt(ternary_list &c,
    forest* cons, forest* arg, forest* trans, forest* res);

  virtual void saturateHelper(node_handle a, unpacked_node& nb) override;
};

// ******************************************************************
// *                                                                *
// *                constrained_saturation_mt  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::constrained_saturation_mt: public operation
{
protected:
  constrained_dfs_mt* parent;

  forest* consF;
  forest* argF;
  forest* resF;

  virtual ~constrained_saturation_mt();

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

  bool checkTerminals(node_handle a, node_handle b, node_handle& c);

  ct_entry_key* findResult(node_handle a, node_handle b, int level, node_handle &c);
  void saveResult(ct_entry_key* Key,
    node_handle a, node_handle b, int level, node_handle c);

public:
  constrained_saturation_mt(constrained_dfs_mt* p,
    forest* cons, forest* arg, forest* res);

  bool matches(const forest* arg1, const forest* arg2,
    const forest* res) const;

  // high-level front-end
  void saturate(node_handle a, node_handle b, node_handle& c);

  void saturate(node_handle a, node_handle b, int level, node_handle& c);
};

// ******************************************************************
// *                                                                *
// *               constrained_bckwd_dfs_evplus class               *
// *                                                                *
// ******************************************************************

class MEDDLY::constrained_bckwd_dfs_evplus: public ternary_operation
{
protected:
  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* minOp;

  dd_edge* splits;

  ct_entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle r, long& dev, node_handle& d);
  void saveResult(ct_entry_key* key,
    long aev, node_handle a, long bev, node_handle b, node_handle r, long dev, node_handle d);

  void splitMxd(const dd_edge& mxd);
  void recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);

public:
  constrained_bckwd_dfs_evplus(ternary_list &c,
    forest* cons, forest* arg, forest* trans, forest* res);

  virtual void compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res);
  void _compute(int aev, node_handle a, int bev, node_handle b, node_handle r, long& cev, node_handle& c);

  void saturateHelper(long aev, node_handle a, unpacked_node& nb);
};

// ******************************************************************
// *                                                                *
// *              constrained_saturation_evplus  class              *
// *                                                                *
// ******************************************************************

class MEDDLY::constrained_saturation_evplus: public operation
{
protected:
  constrained_bckwd_dfs_evplus* parent;

  forest* consF;
  forest* argF;
  forest* resF;

  virtual ~constrained_saturation_evplus();

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  ct_entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c);
  void saveResult(ct_entry_key* Key,
    long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c);

public:
  constrained_saturation_evplus(constrained_bckwd_dfs_evplus* p,
    forest* cons, forest* arg, forest* res);

  bool matches(const forest* arg1, const forest* arg2,
    const forest* res) const;

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
// *              constrained_bckwd_bfs_evplus methods              *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_bckwd_bfs_evplus::constrained_bckwd_bfs_evplus(
        ternary_list &c, forest* cons, forest* arg, forest* trans,
        forest* res) : ternary_operation(c, 0, cons, arg, trans, res)
{
    if (resF->getRangeType() == range_type::INTEGER) {
        plusOp = PLUS(resF, arg1F, resF);
        minOp = UNION(resF, resF, resF);
    } else {
        throw error(error::INVALID_OPERATION);
    }
    imageOp = PRE_IMAGE(arg2F, arg3F, resF);
}

void MEDDLY::constrained_bckwd_bfs_evplus::compute(const dd_edge& a,
        const dd_edge& b, const dd_edge& r, dd_edge& res)
{
    MEDDLY_DCASSERT(res.getForest() == resF);

    iterate(a, b, r, res);
}

//void MEDDLY::constrained_bckwd_bfs_evplus::iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c)
void MEDDLY::constrained_bckwd_bfs_evplus::iterate(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& c)
{
  /*
  cev = bev;
  MEDDLY_DCASSERT(argF == resF);
  c = resF->linkNode(b);

  node_handle prev = 0;
  while (prev != c) {
    resF->unlinkNode(prev);
    prev = c;
    long tev = Inf<long>();
    node_handle t = 0;
    imageOp->computeTemp(cev, c, r, tev, t);
    node_handle oldt = t;
    plusOp->computeTemp(tev, oldt, aev, a, tev, t);
    minOp->computeTemp(cev, c, tev, t, cev, c);
    resF->unlinkNode(oldt);
    resF->unlinkNode(t);
  }
  resF->unlinkNode(prev);
  */
  c = b;
  dd_edge prev(resF), t(resF);
  while (prev != c) {
    prev = c;
    imageOp->computeTemp(c, r, t);
    plusOp->computeTemp(t, a, t);
    minOp->computeTemp(c, t, c);
  }
}

// ******************************************************************
// *                                                                *
// *                   constrained_dfs_mt methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_dfs_mt::constrained_dfs_mt(ternary_list &c,
  forest* cons, forest* arg, forest* trans, forest* res)
  : ternary_operation(c, 1, cons, arg, trans, res)
{
  MEDDLY_DCASSERT(cons->isMultiTerminal() && !cons->isForRelations());
  MEDDLY_DCASSERT(arg->isMultiTerminal() && !arg->isForRelations());
  MEDDLY_DCASSERT(trans->isMultiTerminal() && trans->isForRelations());
  MEDDLY_DCASSERT(res->isMultiTerminal() && !res->isForRelations());

  mxdIntersectionOp = INTERSECTION(arg3F, arg3F, arg3F);
  mxdDifferenceOp = DIFFERENCE(arg3F, arg3F, arg3F);
  unionOp = UNION(resF, resF, resF);

  splits = nullptr;

  ct_entry_type* et = new ct_entry_type(c.getName(), "NNN:N");
  et->setForestForSlot(0, arg1F);
  et->setForestForSlot(1, arg2F);
  et->setForestForSlot(2, arg3F);
  et->setForestForSlot(4, resF);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::ct_entry_key* MEDDLY::constrained_dfs_mt::findResult(
    node_handle a, node_handle b, node_handle r, node_handle &c)
{
  ct_entry_key* key = CT0->useEntryKey(etype[0], 0);
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

void MEDDLY::constrained_dfs_mt::saveResult(ct_entry_key* key,
  node_handle a, node_handle b, node_handle r, node_handle c)
{
  CTresult[0].reset();
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}


// Partition the nsf based on "top level"
void MEDDLY::constrained_dfs_mt::splitMxd(const dd_edge& mxd)
{
  MEDDLY_DCASSERT(arg3F);
  MEDDLY_DCASSERT(splits == nullptr);

  splits = new dd_edge[arg3F->getNumVariables() + 1];

  dd_edge maxDiag(arg3F), root(mxd), Rpdi(arg3F);

  // Build from top down
  for (int level = arg3F->getMaxLevelIndex(); level > 0; level--) {
    splits[level].attach(arg3F);

    if (root.getNode() == 0) {
      // common and easy special case
      continue;
    }

    int mxdLevel = root.getLevel();
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

void MEDDLY::constrained_dfs_mt::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(res.getForest() == resF);

  node_handle c = 0;
  if (a.getNode() == 0) {
    c = smart_cast<forest*>(resF)->linkNode(b.getNode());
  } else {
    // Partition NSF by levels
    splitMxd(r);

    _compute(a.getNode(), b.getNode(), r.getNode(), c);
  }
  res.set(c);
}

void MEDDLY::constrained_dfs_mt::_compute(node_handle a, node_handle b, node_handle r,
  node_handle& c)
{
  // Execute saturation operation
  constrained_saturation_mt* bckwdSatOp = new constrained_saturation_mt(this, consF, argF, resF);
  bckwdSatOp->saturate(a, b, c);

  // Cleanup
//  bckwdSatOp->removeAllComputeTableEntries();
  //delete bckwdSatOp;
//  removeAllComputeTableEntries();
  delete[] splits;
  splits = nullptr;
}

// ******************************************************************
// *                                                                *
// *                constrained_forwd_dfs_mt methods                *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_forwd_dfs_mt::constrained_forwd_dfs_mt(ternary_list &c,
  forest* cons, forest* arg, forest* trans, forest* res)
  : constrained_dfs_mt(c, cons, arg, trans, res)
{
}

void MEDDLY::constrained_forwd_dfs_mt::saturateHelper(node_handle a, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);

  const dd_edge& mxd = splits[nb.getLevel()];
  if (mxd.getNode() == 0) {
    return;
  }

  const int mxdLevel = mxd.getLevel();
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(arg3F, nb.getLevel(), mxd.getNode(), FULL_ONLY)
    : arg3F->newUnpacked(mxd.getNode(), FULL_ONLY);

  unpacked_node* A = isLevelAbove(nb.getLevel(), consF->getNodeLevel(a))
    ? unpacked_node::newRedundant(consF, nb.getLevel(), a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);

  dd_edge nbdj(resF), newst(resF);

  // indices to explore
  std::deque<int> queue;
  std::vector<bool> waiting(nb.getSize(), false);
  for (int i = 0; i < nb.getSize(); i++) {
    if (nb.down(i) != 0 && Ru->down(i) != 0) {
      queue.push_back(i);
      waiting[i] = true;
    }
  }

  unpacked_node** Rps = new unpacked_node*[Ru->getSize()];
  for (int i = 0; i < Ru->getSize(); i++) {
    if (Ru->down(i) == 0) {
      Rps[i] = nullptr;
    }
    else {
      Rps[i] = (arg3F->getNodeLevel(Ru->down(i)) == -nb.getLevel())
        ? arg3F->newUnpacked(Ru->down(i), SPARSE_ONLY)
        : unpacked_node::newIdentity(arg3F, -nb.getLevel(), i, Ru->down(i), SPARSE_ONLY);
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int i = queue.front();
    queue.pop_front();
    waiting[i] = false;

    MEDDLY_DCASSERT(nb.down(i) != 0 && Rps[i] != nullptr);

    for (int jz = 0; jz < Rps[i]->getSize(); jz++) {
      const int j = Rps[i]->index(jz);
      if (A->down(j) == 0) {
        continue;
      }

      if (Rps[i]->down(jz) != 0) {
        node_handle rec = 0;
        recFire(A->down(j), nb.down(i), Rps[i]->down(jz), rec);
        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == 0) {
          continue;
        }

        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == nb.down(j)) {
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (nb.down(j) == 0) {
          nb.setFull(j, rec);
        }
        else {
          nbdj.set(nb.down(j));  // clobber
          newst.set(rec);     // clobber
          unionOp->computeTemp(nbdj, newst, nbdj);
          updated = (nbdj.getNode() != nb.down(j));
          nb.setFull(j, nbdj);
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

  unpacked_node::Recycle(Ru);
  for (int i = 0; i < Ru->getSize(); i++) {
    if (Rps[i] != nullptr) {
      unpacked_node::Recycle(Rps[i]);
    }
  }
  delete[] Rps;
  unpacked_node::Recycle(A);
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
  ct_entry_key* key = findResult(a, b, r, c);
  if (key == 0) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);
  const int rLevel = arg3F->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), bLevel), aLevel);
  const int size = resF->getLevelSize(level);

  dd_edge Tdj(resF), newst(resF);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(argF, level, b, FULL_ONLY)
    : argF->newUnpacked(b, FULL_ONLY);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  if (ABS(rLevel) < level) {
    // Assume identity reduction
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i = 0; i < size; i++) {
      node_handle t = 0;
      recFire(A->down(i), B->down(i), r, t);
      T->setFull(i, t);
    }
  }
  else {
    //
    // Need to process this level in the MXD.

    // clear out result (important!)
    T->clear(0, size);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node* Ru = (rLevel < 0)
      ? unpacked_node::newRedundant(arg3F, -rLevel, r, SPARSE_ONLY)
      : arg3F->newUnpacked(r, SPARSE_ONLY);

    // loop over mxd "rows"
    for (int iz = 0; iz < Ru->getSize(); iz++) {
      const int i = Ru->index(iz);

      unpacked_node* Rp = isLevelAbove(-level, arg3F->getNodeLevel(Ru->down(iz)))
        ? unpacked_node::newIdentity(arg3F, -level, i, Ru->down(iz), SPARSE_ONLY)
        : arg3F->newUnpacked(Ru->down(iz), SPARSE_ONLY);

      // loop over mxd "columns"
      for (int jz = 0; jz < Rp->getSize(); jz++) {
        const int j = Rp->index(jz);
        if (A->down(j) == 0) {
          continue;
        }

        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle n = 0;
        recFire(A->down(j), B->down(i), Rp->down(jz), n);

        if (n == 0) {
          continue;
        }

        if (T->down(j) == 0) {
          T->setFull(j, n);
          continue;
        }

        // there's new states and existing states; union them.
        Tdj.set(T->down(j));
        newst.set(n);
        unionOp->computeTemp(newst, Tdj, Tdj);
        T->setFull(j, Tdj);
      } // for j

      unpacked_node::Recycle(Rp);
    } // for i

    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(B);

  saturateHelper(a, *T);
  c = resF->createReducedNode(-1, T);

  saveResult(key, a, b, r, c);
}

// ******************************************************************
// *                                                                *
// *                constrained_bckwd_dfs_mt methods                *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_bckwd_dfs_mt::constrained_bckwd_dfs_mt(ternary_list &c,
  forest* cons, forest* arg, forest* trans, forest* res)
  : constrained_dfs_mt(c, cons, arg, trans, res)
{
}

void MEDDLY::constrained_bckwd_dfs_mt::saturateHelper(node_handle a, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);

  const dd_edge& mxd = splits[nb.getLevel()];
  if (mxd.getNode() == 0) {
    return;
  }

  const int mxdLevel = mxd.getLevel();
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(arg3F, nb.getLevel(), mxd.getNode(), SPARSE_ONLY)
    : arg3F->newUnpacked(mxd.getNode(), SPARSE_ONLY);

  unpacked_node* A = isLevelAbove(nb.getLevel(), consF->getNodeLevel(a))
    ? unpacked_node::newRedundant(consF, nb.getLevel(), a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);

  dd_edge nbdi(resF), newst(resF);

  // indices to explore
  std::deque<int> queue;
  std::vector<bool> waiting(nb.getSize(), false);
  for (int j = 0; j < nb.getSize(); j++) {
    if (nb.down(j) != 0) {
      queue.push_back(j);
      waiting[j] = true;
    }
  }

  unpacked_node** Rps = new unpacked_node*[Ru->getSize()];
  for (int iz = 0; iz < Ru->getSize(); iz++) {
    const int i = Ru->index(iz);
    if (A->down(i) == 0) {
      Rps[iz] = nullptr;
    }
    else {
      Rps[iz] = (arg3F->getNodeLevel(Ru->down(iz)) == -nb.getLevel())
        ? arg3F->newUnpacked(Ru->down(iz), FULL_ONLY)
        : unpacked_node::newIdentity(arg3F, -nb.getLevel(), i, Ru->down(iz), FULL_ONLY);
    }
  }

  // explore indexes
  while (!queue.empty()) {
    const int j = queue.front();
    queue.pop_front();
    waiting[j] = false;

    MEDDLY_DCASSERT(nb.down(j) != 0);

    for (int iz = 0; iz < Ru->getSize(); iz++) {
      const int i = Ru->index(iz);
      if (A->down(i) == 0) {
        continue;
      }

      MEDDLY_DCASSERT(Rps[iz] != nullptr);

      if (Rps[iz]->down(j) != 0) {
        node_handle rec = 0;
        recFire(A->down(i), nb.down(j), Rps[iz]->down(j), rec);
        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == 0) {
          continue;
        }

        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == nb.down(i)) {
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (nb.down(i) == 0) {
          nb.setFull(i, rec);
        }
        else {
          nbdi.set(nb.down(i));  // clobber
          newst.set(rec);     // clobber
          unionOp->computeTemp(nbdi, newst, nbdi);
          updated = (nbdi.getNode() != nb.down(i));
          nb.setFull(i, nbdi);
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

  unpacked_node::Recycle(Ru);
  for (int iz = 0; iz < Ru->getSize(); iz++) {
    unpacked_node::Recycle(Rps[iz]);
  }
  delete[] Rps;
  unpacked_node::Recycle(A);
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
  ct_entry_key* key = findResult(a, b, r, c);
  if (key == 0) {
    return;
  }

  // check if mxd and evmdd are at the same level
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);
  const int rLevel = arg3F->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), bLevel), aLevel);
  const int size = resF->getLevelSize(level);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(argF, level, b, FULL_ONLY)
    : argF->newUnpacked(b, FULL_ONLY);

  dd_edge Tdi(resF), newst(resF);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  if (ABS(rLevel) < level) {
    // Assume identity reduction
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i = 0; i < size; i++) {
      node_handle t = 0;
      recFire(A->down(i), B->down(i), r, t);
      T->setFull(i, t);
    }
  }
  else {
    //
    // Need to process this level in the MXD.

    // clear out result (important!)
    T->clear(0, size);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node* Ru = (rLevel < 0)
      ? unpacked_node::newRedundant(arg3F, -rLevel, r, SPARSE_ONLY)
      : arg3F->newUnpacked(r, SPARSE_ONLY);

    // loop over mxd "rows"
    for (int iz = 0; iz < Ru->getSize(); iz++) {
      const int i = Ru->index(iz);
      if (A->down(i) == 0) {
        continue;
      }

      unpacked_node* Rp = isLevelAbove(-level, arg3F->getNodeLevel(Ru->down(iz)))
        ? unpacked_node::newIdentity(arg3F, -level, i, Ru->down(iz), SPARSE_ONLY)
        : arg3F->newUnpacked(Ru->down(iz), SPARSE_ONLY);

      // loop over mxd "columns"
      for (int jz = 0; jz < Rp->getSize(); jz++) {
        const int j = Rp->index(jz);

        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        node_handle n = 0;
        recFire(A->down(i), B->down(j), Rp->down(jz), n);

        if (n == 0) {
          continue;
        }

        if (T->down(i) == 0) {
          T->setFull(i, n);
          continue;
        }

        // there's new states and existing states; union them.
        Tdi.set(T->down(i));
        newst.set(n);
        unionOp->computeTemp(newst, Tdi, Tdi);
        T->setFull(i, Tdi);
      } // for j

      unpacked_node::Recycle(Rp);
    } // for i

    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(B);

  saturateHelper(a, *T);
  c = resF->createReducedNode(-1, T);

  saveResult(key, a, b, r, c);
}

// ******************************************************************
// *                                                                *
// *               constrained_saturation_mt  methods               *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_saturation_mt::constrained_saturation_mt(
    constrained_dfs_mt* p, forest* cons, forest* arg, forest* res)
  : operation("constr_sat", 1)
{
    // TBD: throw here instead
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

  ct_entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new ct_entry_type("constrained_saturation_mt", "NNI:N");
    et->setForestForSlot(0, cons);
    et->setForestForSlot(1, arg);
    et->setForestForSlot(4, res);
  } else {
    et = new ct_entry_type("constrained_saturation_mt", "NN:N");
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

MEDDLY::ct_entry_key* MEDDLY::constrained_saturation_mt::findResult(
    node_handle a, node_handle b, int level, node_handle &c)
{
  ct_entry_key* key = CT0->useEntryKey(etype[0], 0);
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

void MEDDLY::constrained_saturation_mt::saveResult(ct_entry_key* key,
  node_handle a, node_handle b, int level, node_handle c)
{
  CTresult[0].reset();
  CTresult[0].writeN(c);
  CT0->addEntry(key, CTresult[0]);
}

void MEDDLY::constrained_saturation_mt::saturate(node_handle a, node_handle b, node_handle& c)
{
  saturate(a, b, argF->getMaxLevelIndex(), c);
}

void MEDDLY::constrained_saturation_mt::saturate(node_handle a, node_handle b, int level, node_handle& c)
{
  if (checkTerminals(a, b, c)) {
    return;
  }

  ct_entry_key* key = findResult(a, b, level, c);
  if (key == 0) {
    return;
  }

  const int sz = argF->getLevelSize(level);
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);

  unpacked_node* A = (aLevel < level)
    ? unpacked_node::newRedundant(consF, level, a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);
  unpacked_node* B = (bLevel < level)
    ? unpacked_node::newRedundant(argF, level, b, FULL_ONLY)
    : argF->newUnpacked(b, FULL_ONLY);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
    if (A->down(i) == 0) {
      MEDDLY_DCASSERT(resF == argF);
      T->setFull(i, resF->linkNode(B->down(i)));
    }
    else {
      node_handle t = 0;
      saturate(A->down(i), B->down(i), level - 1, t);
      T->setFull(i, t);
    }
  }

  // Cleanup
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(B);

  parent->saturateHelper(a, *T);
  c = resF->createReducedNode(-1, T);

  // save in compute table
  saveResult(key, a, b, level, c);
}

// ******************************************************************
// *                                                                *
// *              constrained_bckwd_dfs_evplus methods              *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_bckwd_dfs_evplus::constrained_bckwd_dfs_evplus(
    ternary_list &c, forest* cons, forest* arg, forest* trans, forest* res)
  : common_constrained(c, 1, cons, arg, trans, res)
{
  MEDDLY_DCASSERT(cons->isEVPlus() && !cons->isForRelations());
  MEDDLY_DCASSERT(arg->isEVPlus() && !arg->isForRelations());
  MEDDLY_DCASSERT(trans->isMultiTerminal() && trans->isForRelations());
  MEDDLY_DCASSERT(res->isEVPlus() && !res->isForRelations());

  mxdIntersectionOp = INTERSECTION(arg3F, arg3F, arg3F);
  mxdDifferenceOp = DIFFERENCE(arg3F, arg3F, arg3F);
  minOp = UNION(resF, resF, resF);

  splits = nullptr;

  ct_entry_type* et = new ct_entry_type(c.getName(), "LNNN:LN");
  et->setForestForSlot(1, cons);
  et->setForestForSlot(2, arg);
  et->setForestForSlot(3, trans);
  et->setForestForSlot(6, res);
  registerEntryType(0, et);
  buildCTs();
}

MEDDLY::ct_entry_key* MEDDLY::constrained_bckwd_dfs_evplus::findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle r, long& cev, node_handle &c)
{
  ct_entry_key* key = CT0->useEntryKey(etype[0], 0);
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

void MEDDLY::constrained_bckwd_dfs_evplus::saveResult(ct_entry_key* key,
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
void MEDDLY::constrained_bckwd_dfs_evplus::splitMxd(const dd_edge& mxd)
{
  MEDDLY_DCASSERT(arg3F);
  MEDDLY_DCASSERT(splits == nullptr);

  splits = new dd_edge[arg3F->getNumVariables() + 1];

  // we'll be unlinking later, so...
  dd_edge root(mxd), maxDiag(arg3F), Rpdi(arg3F);

  // Build from top down
  for (int level = arg3F->getMaxLevelIndex(); level > 0; level--) {
    splits[level].attach(arg3F);
    if (root.getNode() == 0) {
      // common and easy special case
      continue;
    }

    int mxdLevel = root.getLevel();
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

void MEDDLY::constrained_bckwd_dfs_evplus::compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res)
{
  MEDDLY_DCASSERT(res.getForest() == resF);

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

void MEDDLY::constrained_bckwd_dfs_evplus::_compute(int aev, node_handle a, int bev, node_handle b, node_handle r,
  long& cev, node_handle& c)
{
  // Execute saturation operation
  constrained_saturation_evplus* bckwdSatOp = new constrained_saturation_evplus(this, consF, argF, resF);
  bckwdSatOp->saturate(aev, a, bev, b, cev, c);

  // Cleanup
//  bckwdSatOp->removeAllComputeTableEntries();
  //delete bckwdSatOp;
//  removeAllComputeTableEntries();
  delete[] splits;
  splits = nullptr;
}

void MEDDLY::constrained_bckwd_dfs_evplus::saturateHelper(long aev, node_handle a, unpacked_node& nb)
{
  MEDDLY_DCASSERT(a != 0);

  const dd_edge& mxd = splits[nb.getLevel()];
  if (mxd.getNode() == 0) {
    return;
  }

  const int mxdLevel = mxd.getLevel();
  MEDDLY_DCASSERT(ABS(mxdLevel) == nb.getLevel());

  // Initialize mxd readers, note we might skip the unprimed level
  unpacked_node* Ru = (mxdLevel < 0)
    ? unpacked_node::newRedundant(arg3F, nb.getLevel(), mxd.getNode(), SPARSE_ONLY)
    : arg3F->newUnpacked(mxd.getNode(), SPARSE_ONLY);

  unpacked_node* A = isLevelAbove(nb.getLevel(), consF->getNodeLevel(a))
    ? unpacked_node::newRedundant(consF, nb.getLevel(), 0L, a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);

  // indices to explore
  std::deque<int> queue;
  std::vector<bool> waiting(nb.getSize(), false);
  for (int j = 0; j < nb.getSize(); j++) {
    if (nb.down(j) != 0) {
      queue.push_back(j);
      waiting[j] = true;
    }
  }

  unpacked_node** Rps = new unpacked_node*[Ru->getSize()];
  for (int iz = 0; iz < Ru->getSize(); iz++) {
    const int i = Ru->index(iz);
    if (A->down(i) == 0) {
      Rps[iz] = nullptr;
    }
    else {
      Rps[iz] = (arg3F->getNodeLevel(Ru->down(iz)) == -nb.getLevel())
        ? arg3F->newUnpacked(Ru->down(iz), FULL_ONLY)
        : unpacked_node::newIdentity(arg3F, -nb.getLevel(), i, Ru->down(iz), FULL_ONLY);
    }
  }

  dd_edge nbdi(resF), newst(resF);

  // explore indexes
  while (!queue.empty()) {
    const int j = queue.front();
    queue.pop_front();
    waiting[j] = false;

    MEDDLY_DCASSERT(nb.down(j) != 0);

    for (int iz = 0; iz < Ru->getSize(); iz++) {
      const unsigned i = Ru->index(iz);
      if (A->down(i) == 0) {
        MEDDLY_DCASSERT(A->edge_long(i) == 0);
        continue;
      }

      MEDDLY_DCASSERT(Rps[iz] != nullptr);

      if (Rps[iz]->down(j) != 0) {
        long recev = 0;
        node_handle rec = 0;
        recFire(aev + A->edge_long(i), A->down(i), nb.edge_long(j), nb.down(j), Rps[iz]->down(j), recev, rec);
        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == 0) {
          MEDDLY_DCASSERT(recev == 0);
          continue;
        }

        MEDDLY_DCASSERT(isLevelAbove(nb.getLevel(), resF->getNodeLevel(rec)));

        if (rec == nb.down(i)) {
          // Compute the minimum
          if (recev < nb.edge_long(i)) {
            nb.setEdgeval(i, recev);
          }
          resF->unlinkNode(rec);
          continue;
        }

        bool updated = true;

        if (nb.down(i) == 0) {
          MEDDLY_DCASSERT(nb.edge_long(i) == 0);
          nb.setFull(i, edge_value(recev), rec);
        }
        else {
          nbdi.set(nb.down(i), nb.edge_long(i));
          newst.set(rec, recev);
          minOp->computeTemp(nbdi, newst, nbdi);
          updated = (nbdi.getNode() != nb.down(i));
          nb.setFull(i, nbdi);
/*
          long accev = Inf<long>();
          node_handle acc = 0;
          minOp->computeTemp(nb.edge_long(i), nb.down(i), recev, rec, accev, acc);
          resF->unlinkNode(rec);
          if (acc != nb.down(i)) {
            resF->unlinkNode(nb.down(i));
            nb.setEdge(i, accev);
            nb.d_ref(i) = acc;
          }
          else {
            MEDDLY_DCASSERT(accev == nb.edge_long(i));
            resF->unlinkNode(acc);
            updated = false;
          }
*/
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
        } // if updated

      }
    }
  }

  unpacked_node::Recycle(Ru);
  for (int iz = 0; iz < Ru->getSize(); iz++) {
    unpacked_node::Recycle(Rps[iz]);
  }
  delete[] Rps;
  unpacked_node::Recycle(A);
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
  ct_entry_key* key = findResult(aev, a, bev, b, r, cev, c);
  if (key == 0) {
    MEDDLY_DCASSERT(cev >= 0);
    return;
  }

  // check if mxd and evmdd are at the same level
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);
  const int rLevel = arg3F->getNodeLevel(r);
  const int level = MAX(MAX(ABS(rLevel), bLevel), aLevel);
  const int size = resF->getLevelSize(level);

  // Initialize evmdd reader
  unpacked_node* A = isLevelAbove(level, aLevel)
    ? unpacked_node::newRedundant(consF, level, 0L, a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);
  unpacked_node* B = isLevelAbove(level, bLevel)
    ? unpacked_node::newRedundant(argF, level, 0L, b, FULL_ONLY)
    : argF->newUnpacked(b, FULL_ONLY);

  unpacked_node* T = unpacked_node::newFull(resF, level, size);
  if (ABS(rLevel) < level) {
    // Assume identity reduction
    // Skipped levels in the MXD,
    // that's an important special case that we can handle quickly.

    for (int i = 0; i < size; i++) {
      long tev = 0;
      node_handle t = 0;
      recFire(aev + A->edge_long(i), A->down(i), bev + B->edge_long(i), B->down(i), r, tev, t);
      T->setFull(i, edge_value(tev), t);
    }
  }
  else {
    //
    // Need to process this level in the MXD.

    // clear out result (important!)
    T->clear(0, size);

    dd_edge Tdi(resF), newst(resF);

    // Initialize mxd readers, note we might skip the unprimed level
    unpacked_node* Ru = (rLevel < 0)
      ? unpacked_node::newRedundant(arg3F, -rLevel, r, SPARSE_ONLY)
      : arg3F->newUnpacked(r, SPARSE_ONLY);

    // loop over mxd "rows"
    for (int iz = 0; iz < Ru->getSize(); iz++) {
      const int i = Ru->index(iz);
      if (A->down(i) == 0) {
        MEDDLY_DCASSERT(A->edge_long(i) == 0);
        continue;
      }

      unpacked_node* Rp = isLevelAbove(-level, arg3F->getNodeLevel(Ru->down(iz)))
        ? unpacked_node::newIdentity(arg3F, -level, i, Ru->down(iz), SPARSE_ONLY)
        : arg3F->newUnpacked(Ru->down(iz), SPARSE_ONLY);

      // loop over mxd "columns"
      for (int jz = 0; jz < Rp->getSize(); jz++) {
        int j = Rp->index(jz);

        // ok, there is an i->j "edge".
        // determine new states to be added (recursively)
        // and add them
        long nev = 0;
        node_handle n = 0;
        recFire(aev + A->edge_long(i), A->down(i), bev + B->edge_long(j), B->down(j), Rp->down(jz), nev, n);

        if (n == 0) {
          MEDDLY_DCASSERT(nev == 0);
          continue;
        }

        if (T->down(i) == 0) {
          MEDDLY_DCASSERT(T->edge_long(i) == 0);
          T->setFull(i, edge_value(nev), n);
          continue;
        }

        // there's new states and existing states; union them.
        /*
        const node_handle oldi = T->down(i);
        long newev = Inf<long>();
        node_handle newstates = 0;
        minOp->computeTemp(nev, n, T->edge_long(i), oldi, newev, newstates);
        T->setEdge(i, newev);
        T->d_ref(i) = newstates;

        resF->unlinkNode(oldi);
        resF->unlinkNode(n);
*/

        // there's new states and existing states; union them.
        Tdi.set(T->down(i), T->edge_long(i));
        newst.set(n, nev);
        minOp->computeTemp(newst, Tdi, Tdi);
        T->setFull(i, Tdi);
      } // for j

      unpacked_node::Recycle(Rp);
    } // for i

    unpacked_node::Recycle(Ru);
  } // else

  // cleanup mdd reader
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(B);

  saturateHelper(aev, a, *T);
  resF->createReducedNode(-1, T, cev, c);
  MEDDLY_DCASSERT(cev >= 0);

  saveResult(key, aev, a, bev, b, r, cev, c);
}


// ******************************************************************
// *                                                                *
// *             constrained_saturation_evplus  methods             *
// *                                                                *
// ******************************************************************

MEDDLY::constrained_saturation_evplus::constrained_saturation_evplus(
    constrained_bckwd_dfs_evplus* p, forest* cons, forest* arg, forest* res)
  : operation("constr_sat_evplus", 1)
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

  ct_entry_type* et;

  if (argF->isFullyReduced()) {
    // CT entry includes level info
    et = new ct_entry_type("constrained_saturation_evplus", "LNNI:LN");
    et->setForestForSlot(1, consF);
    et->setForestForSlot(2, argF);
    et->setForestForSlot(6, resF);
  } else {
    et = new ct_entry_type("constrained_saturation_evplus", "LNN:LN");
    et->setForestForSlot(1, consF);
    et->setForestForSlot(2, arg);
    et->setForestForSlot(5, resF);
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

MEDDLY::ct_entry_key* MEDDLY::constrained_saturation_evplus::findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c)
{
  ct_entry_key* key = CT0->useEntryKey(etype[0], 0);
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

void MEDDLY::constrained_saturation_evplus::saveResult(ct_entry_key* key,
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
  saturate(aev, a, bev, b, argF->getMaxLevelIndex(), cev, c);
}

void MEDDLY::constrained_saturation_evplus::saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c)
{
  if (checkTerminals(aev, a, bev, b, cev, c)) {
    return;
  }

  ct_entry_key* key = findResult(aev, a, bev, b, level, cev, c);
  if (key == 0) {
    return;
  }

  const int sz = argF->getLevelSize(level);
  const int aLevel = consF->getNodeLevel(a);
  const int bLevel = argF->getNodeLevel(b);

  unpacked_node* A = (aLevel < level)
    ? unpacked_node::newRedundant(consF, level, 0L, a, FULL_ONLY)
    : consF->newUnpacked(a, FULL_ONLY);
  unpacked_node* B = (bLevel < level)
    ? unpacked_node::newRedundant(argF, level, 0L, b, FULL_ONLY)
    : argF->newUnpacked(b, FULL_ONLY);

  // Do computation
  unpacked_node* T = unpacked_node::newFull(resF, level, sz);
  for (int i = 0; i < sz; i++) {
    if (A->down(i) == 0) {
      MEDDLY_DCASSERT(resF == argF);
      T->setFull(i, B->edgeval(i), resF->linkNode(B->down(i)));
    }
    else {
      long tev = Inf<long>();
      node_handle t = 0;
      saturate(aev + A->edge_long(i), A->down(i), B->edge_long(i), B->down(i), level - 1, tev, t);
      T->setFull(i, edge_value(tev), t);
    }
  }

  // Cleanup
  unpacked_node::Recycle(A);
  unpacked_node::Recycle(B);

  parent->saturateHelper(aev, a, *T);
  resF->createReducedNode(-1, T, cev, c);
  cev += bev;

  // save in compute table
  saveResult(key, aev, a, bev, b, level, cev, c);
}


// ******************************************************************

/*
class MEDDLY::common_constrained: public specialized_operation
{
protected:
  forest* consF;
  forest* argF;
  forest* transF;
  forest* resF;

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

public:
  common_constrained(constrained_opname* code, unsigned slots,
    forest* cons, forest* arg, forest* trans, forest* res);
  ~common_constrained();
};

class MEDDLY::constrained_bfs_opname : public constrained_opname {
protected:
  bool forward;

public:
  constrained_bfs_opname(bool fwd);

  virtual specialized_operation* buildOperation(arguments* a);
};
*/


// ******************************************************************
// *                                                                *
// *                     constraint_bfs_opname                      *
// *                                                                *
// ******************************************************************

/*
class MEDDLY::constrained_dfs_opname : public constrained_opname {
protected:
  bool forward;

public:
  constrained_dfs_opname(bool fwd);

  virtual specialized_operation* buildOperation(arguments* a);
};

MEDDLY::constrained_bfs_opname::constrained_bfs_opname(bool fwd)
 : constrained_opname(fwd ? "PostImage" : "PreImage")
{
  forward = fwd;
}

MEDDLY::specialized_operation* MEDDLY::constrained_bfs_opname::buildOperation(arguments* a)
{
  constrained_opname::constrained_args* args = dynamic_cast<constrained_opname::constrained_args*>(a);
  specialized_operation* op = 0;
  if (forward) {
    throw error(error::NOT_IMPLEMENTED);
  }
  else {
    op = new constrained_bckwd_bfs_evplus(this,
      args->consForest,
      args->inForest,
      args->relForest,
      args->outForest);
  }
  return op;
}
*/


// ******************************************************************
// *                                                                *
// *                     constraint_dfs_opname                      *
// *                                                                *
// ******************************************************************

/*
MEDDLY::constrained_dfs_opname::constrained_dfs_opname(bool fwd)
 : constrained_opname(fwd ? "PostImage" : "PreImage")
{
  forward = fwd;
}

MEDDLY::specialized_operation* MEDDLY::constrained_dfs_opname::buildOperation(arguments* a)
{
  constrained_opname::constrained_args* args = dynamic_cast<constrained_opname::constrained_args*>(a);
  specialized_operation* op = 0;
  if (forward) {
    if (args->consForest->isMultiTerminal() && args->inForest->isMultiTerminal()
      && args->relForest->isMultiTerminal() && args->outForest->isMultiTerminal()) {
      op = new constrained_forwd_dfs_mt(this,
        args->consForest,
        args->inForest,
        args->relForest,
        args->outForest);
    }
    else {
      throw error(error::NOT_IMPLEMENTED);
    }
  }
  else {
    if (args->consForest->isMultiTerminal() && args->inForest->isMultiTerminal()
      && args->relForest->isMultiTerminal() && args->outForest->isMultiTerminal()) {
      op = new constrained_bckwd_dfs_mt(this,
        args->consForest,
        args->inForest,
        args->relForest,
        args->outForest);
    }
    else if (args->consForest->isEVPlus() && args->inForest->isEVPlus()
      && args->relForest->isMultiTerminal() && args->outForest->isEVPlus()) {
      op = new constrained_bckwd_dfs_evplus(this,
        args->consForest,
        args->inForest,
        args->relForest,
        args->outForest);
    }
    else {
      throw error(error::NOT_IMPLEMENTED);
    }
  }
  return op;
}
*/




// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

/*
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
*/

