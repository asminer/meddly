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

#ifndef MEDDLY_EVMXD_H
#define MEDDLY_EVMXD_H

#include "ev.h"
#include "../oper_binary.h"

namespace MEDDLY {
  class evmxd_forest;
};

class MEDDLY::evmxd_forest : public ev_forest {
  public:
    evmxd_forest(domain* d, range_type t, edge_labeling ev, const policies &p);

    virtual void reorderVariables(const int* level2var);
    virtual void swapAdjacentVariables(int level);
    virtual void moveDownVariable(int high, int low);
    virtual void moveUpVariable(int low, int high);

#ifdef ALLOW_DEPRECATED_0_17_7
  protected:
    template <class OPERATION, typename TYPE>
    inline void evaluateT(const dd_edge &f, const int* vlist,
      const int* vplist, TYPE &val) const
    {
      if (!f.isAttachedTo(this))
          throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
      if (vlist == 0)
          throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);
      if (vplist == 0)
          throw error(error::INVALID_VARIABLE, __FILE__, __LINE__);

      // Assumptions:
      // (1) vlist and vplist do not contain special values (-1, -2, etc).
      // (2) vlist and vplist contains a single element.
      node_handle node = f.getNode();
      f.getEdgeValue(val);

      while (!isTerminalNode(node)) {
        TYPE ev;
        int k = getNodeLevel(node);
        getDownPtr(node, ((k > 0) ? vlist[k] : vplist[-k]), ev, node);
        val = (node) ? OPERATION::apply(val, ev) : ev;
      }
    }

  public:
    /// Special case for createEdge(), with only one minterm.
    template <class OPERATION, typename TYPE>
    inline void
    createEdgePath(int k, const int* vlist, const int* vplist,
      TYPE &ev, node_handle &ed)
    {
      if (0==ed) return;

      for (int i=1; i<=k; i++) {
        if (DONT_CHANGE == vplist[i]) {
          //
          // Identity node
          //
          MEDDLY_DCASSERT(DONT_CARE == vlist[i]);
          if (isIdentityReduced()) continue;
          // Build an identity node by hand
          unsigned sz = unsigned(getLevelSize(i));
          unpacked_node* nb = unpacked_node::newFull(this, i, sz);
          for (unsigned v=0; v<sz; v++) {
            unpacked_node* nbp = unpacked_node::newSparse(this, -i, 1);
            nbp->setSparse(0, v, ev, linkNode(ed));
            TYPE pev;
            node_handle pd;
            createReducedNode(v, nbp, pev, pd);
            nb->setFull(v, pev, pd);
          }
          unlinkNode(ed);
          createReducedNode(-1, nb, ev, ed);
          continue;
        }
        //
        // process primed level
        //
        node_handle edpr;
        TYPE evpr;
        if (DONT_CARE == vplist[i]) {
          if (isFullyReduced()) {
            // DO NOTHING
            edpr = ed;
            evpr = ev;
          } else {
            // build redundant node
            int sz = getLevelSize(-i);
            unpacked_node* nb = unpacked_node::newFull(this, -i, sz);
            for (int v=0; v<sz; v++) {
              nb->setFull(v, ev, linkNode(ed));
            }
            unlinkNode(ed);
            createReducedNode(-1, nb, evpr, edpr);
          }
        } else {
          // sane value
          unpacked_node* nb = unpacked_node::newSparse(this, -i, 1);
          nb->setSparse(0, vplist[i], ev, ed);
          createReducedNode(vlist[i], nb, evpr, edpr);
        }
        //
        // process unprimed level
        //
        if (DONT_CARE == vlist[i]) {
          if (isFullyReduced()) {
            ed = edpr;
            ev = evpr;
            continue;
          }
          // build redundant node
          int sz = getLevelSize(i);
          unpacked_node *nb = unpacked_node::newFull(this, i, sz);
          if (isIdentityReduced()) {
            // Below is likely a singleton, so check for identity reduction
            // on the appropriate v value
            for (int v=0; v<sz; v++) {
              node_handle dpr = (v == vplist[i]) ? ed : edpr;
              nb->setFull(v, evpr, linkNode(dpr));
            }
          } else {
            // Doesn't matter what happened below
            for (int v=0; v<sz; v++) {
              nb->setFull(v, evpr, linkNode(edpr));
            }
          }
          unlinkNode(edpr);
          createReducedNode(-1, nb, ev, ed);
        } else {
          // sane value
          unpacked_node* nb = unpacked_node::newSparse(this, i, 1);
          nb->setSparse(0, vlist[i], evpr, edpr);
          createReducedNode(-1, nb, ev, ed);
        }
      } // for i
    }

#endif  // ALLOW_DEPRECATED_0_17_7

};


//
// Helper class for createEdge
//

#ifdef ALLOW_DEPRECATED_0_17_7

namespace MEDDLY {

  template <class OPERATION, typename T>
  class evmxd_edgemaker {
      evmxd_forest* F;
      const int* const* vulist;
      const int* const* vplist;
      const T* values;
      int* order;
      int N;
      int K;
      binary_operation* unionOp;
    public:
      evmxd_edgemaker(evmxd_forest* f,
        const int* const* mt, const int* const* mp, const T* v,
        int* o, int n, int k, binary_operation* unOp)
      {
        F = f;
        vulist = mt;
        vplist = mp;
        values = v;
        order = o;
        N = n;
        K = k;
        unionOp = unOp;
      }

      inline const int* unprimed(int i) const {
        MEDDLY_CHECK_RANGE(0, i, N);
        return vulist[order[i]];
      }
      inline int unprimed(int i, int k) const {
        MEDDLY_CHECK_RANGE(0, i, N);
        MEDDLY_CHECK_RANGE(1, k, K+1);
        return vulist[order[i]][k];
      }
      inline const int* primed(int i) const {
        MEDDLY_CHECK_RANGE(0, i, N);
        return vplist[order[i]];
      }
      inline int primed(int i, int k) const {
        MEDDLY_CHECK_RANGE(0, i, N);
        MEDDLY_CHECK_RANGE(1, k, K+1);
        return vplist[order[i]][k];
      }
      inline T term(int i) const {
        MEDDLY_CHECK_RANGE(0, i, N);
        return values ? values[order[i]]: 1;
      }
      inline void swap(int i, int j) {
        MEDDLY_CHECK_RANGE(0, i, N);
        MEDDLY_CHECK_RANGE(0, j, N);
        MEDDLY::SWAP(order[i], order[j]);
      }

      inline void createEdge(T &ev, node_handle &ed) {
        createEdgeUn(K, 0, N, ev, ed);
      }

    protected:

      /**
          Recursive implementation of createEdge(),
          unprimed levels, for use by evmxd_forest descendants.
      */
      void createEdgeUn(int k, int start, int stop, T &ev, node_handle &ed) {
        MEDDLY_DCASSERT(k>=0);
        MEDDLY_DCASSERT(stop > start);
        //
        // Fast special case
        //
        if (1==stop-start) {
          ev = term(start);
          terminal t;
          t.setOmega();
          ed = t.getHandle();
          F->createEdgePath<OPERATION, T>(k, unprimed(start), primed(start),
            ev, ed);
          return;
        }
        //
        // Check terminal case
        //
        if (0==k) {
          ev = term(start);
          for (int i=start+1; i<stop; i++) {
            OPERATION::unionEq(ev, term(i));
          }
          terminal t;
          t.setOmega();
          ed = t.getHandle();
          return;
        }

        // size of variables at level k
        unsigned lastV = unsigned(F->getLevelSize(k));
        // index of end of current batch
        int batchP = start;

        //
        // Move any "don't cares" to the front, and process them
        //
        unsigned nextV = lastV;
        for (int i=start; i<stop; i++) {
          if (DONT_CARE == unprimed(i, k)) {
            if (batchP != i) {
              swap(batchP, i);
            }
            batchP++;
          } else {
            MEDDLY_DCASSERT(unprimed(i, k) >= 0);
            nextV = MIN(nextV, unsigned(unprimed(i, k)));
          }
        }

        //
        // Move any "don't changes" below the "don't cares", to the front,
        // and process them to construct a new level-k node.
        int dch = start;
        for (int i=start; i<batchP; i++) {
          if (DONT_CHANGE == primed(i, k)) {
            if (dch != i) {
              swap(dch, i);
            }
            dch++;
          }
        }

        //
        // Process "don't care, don't change" pairs, if any
        //
        dd_edge dontcare(F);
        F->getTransparentEdge(dontcare);

        if (dch > start) {
          T dc_ev;
          node_handle dc;
          createEdgeUn(k-1, start, dch, dc_ev, dc);
          makeIdentityEdge(k, dc_ev, dc);
          dontcare.set(dc_ev, dc);
          // done with those
          start = dch;
        }

        //
        // Process "don't care, ordinary" pairs, if any
        // (producing a level-k node)
        //
        if (batchP > start) {
          T dcnormal_ev;
          node_handle dcnormal_nh;
          createEdgePr(-1, -k, start, batchP, dcnormal_ev, dcnormal_nh);
          F->makeNodeAtLevel<OPERATION, T>(k, dcnormal_ev, dcnormal_nh);
          dd_edge dcnormal(F);
          dcnormal.set(dcnormal_ev, dcnormal_nh);

          MEDDLY_DCASSERT(unionOp);
          unionOp->computeTemp(dontcare, dcnormal, dontcare);
        }

        //
        // Start new node at level k
        //
        unpacked_node* nb = unpacked_node::newSparse(F, k, lastV);
        unsigned z = 0; // number of nonzero edges in our sparse node

        //
        // For each value v,
        //  (1) move those values to the front
        //  (2) process them, if any
        // Then when we are done, union with any don't cares
        //
        for (unsigned v=nextV; v<lastV; v=nextV) {
          nextV = lastV;
          //
          // neat trick!
          // shift the array over, because we're done with the previous batch
          //
          start = batchP;

          //
          // (1) move anything with value v, to the "new" front
          //
          for (int i=start; i<stop; i++) {
            if (v == unprimed(i, k)) {
              if (batchP != i) {
                swap(batchP, i);
              }
              batchP++;
            } else {
              nextV = MIN(nextV, unsigned(unprimed(i, k)));
            }
          }

          //
          // (2) recurse if necessary
          //
          if (0==batchP) continue;
          T these_ev;
          node_handle these_ptr;
          createEdgePr(v, -k, start, batchP, these_ev, these_ptr);
          nb->setSparse(z, v, these_ev, these_ptr);
          z++;
        } // for v

        //
        // Union with don't cares
        //
        MEDDLY_DCASSERT(unionOp);
        node_handle built_nh;
        T built_ev;
        nb->shrink(z);
        F->createReducedNode(-1, nb, built_ev, built_nh);
        dd_edge built(F);
        built.set(built_ev, built_nh);

        unionOp->computeTemp(dontcare, built, built);
        ed = F->linkNode(built);
        built.getEdgeValue(ev);
      };

    protected:

      /**
          Recursive implementation of createEdge(),
          primed levels
      */
      void createEdgePr(int in, int k, int start, int stop,
        T &ev, node_handle &ed) {
        MEDDLY_DCASSERT(k<0);
        MEDDLY_DCASSERT(stop > start);

        //
        // Don't need to check for terminals
        //

        // size of variables at level k
        unsigned lastV = unsigned(F->getLevelSize(k));
        // current batch size
        int batchP = start;

        //
        // Move any "don't cares" to the front, and process them
        //
        unsigned nextV = lastV;
        for (int i=start; i<stop; i++) {
          if (DONT_CARE == primed(i, -k)) {
            if (batchP != i) {
              swap(batchP, i);
            }
            batchP++;
          } else {
            nextV = MIN(nextV, unsigned(primed(i, -k)));
          }
        }

        dd_edge dontcare(F);

        if (batchP > start) {
          node_handle dc_ptr;
          T dc_val;
          createEdgeUn(-k-1, start, batchP, dc_val, dc_ptr);
          dontcare.set(dc_val, dc_ptr);
        } else {
          F->getTransparentEdge(dontcare);
        }

        //
        // Start new node at level k
        //
        unpacked_node* nb = unpacked_node::newSparse(F, k, lastV);
        unsigned z = 0; // number of nonzero edges in our sparse node

        //
        // For each value v,
        //  (1) move those values to the front
        //  (2) process them, if any
        //  (3) union with don't cares
        //
        dd_edge these(F);

        for (unsigned v = (dontcare.getNode()) ? 0 : nextV;
             v<lastV;
             v = (dontcare.getNode()) ? v+1 : nextV)
        {
          nextV = lastV;
          //
          // neat trick!
          // shift the array over, because we're done with the previous batch
          //
          start = batchP;

          //
          // (1) move anything with value v, or don't change if v=in,
          //     to the "new" front
          //
          bool veqin = (v==in);
          for (int i=start; i<stop; i++) {
            if (v == primed(i, -k) || (veqin && DONT_CHANGE==primed(i, -k))) {
              if (batchP != i) {
                swap(batchP, i);
              }
              batchP++;
            } else {
              nextV = MIN(nextV, unsigned(primed(i, -k)));
            }
          }

          //
          // (2) recurse if necessary
          //
          if (batchP > start) {
            node_handle these_ptr;
            T these_val;
            createEdgeUn(-k-1, start, batchP, these_val, these_ptr);
            these.set(these_val, these_ptr);
          } else {
            F->getTransparentEdge(dontcare);
          }

          //
          // (3) union with don't cares
          //
          MEDDLY_DCASSERT(unionOp);
          unionOp->computeTemp(dontcare, these, these);

          //
          // add to sparse node, unless empty
          //
          if (0==these.getNode()) continue;
          nb->setSparse(z, v, these.getEdgeValue(), F->linkNode(these));
          z++;
        } // for v

        //
        // Cleanup
        //
        nb->shrink(z);
        F->createReducedNode(in, nb, ev, ed);
      };


      //
      //
      // Helper for createEdge
      //
      inline void makeIdentityEdge(int k, T& pev, node_handle& p) {
        if (F->isIdentityReduced()) return;

        // build an identity node by hand
        unsigned lastV = unsigned(F->getLevelSize(k));
        unpacked_node* nb = unpacked_node::newFull(F, k, lastV);
        for (unsigned v=0; v<lastV; v++) {
          unpacked_node* nbp = unpacked_node::newSparse(F, -k, 1);
          nbp->setSparse(0, v, pev, F->linkNode(p));
          node_handle pr;
          T pr_ev;
          F->createReducedNode(v, nbp, pr_ev, pr);
          nb->setFull(v, pr_ev, pr);
        } // for v
        F->unlinkNode(p);
        F->createReducedNode(-1, nb, pev, p);
      }



  }; // class evmxd_edgemaker

}; // namespace MEDDLY

#endif  // ALLOW_DEPRECATED_0_17_7

#endif
