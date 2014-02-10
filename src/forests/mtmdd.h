
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

#ifndef MTMDD_H
#define MTMDD_H

#include "mt.h"

namespace MEDDLY {
  class mtmdd_forest;
};

class MEDDLY::mtmdd_forest : public mt_forest {
  public:
    mtmdd_forest(int dsl, domain* d, range_type t, const policies &p);

  protected:
    inline node_handle evaluateRaw(const dd_edge &f, const int* vlist) const {
      node_handle p = f.getNode();
      while (!isTerminalNode(p)) {
        p = getDownPtr(p, vlist[getNodeHeight(p)]);
      }
      return p;
    }
    
  public:
    /// Special case for createEdge(), with only one minterm.
    inline node_handle 
    createEdgePath(int k, const int* vlist, node_handle bottom) 
    {
        if (0==bottom) return bottom;
        for (int i=1; i<=k; i++) {
          if (DONT_CARE == vlist[i]) {
            // make a redundant node
            if (isFullyReduced()) continue; 
            int sz = getLevelSize(i);
            node_builder& nb = useNodeBuilder(i, sz);
            nb.d(0) = bottom;
            for (int v=1; v<sz; v++) {
              nb.d(v) = linkNode(bottom);
            }
            bottom = createReducedNode(-1, nb);
          } else {
            // make a singleton node
            node_builder& nb = useSparseBuilder(i, 1);
            nb.i(0) = vlist[i];
            nb.d(0) = bottom;
            bottom = createReducedNode(-1, nb);
          }
        } // for i
        return bottom;
    }

};


//
// Helper class for createEdge
//

namespace MEDDLY {

  template <class ENCODER, typename T>
  class mtmdd_edgemaker {
      mtmdd_forest* F;
      const int* const* vlist;
      const T* values;
      int* order;
      int N;
      int K;
      binary_operation* unionOp;
    public:
      mtmdd_edgemaker(mtmdd_forest* f, const int* const* mt, const T* v, 
        int* o, int n, int k, binary_operation* unOp) 
      {
        F = f;
        vlist = mt;
        values = v;
        order = o;
        N = n;
        K = k;
        unionOp = unOp;
      }

      inline const int* unprimed(int i) const {
        MEDDLY_CHECK_RANGE(0, i, N);
        return vlist[order[i]];
      }
      inline int unprimed(int i, int k) const {
        MEDDLY_CHECK_RANGE(0, i, N);
        MEDDLY_CHECK_RANGE(1, k, K+1);
        return vlist[order[i]][k];
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

      inline node_handle createEdge() {
        return createEdge(K, 0, N);
      }

      /**
          Recursive implementation of createEdge(),
          for use by mtmdd_forest descendants.
      */
      node_handle createEdge(int k, int start, int stop) {
        MEDDLY_DCASSERT(k>=0);
        MEDDLY_DCASSERT(stop > start);
        // 
        // Fast special case
        //
        if (1==stop-start) {
          return F->createEdgePath(k, unprimed(start),
            ENCODER::value2handle(term(start))
          );
        }
        //
        // Check terminal case
        //
        if (0==k) {
          T accumulate = term(start);
          for (int i=start; i<stop; i++) {
            accumulate += term(i);
          }
          return ENCODER::value2handle(accumulate);
        }

        // size of variables at level k
        int lastV = F->getLevelSize(k);
        // index of end of current batch
        int batchP = start;

        //
        // Move any "don't cares" to the front, and process them
        //
        int nextV = lastV;
        for (int i=start; i<stop; i++) {
          if (DONT_CARE == unprimed(i, k)) {
            if (batchP != i) {
              swap(batchP, i);
            }
            batchP++;
          } else {
            MEDDLY_DCASSERT(unprimed(i, k) >= 0);
            nextV = MIN(nextV, unprimed(i, k));
          }
        }
        node_handle dontcares;
        if (batchP > start) {
          dontcares = createEdge(k-1, start, batchP);
        } else {
          dontcares = 0;
        }

        //
        // Start new node at level k
        //
        node_builder& nb = F->useSparseBuilder(k, lastV);
        int z = 0; // number of nonzero edges in our sparse node

        //
        // For each value v, 
        //  (1) move those values to the front
        //  (2) process them, if any
        //  (3) union with don't cares
        //
        for (int v = (dontcares) ? 0 : nextV; 
             v<lastV; 
             v = (dontcares) ? v+1 : nextV) 
        {
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
              nextV = MIN(nextV, unprimed(i, k));
            }
          }

          //
          // (2) recurse if necessary
          //
          node_handle these;
          if (batchP > start) {
            these = createEdge(k-1, start, batchP);
          } else {
            these = 0;
          }

          //
          // (3) union with don't cares
          //
          MEDDLY_DCASSERT(unionOp);
          node_handle total = unionOp->compute(dontcares, these);
          F->unlinkNode(these);

          //
          // add to sparse node, unless empty
          //
          if (0==total) continue;
          nb.i(z) = v;
          nb.d(z) = total;
          z++;
        } // for v

        //
        // Cleanup
        //
        F->unlinkNode(dontcares);
        nb.shrinkSparse(z);
        return F->createReducedNode(-1, nb);
      };

  }; // class mtmdd_edgemaker

}; // namespace MEDDLY

#endif
