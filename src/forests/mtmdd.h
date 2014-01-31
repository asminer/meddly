
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



// TODO: Testing
// TODO: mdd_node_manager
// TODO: mtmdd_node_manager

// TODO: inPlaceSortBuild() must be modified to deal with don't care and
//       don't change while building the node instead of deal with them
//       separately (before the call to inPlaceSortBuild()).
//       For this purpose, verify that compute_manager::UNION and PLUS
//       work with nodes that are not at the top-level correctly.


/* 
  TODO: ensure this rule
  All extensions must over-ride either reduceNode() or normalizeAndReduceNode().
  Normalizing is only for edge-valued decision diagrams.
*/
  
#ifndef MTMDD_H
#define MTMDD_H

#include "mt.h"

namespace MEDDLY {

  /**
      Base class for all multi-terminal MDDs.
      I.e., everything multi-terminal and not for relations.
  */
  template <class ENCODER>
  class mtmdd_forest : public mt_forest<ENCODER> {
    protected:
      mtmdd_forest(int dsl, domain* d, forest::range_type t, 
        const forest::policies &p) : mt_forest<ENCODER>(dsl, d, false, t, p)
      {
        // nothing to construct
      }

    protected:
      /**
          Template implementation of evaluate().
          Derived classes should call this.
      */
      template <typename T>
      inline void evaluateTempl(const dd_edge &f, const int* vlist, T &term) const
      {
        node_handle p = f.getNode();
        while (!mt_forest<ENCODER>::isTerminalNode(p)) {
          int i = vlist[mt_forest<ENCODER>::getNodeHeight(p)];
          p = mt_forest<ENCODER>::getDownPtr(p, i);
        }
        term = ENCODER::handle2value(p);
      }

      /**
          Recursive template implementation of createEdge(),
          the one that uses an array of minterms.
          Derived classes should call this.
          Note: null array of return values
          corresponds to "all 1s".
      */
      template <typename T>
      inline node_handle createEdgeRT(int k, int** vlist, T* terms, int N) {
        MEDDLY_DCASSERT(k>=0);
        // 
        // Fast special case
        //
        if (1==N) {
          return createEdgePath(k, vlist[0], 
            ENCODER::value2handle(terms ? terms[0] : 1)
          );
        }
        //
        // Check terminal case
        //
        if (0==k) {
          T accumulate;
          if (terms) {
            accumulate = terms[0];
            for (int i=1; i<N; i++) {
              accumulate += terms[i];
            }
          } else {
            accumulate = true;
          }
          return ENCODER::value2handle(accumulate);
        }

        // size of variables at level k
        int lastV = mt_forest<ENCODER>::getLevelSize(k);
        // current batch size
        int batchP = 0;

        //
        // Move any "don't cares" to the front, and process them
        //
        int nextV = lastV;
        for (int i=0; i<N; i++) {
          if (DONT_CARE == vlist[i][k]) {
            if (batchP != i) {
              SWAP(vlist[batchP], vlist[i]);
              if (terms) SWAP(terms[batchP], terms[i]);
            }
            batchP++;
          } else {
            MEDDLY_DCASSERT(vlist[i][k] >= 0);
            nextV = MIN(nextV, vlist[i][k]);
          }
        }
        node_handle dontcares;
        if (batchP) {
          dontcares = createEdgeRT(k-1, vlist, terms, batchP);
        } else {
          dontcares = 0;
        }

        //
        // Start new node at level k
        //
        node_builder& nb = mt_forest<ENCODER>::useSparseBuilder(k, lastV);
        int z = 0; // number of nonzero edges in our sparse node

        //
        // For each value v, 
        //  (1) move those values to the front
        //  (2) process them, if any
        //  (3) union with don't cares
        //
        int v = (dontcares) ? 0 : nextV;
        for (int v=0; v<lastV; v = (dontcares) ? v+1 : nextV) {
          nextV = lastV;
          //
          // neat trick!
          // shift the array over, because we're done with the previous batch
          //
          vlist += batchP;
          if (terms) terms += batchP;
          N -= batchP;
          batchP = 0;

          //
          // (1) move anything with value v, to the "new" front
          //
          for (int i=0; i<N; i++) {
            if (v == vlist[i][k]) {
              if (batchP != i) {
                SWAP(vlist[batchP], vlist[i]);
                if (terms) SWAP(terms[batchP], terms[i]);
              }
              batchP++;
            } else {
              nextV = MIN(nextV, vlist[i][k]);
            }
          }

          //
          // (2) recurse if necessary
          //
          node_handle these;
          if (batchP) {
            these = createEdgeRT(k-1, vlist, terms, batchP);
          } else {
            these = 0;
          }

          //
          // (3) union with don't cares
          //
          MEDDLY_DCASSERT(mt_forest<ENCODER>::unionOp);
          node_handle total = mt_forest<ENCODER>::unionOp->compute(dontcares, these);
          mt_forest<ENCODER>::unlinkNode(these);

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
        mt_forest<ENCODER>::unlinkNode(dontcares);
        nb.shrinkSparse(z);
        return mt_forest<ENCODER>::createReducedNode(-1, nb);
      };
  }; // class

};  // namespace

#endif

