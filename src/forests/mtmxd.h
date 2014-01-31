
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
// TODO: mxd_node_manager
// TODO: mtmxd_node_manager (??)

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
  
#ifndef MTMXD_H
#define MTMXD_H

#include "mt.h"

namespace MEDDLY {

  /**
      Base class for all multi-terminal MxDs.
      I.e., everything multi-terminal and for relations.
  */
  template <class ENCODER>
  class mtmxd_forest : public mt_forest<ENCODER> {
    protected:
      mtmxd_forest(int dsl, domain* d, forest::range_type t, 
        const forest::policies &p) : mt_forest<ENCODER>(dsl, d, true, t, p)
      {
        // nothing to construct
      }

    protected:
      /**
          Template implementation of evaluate().
          Derived classes should call this.
      */
      template <typename T>
      inline void evaluateTempl(const dd_edge &f, const int* vlist, 
        const int* vplist, T &term) const
      {
        node_handle p = f.getNode();
        while (!mt_forest<ENCODER>::isTerminalNode(p)) {
          int k = mt_forest<ENCODER>::getNodeLevel(p);
          int i = (k<0) ? vplist[-k] : vlist[k];
          p = mt_forest<ENCODER>::getDownPtr(p, i);
        } 
        term = ENCODER::handle2value(p);
      }

      /**
          Recursive template implementation of createEdge(),
          the one that uses an array of minterms.
          Derived classes should call this.
            @param  k       Level, should be non-negative
            @param  vlist   Array of "from" minterms
            @param  vplist  Array of "to" minterms
            @param  terms   Array of terminal values, or null
                            to indicate "all true"
            @param  N       Dimension of vlist/vplist/terms arrays.
            @return         New node handle
      */
      template <typename T>
      inline node_handle createEdgeRT(int k, int** vlist, 
        int** vplist, T* terms, int N) 
      {
        MEDDLY_DCASSERT(k>=0);
        // 
        // Fast special case
        //
        if (1==N) {
          return createEdgePath(k, vlist[0], vplist[0], 
            ENCODER::value2handle(terms ? terms[0] : true)
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
              SWAP(vplist[batchP], vplist[i]);
              if (terms) SWAP(terms[batchP], terms[i]);
            }
            batchP++;
          } else {
            MEDDLY_DCASSERT(vlist[i][k] >= 0);
            nextV = MIN(nextV, vlist[i][k]);
          }
        }
        node_handle dontcares = 0;
        //
        // Move any "don't changes" below the "don't cares", to the front,
        // and process them to construct a new level-k node.
        int dch = 0;
        for (int i=0; i<batchP; i++) {
          if (DONT_CHANGE == vplist[i][k]) {
            if (dch != i) {
              SWAP(vlist[dch], vlist[i]);
              SWAP(vplist[dch], vplist[i]);
              if (terms) SWAP(terms[dch], terms[i]);
            }
            dch++;
          }
        } 
        // process "don't care, don't change" pairs
        if (dch) {
          node_handle below = createEdgeRT(k-1, vlist, vplist, terms, dch);
          dontcares = makeIdentityEdge(k, below);
          // done with those
          vlist += dch;
          vplist += dch;
          if (terms) terms += dch;
          N -= dch;
          batchP -= dch;
        }
        


        //
        // process the don't care, ordinary pairs
        // (again, producing a level-k node.)
        //
        if (batchP) {
          node_handle dcnormal = 
            mt_forest<ENCODER>::makeNodeAtLevel(k,
              createEdgeRTpr(-1, -k, vlist, vplist, terms, batchP)
            );

          MEDDLY_DCASSERT(mt_forest<ENCODER>::unionOp);
          node_handle total 
          = mt_forest<ENCODER>::unionOp->compute(dontcares, dcnormal);
          mt_forest<ENCODER>::unlinkNode(dcnormal);
          mt_forest<ENCODER>::unlinkNode(dontcares);
          dontcares = total;
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
        // Then when we are done, union with any don't cares
        //
        for (int v=nextV; v<lastV; v=nextV) {
          nextV = lastV;
          //
          // neat trick!
          // shift the array over, because we're done with the previous batch
          //
          vlist += batchP;
          vplist += batchP;
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
                SWAP(vplist[batchP], vplist[i]);
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
          if (0==batchP) continue;
          nb.i(z) = v;
          nb.d(z) = createEdgeRTpr(v, -k, vlist, vplist, terms, batchP);
          z++;
        } // for v

        //
        // Union with don't cares
        //
        MEDDLY_DCASSERT(mt_forest<ENCODER>::unionOp);
        nb.shrinkSparse(z);
        node_handle built = mt_forest<ENCODER>::createReducedNode(-1, nb);
        node_handle total = mt_forest<ENCODER>::unionOp->compute(dontcares, built);
        mt_forest<ENCODER>::unlinkNode(dontcares);
        mt_forest<ENCODER>::unlinkNode(built);
        return total; 
      } // createEdgeRT

    private:
      /**
          Recursive template implementation of createEdge(),
          the one that uses an array of minterms.
          Helper function - for primed levels.
            @param  in      Index of incoming edge.  
                            Only matters for primed levels.
            @param  k       Level, should be negative
            @param  vlist   Array of "from" minterms
            @param  vplist  Array of "to" minterms
            @param  terms   Array of terminal values, or null
                            to indicate "all true"
            @param  N       Dimension of vlist/vplist/terms arrays.
            @return         New node handle
      */
      template <typename T>
      inline node_handle createEdgeRTpr(int in, int k, int** vlist, 
        int** vplist, T* terms, int N) 
      {
        MEDDLY_DCASSERT(k<0);

        //
        // Don't need to check for terminals
        //

        // size of variables at level k
        int lastV = mt_forest<ENCODER>::getLevelSize(k);
        // current batch size
        int batchP = 0;

        //
        // Move any "don't cares" to the front, and process them
        //
        int nextV = lastV;
        for (int i=0; i<N; i++) {
          if (DONT_CARE == vplist[i][-k]) {
            if (batchP != i) {
              SWAP(vlist[batchP], vlist[i]);
              SWAP(vplist[batchP], vplist[i]);
              if (terms) SWAP(terms[batchP], terms[i]);
            }
            batchP++;
          } else {
            nextV = MIN(nextV, vplist[i][-k]);
          }
        }
        node_handle dontcares;
        if (batchP) {
          dontcares = createEdgeRT(-k-1, vlist, vplist, terms, batchP);
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
          vplist += batchP;
          if (terms) terms += batchP;
          N -= batchP;
          batchP = 0;

          //
          // (1) move anything with value v, or don't change if v=in,
          //     to the "new" front
          //
          bool veqin = (v==in);
          for (int i=0; i<N; i++) {
            if (v == vplist[i][-k] || (veqin && DONT_CHANGE==vplist[i][-k])) {
              if (batchP != i) {
                SWAP(vlist[batchP], vlist[i]);
                SWAP(vplist[batchP], vplist[i]);
                if (terms) SWAP(terms[batchP], terms[i]);
              }
              batchP++;
            } else {
              nextV = MIN(nextV, vplist[i][-k]);
            }
          }

          //
          // (2) recurse if necessary
          //
          node_handle these;
          if (batchP) {
            these = createEdgeRT(-k-1, vlist, vplist, terms, batchP);
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
        return mt_forest<ENCODER>::createReducedNode(in, nb);
      };


      //
      //
      // Helper for createEdgeRT
      //
      inline node_handle makeIdentityEdge(int k, node_handle p) {
        if (mt_forest<ENCODER>::isIdentityReduced()) {
          return p;
        }
        // build an identity node by hand
        int lastV = mt_forest<ENCODER>::getLevelSize(k);
        node_builder& nb = mt_forest<ENCODER>::useNodeBuilder(k, lastV);
        for (int v=0; v<lastV; v++) {
          node_builder& nbp = mt_forest<ENCODER>::useSparseBuilder(-k, 1);
          nbp.i(0) = v;
          nbp.d(0) = mt_forest<ENCODER>::linkNode(p);
          nb.d(v) = mt_forest<ENCODER>::createReducedNode(v, nbp);
        } // for v
        mt_forest<ENCODER>::unlinkNode(p);
        return mt_forest<ENCODER>::createReducedNode(-1, nb);
      }


  }; // class

}; // namespace

#endif

