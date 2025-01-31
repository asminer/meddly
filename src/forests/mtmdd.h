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

#ifndef MEDDLY_MTMDD_H
#define MEDDLY_MTMDD_H

#include "mt.h"

namespace MEDDLY {
  class mtmdd_forest;
};

class MEDDLY::mtmdd_forest : public mt_forest {
  public:
    mtmdd_forest(domain* d, range_type t, const policies &p, int* level_reduction_rule=NULL);

    virtual void swapAdjacentVariables(int level);
    virtual void moveDownVariable(int high, int low);
    virtual void moveUpVariable(int low, int high);

    virtual void dynamicReorderVariables(int top, int bottom);

#ifdef ALLOW_DEPRECATED_0_17_7
    virtual enumerator::iterator* makeFullIter() const
    {
      return new mtmdd_iterator(this);
    }
#endif

  protected:
#ifdef ALLOW_DEPRECATED_0_17_7
    inline node_handle evaluateRaw(const dd_edge &f, const int* vlist) const {
      node_handle p = f.getNode();
      while (!isTerminalNode(p)) {
    	int level = getNodeLevel(p);
        p = getDownPtr(p, vlist[getVarByLevel(level)]);
      }
      return p;
    }

    class mtmdd_iterator : public mt_iterator {
      public:
        mtmdd_iterator(const forest* F);
        virtual ~mtmdd_iterator();
        virtual bool start(const dd_edge &e);
        virtual bool next();
      private:
        bool first(int k, node_handle p);
    };
#endif

    // Move the variable to the optimal level between top and bottom
    void sifting(int var, int top, int bottom);

};


//
// Helper class for createEdge
//

#ifdef ALLOW_DEPRECATED_0_17_7

namespace MEDDLY {

  template <typename T>
  class mtmdd_edgemaker {
      mtmdd_forest* F;
      const int* const* vlist;
      const T* values;
      int* order;
      int N;
      int K;
      binary_operation* unionOp;
      node_handle zero_terminal;
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
        terminal t(T(0));
        zero_terminal = t.getHandle();
      }

      inline const int* unprimed(int i) const {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        return vlist[order[i]];
      }
      inline int unprimed(int i, int k) const {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, k, K+1);
        return vlist[order[i]][k];
      }
      inline T term(int i) const {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        return values ? values[order[i]]: 1;
      }
      inline void swap(int i, int j) {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, N);
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, j, N);
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
          terminal t(term(start));
          return createEdgePath(k, unprimed(start),
            t.getHandle()
          );
        }
        //
        // Check terminal case
        //
        if (0==k) {
          T accumulate = term(start);
          for (int i=start+1; i<stop; i++) {
            accumulate += term(i);
          }
          terminal t(accumulate);
          return t.getHandle();
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

        node_handle dontcares = (batchP > start) ? createEdge(k-1, start, batchP) : 0;
        if(dontcares && F->isQuasiReduced()){
          // Add the redundant node at level k
          unpacked_node* nb = unpacked_node::newFull(F, k, lastV);
          for (unsigned v = 0; v<lastV; v++) {
              nb->setFull(v, F->linkNode(dontcares));
          }
#ifdef ALLOW_EXTENSIBLE
          if (F->isExtensibleLevel(k)) nb->markAsExtensible();
#endif
          node_handle built=F->createReducedNode(-1, nb);
          F->unlinkNode(dontcares);
          dontcares=built;
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
		    if(F->isQuasiReduced() && F->getTransparentNode()!=zero_terminal){
		      // TODO: Can be simpler
		      node_handle zero=makeOpaqueZeroNodeAtLevel(k-1);

		      while(z<nextV) {
                  nb->setSparse(z, z, F->linkNode(zero));
			      z++;
		      }
		      unsigned v = nextV;
		      while(v<lastV) {
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
			        if (v == unsigned(unprimed(i, k))) {
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
			      MEDDLY_DCASSERT(batchP > start);
			      node_handle total = createEdge(k-1, start, batchP);

			      //
			      // add to sparse node, unless transparent
			      //
			      if (total!=F->getTransparentNode()) {
                      nb->setSparse(z, v, total);
			          z++;
			      }

			      v++;
			      while(v<nextV) {
                    nb->setSparse(z, v, F->linkNode(zero));
			        z++;
			        v++;
			      }
		      }

#ifdef ALLOW_EXTENSIBLE
          if (F->isExtensibleLevel(k)) {
            nb->resize(lastV+1);
            nb->setSparse(z, v, F->linkNode(zero));
            z++;
            v++;
            nb->markAsExtensible();
          }
#endif

		      F->unlinkNode(zero);
		    }
		    else{
		      for (unsigned v = nextV; v<lastV; v = nextV) {
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
			      if (v == unsigned(unprimed(i, k))) {
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
			    MEDDLY_DCASSERT(batchP > start);
			    node_handle total = createEdge(k-1, start, batchP);

			    //
			    // add to sparse node, unless transparent
			    //
			    if (total==F->getTransparentNode()) continue;
                nb->setSparse(z, v, total);
			    z++;
		      }
		    }

        //
        // Cleanup
        //
        nb->shrink(z);

        MEDDLY_DCASSERT(unionOp);
        dd_edge dontcaresE(F), built(F);
        dontcaresE.set(dontcares);
        built.set( F->createReducedNode(-1, nb) );
        unionOp->computeTemp(dontcaresE, built, built);
        return F->linkNode(built);
      }

    protected:
      /// Special case for createEdge(), with only one minterm.
      inline node_handle
      createEdgePath(int k, const int* _vlist, node_handle bottom)
      {
          if (bottom==0 && (!F->isQuasiReduced() || F->getTransparentNode()==zero_terminal)) {
            return bottom;
          }

          for (int i=1; i<=k; i++) {
            if (DONT_CARE == _vlist[i]) {
              // make a redundant node
              if (F->isFullyReduced()) continue;
              int sz = F->getLevelSize(i);
              unpacked_node* nb = unpacked_node::newFull(F, i, sz);
              nb->setFull(0, bottom);
              for (int v=1; v<sz; v++) {
                nb->setFull(v, F->linkNode(bottom));
              }
#ifdef ALLOW_EXTENSIBLE
              if (F->isExtensibleLevel(i)) nb->markAsExtensible();
#endif
              bottom = F->createReducedNode(-1, nb);
            } else {
              if(F->isQuasiReduced() && F->getTransparentNode()!=zero_terminal){
                int sz = F->getLevelSize(i);
                if (F->isExtensibleLevel(i) && (sz == _vlist[i]+1)) sz++;
                unpacked_node* nb = unpacked_node::newFull(F, i, sz);
                node_handle zero=makeOpaqueZeroNodeAtLevel(i-1);
                // add opaque zero nodes
                for(int v=0; v<sz; v++) {
                    nb->setFull(v, (v==_vlist[i] ? bottom : F->linkNode(zero)));
                }
                F->unlinkNode(zero);
#ifdef ALLOW_EXTENSIBLE
                if (F->isExtensibleLevel(i)) nb->markAsExtensible();
#endif
                bottom=F->createReducedNode(-1, nb);
              }
              else{
                // make a singleton node
                unpacked_node* nb = unpacked_node::newSparse(F, i, 1);
                MEDDLY_DCASSERT(_vlist[i] >= 0);
                nb->setSparse(0, unsigned(_vlist[i]), bottom);
                bottom = F->createReducedNode(-1, nb);
              }
            }
          } // for i
          return bottom;
      }

      // Make zero nodes recursively when:
      // 1. quasi reduction
      // 2. the transparent value is not zero
      node_handle makeOpaqueZeroNodeAtLevel(int k)
      {
  	    MEDDLY_DCASSERT(F->isQuasiReduced());
  	    MEDDLY_DCASSERT(F->getTransparentNode()!=zero_terminal);

  	    return F->makeNodeAtLevel(k, zero_terminal);
      }
  }; // class mtmdd_edgemaker

}; // namespace MEDDLY

#endif // ALLOW_DEPRECATED_0_17_7

#endif
