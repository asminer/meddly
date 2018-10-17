
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

#ifndef CBDD_H
#define CBDD_H

#include "evmdd.h"

namespace MEDDLY {
  class cbdd;
};


class MEDDLY::cbdd : public evmdd_forest {
  public:
    class OP : public EVencoder<long> {
      public:
        static inline void setEdge(void* ptr, long v) {
          writeValue(ptr, v);
        }
        static inline bool isIdentityEdge(const void* p) {
          long ev;
          readValue(p, ev);
          return (0 == ev);
        }
        static inline bool isTransparentEdge(const void* p) {
          return isIdentityEdge(p);
        }
        static inline long getRedundantEdge() {
          return 0;
        }
        static inline long apply(long a, long b) {
          return a + b;
        }
        static inline void makeEmptyEdge(long &ev, node_handle &ep) {
          ev = 0;
          ep = 0;
        }
        static inline void unionEq(long &a, long b) {
          if (b < a) {
            a = b;
          }
        }
    };

  public:
    cbdd(int dsl, domain *d, const policies &p, int* level_reduction_rule=NULL);
    ~cbdd();

    virtual void createEdge(long reduction, dd_edge &e);
//    virtual void createEdge(const int* const* vlist, const long* terms, int N, dd_edge &e);
//    virtual void createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a);
    virtual void evaluate(const dd_edge &f, const int* vlist, long &term) const;

    virtual bool isTransparentEdge(node_handle p, const void* v) const;
    virtual void getTransparentEdge(node_handle &p, void* v) const;
    virtual bool areEdgeValuesEqual(const void* eva, const void* evb) const;
    virtual bool isRedundant(const unpacked_node &nb) const;
    virtual bool isIdentityEdge(const unpacked_node &nb, int i) const;
    virtual bool isZeroSuppressed(const unpacked_node &nb) const;
//    virtual bool isOneSuppressed(const unpacked_node &nb) const;

    virtual enumerator::iterator* makeFullIter() const {
      return new evpimdd_iterator(this);
    }

  protected:
    virtual void normalize(unpacked_node &nb, long& ev) const;
    virtual void showEdgeValue(output &s, const void* edge) const;
    virtual void writeEdgeValue(output &s, const void* edge) const;
    virtual void readEdgeValue(input &s, void* edge);
    virtual void showUnhashedHeader(output &s, const void* uh) const;
    virtual void writeUnhashedHeader(output &s, const void* uh) const;
    virtual void readUnhashedHeader(input &s, unpacked_node &nb) const;
    virtual const char* codeChars() const;

    void createReducedHelper(int in, unpacked_node &nb, long& r, node_handle& node);

  public:
    void createReducedNode(int in, MEDDLY::unpacked_node *un, long& ev,
        node_handle& node) {
      MEDDLY_DCASSERT(un);
      normalize(*un, ev);
      //un->computeHash();
      createReducedHelper(in, *un, ev, node);
      ev = -1;
    #ifdef TRACK_DELETIONS
      printf("Created node %d\n", node);
    #endif
      unpacked_node::recycle(un);
    }

  protected:
    class evpimdd_iterator : public enumerator::iterator {
      public:
        evpimdd_iterator(const expert_forest* F);
        virtual ~evpimdd_iterator();

        virtual void getValue(long &termVal) const;
        virtual bool start(const dd_edge &e);
        virtual bool next();
      private:
        bool first(int k, node_handle p);

      protected:
        long* acc_evs;  // for accumulating edge values
    };
};

#endif

