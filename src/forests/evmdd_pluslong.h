
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

#ifndef EVMDD_PLUSLONG_H
#define EVMDD_PLUSLONG_H

#include "evmdd.h"

namespace MEDDLY {
  class evmdd_pluslong;
  class evmdd_index_set_long;
};


class MEDDLY::evmdd_pluslong : public evmdd_forest {
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
        static inline void makeEmptyEdge(dd_edge& e) {
          e.set(0, long(0));
        }
        static inline void unionEq(long &a, long b) {
          if (b < a) {
            a = b;
          }
        }
    };

  public:
    evmdd_pluslong(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule=NULL, bool index_set=false);
    ~evmdd_pluslong();

    virtual void createEdge(long val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const long* terms, int N, dd_edge &e);
    virtual void createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a);
    virtual void evaluate(const dd_edge &f, const int* vlist, long &term) const;

    virtual bool isTransparentEdge(node_handle p, const void* v) const;
    virtual void getTransparentEdge(node_handle &p, void* v) const;
    virtual bool areEdgeValuesEqual(const void* eva, const void* evb) const;
    virtual bool isRedundant(const unpacked_node &nb) const;
    virtual bool isIdentityEdge(const unpacked_node &nb, int i) const;

    virtual enumerator::iterator* makeFullIter() const {
      return new evpimdd_iterator(this);
    }

    virtual void swapAdjacentVariables(int level);

  protected:
    virtual void normalize(unpacked_node &nb, long& ev) const;
    virtual void showEdgeValue(output &s, const void* edge) const;
    virtual void writeEdgeValue(output &s, const void* edge) const;
    virtual void readEdgeValue(input &s, void* edge);
    virtual void showUnhashedHeader(output &s, const void* uh) const;
    virtual void writeUnhashedHeader(output &s, const void* uh) const;
    virtual void readUnhashedHeader(input &s, unpacked_node &nb) const;
    virtual const char* codeChars() const;

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


class MEDDLY::evmdd_index_set_long : public evmdd_pluslong {
  public:
    evmdd_index_set_long(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule);
    ~evmdd_index_set_long();

    virtual void getElement(const dd_edge& a, int index, int* e);
    virtual void getElement(const dd_edge& a, long index, int* e);
};

#endif

