
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

#ifndef EVMXD_PLUSLONG_H
#define EVMXD_PLUSLONG_H

#include "evmxd.h"

namespace MEDDLY {
  class evmxd_pluslong;
};


class MEDDLY::evmxd_pluslong : public evmxd_forest {
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
        static inline long apply(double a, double b) {
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
    evmxd_pluslong(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule);
    ~evmxd_pluslong();

    // using evmxd_forest::createEdgeForVar;

    virtual void createEdge(long val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
      const long* terms, int N, dd_edge &e);
    virtual void createEdgeForVar(int vh, bool vp, const long* terms,
      dd_edge& a);
    virtual void evaluate(const dd_edge &f, const int* vlist,
      const int* vplist, long &term) const;


    virtual bool isTransparentEdge(node_handle p, const void* v) const;
    virtual void getTransparentEdge(node_handle &p, void* v) const;
    virtual bool areEdgeValuesEqual(const void* eva, const void* evb) const;
    virtual bool isRedundant(const unpacked_node &nb) const;
    virtual bool isIdentityEdge(const unpacked_node &nb, int i) const;

    virtual enumerator::iterator* makeFullIter() const {
      return new evtrmxd_iterator(this);
    }

    virtual enumerator::iterator* makeFixedRowIter() const 
    {
      return new evtrmxd_fixedrow_iter(this);
    }

    virtual enumerator::iterator* makeFixedColumnIter() const 
    {
      return new evtrmxd_fixedcol_iter(this);
    }


  protected:
    virtual void normalize(unpacked_node &nb, long& ev) const;
    virtual void showEdgeValue(output &s, const void* edge) const;
    virtual void writeEdgeValue(output &s, const void* edge) const;
    virtual void readEdgeValue(input &s, void* edge);
    virtual const char* codeChars() const;

  protected:
    class evtrmxd_baseiter : public enumerator::iterator {
      public:
        evtrmxd_baseiter(const expert_forest* F);
        virtual ~evtrmxd_baseiter();
        virtual void getValue(long &termVal) const;
      protected:
        long* acc_evs;  // for accumulating edge values
      private:
        long* raw_acc_evs;
    };

    class evtrmxd_iterator : public evtrmxd_baseiter {
      public:
        evtrmxd_iterator(const expert_forest* F) : evtrmxd_baseiter(F) { }
        virtual ~evtrmxd_iterator() { }

        virtual bool start(const dd_edge &e);
        virtual bool next();
      private:
        bool first(int k, node_handle p);
    };

    class evtrmxd_fixedrow_iter : public evtrmxd_baseiter {
      public:
        evtrmxd_fixedrow_iter(const expert_forest* F) : evtrmxd_baseiter(F) {}
        virtual ~evtrmxd_fixedrow_iter() { }

        virtual bool start(const dd_edge &e, const int*);
        virtual bool next();
      private:
        bool first(int k, node_handle p);
    };

    class evtrmxd_fixedcol_iter : public evtrmxd_baseiter {
      public:
        evtrmxd_fixedcol_iter(const expert_forest* F) : evtrmxd_baseiter(F) {}
        virtual ~evtrmxd_fixedcol_iter() { }

        virtual bool start(const dd_edge &e, const int*);
        virtual bool next();
      private:
        bool first(int k, node_handle p);
    };

};

#endif

