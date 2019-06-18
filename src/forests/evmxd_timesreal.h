
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

#ifndef EVMXD_TIMESREAL_H
#define EVMXD_TIMESREAL_H

#include "evmxd.h"

namespace MEDDLY {
  class evmxd_timesreal;
};


class MEDDLY::evmxd_timesreal : public evmxd_forest {
  public:
    class OP : public EVencoder<float> {
      public:
        static inline void setEdge(void* ptr, float v) {
          writeValue(ptr, v);
        }
        static inline bool isIdentityEdge(const void* p) {
          float ev;
          readValue(p, ev);
          return !notClose(ev, 1.0);
        }
        static inline bool isTransparentEdge(const void* p) {
          float ev;
          readValue(p, ev);
          return (0.0 == ev);
        }
        static inline double getRedundantEdge() {
          return 1.0f;
        }
        static inline double apply(double a, double b) {
          return a*b;
        }
        static inline void makeEmptyEdge(float &ev, node_handle &ep) {
          ev = 0;
          ep = 0;
        }
        static inline void makeEmptyEdge(node_handle &ep, void* ev) {
          ep = 0;
          writeValue(ev, 0);
        }
        static inline void unionEq(float &a, float b) {
          a += b;
        }
        // bonus
        static inline bool notClose(float a, float b) {
          if (a) {
            double diff = a-b;
            return ABS(diff/a) > 1e-6;
          } else {
            return ABS(b) > 1e-10;
          }
        }
    };

  public:
    evmxd_timesreal(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule=NULL);
    ~evmxd_timesreal();

    virtual void createEdge(float val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
      const float* terms, int N, dd_edge &e);
    virtual void createEdgeForVar(int vh, bool vp, const float* terms,
      dd_edge& a);
    virtual void evaluate(const dd_edge &f, const int* vlist,
      const int* vplist, float &term) const;

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
    virtual void normalize(unpacked_node &nb, float& ev) const;
    virtual void showEdgeValue(output &s, const void* edge) const;
    virtual void writeEdgeValue(output &s, const void* edge) const;
    virtual void readEdgeValue(input &s, void* edge);
    virtual const char* codeChars() const;

  protected:
    class evtrmxd_baseiter : public enumerator::iterator {
      public:
        evtrmxd_baseiter(const expert_forest* F);
        virtual ~evtrmxd_baseiter();
        virtual void getValue(float &termVal) const;
      protected:
        double* acc_evs;  // for accumulating edge values
      private:
        double* raw_acc_evs;  
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

