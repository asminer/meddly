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

#ifndef MEDDLY_EVMXD_TIMESREAL_H
#define MEDDLY_EVMXD_TIMESREAL_H

#include "evmxd.h"

namespace MEDDLY {
  class evmxd_timesreal;
};


class MEDDLY::evmxd_timesreal : public evmxd_forest {
  public:
    class OP {
      public:
        static inline double apply(double a, double b) {
          return a*b;
        }
        static inline void unionEq(float &a, float b) {
          a += b;
        }
    };

  public:
    evmxd_timesreal(domain *d, const policies &p, int* level_reduction_rule=NULL);
    ~evmxd_timesreal();

    virtual void createEdgeForVar(int vh, bool vp, const float* terms,
      dd_edge& a);
#ifdef ALLOW_DEPRECATED_0_17_7
    virtual void createEdge(float val, dd_edge &e);
    virtual void createEdge(double val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
      const float* terms, int N, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist,
      const int* vplist, float &term) const;

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
#endif

    virtual void showEdge(output &s, const edge_value &ev, node_handle d) const;

#ifdef ALLOW_DEPRECATED_0_17_7
  protected:
    class evtrmxd_baseiter : public enumerator::iterator {
      public:
        evtrmxd_baseiter(const forest* F);
        virtual ~evtrmxd_baseiter();
        virtual void getValue(float &termVal) const;
      protected:
        double* acc_evs;  // for accumulating edge values
      private:
        double* raw_acc_evs;
    };

    class evtrmxd_iterator : public evtrmxd_baseiter {
      public:
        evtrmxd_iterator(const forest* F) : evtrmxd_baseiter(F) { }
        virtual ~evtrmxd_iterator() { }

        virtual bool start(const dd_edge &e);
        virtual bool next();
      private:
        bool first(int k, node_handle p);
    };

    class evtrmxd_fixedrow_iter : public evtrmxd_baseiter {
      public:
        evtrmxd_fixedrow_iter(const forest* F) : evtrmxd_baseiter(F) {}
        virtual ~evtrmxd_fixedrow_iter() { }

        virtual bool start(const dd_edge &e, const int*);
        virtual bool next();
      private:
        bool first(int k, node_handle p);
    };

    class evtrmxd_fixedcol_iter : public evtrmxd_baseiter {
      public:
        evtrmxd_fixedcol_iter(const forest* F) : evtrmxd_baseiter(F) {}
        virtual ~evtrmxd_fixedcol_iter() { }

        virtual bool start(const dd_edge &e, const int*);
        virtual bool next();
      private:
        bool first(int k, node_handle p);
    };
#endif
};

#endif

