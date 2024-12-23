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

// This hasn't been tested very much yet...

#ifndef MEDDLY_EVMDD_TIMESREAL_H
#define MEDDLY_EVMDD_TIMESREAL_H

#include "evmdd.h"

namespace MEDDLY {
  class evmdd_timesreal;
};


class MEDDLY::evmdd_timesreal : public evmdd_forest {
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
    evmdd_timesreal(domain *d, const policies &p, int* level_reduction_rule=NULL);
    ~evmdd_timesreal();

    virtual void createEdge(float val, dd_edge &e);
    virtual void createEdge(double val, dd_edge &e);
    virtual void createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a);

#ifdef ALLOW_DEPRECATED_0_17_7
    virtual void createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, float &term) const;

    virtual enumerator::iterator* makeFullIter() const {
      return new evtrmdd_iterator(this);
    }
#endif

    virtual void showEdge(output &s, const edge_value &ev, node_handle d) const;

  protected:
#ifdef ALLOW_DEPRECATED_0_17_3
    virtual const char* codeChars() const;
#endif

#ifdef ALLOW_DEPRECATED_0_17_7
  protected:
    class evtrmdd_iterator : public enumerator::iterator {
      public:
        evtrmdd_iterator(const forest* F);
        virtual ~evtrmdd_iterator();

        virtual void getValue(float &termVal) const;
        virtual bool start(const dd_edge &e);
        virtual bool next();
      private:
        bool first(int k, node_handle p);

      protected:
        double* acc_evs;  // for accumulating edge values
    };
#endif
};

#endif

