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

#ifndef MEDDLY_EVMDD_PLUSLONG_H
#define MEDDLY_EVMDD_PLUSLONG_H

#include "evmdd.h"

namespace MEDDLY {
  class evmdd_pluslong;
  class evmdd_index_set_long;
};


class MEDDLY::evmdd_pluslong : public evmdd_forest {
  public:
    class OP {
      public:
        static inline long apply(long a, long b) {
          return a + b;
        }
        static inline void unionEq(long &a, long b) {
          if (b < a) {
            a = b;
          }
        }
    };

  public:
    evmdd_pluslong(domain *d, const policies &p, bool index_set=false);
    ~evmdd_pluslong();

#ifdef ALLOW_DEPRECATED_0_18_0
    virtual void createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a);
#endif

    virtual void swapAdjacentVariables(int level);

};


class MEDDLY::evmdd_index_set_long : public evmdd_pluslong {
  public:
    evmdd_index_set_long(domain *d, const policies &p);
    virtual ~evmdd_index_set_long();
};

#endif

