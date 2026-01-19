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

#ifndef MEDDLY_MTMDDBOOL_H
#define MEDDLY_MTMDDBOOL_H

#include "mtmdd.h"

namespace MEDDLY {
  class mt_mdd_bool;
};

// ******************************************************************

/**
    Forest for multi-terminal, mdd, boolean range.
*/
class MEDDLY::mt_mdd_bool : public mtmdd_forest {
  public:

    mt_mdd_bool(domain *d, const policies &p);
    ~mt_mdd_bool();

#ifdef ALLOW_DEPRECATED_0_18_0
    virtual void createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a);
#endif
#ifdef ALLOW_DEPRECATED_0_17_7
    virtual void createEdge(bool val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;
#endif

};

#endif

