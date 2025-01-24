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

#ifndef MEDDLY_MTMDDREAL_H
#define MEDDLY_MTMDDREAL_H

#include "mtmdd.h"

namespace MEDDLY {
  class mt_mdd_real;
};

// ******************************************************************

/**
    Forest for multi-terminal, mdd, real range.
*/
class MEDDLY::mt_mdd_real : public mtmdd_forest {
  public:

    mt_mdd_real(domain *d, const policies &p, int* level_reduction_rule=NULL, float tv=0);
    ~mt_mdd_real();

    virtual void createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a);
#ifdef ALLOW_DEPRECATED_0_17_7
    virtual void createEdge(float val, dd_edge &e);
    virtual void createEdge(double val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const float* terms, int N, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, float &term)
      const;
#endif

    virtual void showEdge(output &s, const edge_value &ev, node_handle d) const;

};

#endif

