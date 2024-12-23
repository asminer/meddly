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

#ifndef MEDDLY_MTMXDREAL_H
#define MEDDLY_MTMXDREAL_H

#include "mtmxd.h"

namespace MEDDLY {
  class mt_mxd_real;
};

// ******************************************************************

/**
    Forest for multi-terminal, mxd, real range.
*/
class MEDDLY::mt_mxd_real : public mtmxd_forest {
  public:

    mt_mxd_real(domain *d, const policies &p, int* level_reduction_rule=NULL, float tv=0);
    ~mt_mxd_real();

    virtual void createEdge(float val, dd_edge &e);
    virtual void createEdge(double val, dd_edge &e);

    virtual void createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a);

#ifdef ALLOW_DEPRECATED_0_17_7
    void createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge& e);

    void evaluate(const dd_edge& f, const int* vlist, const int* vplist,
        float &term) const;
#endif

    virtual void showEdge(output &s, const edge_value &ev, node_handle d) const;

#ifdef ALLOW_DEPRECATED_0_17_3
  protected:
    virtual const char* codeChars() const;
#endif
};

#endif

