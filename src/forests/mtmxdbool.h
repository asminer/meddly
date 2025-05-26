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

#ifndef MEDDLY_MTMXDBOOL_H
#define MEDDLY_MTMXDBOOL_H

#include "mtmxd.h"

namespace MEDDLY {
  class mt_mxd_bool;
};

// ******************************************************************

/**
    Forest for multi-terminal, mxd, boolean range.
*/
class MEDDLY::mt_mxd_bool : public mtmxd_forest {
  public:

    mt_mxd_bool(domain *d, const policies &p, bool tv=false);
    ~mt_mxd_bool();

#ifdef ALLOW_DEPRECATED_0_17_9
    virtual void createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a);
#endif
#ifdef ALLOW_DEPRECATED_0_17_7
    virtual void createEdge(bool val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge& e);
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;
#endif

    virtual node_handle unionOneMinterm(node_handle a,  int* from,  int* to, int level);

#ifdef VIRTUAL_IO_METHODS
    virtual void showEdge(output &s, const edge_value &ev, node_handle d) const;
#endif

  protected:
    node_handle unionOneMinterm_r(int index, int vh, node_handle nh,  int* from,  int* to);
    bool checkTerminalMinterm(node_handle a,  int* from,  int* to,int levelRead, node_handle& c);

};

#endif

