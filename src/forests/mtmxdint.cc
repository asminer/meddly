
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


#include "mtmxdint.h"


MEDDLY::mt_mxd_int::mt_mxd_int(int dsl, domain *d, const policies &p)
: MEDDLY::mtmxd_forest(dsl, d, true, INTEGER, MULTI_TERMINAL, p)
{ }


MEDDLY::mt_mxd_int::~mt_mxd_int()
{ }

void MEDDLY::mt_mxd_int::createEdge(int val, dd_edge &e)
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);

  int node = createEdgeTo(getTerminalNode(val));
  e.set(node, 0, getNodeLevel(node));
}


void MEDDLY::mt_mxd_int::createEdge(const int* const* vlist,
    const int* const* vplist, const int* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0 || terms == 0 || N <= 0)
    throw error(error::INVALID_VARIABLE);

  createEdgeInternal(vlist, vplist, terms, N, e);
}


void MEDDLY::mt_mxd_int::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, int &term) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);
  if (f.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0) 
    throw error(error::INVALID_VARIABLE);

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getInteger(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
}

