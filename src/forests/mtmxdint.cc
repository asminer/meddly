
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
: mtmxd_forest<expert_forest::int_encoder>(dsl, d, INTEGER, p)
{ 
  initializeForest();
}

MEDDLY::mt_mxd_int::~mt_mxd_int()
{ }

void MEDDLY::mt_mxd_int::createEdge(int term, dd_edge& e)
{
  createEdgeTempl(term, e);
}

void MEDDLY::mt_mxd_int
::createEdge(int** vlist, int** vplist, int* terms, int N, dd_edge &e)
{
  unionOp = getOperation(PLUS, this, this, this);
  enlargeVariables(vlist, N, false);
  enlargeVariables(vplist, N, true);
  e.set(createEdgeRT(getDomain()->getNumVariables(), vlist, vplist, terms, N), 0);
}

void MEDDLY::mt_mxd_int::
createEdgeForVar(int vh, bool vp, const int* terms, dd_edge& a)
{
  createEdgeForVarTempl(vh, vp, terms, a);
}

void MEDDLY::mt_mxd_int::evaluate(const dd_edge &f, const int* vlist, 
  const int* vplist, int &term) const
{
  evaluateTempl(f, vlist, vplist, term);
}

const char* MEDDLY::mt_mxd_int::codeChars() const
{
  return "dd_txi";
}

