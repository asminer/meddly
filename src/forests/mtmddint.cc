
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


#include "mtmddint.h"

MEDDLY::mt_mdd_int::mt_mdd_int(int dsl, domain *d, const policies &p)
: mtmdd_forest(dsl, d, INTEGER, p)
{ 
  initializeForest();
}

MEDDLY::mt_mdd_int::~mt_mdd_int()
{ }

void MEDDLY::mt_mdd_int::createEdge(int term, dd_edge& e)
{
  createEdgeTempl<int_encoder, int>(term, e);
}

void MEDDLY::mt_mdd_int::createEdge(const int* const* vlist, const int* terms, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(PLUS, this, this, this);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);

  mtmdd_edgemaker<int_encoder, int>
  EM(this, vlist, terms, order, N, getDomain()->getNumVariables(), unionOp);

  e.set(EM.createEdge(), 0);
}

void MEDDLY::mt_mdd_int::
createEdgeForVar(int vh, bool vp, const int* terms, dd_edge& a)
{
  createEdgeForVarTempl<int_encoder, int>(vh, vp, terms, a);
}

void MEDDLY::mt_mdd_int
::evaluate(const dd_edge &f, const int* vlist, int &term) const
{
  term = int_encoder::handle2value(evaluateRaw(f, vlist));
}

void MEDDLY::mt_mdd_int::showTerminal(FILE* s, node_handle tnode) const
{
  int_encoder::show(s, tnode);
}

void MEDDLY::mt_mdd_int::writeTerminal(FILE* s, node_handle tnode) const
{
  int_encoder::write(s, tnode);
}

MEDDLY::node_handle MEDDLY::mt_mdd_int::readTerminal(FILE* s)
{
  return int_encoder::read(s);
}

const char* MEDDLY::mt_mdd_int::codeChars() const
{
  return "dd_tvi";
}

