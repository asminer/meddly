
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


#include "mtmxdbool.h"

MEDDLY::mt_mxd_bool::mt_mxd_bool(int dsl, domain *d, const policies &p)
: mtmxd_forest(dsl, d, BOOLEAN, p)
{ 
  initializeForest();
}

MEDDLY::mt_mxd_bool::~mt_mxd_bool()
{ }

void MEDDLY::mt_mxd_bool::createEdge(bool term, dd_edge& e)
{
  createEdgeTempl<bool_encoder, bool>(term, e);
}

void MEDDLY::mt_mxd_bool
::createEdge(const int* const* vlist, const int* const* vplist, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(UNION, this, this, this);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);
  enlargeVariables(vplist, N, true);
  mtmxd_edgemaker<bool_encoder, bool>
  EM(this, vlist, vplist, 0, order, N, getDomain()->getNumVariables(), unionOp);

  e.set(EM.createEdge(), 0);
}

void MEDDLY::mt_mxd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  createEdgeForVarTempl<bool_encoder, bool>(vh, vp, terms, a);
}

void MEDDLY::mt_mxd_bool::evaluate(const dd_edge &f, const int* vlist, 
  const int* vplist, bool &term) const
{
  term = bool_encoder::handle2value(evaluateRaw(f, vlist, vplist));
}

void MEDDLY::mt_mxd_bool::showTerminal(FILE* s, node_handle tnode) const
{
  bool_encoder::show(s, tnode);
}

void MEDDLY::mt_mxd_bool::writeTerminal(FILE* s, node_handle tnode) const
{
  bool_encoder::write(s, tnode);
}

MEDDLY::node_handle MEDDLY::mt_mxd_bool::readTerminal(FILE* s)
{
  return bool_encoder::read(s);
}

const char* MEDDLY::mt_mxd_bool::codeChars() const
{
  return "dd_txb";
}

