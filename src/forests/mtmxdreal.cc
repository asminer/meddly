
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


#include "mtmxdreal.h"

MEDDLY::mt_mxd_real::mt_mxd_real(int dsl, domain *d, const policies &p, float tv)
: mtmxd_forest(dsl, d, REAL, p)
{ 
  initializeForest();

  transparent=float_Tencoder::value2handle(tv);
}

MEDDLY::mt_mxd_real::~mt_mxd_real()
{ }

void MEDDLY::mt_mxd_real::createEdge(float term, dd_edge& e)
{
  createEdgeTempl<float_Tencoder, float>(term, e);
}

void MEDDLY::mt_mxd_real
::createEdge(const int* const* vlist, const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(PLUS, this, this, this);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);
  enlargeVariables(vplist, N, true);

  mtmxd_edgemaker<float_Tencoder, float>
  EM(this, vlist, vplist, terms, order, N, getDomain()->getNumVariables(), unionOp);

  e.set(EM.createEdge());
}

void MEDDLY::mt_mxd_real::
createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a)
{
  createEdgeForVarTempl<float_Tencoder, float>(vh, vp, terms, a);
}

void MEDDLY::mt_mxd_real::evaluate(const dd_edge &f, const int* vlist, 
  const int* vplist, float &term) const
{
  term = float_Tencoder::handle2value(evaluateRaw(f, vlist, vplist));
}

void MEDDLY::mt_mxd_real::showTerminal(FILE* s, node_handle tnode) const
{
  float_Tencoder::show(s, tnode);
}

void MEDDLY::mt_mxd_real::writeTerminal(FILE* s, node_handle tnode) const
{
  float_Tencoder::write(s, tnode);
}

MEDDLY::node_handle MEDDLY::mt_mxd_real::readTerminal(FILE* s)
{
  return float_Tencoder::read(s);
}

const char* MEDDLY::mt_mxd_real::codeChars() const
{
  return "dd_txr";
}

