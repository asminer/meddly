
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

MEDDLY::mt_mxd_bool::mt_mxd_bool(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule, bool tv)
: mtmxd_forest(dsl, d, BOOLEAN, p, level_reduction_rule)
{ 
  initializeForest();

  transparent=bool_Tencoder::value2handle(tv);
}

MEDDLY::mt_mxd_bool::~mt_mxd_bool()
{ }

void MEDDLY::mt_mxd_bool::createEdge(bool term, dd_edge& e)
{
  createEdgeTempl<bool_Tencoder, bool>(term, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_bool
::createEdge(const int* const* vlist, const int* const* vplist, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(UNION, this, this, this);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);
  enlargeVariables(vplist, N, true);

  int num_vars=getNumVariables();

  // Create vlist and vplist following the mapping between variable and level
  int** ordered_vlist=static_cast<int**>(malloc(N*sizeof(int*)+(num_vars+1)*N*sizeof(int)));
  if(ordered_vlist==0){
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  int** ordered_vplist=static_cast<int**>(malloc(N*sizeof(int*)+(num_vars+1)*N*sizeof(int)));
  if(ordered_vplist==0){
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }

  ordered_vlist[0]=reinterpret_cast<int*>(&ordered_vlist[N]);
  ordered_vplist[0]=reinterpret_cast<int*>(&ordered_vplist[N]);
  for(int i=1; i<N; i++) {
	  ordered_vlist[i]=(ordered_vlist[i-1]+num_vars+1);
	  ordered_vplist[i]=(ordered_vplist[i-1]+num_vars+1);
  }
  for(int i=0; i<=num_vars; i++) {
	  int level=getLevelByVar(i);
	  for(int j=0; j<N; j++) {
		  ordered_vlist[j][level]=vlist[j][i];
		  ordered_vplist[j][level]=vplist[j][i];
	  }
  }

  mtmxd_edgemaker<bool_Tencoder, bool>
  EM(this, ordered_vlist, ordered_vplist, 0, order, N, num_vars, unionOp);

  e.set(EM.createEdge());

  free(ordered_vlist);
  free(ordered_vplist);

#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_bool::evaluate(const dd_edge &f, const int* vlist, 
  const int* vplist, bool &term) const
{
  term = bool_Tencoder::handle2value(evaluateRaw(f, vlist, vplist));
}

void MEDDLY::mt_mxd_bool::showTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::show(s, tnode);
}

void MEDDLY::mt_mxd_bool::writeTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::write(s, tnode);
}

MEDDLY::node_handle MEDDLY::mt_mxd_bool::readTerminal(input &s)
{
  return bool_Tencoder::read(s);
}

const char* MEDDLY::mt_mxd_bool::codeChars() const
{
  return "dd_txb";
}

