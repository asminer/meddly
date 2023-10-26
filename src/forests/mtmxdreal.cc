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

MEDDLY::mt_mxd_real::mt_mxd_real(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule, float tv)
: mtmxd_forest(dsl, d, range_type::REAL, p, level_reduction_rule)
{
    terminal t(tv);
    setTransparentEdge(t.getHandle());
    // setTransparentEdge(float_Tencoder::value2handle(tv));
    initializeForest();
}

MEDDLY::mt_mxd_real::~mt_mxd_real()
{ }

void MEDDLY::mt_mxd_real::createEdge(float term, dd_edge& e)
{
    // createEdgeTempl<float_Tencoder, float>(term, e);
    createEdgeTempl<float>(term, e);
#ifdef DEVELOPMENT_CODE
    validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_real
::createEdge(const int* const* vlist, const int* const* vplist, const float* terms, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(PLUS, this, this, this);
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

  // mtmxd_edgemaker<float_Tencoder, float>
  mtmxd_edgemaker<float>
  EM(this, ordered_vlist, ordered_vplist, terms, order, N, num_vars, unionOp);

  e.set(EM.createEdge());

  free(ordered_vlist);
  free(ordered_vplist);

#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_real::
createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a)
{
    // createEdgeForVarTempl<float_Tencoder, float>(vh, vp, terms, a);
    createEdgeForVarTempl<float>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
    validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_real::evaluate(const dd_edge &f, const int* vlist,
  const int* vplist, float &term) const
{
    terminal t;
    t.setFromHandle(terminal_type::REAL, evaluateRaw(f, vlist, vplist));
    term = t.getReal();
    // term = float_Tencoder::handle2value(evaluateRaw(f, vlist, vplist));
}

void MEDDLY::mt_mxd_real::showEdge(output &s, const edge_value &ev,
        node_handle d) const
{
    if (d>0) {
        s.put('#');
        s.put(d);
    } else {
        terminal t;
        t.setFromHandle(terminal_type::REAL, d);
        s.put('t');
        s.put( t.getReal() );
    }
}


const char* MEDDLY::mt_mxd_real::codeChars() const
{
  return "dd_txr";
}

