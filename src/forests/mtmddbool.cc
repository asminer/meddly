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


#include "mtmddbool.h"
#include "../terminal.h"
#include "../opname.h"

MEDDLY::mt_mdd_bool::mt_mdd_bool(domain *d, const policies &p, int* level_reduction_rule, bool tv)
: mtmdd_forest(d, range_type::BOOLEAN, p, level_reduction_rule)
{
    terminal t(tv);
    setTransparentEdge(t.getHandle());
    // setTransparentEdge(bool_Tencoder::value2handle(tv));
    initializeStorage();
}

MEDDLY::mt_mdd_bool::~mt_mdd_bool()
{ }

void MEDDLY::mt_mdd_bool::createEdge(bool term, dd_edge& e)
{
  // createEdgeTempl<bool_Tencoder, bool>(term, e);
  createEdgeTempl<bool>(term, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool::createEdge(const int* const* vlist, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(UNION, e, e, e);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);

  int num_vars=getNumVariables();

  // Create vlist following the mapping between variable and level
  int** ordered_vlist=static_cast<int**>(malloc(N*sizeof(int*)+(num_vars+1)*N*sizeof(int)));
  if(ordered_vlist==0){
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }

  ordered_vlist[0]=reinterpret_cast<int*>(&ordered_vlist[N]);
  for(int i=1; i<N; i++) {
	  ordered_vlist[i]=(ordered_vlist[i-1]+num_vars+1);
  }
  for(int i=0; i<=num_vars; i++) {
	  int level=getLevelByVar(i);
	  for(int j=0; j<N; j++) {
		  ordered_vlist[j][level]=vlist[j][i];
	  }
  }

  // mtmdd_edgemaker<bool_Tencoder, bool>
  mtmdd_edgemaker<bool>
    EM(this, ordered_vlist, 0, order, N, getDomain()->getNumVariables(), unionOp);

  e.set(EM.createEdge());

  free(ordered_vlist);

#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  //createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
  createEdgeForVarTempl<bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool
::evaluate(const dd_edge &f, const int* vlist, bool &term) const
{
    terminal t(terminal_type::BOOLEAN, evaluateRaw(f, vlist));
    term = t.getBoolean();
    // term = bool_Tencoder::handle2value(evaluateRaw(f, vlist));
}

void MEDDLY::mt_mdd_bool::showEdge(output &s, const edge_value &ev,
        node_handle d) const
{
    if (d>0) {
        s.put('#');
        s.put(d);
    } else {
        terminal t(terminal_type::BOOLEAN, d);
        s.put( t.getBoolean() ? 'T' : 'F' );
    }
}

#ifdef ALLOW_DEPRECATED_0_17_3

const char* MEDDLY::mt_mdd_bool::codeChars() const
{
  return "dd_tvb";
}

#endif
