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


#include "mtmddreal.h"

MEDDLY::mt_mdd_real::mt_mdd_real(domain *d, const policies &p)
: mtmdd_forest(d, range_type::REAL, p)
{
}

MEDDLY::mt_mdd_real::~mt_mdd_real()
{ }

#ifdef ALLOW_DEPRECATED_0_18_0

void MEDDLY::mt_mdd_real::
createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a)
{
  // createEdgeForVarTempl<float_Tencoder, float>(vh, vp, terms, a);
  createEdgeForVarTempl<float>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}

#endif

