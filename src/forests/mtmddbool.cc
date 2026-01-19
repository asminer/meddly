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

MEDDLY::mt_mdd_bool::mt_mdd_bool(domain *d, const policies &p)
: mtmdd_forest(d, range_type::BOOLEAN, p)
{
}

MEDDLY::mt_mdd_bool::~mt_mdd_bool()
{ }

#ifdef ALLOW_DEPRECATED_0_18_0

void MEDDLY::mt_mdd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  //createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
  createEdgeForVarTempl<bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}

#endif

