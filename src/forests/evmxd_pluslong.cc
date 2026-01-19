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

#include "evmxd_pluslong.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    evmxd_pluslong  methods                    *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::evmxd_pluslong::evmxd_pluslong(domain *d, const policies &p)
 : evmxd_forest(d, range_type::INTEGER, edge_labeling::EVPLUS, p)
{
}

MEDDLY::evmxd_pluslong::~evmxd_pluslong()
{ }

#ifdef ALLOW_DEPRECATED_0_18_0

void MEDDLY::evmxd_pluslong
::createEdgeForVar(int vh, bool vp, const long* terms, dd_edge& a)
{
  createEdgeForVarTempl<OP, long>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}

#endif

