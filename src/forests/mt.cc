
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


// TODO: implement changes made to mt_forest interface
//
// TODO: HERE: go through every function in mdds.h and mdds.cc


#include "mt.h"
#include "../unique_table.h"

#define ENABLE_CACHE_COUNTING 0
#define ENABLE_IN_COUNTING 0

// ******************************************************************
// *                                                                *
// *                     mt_base_forest methods                     *
// *                                                                *
// ******************************************************************


MEDDLY::mt_base_forest::mt_base_forest(int dsl, domain *d, bool rel,
  range_type t, const policies &p)
: expert_forest(dsl, d, rel, t, MULTI_TERMINAL, p)
{
  unionOp = 0;
}

bool MEDDLY::mt_base_forest::isRedundant(const node_builder &nb) const
{
  if (isQuasiReduced()) return false;
  if (nb.getLevel() < 0 && isIdentityReduced()) return false;
  int common = nb.d(0);
  for (int i=1; i<nb.rawSize(); i++) 
    if (nb.d(i) != common) return false;
  return true;
}

bool MEDDLY::mt_base_forest::isIdentityEdge(const node_builder &nb, int i) const
{
  if (nb.getLevel() > 0) return false;
  if (!isIdentityReduced()) return false;
  if (i<0) return false;
  return nb.d(i) != 0;
}

