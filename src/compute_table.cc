
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


#include "defines.h"

// **********************************************************************
// *                                                                    *
// *                         ct_object  methods                         *
// *                                                                    *
// **********************************************************************

MEDDLY::ct_object::ct_object()
{
}

MEDDLY::ct_object::~ct_object()
{
}

// **********************************************************************
// *                                                                    *
// *                    compute_table_style  methods                    *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table_style::compute_table_style()
{
}

MEDDLY::compute_table_style::~compute_table_style()
{
}

MEDDLY::compute_table* 
MEDDLY::compute_table_style::create(const settings::computeTableSettings &s)
      const
{
  throw error(error::TYPE_MISMATCH);
}


MEDDLY::compute_table* 
MEDDLY::compute_table_style::create(const settings::computeTableSettings &s, 
      operation* op) const
{
  throw error(error::TYPE_MISMATCH);
}

// **********************************************************************
// *                                                                    *
// *                       compute_table  methods                       *
// *                                                                    *
// **********************************************************************

MEDDLY::compute_table::compute_table(const settings::computeTableSettings &s)
{
  maxSize = s.maxSize;
  if (0==maxSize)
    throw error(error::INVALID_ASSIGNMENT);

  switch (s.staleRemoval) {
    case settings::computeTableSettings::Aggressive:
            checkStalesOnFind = true;
            checkStalesOnResize = true;
            break;
    case settings::computeTableSettings::Moderate:
            checkStalesOnFind = false;
            checkStalesOnResize = true;
            break;
    case settings::computeTableSettings::Lazy:
            checkStalesOnFind = false;
            checkStalesOnResize = false;
            break;
  }

  perf.numEntries = 0;
  perf.hits = 0;
  perf.pings = 0;
  perf.numLargeSearches = 0;
  perf.maxSearchLength = 0;
  for (int i=0; i<perf.searchHistogramSize; i++)
    perf.searchHistogram[i] = 0;
}

MEDDLY::compute_table::~compute_table()
{
}

// **********************************************************************

MEDDLY::compute_table::search_key::search_key(operation* _op)
{
  op = _op;
}

MEDDLY::compute_table::search_key::~search_key()
{
}

// **********************************************************************

MEDDLY::compute_table::search_result::search_result()
{

}
MEDDLY::compute_table::search_result::~search_result()
{
}

// **********************************************************************

MEDDLY::compute_table::entry_builder::entry_builder()
{

}
MEDDLY::compute_table::entry_builder::~entry_builder()
{
}

