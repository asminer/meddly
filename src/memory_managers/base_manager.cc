
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

// ******************************************************************
// *                                                                *
// *                  memory_manager_style methods                  *
// *                                                                *
// ******************************************************************

MEDDLY::memory_manager_style::memory_manager_style(const char* n)
{
  name = n;
}

MEDDLY::memory_manager_style::~memory_manager_style()
{
}

// ******************************************************************
// *                                                                *
// *                     memory_manager methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::memory_manager::memory_manager(const char* n, memstats &stats)
 : my_mem(stats)
{
  style_name = n;
  chunk_base = 0;
  chunk_multiplier = 0;
}

MEDDLY::memory_manager::~memory_manager()
{
}

void* MEDDLY::memory_manager::slowChunkAddress(node_address h) const
{
  MEDDLY_DCASSERT(0);
  return 0;
}

