
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "init_managers.h"

#include "orig_grid.h"
#include "array_grid.h"
#include "malloc_style.h"
#include "heap_manager.h"
#include "freelists.h"

namespace MEDDLY {
  const memory_manager_style* ORIGINAL_GRID = 0;
  const memory_manager_style* ARRAY_PLUS_GRID = 0;
  const memory_manager_style* MALLOC_MANAGER = 0;
  const memory_manager_style* HEAP_MANAGER = 0;
  const memory_manager_style* FREELISTS = 0;
};




MEDDLY::memman_initializer::memman_initializer(initializer_list *p)
 : initializer_list(p)
{
  original_grid = 0;
  array_plus_grid = 0;
  malloc_manager = 0;
  heap_manager = 0;
  freelists = 0;
}

void MEDDLY::memman_initializer::setup()
{
  ORIGINAL_GRID = (original_grid  = new orig_grid_style("ORIGINAL_GRID"));
  ARRAY_PLUS_GRID = (array_plus_grid  = new array_grid_style("ARRAY_PLUS_GRID"));
  MALLOC_MANAGER = (malloc_manager = new malloc_style("MALLOC_MANAGER"));
  HEAP_MANAGER = (heap_manager = new heap_style("HEAP_MANAGER"));
  FREELISTS = (freelists = new freelist_style("FREELISTS"));
}

void MEDDLY::memman_initializer::cleanup()
{
  delete original_grid;
  delete array_plus_grid;
  delete malloc_manager;
  delete heap_manager;
  delete freelists;
  ORIGINAL_GRID = (original_grid  = 0);
  ARRAY_PLUS_GRID = (array_plus_grid  = 0);
  MALLOC_MANAGER = (malloc_manager  = 0);
  HEAP_MANAGER = (heap_manager  = 0);
  FREELISTS = (freelists = 0);
}

