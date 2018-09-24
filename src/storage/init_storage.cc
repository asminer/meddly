
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
#include "init_storage.h"

#include "simple.h"
#include "pattern.h"
#include "best.h"


namespace MEDDLY {
  const node_storage_style* SIMPLE_STORAGE = 0;
  const node_storage_style* PATTERN_STORAGE = 0;
  const node_storage_style* BEST_STORAGE = 0;
  
};

MEDDLY::storage_initializer::storage_initializer(initializer_list *p)
: initializer_list(p)
{
  simple = 0;
  pattern = 0;
  best = 0;
}

void MEDDLY::storage_initializer::setup()
{
  SIMPLE_STORAGE = (simple = new simple_separated_style("SIMPLE_STORAGE")); 
  PATTERN_STORAGE = (pattern = new pattern_storage_style("PATTERN_STORAGE")); 
  BEST_STORAGE = (best = new best_storage_style("BEST_STORAGE"));
}

void MEDDLY::storage_initializer::cleanup()
{
  delete simple;
  SIMPLE_STORAGE = (simple = 0);
  
  delete pattern;
  PATTERN_STORAGE = (pattern = 0);
  
  delete best;
  BEST_STORAGE = (best = 0);
}

