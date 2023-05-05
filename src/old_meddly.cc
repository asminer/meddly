
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

/*
    Implementation of the "global" functions given in
    meddly.h and meddly_expert.h.

*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "defines.h"

#include "variable.h"
#include "domain.h"
#include "unpacked_node.h"

#include "old_meddly.h"
// #include "compute_table.h"
#include "memory_managers/init_managers.h"
#include "forests/init_forests.h"
#include "storage/init_storage.h"

#include "ct_initializer.h"
#include "compute_table.h"

#include "ops_builtin.h"

#include "oper.h"
#include "oper_unary.h"
#include "oper_binary.h"
#include "oper_special.h"

// #define STATS_ON_DESTROY

namespace MEDDLY {
  // "global" variables


  //
  // List of all domains
  //
  domain** domain::dom_list = 0;
  int* domain::dom_free = 0;
  int domain::dom_list_size = 0;
  int domain::free_list = -1;

  //
  // List of free unpacked nodes
  unpacked_node* unpacked_node::freeList = 0;
  unpacked_node* unpacked_node::buildList = 0;

};


