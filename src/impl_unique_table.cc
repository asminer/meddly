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

#include "impl_unique_table.h"
#include "relation_node.h"

#define TERMINAL_NODE 1

MEDDLY::impl_unique_table::impl_unique_table(forest* ef)
: parent(ef)
{
  size = 0;
  num_entries = 0;
  last_handle = 0;
}

MEDDLY::impl_unique_table::~impl_unique_table()
{
  table.clear();
}

MEDDLY::node_handle
MEDDLY::impl_unique_table::add(node_handle rnh, relation_node *rnb)
{
      std::pair<node_handle, relation_node*> newNode(rnh,rnb);
      table.insert(newNode);
      return rnh;
}

MEDDLY::relation_node*
MEDDLY::impl_unique_table::getNode(node_handle rnh)
{
  std::unordered_map<node_handle, relation_node*>::iterator it = table.find(rnh);
  if(it!=table.end())
    return it->second;
  else
    return NULL;
}

MEDDLY::node_handle
MEDDLY::impl_unique_table::isDuplicate(relation_node *rnb)
{
  std::unordered_map<node_handle, relation_node*>::iterator it = table.begin();
  while(it != table.end())
  {
    if((it->second)->equals(rnb))
      return (it->second)->getID();
    ++it;
  }
  return 0;
}
