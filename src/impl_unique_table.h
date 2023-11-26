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

#ifndef MEDDLY_IMPL_UNIQUE_TABLE_H
#define MEDDLY_IMPL_UNIQUE_TABLE_H

#include "defines.h"
#include "forest.h"

#include <unordered_map>

namespace MEDDLY {
  class impl_unique_table;
};

/** Unique table for discovering duplicate implcit nodes.
    Used for defining events in dicsrete-state systems

 */
class MEDDLY::impl_unique_table {
  public:
  impl_unique_table(forest *ef);
  ~impl_unique_table();

  inline unsigned getNumEntries() { return table.size(); }
  inline unsigned getSize() { return size; }
  inline unsigned getLastHandle() { return last_handle; }

  /** Add a node to the unique table.
      Returns handle to the new node.
      Otherwise, return the handle to the existing node.
   */
  node_handle add(node_handle rnh, relation_node* rnb);


  /** Get the node asscoiated with the handle
   */
  relation_node* getNode(node_handle rnh);

  /** Check if node already exist in table
      @param rnb  The relation node.
      @return     If unique, 0
                  Else, existing node handle
   */
  node_handle isDuplicate(relation_node* rnb);

  private:
    // To which forest does this unique table belong to
    forest *parent;

    // The contents of the table
    std::unordered_map<node_handle, relation_node*> table;

    //Top nodes of each event
    std::vector<std::vector<node_handle>> levelTopEvent;

    // Unique Table Stats
    unsigned num_entries;
    unsigned size;
    unsigned last_handle;
  };

#endif // #include guard
