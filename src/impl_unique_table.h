#ifndef IMPL_UNIQUE_TABLE_H
#define IMPL_UNIQUE_TABLE_H

#include "defines.h"

namespace MEDDLY {
  class impl_unique_table;
};

/** Unique table for discovering duplicate implcit nodes.
    Used for defining events in dicsrete-state systems

 */
class MEDDLY::impl_unique_table {
  public:
  impl_unique_table(expert_forest *ef);
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
    expert_forest *parent;
  
    // The contents of the table
    std::unordered_map<node_handle, relation_node*> table;
  
    //Top nodes of each event
    std::vector<std::vector<node_handle>> levelTopEvent;
  
    // Unique Table Stats
    unsigned num_entries;
    unsigned size;
    unsigned last_handle;
  };

#endif
