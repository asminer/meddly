#include "impl_unique_table.h"
#include "relation_node.h"

#define TERMINAL_NODE 1

MEDDLY::impl_unique_table::impl_unique_table(expert_forest* ef)
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
