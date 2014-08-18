
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

#include <vector>
#include <algorithm>

#include "mtmdd.h"
#include "../unique_table.h"

MEDDLY::mtmdd_forest
::mtmdd_forest(int dsl, domain* d, range_type t, const policies &p)
 : mt_forest(dsl, d, false, t, p)
{
  // anything to construct?
}

void MEDDLY::mtmdd_forest::swapAdjacentVariables(int level)
{
	MEDDLY_DCASSERT(level>=1);
	MEDDLY_DCASSERT(level<getNumVariables());

	removeAllComputeTableEntries();

	// Domain has already update the variable-level mapping
	int high_var=getVarByLevel(level);
	int high_node_size=unique->getNumEntries(high_var);
	node_handle* high_nodes=static_cast<node_handle*>(malloc(high_node_size*sizeof(node_handle)));
	unique->getItems(high_var, high_nodes, high_node_size);
	unique->clear(high_var);

	int low_var=getVarByLevel(level+1);
	int low_node_size=unique->getNumEntries(low_var);
	node_handle* low_nodes=static_cast<node_handle*>(malloc(low_node_size*sizeof(node_handle)));
	unique->getItems(low_var, low_nodes, low_node_size);
	unique->clear(low_var);

	int high_size=getLevelSize(level);
	int low_size=getLevelSize(level+1);

	// Identify all nodes at level+1 whose children are lower than level
	int i=0, j=0;
	for(i=0; i<high_node_size; i++) {
		node_reader* nr = initNodeReader(high_nodes[i], true);
		MEDDLY_DCASSERT(nr->getLevel()==level+1);
		MEDDLY_DCASSERT(nr->getSize()==high_size);

		bool flag=true;

		for(int k=0; k<high_size; k++){
			if(getNodeLevel(nr->d(k))>=level){
				flag=false;
				break;
			}
		}

		if(flag) {
			// Move to level
			node_builder& nb = useNodeBuilder(level, high_size);
			for(int k=0; k<high_size; k++){
				nb.d(k)=nr->d(k);
			}
			modifyReducedNodeInPlace(nb, high_nodes[i]);
		}
		else{
			// Remove the identified nodes from high_nodes
			high_nodes[j++]=high_nodes[i];
		}
	}
	high_node_size=j;

	std::vector<int> unlink;
	// Process the rest of nodes at level+1
	for(i=0; i<high_node_size; i++) {
		node_reader* high_nr = initNodeReader(high_nodes[i], true);
		node_builder& high_nb = useNodeBuilder(level+1, low_size);
		for(int k=0; k<low_size; k++) {
			node_builder& low_nb = useNodeBuilder(level, high_size);
			for(int l=0; l<high_size; l++) {
				if(getNodeLevel(high_nr->d(l))<level){
					low_nb.d(l)=linkNode(high_nr->d(l));
				}
				else{
					node_reader* low_nr = initNodeReader(high_nr->d(l), true);
					MEDDLY_DCASSERT(low_nr->getSize()==low_size);
					low_nb.d(l)=linkNode(low_nr->d(k));

					if(k==low_size-1){
						unlink.push_back(high_nr->d(l));
					}
				}
			}
			high_nb.d(k)=createReducedNode(-1, low_nb);
		}
		node_reader::recycle(high_nr);
		modifyReducedNodeInPlace(high_nb, high_nodes[i]);
	}

	for(std::vector<int>::iterator itr=unlink.begin(); itr!=unlink.end(); itr++) {
		unlinkNode(*itr);
	}

	// Identify nodes at level which are pointed from a node above level+1
	for(i=0; i<low_node_size; i++) {
		if(isActiveNode(low_nodes[i]) && getInCount(low_nodes[i])>0) {
			node_reader* nr = initNodeReader(low_nodes[i], true);
			node_builder& nb = useNodeBuilder(level+1, low_size);
			for(int k=0; k<low_size; k++) {
				// FIXME: May not increase the reference number
				nb.d(k)=linkNode(nr->d(k));
			}
			node_reader::recycle(nr);

			modifyReducedNodeInPlace(nb, low_nodes[i]);
		}
	}

	free(high_nodes);
	free(low_nodes);
}

typedef struct NodeOrder
{
	MEDDLY::node_handle node;
	int weight;
}NodeOrder;

typedef struct NodeOrderLt
{
	bool operator() (const NodeOrder& no1, const NodeOrder& no2)
	{
		return no1.weight<no2.weight;
	}
}NodeOrderLt;

void MEDDLY::mtmdd_forest::moveDownVariable(int high, int low)
{
	// Pre-condition: The compute table has been cleared

	MEDDLY_DCASSERT(low<high);
	MEDDLY_DCASSERT(low>=1);
	MEDDLY_DCASSERT(high<=getNumVariables());

	// Modify the level of nodes from level high-1 to level low
	for(int level=high-1; level>=low; level--) {
		int var = getVarByLevel(level+1);
		int size = unique->getNumEntries(var);
		node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
		unique->getItems(var, nodes, size);
		for(int i=0; i<size; i++) {
			MEDDLY_DCASSERT(isActiveNode(nodes[i]));
			MEDDLY_DCASSERT(getNodeLevel(nodes[i])==level);
			setNodeLevel(nodes[i], level+1);
		}
		free(nodes);
	}

	// Process the nodes at level high (before)
	int var = getVarByLevel(low);
	int size = unique->getNumEntries(var);
	node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
	unique->getItems(var, nodes, size);
	unique->clear(var);
	int varSize = getLevelSize(low);

	// Sort the nodes at level high (before) by the maximum level of their children
	std::vector<NodeOrder> order;
	for(int i=0; i<size; i++) {
		node_handle node = nodes[i];
		node_reader* nr = initNodeReader(node, true);
		int weight = getNodeLevel(nr->d(0));
		for(int val=1; val<varSize; val++) {
			if(weight<getNodeLevel(nr->d(val))) {
				weight = getNodeLevel(nr->d(val));
			}
		}
		node_reader::recycle(nr);

		order.push_back(NodeOrder{node, weight});
	}
	std::sort(order.begin(), order.end(), NodeOrderLt());
	for(int i=0; i<size; i++) {
		nodes[i] = order[i].node;
	}

	binary_operation* op = getOperation(UNION, this, this, this);
	for(int i=0; i<size; i++) {
		node_handle node = nodes[i];
		node_reader* nr = initNodeReader(node, true);

		node_handle h1 = recursiveReduceDown(nr->d(0), low, 0);
		unlinkNode(nr->d(0));
		for(int val=1; val<varSize; val++) {
			node_handle h2 = recursiveReduceDown(nr->d(val), low, val);
			unlinkNode(nr->d(val));
			node_handle result = op->compute(h1, h2);
			unlinkNode(h1);
			unlinkNode(h2);
			h1 = result;
		}

		node_reader::recycle(nr);
		removeAllComputeTableEntries();

		nr = initNodeReader(h1, true);
		int level = getNodeLevel(h1);
		node_builder& nb = useNodeBuilder(level, getLevelSize(level));
		for(int j=0; j<getLevelSize(level); j++) {
			nb.d(j) = linkNode(nr->d(j));
		}
		node_reader::recycle(nr);

		MEDDLY_DCASSERT(getInCount(h1)==1);
		unlinkNode(h1);

		modifyReducedNodeInPlace(nb, node);
	}

	free(nodes);
}

MEDDLY::node_handle MEDDLY::mtmdd_forest::recursiveReduceDown(node_handle node, int low, int val)
{
	int level = getNodeLevel(node);
	if(level < low) {
		if(node!=getTransparentNode()) {
			node_builder& nb = useSparseBuilder(low, 1);
			nb.i(0) = val;
			nb.d(0) = linkNode(node);
			return createReducedNode(-1, nb);
		}
		else {
			return linkNode(node);
		}
	}

	MEDDLY_DCASSERT(level>=low+1);

	node_builder& nb = useNodeBuilder(level, getLevelSize(level));
	node_reader* nr = initNodeReader(node, true);
	int varSize = getLevelSize(level);
	for(int v=0; v<varSize; v++) {
		nb.d(v) = recursiveReduceDown(nr->d(v), low, val);
	}
	node_reader::recycle(nr);

	return createReducedNode(-1, nb);
}

void MEDDLY::mtmdd_forest::moveUpVariable(int low, int high)
{
	// Pre-condition:
	// - The compute table has been cleared
	// - The change in variable-level mapping in domain is done

	MEDDLY_DCASSERT(low<high);
	MEDDLY_DCASSERT(low>=1);
	MEDDLY_DCASSERT(high<=getNumVariables());

	// Mark the nodes in the functions depending on the variable to be moved up
	for(int level=low+1; level<=high; level++) {
		int var = getVarByLevel(level-1);
		int size = unique->getNumEntries(var);
		node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
		unique->getItems(var, nodes, size);
		for(int i=0; i<size; i++) {
			MEDDLY_DCASSERT(isActiveNode(nodes[i]));
			MEDDLY_DCASSERT(getNodeLevel(nodes[i])==level);

			node_reader* nr = initNodeReader(nodes[i], true);
			bool flag = false;
			for(int v=0; v<getLevelSize(level-1); v++) {
				int childLevel=getNodeLevel(nr->d(v));
				if(childLevel==low || (childLevel>low && getNode(nr->d(v)).isMarked())) {
					flag = true;
					break;
				}
			}
			if(flag) {
				getNode(nodes[i]).mark();
			}
		}
		free(nodes);
	}

	// Modify the level of nodes
	for(int level=low+1; level<=high; level++) {
		int var = getVarByLevel(level-1);
		int size = unique->getNumEntries(var);
		node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
		unique->getItems(var, nodes, size);
		for(int i=0; i<size; i++) {
			MEDDLY_DCASSERT(isActiveNode(nodes[i]));
			MEDDLY_DCASSERT(getNodeLevel(nodes[i])==level);

			setNodeLevel(nodes[i], level-1);
		}
		free(nodes);
	}
	int var = getVarByLevel(high);
	int size = unique->getNumEntries(var);
	node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
	unique->getItems(var, nodes, size);
	for(int i=0; i<size; i++) {
		MEDDLY_DCASSERT(isActiveNode(nodes[i]));
		MEDDLY_DCASSERT(getNodeLevel(nodes[i])==low);

		setNodeLevel(nodes[i], high);
	}
	free(nodes);

	// Process the nodes from level high to level low+1 (after)
	int varSize = getLevelSize(high);
	for(int level=high; level>low; level--) {
		int var = getVarByLevel(level-1);
		int size = unique->getNumEntries(var);
		node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
		unique->getItems(var, nodes, size);

		for(int i=0; i<size; i++) {
			node_handle node = nodes[i];
			if(getNode(node).isMarked()) {
				node_builder& nb = useNodeBuilder(high, getLevelSize(high));
				for(int val=0; val<varSize; val++) {
					nb.d(val) = recursiveReduceUp(node, low, high, val);
				}
				node_reader* nr = initNodeReader(node, true);
				for(int val=0; val<getLevelSize(level-1); val++) {
					unlinkNode(nr->d(val));
				}
				node_reader::recycle(nr);

				// h may not be at level high
				node_handle h = createReducedNode(-1, nb);
				nr = initNodeReader(h, true);
				int hLevel = getNodeLevel(h);
				nb = useNodeBuilder(hLevel, getLevelSize(hLevel));
				for(int val=0; val<getLevelSize(hLevel); val++) {
					nb.d(val) = linkNode(nr->d(val));
				}
				node_reader::recycle(nr);

				MEDDLY_DCASSERT(getInCount(h)==1);
				unlinkNode(h);

				unique->remove(hashNode(node), node);
				modifyReducedNodeInPlace(nb, node);
			}
		}
		free(nodes);
	}
}

MEDDLY::node_handle MEDDLY::mtmdd_forest::recursiveReduceUp(node_handle node, int low, int high, int val)
{
	int level = getNodeLevel(node);
	if(level < low || (level!=high && !getNode(node).isMarked())) {
		return linkNode(node);
	}

	node_reader* nr = initNodeReader(node, true);
	if(level == high) {
		return linkNode(nr->d(val));
	}

	node_builder& nb = useNodeBuilder(level, getLevelSize(level));
	int varSize = getLevelSize(level);
	for(int v=0; v<varSize; v++) {
		nb.d(v) = recursiveReduceUp(nr->d(v), low, high, val);
	}
	node_reader::recycle(nr);

	return createReducedNode(-1, nb);
}

// ******************************************************************
// *                                                                *
// *              mtmdd_forest::mtmdd_iterator methods              *
// *                                                                *
// ******************************************************************

MEDDLY::mtmdd_forest::mtmdd_iterator::mtmdd_iterator(const expert_forest *F)
 : mt_iterator(F)
{
}

MEDDLY::mtmdd_forest::mtmdd_iterator::~mtmdd_iterator()
{
}

bool MEDDLY::mtmdd_forest::mtmdd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmdd_forest::mtmdd_iterator::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  int k;
  node_handle down = 0;
  for (k=1; k<=maxLevel; k++) {
    nzp[k]++;
    if (nzp[k] < path[k].getNNZs()) {
      int var = F->getVarByLevel(k);
      index[var] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      break;
    }
  }
  level_change = k;
  if (k>maxLevel) {
    return false;
  }

  return first(k-1, down);
}

bool MEDDLY::mtmdd_forest::mtmdd_iterator::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(!F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==down) return false;

  for ( ; k; k--) {
    MEDDLY_DCASSERT(down);
    int kdn = F->getNodeLevel(down);
    MEDDLY_DCASSERT(kdn <= k);
    if (kdn < k)  F->initRedundantReader(path[k], k, down, false);
    else          F->initNodeReader(path[k], down, false);

    nzp[k] = 0;

    int var = F->getVarByLevel(k);
    index[var] = path[k].i(0);
    down = path[k].d(0);
  }
  // save the terminal value
  index[0] = down;
  return true;
}
