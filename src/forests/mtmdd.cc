
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
#include <ext/hash_map>

#include "mtmdd.h"
#include "../unique_table.h"
#include "../hash_stream.h"
#include "../heap.h"

MEDDLY::mtmdd_forest
::mtmdd_forest(int dsl, domain* d, range_type t, const policies &p)
 : mt_forest(dsl, d, false, t, p)
{
  // anything to construct?
}

void MEDDLY::mtmdd_forest::reorderVariables(const int* order)
{
	removeAllComputeTableEntries();

	int size=getDomain()->getNumVariables();
	for(int i=1; i<=size; i++) {
		printf("Lv %d: %d\n", i, unique->getNumEntries(getVarByLevel(i)));
	}
	printf("#Node: %d\n", getCurrentNumNodes());

	if(isLowestInversion()) {
		reorderVariablesLowestInversion(order);
	}
	else if(isHighestInversion()) {
		reorderVariablesHighestInversion(order);
	}
	else if(isBubbleDown()) {
		reorderVariablesBubbleDown(order);
	}
	else if(isBubbleUp()) {
		reorderVariablesBubbleUp(order);
	}
	else if(isLowestCost()) {
		reorderVariablesLowestCost(order);
	}

//	for(int i=1; i<=size; i++) {
//		printf("Lv %d: %d\n", i, unique->getNumEntries(getVarByLevel(i)));
//	}
	printf("#Node: %d\n", getCurrentNumNodes());
}

void MEDDLY::mtmdd_forest::reorderVariablesLowestInversion(const int* order)
{
	int size = getDomain()->getNumVariables();

	// The variables below ordered_level are ordered
	int ordered_level = 1;
	int level = 1;
	int swap = 0;
	while (ordered_level<size) {
		level = ordered_level;
		while(level>0 && (order[getVarByLevel(level)] > order[getVarByLevel(level+1)])) {
			swapAdjacentVariables(level);
			level--;
			swap++;
		}
		ordered_level++;
	}

	printf("Total Swap: %d\n", swap);
}

void MEDDLY::mtmdd_forest::reorderVariablesHighestInversion(const int* order)
{
	int size = getDomain()->getNumVariables();

	// The variables above ordered_level are ordered
	int ordered_level = size-1;
	int level = size-1;
	int swap = 0;
	while (ordered_level>0) {
		level = ordered_level;
		while(level<size && (order[getVarByLevel(level)] > order[getVarByLevel(level+1)])) {
			swapAdjacentVariables(level);
			level++;
			swap++;
		}
		ordered_level--;
	}

	printf("Total Swap: %d\n", swap);
}

void MEDDLY::mtmdd_forest::reorderVariablesBubbleDown(const int* order)
{
	int size = getDomain()->getNumVariables();

	// Construct the mapping from level to variable
	int* level_to_var = static_cast<int*>(calloc(sizeof(int), size+1));
	for(int i=1; i<size+1; i++) {
		level_to_var[order[i]] = i;
	}

	int swap=0;
	for(int i=1; i<size; i++) {
		int level = getLevelByVar(level_to_var[i]);
		while(level>i) {
			swapAdjacentVariables(level-1);
			level--;
			swap++;
		}
	}

	free(level_to_var);
	printf("Total Swap: %d\n", swap);
}

void MEDDLY::mtmdd_forest::reorderVariablesBubbleUp(const int* order)
{
	int size = getDomain()->getNumVariables();

	// Construct the mapping from level to variable
	int* level_to_var = static_cast<int*>(calloc(sizeof(int), size+1));
	for(int i=1; i<size+1; i++) {
		level_to_var[order[i]] = i;
	}

	int swap=0;
	for(int i=size; i>1; i--) {
		int level = getLevelByVar(level_to_var[i]);
		while(level<i) {
			swapAdjacentVariables(level);
			level++;
			swap++;
		}
	}

	free(level_to_var);
	printf("Total Swap: %d\n", swap);
}

void MEDDLY::mtmdd_forest::reorderVariablesLowestCost(const int* order)
{
	int size = getDomain()->getNumVariables();
	InversionHeap heap(size);

	for(int i=1; i<size; i++) {
		if(order[getVarByLevel(i)] > order[getVarByLevel(i+1)]) {
			long cost = calculate_swap_cost(i);
			heap.push(i, cost);
		}
	}

	int swap = 0;
	while(!heap.empty()) {
		int level = heap.top();
		swapAdjacentVariables(level);
		swap++;
		heap.pop();

		if(level<size-1 && (order[getVarByLevel(level+1)] > order[getVarByLevel(level+2)])) {
			long cost = calculate_swap_cost(level+1);
			heap.push(level+1, cost);
		}
		if(level>1 && (order[getVarByLevel(level-1)] > order[getVarByLevel(level)])) {
			long cost = calculate_swap_cost(level-1);
			heap.push(level-1, cost);
		}
	}

	printf("Total Swap: %d\n", swap);
}

long MEDDLY::mtmdd_forest::calculate_swap_cost(int level)
{
//	int low_var=getVarByLevel(level);
//	int high_var=getVarByLevel(level+1);
//	return unique->getSize(high_var)*getVariableSize(high_var)*getVariableSize(low_var)
//			+ unique->getSize(low_var);
	int low_var=getVarByLevel(level);
	int high_var=getVarByLevel(level+1);
	return static_cast<long>(unique->getNumEntries(high_var))*getVariableSize(high_var)*getVariableSize(low_var)
			+ unique->getNumEntries(low_var);
//	return unique->getSize(high_var) + unique->getSize(low_var);
}

void MEDDLY::mtmdd_forest::swapAdjacentVariables(int level)
{
	MEDDLY_DCASSERT(level>=1);
	MEDDLY_DCASSERT(level<getNumVariables());

	removeAllComputeTableEntries();

	int high_var=getVarByLevel(level+1);
	int low_var=getVarByLevel(level);
	int high_size=getVariableSize(high_var);
	int low_size=getVariableSize(low_var);

	int high_node_size=unique->getNumEntries(high_var);
	node_handle* high_nodes=static_cast<node_handle*>(malloc(high_node_size*sizeof(node_handle)));
	unique->getItems(high_var, high_nodes, high_node_size);

	int low_node_size=unique->getNumEntries(low_var);
	node_handle* low_nodes=static_cast<node_handle*>(malloc(low_node_size*sizeof(node_handle)));
	unique->getItems(low_var, low_nodes, low_node_size);

	printf("Before: Level %d : %d, Level %d : %d\n",
			level+1, high_node_size,
			level, low_node_size);

	int i=0, j=0;
	// Renumber the level of nodes for the variable to be moved down
	for(i=0; i<high_node_size; i++) {
		node_reader* nr = initNodeReader(high_nodes[i], true);
		MEDDLY_DCASSERT(nr->getLevel()==level+1);
		MEDDLY_DCASSERT(nr->getSize()==high_size);

		for(int k=0; k<high_size; k++){
			if(getNodeLevel(nr->d(k))>=level){
				// Remove the nodes corresponding to functions that
				// are independent of the variable to be moved up
				high_nodes[j++]=high_nodes[i];
				break;
			}
		}
		node_reader::recycle(nr);

		setNodeLevel(high_nodes[i], level);
	}
	high_node_size=j;

	// Renumber the level of nodes for the variable to be moved up
	for(i=0; i<low_node_size; i++) {
		setNodeLevel(low_nodes[i], level+1);
	}

	// Update the variable order
	order_var[high_var] = level;
	order_var[low_var] = level+1;
	order_level[level+1] = low_var;
	order_level[level] = high_var;

	// Process the rest of nodes for the variable to be moved down
	for(i=0; i<high_node_size; i++) {
		node_reader* high_nr = initNodeReader(high_nodes[i], true);
		node_builder& high_nb = useNodeBuilder(level+1, low_size);
		for(j=0; j<low_size; j++) {
			node_builder& low_nb = useNodeBuilder(level, high_size);
			for(int k=0; k<high_size; k++) {
				if(getNodeLevel(high_nr->d(k))<level){
					low_nb.d(k)=linkNode(high_nr->d(k));
				}
				else{
					node_reader* low_nr = initNodeReader(high_nr->d(k), true);
					MEDDLY_DCASSERT(low_nr->getSize()==low_size);
					low_nb.d(k)=linkNode(low_nr->d(j));
					node_reader::recycle(low_nr);
				}
			}
			high_nb.d(j)=createReducedNode(-1, low_nb);
		}

//		for(int k=0; k<high_size; k++) {
//			unlinkNode(high_nr->d(k));
//		}
		node_reader::recycle(high_nr);

		unique->remove(hashNode(high_nodes[i]), high_nodes[i]);
		// The reduced node of high_nb must be at level+1
		// Assume the reduced node is at level
		// Then high_nodes[i] corresponds to a function that
		// is independent of the variable to be moved up
		// This is a contradiction
		modifyReducedNodeInPlace(high_nb, high_nodes[i]);
	}

	free(high_nodes);
	free(low_nodes);

	printf("After: Level %d : %d, Level %d : %d\n",
			level+1, unique->getNumEntries(low_var),
			level, unique->getNumEntries(high_var));
	printf("#Node: %d\n", getCurrentNumNodes());
}

typedef struct ProcessedItem
{
	MEDDLY::node_handle node;
	int index;

	bool operator==(const ProcessedItem& item) const
	{
		return (node==item.node) && (index==item.index);
	}
}ProcessedItem;

typedef struct HashProcessedItem
{
    size_t operator()(const ProcessedItem& item) const
    {
    	MEDDLY::hash_stream s;
    	s.start(0);
    	s.push(item.node, item.index);
    	return s.finish();
    }
}HashProcessedItem;

// Temporary solution
// Need to redesign in future
typedef __gnu_cxx::hash_map<ProcessedItem, MEDDLY::node_handle, HashProcessedItem> ProcessedHashMap;
static ProcessedHashMap processed;

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

	// Renumber the level of nodes
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

	int var = getVarByLevel(low);
	int size = unique->getNumEntries(var);
	node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
	unique->getItems(var, nodes, size);
	unique->clear(var);
	int varSize = getLevelSize(low);

	// Sort the nodes for the variable to be moved down
	// by the maximum level of their children in ascending order
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
		// If not, h1 will not be recycled and will stay in the unique table
		// Will cause problems when updating node in place
		// because node and h1 are the same
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

//#ifdef __USE_CACHE
//	processed.clear();
//#endif
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

//#ifdef __USE_CACHE
//	ProcessedHashMap::const_iterator itr=processed.find(ProcessedItem{node, val});
//	if(itr!=processed.end()) {
//		MEDDLY_DCASSERT(isActiveNode(itr->second));
//		return linkNode(itr->second);
//	}
//#endif

	node_builder& nb = useNodeBuilder(level, getLevelSize(level));
	node_reader* nr = initNodeReader(node, true);
	int varSize = getLevelSize(level);
	for(int v=0; v<varSize; v++) {
		nb.d(v) = recursiveReduceDown(nr->d(v), low, val);
	}
	node_reader::recycle(nr);

	MEDDLY::node_handle h = createReducedNode(-1, nb);
//#ifdef __USE_CACHE
//	processed.insert(__gnu_cxx::pair<ProcessedItem, node_handle>(ProcessedItem{node, val}, h));
//#endif
	return h;
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

		int varSize = getLevelSize(level-1);
		for(int i=0; i<size; i++) {
			MEDDLY_DCASSERT(isActiveNode(nodes[i]));
			MEDDLY_DCASSERT(getNodeLevel(nodes[i])==level);

			node_reader* nr = initNodeReader(nodes[i], true);
			for(int v=0; v<varSize; v++) {
				int childLevel=getNodeLevel(nr->d(v));
				if(childLevel==low || (childLevel>low && getNode(nr->d(v)).isMarked())) {
					getNode(nodes[i]).mark();
					break;
				}
			}
			node_reader::recycle(nr);
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

	// Process the nodes corresponding to functions
	// depending on the variable to be moved up
	int varSize = getLevelSize(high);
	for(int level=high; level>low; level--) {
		int var = getVarByLevel(level-1);
		int size = unique->getNumEntries(var);
		node_handle* nodes = static_cast<node_handle*>(malloc(size*sizeof(node_handle)));
		unique->getItems(var, nodes, size);

		for(int i=0; i<size; i++) {
			node_handle node = nodes[i];
			if(getNode(node).isMarked()) {
				node_builder& nb = useNodeBuilder(high, varSize);
				for(int val=0; val<varSize; val++) {
					nb.d(val) = recursiveReduceUp(node, low, high, val);
				}
				node_reader* nr = initNodeReader(node, true);
				for(int val=0; val<getLevelSize(level-1); val++) {
					unlinkNode(nr->d(val));
				}
				node_reader::recycle(nr);

				// The reduced node of nb must be at high
				// Assume the reduced node is below high and above low
				// Then node corresponds to a function that
				// is independent of the variable to be moved up
				// This is a contradiction
				unique->remove(hashNode(node), node);
				modifyReducedNodeInPlace(nb, node);
			}
		}
		free(nodes);
	}

	processed.clear();
}

MEDDLY::node_handle MEDDLY::mtmdd_forest::recursiveReduceUp(node_handle node, int low, int high, int val)
{
	int level = getNodeLevel(node);
	if(level < low || (level!=high && !getNode(node).isMarked())) {
		return linkNode(node);
	}

	if(level == high) {
		node_reader* nr = initNodeReader(node, true);
		node_handle h = linkNode(nr->d(val));
		node_reader::recycle(nr);
		return h;
	}

	ProcessedHashMap::const_iterator itr=processed.find(ProcessedItem{node, val});
	if(itr!=processed.end()) {
		MEDDLY_DCASSERT(isActiveNode(itr->second));
		return linkNode(itr->second);
	}

	node_reader* nr = initNodeReader(node, true);
	int varSize = getLevelSize(level);
	node_builder& nb = useNodeBuilder(level, varSize);
	for(int v=0; v<varSize; v++) {
		nb.d(v) = recursiveReduceUp(nr->d(v), low, high, val);
	}
	node_reader::recycle(nr);

	node_handle h = createReducedNode(-1, nb);
	if(getInCount(node)>1) {
		processed.insert(__gnu_cxx::pair<ProcessedItem, node_handle>(ProcessedItem{node, val}, h));
	}

	return h;
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
