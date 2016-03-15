
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

#include <cstdlib>
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

//	int size=getDomain()->getNumVariables();
//	for(int i=size; i>=1; i--) {
//		printf("Lv %d Var %d: %d\n", i, getVarByLevel(i), unique->getNumEntries(getVarByLevel(i)));
//	}
	printf("#Node: %d\n", getCurrentNumNodes());

	resetPeakNumNodes();
	resetPeakMemoryUsed();

	if(isLowestInversion()) {
		reorderVariablesLowestInversion(order);
	}
	else if(isHighestInversion()) {
		reorderVariablesHighestInversion(order);
	}
	else if(isSinkDown()) {
		reorderVariablesSinkDown(order);
	}
	else if(isBubbleUp()) {
		reorderVariablesBubbleUp(order);
	}
	else if(isLowestCost()) {
		reorderVariablesLowestCost(order);
	}
	else if(isLowestMemory()) {
		reorderVariablesLowestMemory(order);
	}
	else if(isRandom()) {
		reorderVariablesRandom(order);
	}
	else if(isLARC()) {
		reorderVariablesLARC(order);
	}

//	for(int i=size; i>=1; i--) {
//		printf("Lv %d Var %d: %d\n", i, getVarByLevel(i), unique->getNumEntries(getVarByLevel(i)));
//	}
	printf("#Node: %d\n", getCurrentNumNodes());
	printf("Peak #Node: %d\n", getPeakNumNodes());
	printf("Peak Memory: %ld\n", getPeakMemoryUsed());
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

void MEDDLY::mtmdd_forest::reorderVariablesSinkDown(const int* order)
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
	IndexedHeap<long, less<long> > heap(size);

	for(int i=1; i<size; i++) {
		if(order[getVarByLevel(i)] > order[getVarByLevel(i+1)]) {
			long weight = calculate_swap_cost(i);
			heap.push(i, weight);
		}
	}

	int swap = 0;
	while(!heap.empty()) {
		int level = heap.top_key();
		swapAdjacentVariables(level);
		swap++;
		heap.pop();

		if(level<size-1 && (order[getVarByLevel(level+1)] > order[getVarByLevel(level+2)])) {
			long weight = calculate_swap_cost(level+1);
			heap.push(level+1, weight);
		}
		if(level>1 && (order[getVarByLevel(level-1)] > order[getVarByLevel(level)])) {
			long weight = calculate_swap_cost(level-1);
			heap.push(level-1, weight);
		}
	}

	printf("Total Swap: %d\n", swap);
}

long MEDDLY::mtmdd_forest::calculate_swap_cost(int level)
{
	int lvar=getVarByLevel(level);
	int hvar=getVarByLevel(level+1);
	return static_cast<long>(unique->getNumEntries(hvar))*getVariableSize(hvar);
//	return unique->getSize(high_var) + unique->getSize(low_var);
}

void MEDDLY::mtmdd_forest::reorderVariablesLARC(const int* order)
{
	int size = getDomain()->getNumVariables();
	IndexedHeap<long, less<double> > heap(size);

	for(int i=1; i<size; i++) {
		if(order[getVarByLevel(i)] > order[getVarByLevel(i+1)]) {
			double weight = calculate_avg_ref_count(i);
			heap.push(i, weight);
		}
	}

	int swap = 0;
	while(!heap.empty()) {
		int level = heap.top_key();
		swapAdjacentVariables(level);
		swap++;
		heap.pop();

		if(level<size-1 && (order[getVarByLevel(level+1)] > order[getVarByLevel(level+2)])) {
			double weight = calculate_avg_ref_count(level+1);
			heap.push(level+1, weight);
		}
		if(level>1 && (order[getVarByLevel(level-1)] > order[getVarByLevel(level)])) {
			double weight = calculate_avg_ref_count(level-1);
			heap.push(level-1, weight);
		}
	}

	printf("Total Swap: %d\n", swap);
}

double MEDDLY::mtmdd_forest::calculate_avg_ref_count(int level)
{
	int lvar=getVarByLevel(level);
	int lnum=unique->getNumEntries(lvar);
	if(lnum==0) {
		return 0;
	}
	else {
		node_handle* lnodes=static_cast<node_handle*>(malloc(lnum*sizeof(node_handle)));
		unique->getItems(lvar, lnodes, lnum);

		int edges=0;
		for(int i=0; i<lnum; i++){
			edges+=getInCount(lnodes[i]);
		}

		free(lnodes);
		return (double)edges/lnum;
	}
}

void MEDDLY::mtmdd_forest::reorderVariablesLowestMemory(const int* order)
{
	int size = getDomain()->getNumVariables();
	IndexedHeap<long, less<long> > heap(size);

	for(int i=1; i<size; i++) {
		if(order[getVarByLevel(i)] > order[getVarByLevel(i+1)]) {
			long cost = calculate_swap_memory_cost(i);
			heap.push(i, cost);
		}
	}

	int swap = 0;
	bool swapped[] = {false, false};
	int saved = 0;
	while(!heap.empty()) {
		int level = heap.top_key();
		if(swapped[0] || swapped[1]) {
			saved++;
			swapped[0] = false;
			swapped[1] = false;
		}
		else {
			swapAdjacentVariables(level);
		}
		swap++;
		heap.pop();

		if(level<size-1 && (order[getVarByLevel(level+1)] > order[getVarByLevel(level+2)])) {
			int lvar=getVarByLevel(level+1);
			int hvar=getVarByLevel(level+2);
			int before=unique->getNumEntries(hvar)+unique->getNumEntries(lvar);

			swapAdjacentVariables(level+1);
			int after=unique->getNumEntries(hvar)+unique->getNumEntries(lvar);

			long cost = static_cast<long>(after-before);
			heap.push(level+1, cost);

			swapped[0]=true;
		}
		if(level>1 && (order[getVarByLevel(level-1)] > order[getVarByLevel(level)])) {
			int lvar=getVarByLevel(level-1);
			int hvar=getVarByLevel(level);
			int before=unique->getNumEntries(hvar)+unique->getNumEntries(lvar);

			swapAdjacentVariables(level-1);
			int after=unique->getNumEntries(hvar)+unique->getNumEntries(lvar);

			long cost = static_cast<long>(after-before);
			heap.push(level-1, cost);

			swapped[1]=true;
		}

		if(swapped[0] && heap.top_key()!=level+1) {
			// Undo
			swapAdjacentVariables(level+1);
			swapped[0]=false;
		}

		if(swapped[1] && heap.top_key()!=level-1) {
			// Undo
			swapAdjacentVariables(level-1);
			swapped[1]=false;
		}
	}

	printf("Total Swap: %d\n", swap);
	printf("Saved Swap: %d\n", saved);
}

long MEDDLY::mtmdd_forest::calculate_swap_memory_cost(int level)
{
	int lvar=getVarByLevel(level);
	int hvar=getVarByLevel(level+1);
	int before=unique->getNumEntries(hvar)+unique->getNumEntries(lvar);

	swapAdjacentVariables(level);
	int after=unique->getNumEntries(hvar)+unique->getNumEntries(lvar);
	swapAdjacentVariables(level);

	return static_cast<long>(after-before);
}

void MEDDLY::mtmdd_forest::reorderVariablesRandom(const int* order)
{
	int size = getDomain()->getNumVariables();
	std::vector<bool> inversions(size+1, false);
	std::vector<int> levels;

	for(int i=1; i<size; i++){
		if(order[getVarByLevel(i)] > order[getVarByLevel(i+1)]) {
			inversions[i]=true;
			levels.push_back(i);
		}
	}

	srand(time(NULL));
	int seed = rand();
	printf("Seed: %d\n", seed);
	srand(seed);

	int swap = 0;
	while(!levels.empty()){
		int index = rand()%levels.size();
		int level = levels[index];
		MEDDLY_DCASSERT(inversions[level]);
		swapAdjacentVariables(level);
		swap++;

		inversions[level] = false;
		if(level>1){
			if(!inversions[level-1] && (order[getVarByLevel(level-1)] > order[getVarByLevel(level)])){
				// New inversion at lower level
				inversions[level-1] = true;
				levels.push_back(level-1);
			}
		}
		if(level<size-1){
			if(!inversions[level+1] && (order[getVarByLevel(level+1)] > order[getVarByLevel(level+2)])){
				// New inversion at upper level
				inversions[level+1] = true;
				levels.push_back(level+1);
			}
		}
		levels[index]=levels.back();
		levels.pop_back();
	}

	printf("Total Swap: %d\n", swap);
}

void MEDDLY::mtmdd_forest::swapAdjacentVariables(int level)
{
	MEDDLY_DCASSERT(level>=1);
	MEDDLY_DCASSERT(level<getNumVariables());

	int hvar = getVarByLevel(level+1);  // The variable at the higher level
	int lvar = getVarByLevel(level);    // The variable at the lower level
	int hsize = getVariableSize(hvar);  // The size of the variable at the higher level
	int lsize = getVariableSize(lvar);  // The size of the variable at the lower level

	int hnum = unique->getNumEntries(hvar); // The number of nodes associated with the variable at the higher level
	node_handle* hnodes = static_cast<node_handle*>(malloc(hnum*sizeof(node_handle)));
	unique->getItems(hvar, hnodes, hnum);

	int lnum = unique->getNumEntries(lvar); // The nubmer of nodes associated with the variable at the lower level
	node_handle* lnodes = static_cast<node_handle*>(malloc(lnum*sizeof(node_handle)));
	unique->getItems(lvar, lnodes, lnum);

//	printf("Before: Level %d : %d, Level %d : %d\n",
//			level+1, hnum,
//			level, lnum);

	int num = 0;
	// Renumber the level of nodes for the variable to be moved down
	for(int i=0; i<hnum; i++) {
		node_reader* nr = initNodeReader(hnodes[i], true);
		MEDDLY_DCASSERT(nr->getLevel() == level+1);
		MEDDLY_DCASSERT(nr->getSize() == hsize);

		for(int j=0; j<hsize; j++){
			if(isLevelAbove(getNodeLevel(nr->d(j)), level-1)){
				// Remove the nodes corresponding to functions that
				// are independent of the variable to be moved up
				hnodes[num++] = hnodes[i];
				break;
			}
		}
		node_reader::recycle(nr);

		setNodeLevel(hnodes[i], level);
	}
	hnum = num;

	// Renumber the level of nodes for the variable to be moved up
	for(int i=0; i<lnum; i++) {
		setNodeLevel(lnodes[i], level+1);
	}

	// Update the variable order
	order_var[hvar] = level;
	order_var[lvar] = level+1;
	order_level[level+1] = lvar;
	order_level[level] = hvar;

	node_handle** children = static_cast<node_handle**>(malloc(hsize*sizeof(node_handle*)));
	for(int i=0; i<hsize; i++) {
		children[i] = static_cast<node_handle*>(malloc(lsize*sizeof(node_handle)));
	}

	// Process the rest of nodes for the variable to be moved down
	for(int i=0; i<hnum; i++) {
		node_reader* high_nr = initNodeReader(hnodes[i], true);
		node_builder& high_nb = useNodeBuilder(level+1, lsize);

		for(int j=0; j<hsize; j++) {
			if(isLevelAbove(level, getNodeLevel(high_nr->d(j)))) {
				for(int k=0; k<lsize; k++) {
					children[j][k] = high_nr->d(j);
				}
			}
			else {
				node_reader* nr = initNodeReader(high_nr->d(j), true);
				MEDDLY_DCASSERT(nr->getSize()==lsize);
				for(int k=0; k<lsize; k++) {
					children[j][k] = nr->d(k);
				}
				node_reader::recycle(nr);
			}
		}

		for(int j=0; j<lsize; j++) {
			node_builder& low_nb = useNodeBuilder(level, hsize);
			for(int k=0; k<hsize; k++) {
				low_nb.d(k) = linkNode(children[k][j]);
			}
			high_nb.d(j) = createReducedNode(-1, low_nb);
		}

		node_reader::recycle(high_nr);

		// The reduced node of high_nb must be at level+1
		// Assume the reduced node is at level
		// Then high_nodes[i] corresponds to a function that
		// is independent of the variable to be moved up
		// This is a contradiction
		modifyReducedNodeInPlace(high_nb, hnodes[i]);
	}

	for(int i=0; i<hsize; i++) {
		free(children[i]);
	}
	free(children);

	free(hnodes);
	free(lnodes);

//	printf("After: Level %d : %d, Level %d : %d\n",
//			level+1, unique->getNumEntries(lvar),
//			level, unique->getNumEntries(hvar));
//	printf("#Node: %d\n", getCurrentNumNodes());
}

void MEDDLY::mtmdd_forest::moveDownVariable(int high, int low)
{
	MEDDLY_DCASSERT(low<high);
	MEDDLY_DCASSERT(low>=1);
	MEDDLY_DCASSERT(high<=getNumVariables());

	removeAllComputeTableEntries();

	for(int level=high-1; level>=low; level--) {
		swapAdjacentVariables(level);
	}
}

void MEDDLY::mtmdd_forest::moveUpVariable(int low, int high)
{
	MEDDLY_DCASSERT(low<high);
	MEDDLY_DCASSERT(low>=1);
	MEDDLY_DCASSERT(high<=getNumVariables());

	removeAllComputeTableEntries();

	for(int level=low; level<high; level++) {
		swapAdjacentVariables(level);
	}
}

void MEDDLY::mtmdd_forest::dynamicReorderVariables(int top, int bottom)
{
	MEDDLY_DCASSERT(top > bottom);
	MEDDLY_DCASSERT(top <= getNumVariables());
	MEDDLY_DCASSERT(bottom >= 1);

	removeAllComputeTableEntries();

	vector<int> vars(top - bottom + 1);
	for(int i = 0; i < vars.size(); i++){
		vars[i] = getVarByLevel(bottom + i);
	}

	for(int i = 0; i < vars.size() / 4; i++){
		int max = i;
		size_t max_num = unique->getNumEntries(vars[max]);
		for(int j = i+1; j < vars.size(); j++){
			if(unique->getNumEntries(vars[j]) > max_num){
				max = j;
				max_num = unique->getNumEntries(vars[j]);
			}
		}

		int temp = vars[max];
		vars[max] = vars[i];
		vars[i] = temp;

		sifting(vars[i], top, bottom);
//		printf("%d : %d\n", vars[i], getCurrentNumNodes());
	}
}

void MEDDLY::mtmdd_forest::sifting(int var, int top, int bottom)
{
	int level = getLevelByVar(var);

	MEDDLY_DCASSERT(level <= top && level >= bottom);

	int num = getCurrentNumNodes();
	if(level <= (top + bottom) / 2) {
		// Move to the bottom
		while(level > bottom) {
//			int low_var = getVarByLevel(level - 1);
//			int old_sum = unique->getNumEntries(var) + unique->getNumEntries(low_var);
			swapAdjacentVariables(level - 1);
//			int new_sum = unique->getNumEntries(var) + unique->getNumEntries(low_var);
//			change += (new_sum - old_sum);
			level--;
//
//			if(change < min) {
//				min_level = level;
//				min = change;
//			}
		}

		int change = 0;
		int min_level = bottom;

		MEDDLY_DCASSERT(level == bottom);
		// Move to the top
		while(level < top) {
			int high_var = getVarByLevel(level + 1);
			size_t old_sum = unique->getNumEntries(var) + unique->getNumEntries(high_var);
			swapAdjacentVariables(level);
			size_t new_sum = unique->getNumEntries(var) + unique->getNumEntries(high_var);
			change += (new_sum - old_sum);
			level++;

			if(change <= 0) {
				min_level = level;
				change = 0;
			}
		}

		MEDDLY_DCASSERT(level == top);
		while(level > min_level) {
			swapAdjacentVariables(level - 1);
			level--;
		}
	}
	else {
		// Move to the top
		while(level < top) {
//			int high_var = getVarByLevel(level + 1);
//			int old_sum = unique->getNumEntries(var) + unique->getNumEntries(high_var);
			swapAdjacentVariables(level);
//			int new_sum = unique->getNumEntries(var) + unique->getNumEntries(high_var);
//			change += (new_sum - old_sum);
			level++;
//
//			if(change < min) {
//				min_level = level;
//				min = change;
//			}
		}

		int change = 0;
		int min_level = top;
		int min = change;

		MEDDLY_DCASSERT(level == top);
		// Move to the bottom
		while(level > bottom) {
			int low_var = getVarByLevel(level - 1);
			size_t old_sum = unique->getNumEntries(var) + unique->getNumEntries(low_var);
			swapAdjacentVariables(level - 1);
			size_t new_sum = unique->getNumEntries(var) + unique->getNumEntries(low_var);
			change += (new_sum - old_sum);
			level--;

			if(change <= min) {
				min_level = level;
				min = change;
			}
		}

		MEDDLY_DCASSERT(level == bottom);
		MEDDLY_DCASSERT(min <= 0);
		while(level < min_level) {
			swapAdjacentVariables(level);
			level++;
		}
	}

	MEDDLY_DCASSERT(getCurrentNumNodes() <= num);
	if(getCurrentNumNodes() > num) {
		printf("Error: %d > %d\n", getCurrentNumNodes(), num);
	}
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
      index[k] = path[k].i(nzp[k]);
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
