
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
#include <map>

#include "mtmxd.h"
#include "../unique_table.h"

MEDDLY::mtmxd_forest
::mtmxd_forest(int dsl, domain* d, range_type t, const policies &p)
 : mt_forest(dsl, d, true, t, p)
{
  // anything to construct?
}

// ******************************************************************
// *                                                                *
// *              mtmxd_forest::mtmxd_iterator methods              *
// *                                                                *
// ******************************************************************

MEDDLY::mtmxd_forest::mtmxd_iterator::mtmxd_iterator(const expert_forest *F)
 : mt_iterator(F)
{
}

void MEDDLY::mtmxd_forest::reorderVariables(const int* order)
{
//	removeAllComputeTableEntries();

//	int size=getDomain()->getNumVariables();
//	for(int i=1; i<=size; i++) {
//		printf("Lv %d: %d\n", i, unique->getNumEntries(getVarByLevel(i)));
//		printf("Lv %d: %d\n", -i, unique->getNumEntries(-getVarByLevel(-i)));
//	}
//	printf("#Node: %d\n", getCurrentNumNodes());

	resetPeakNumNodes();
	resetPeakMemoryUsed();

	reorderVariablesHighestInversion(order);

//	for(int i=1; i<=size; i++) {
//		printf("Lv %d: %d\n", i, unique->getNumEntries(getVarByLevel(i)));
//	}
	printf("#Node: %d\n", getCurrentNumNodes());
	printf("Peak #Node: %d\n", getPeakNumNodes());
	printf("Peak Memory: %ld\n", getPeakMemoryUsed());
}

void MEDDLY::mtmxd_forest::reorderVariablesHighestInversion(const int* order)
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

//
// Complete adjacent variable swap by swapping two levels 4 times
// Work for fully-fully reduction only
//
void MEDDLY::mtmxd_forest::swapAdjacentVariables(int level)
{
	if(isVarSwap()){
		swapAdjacentVariablesByVarSwap(level);
	}
	else if(isLevelSwap()){
		swapAdjacentVariablesByLevelSwap(level);
	}
}

void MEDDLY::mtmxd_forest::swapAdjacentVariablesByVarSwap(int level)
{
	// Swap VarHigh and VarLow

	MEDDLY_DCASSERT(level>=1);
	MEDDLY_DCASSERT(level<getNumVariables());

	removeAllComputeTableEntries();

	int var_high=getVarByLevel(level+1);
	int var_low=getVarByLevel(level);
	int size_high=getVariableSize(var_high);
	int size_low=getVariableSize(var_low);

	// Renumber the level of nodes for VarHigh
	int num_high=unique->getNumEntries(var_high);
	node_handle* high_nodes=static_cast<node_handle*>(malloc(num_high*sizeof(node_handle)));
	unique->getItems(var_high, high_nodes, num_high);
	for(int i=0; i<num_high; i++) {
		setNodeLevel(high_nodes[i], level);
	}

	// Renumber the level of nodes for VarHigh'
	int num_phigh=unique->getNumEntries(-var_high);
	node_handle* phigh_nodes=static_cast<node_handle*>(malloc(num_phigh*sizeof(node_handle)));
	unique->getItems(-var_high, phigh_nodes, num_phigh);
	for(int i=0; i<num_phigh; i++) {
		setNodeLevel(phigh_nodes[i], -level);
		getNode(phigh_nodes[i]).mark();
	}

	// Renumber the level of nodes for VarLow
	int num_low=unique->getNumEntries(var_low);
	node_handle* low_nodes=static_cast<node_handle*>(malloc(num_low*sizeof(node_handle)));
	unique->getItems(var_low, low_nodes, num_low);
	for(int i=0; i<num_low; i++) {
		setNodeLevel(low_nodes[i], level+1);
	}
	free(low_nodes);

	// Renumber the level of nodes for VarLow'
	int num_plow=unique->getNumEntries(-var_low);
	node_handle* plow_nodes=static_cast<node_handle*>(malloc(num_plow*sizeof(node_handle)));
	unique->getItems(-var_low, plow_nodes, num_plow);
	for(int i=0; i<num_plow; i++) {
		setNodeLevel(plow_nodes[i], -(level+1));
	}
	free(plow_nodes);

//	printf("Before: Level %d : %d, Level %d : %d, Level %d : %d, Level %d : %d\n",
//			level+1, num_high, -(level+1), num_phigh,
//			level, num_low, -level, num_plow);

	// Update the variable order
	order_var[var_high] = level;
	order_var[-var_high] = -level;
	order_var[var_low] = level+1;
	order_var[-var_low] = -(level+1);
	order_level[level+1] = var_low;
	order_level[-(level+1)] = -var_low;
	order_level[level] = var_high;
	order_level[-level] = -var_high;

	std::vector<node_handle> t;
	std::map<node_handle, node_handle> m;
	std::map<node_handle, int> p;
	// Reconstruct nodes for VarHigh
	for(int i=0; i<num_high; i++) {
		node_handle node=swapAdjacentVariablesInMxD(high_nodes[i]);
		if(high_nodes[i]==node){
			// VarLow is DONT_CHANGE in the MxD
			unlinkNode(node);
		}
		else if(getInCount(node)>1){
			assert(getNodeLevel(node)==-(level+1));

			// Duplication conflict
			node_reader* nr = initNodeReader(high_nodes[i], true);
			for(int j=0; j<size_high; j++){
				if(getNodeLevel(nr->d(j))==-level){
					if(p.find(nr->d(j))==p.end()){
						p[nr->d(j)]=1;
					}
					else{
						p[nr->d(j)]++;
					}
				}
			}
			node_reader::recycle(nr);

			m.insert(std::pair<node_handle, node_handle>(node, high_nodes[i]));
//			printf("UPDATE: %d -> %d\n", node, high_nodes[i]);
		}
		else{
			// Newly created node
			swapNodes(high_nodes[i], node);
			unlinkNode(node);
			if(getNodeLevel(high_nodes[i])==(level+1)){
				t.push_back(high_nodes[i]);
			}
		}
	}

	// Reconstruct nodes for VarHigh'
	for(int i=0; i<num_phigh; i++) {
		if(!isActiveNode(phigh_nodes[i]) || !getNode(phigh_nodes[i]).isMarked()){
			continue;
		}
		getNode(phigh_nodes[i]).unmark();
		if(p.find(phigh_nodes[i])!=p.end() && p[phigh_nodes[i]]==getInCount(phigh_nodes[i])){
			continue;
		}

		// VarLow is DONT_CHANGE in the MxD
		node_reader* nr = initNodeReader(phigh_nodes[i], true);
		bool skip=true;
		for(int j=0; j<size_high; j++){
			if(!isLevelAbove(-level, getNodeLevel(nr->d(j)))){
				skip=false;
				break;
			}
		}
		node_reader::recycle(nr);

		if(!skip){
			node_handle node=swapAdjacentVariablesInMxD(phigh_nodes[i]);
			assert(phigh_nodes[i]!=node);

			if(getInCount(node)>1){
				assert(getNodeLevel(node)==-(level+1));

				// Duplication conflict
				m.insert(std::pair<node_handle, node_handle>(node, phigh_nodes[i]));
//				printf("UPDATE: %d -> %d\n", node, phigh_nodes[i]);
			}
			else{
				// Newly created node
				swapNodes(phigh_nodes[i], node);
				unlinkNode(node);
				if(getNodeLevel(phigh_nodes[i])==(level+1)){
					t.push_back(phigh_nodes[i]);
				}
			}
		}
	}

	if(!m.empty()){
		for(std::vector<node_handle>::iterator itr=t.begin(); itr!=t.end(); itr++){
			assert(getNodeLevel(*itr)==(level+1));
			node_reader* nr = initNodeReader(*itr, true);
			bool update=false;
			for(int i=0; i<size_high; i++){
				if(m.find(nr->d(i))!=m.end()){
					update=true;
					break;
				}
			}
			if(update){
				node_builder& nb = useNodeBuilder(level+1, size_high);
				for(int i=0; i<size_high; i++){
					if(m.find(nr->d(i))==m.end()){
						nb.d(i)=linkNode(nr->d(i));
					}
					else{
						nb.d(i)=linkNode(m[nr->d(i)]);
					}
				}
				node_handle node = createReducedNode(-1, nb);
				assert(getInCount(node)==1 && getNodeLevel(node)==level+1);
				swapNodes(*itr, node);
				unlinkNode(node);
			}
			node_reader::recycle(nr);
		}

		for(std::map<node_handle, node_handle>::iterator itr=m.begin(); itr!=m.end(); itr++){
			assert(getInCount(itr->first)==1);
			swapNodes(itr->first, itr->second);
			unlinkNode(itr->first);
		}
	}

	free(high_nodes);
	free(phigh_nodes);

//	printf("After: Level %d : %d,  Level %d : %d, Level %d : %d, Level %d : %d\n",
//			level+1, unique->getNumEntries(var_low),
//			-(level+1), unique->getNumEntries(-var_low),
//			level, unique->getNumEntries(var_high),
//			-level, unique->getNumEntries(-var_high));
//	printf("#Node: %d\n", getCurrentNumNodes());
}

MEDDLY::node_handle MEDDLY::mtmxd_forest::swapAdjacentVariablesInMxD(node_handle node)
{
	int level=ABS(getNodeLevel(node));
	int hi_var=getVarByLevel(level);
	int lo_var=getVarByLevel(level+1);
	int hi_size=getVariableSize(hi_var);
	int lo_size=getVariableSize(lo_var);

	node_builder& hi_nb = useNodeBuilder(level+1, lo_size);

	if(isFullyReduced() || isQuasiReduced()){
		for(int m=0; m<lo_size; m++) {
			node_builder& pr_hi_nb = useNodeBuilder(-(level+1), lo_size);
			for(int n=0; n<lo_size; n++) {
				node_builder& lo_nb = useNodeBuilder(level, hi_size);
				for(int p=0; p<hi_size; p++){
					node_builder& pr_lo_nb = useNodeBuilder(-level, hi_size);
					for(int q=0; q<hi_size; q++){
						node_handle node_p=(getNodeLevel(node)==level ? getDownPtr(node, p) : node);
						node_handle node_pq=(getNodeLevel(node_p)==-(level) ? getDownPtr(node_p, q) : node_p);
						node_handle node_pqm=(getNodeLevel(node_pq)==(level+1) ? getDownPtr(node_pq, m) : node_pq);
						pr_lo_nb.d(q)=linkNode(getNodeLevel(node_pqm)==-(level+1) ? getDownPtr(node_pqm, n) : node_pqm);
					}
					lo_nb.d(p)=createReducedNode(p, pr_lo_nb);
				}
				pr_hi_nb.d(n)=createReducedNode(-1, lo_nb);
			}
			hi_nb.d(m)=createReducedNode(m, pr_hi_nb);
		}
	}
	else if(isIdentityReduced()){
		for(int m=0; m<lo_size; m++) {
			node_builder& pr_hi_nb = useNodeBuilder(-(level+1), lo_size);
			for(int n=0; n<lo_size; n++) {
				node_builder& lo_nb = useNodeBuilder(level, hi_size);
				for(int p=0; p<hi_size; p++){
					node_builder& pr_lo_nb = useNodeBuilder(-level, hi_size);
					for(int q=0; q<hi_size; q++){
						node_handle node_p=(getNodeLevel(node)==level ? getDownPtr(node, p) : node);
						if(getNodeLevel(node_p)!=-(level) && q!=p){
							pr_lo_nb.d(q)=linkNode(getTransparentNode());
						}
						else{
							node_handle node_pq=(getNodeLevel(node_p)==-(level) ? getDownPtr(node_p, q) : node_p);
							node_handle node_pqm=(getNodeLevel(node_pq)==(level+1) ? getDownPtr(node_pq, m) : node_pq);
							if(getNodeLevel(node_pqm)!=-(level+1) && n!=m){
								pr_lo_nb.d(q)=linkNode(getTransparentNode());
							}
							else{
								pr_lo_nb.d(q)=linkNode(getNodeLevel(node_pqm)==-(level+1) ? getDownPtr(node_pqm, n) : node_pqm);
							}
						}
					}
					lo_nb.d(p)=createReducedNode(p, pr_lo_nb);
				}
				pr_hi_nb.d(n)=createReducedNode(-1, lo_nb);
			}
			hi_nb.d(m)=createReducedNode(m, pr_hi_nb);
		}
	}
	else{
		throw error(error::NOT_IMPLEMENTED);
	}

	return createReducedNode(-1, hi_nb);
}

void MEDDLY::mtmxd_forest::swapAdjacentVariablesByLevelSwap(int level)
{
	// Swap VarHigh and VarLow
	MEDDLY_DCASSERT(level>=1);
	MEDDLY_DCASSERT(level<getNumVariables());

	if(!isFullyReduced() && !isQuasiReduced()){
		throw error(error::INVALID_OPERATION);
	}

	removeAllComputeTableEntries();

	// x > x' > y > y'
	swapAdjacentLevels(level);
	// x > y > x' > y'
	swapAdjacentLevels(-(level+1));
	// y > x > x' > y'
	swapAdjacentLevels(-level);
	// y > x > y' > x'
	swapAdjacentLevels(level);
	// y > y' > x > x'
}

void MEDDLY::mtmxd_forest::swapAdjacentLevels(int level)
{
	MEDDLY_DCASSERT(ABS(level)>=1);
	MEDDLY_DCASSERT(ABS(level)<=getNumVariables());

	int high_level=(level<0 ? -level : (-level-1));
	int high_var=getVarByLevel(high_level);
	int low_var=getVarByLevel(level);
	int high_size=getVariableSize(ABS(high_var));
	int low_size=getVariableSize(ABS(low_var));

	int high_node_size=unique->getNumEntries(high_var);
	node_handle* high_nodes=static_cast<node_handle*>(malloc(high_node_size*sizeof(node_handle)));
	unique->getItems(high_var, high_nodes, high_node_size);

	int low_node_size=unique->getNumEntries(low_var);
	node_handle* low_nodes=static_cast<node_handle*>(malloc(low_node_size*sizeof(node_handle)));
	unique->getItems(low_var, low_nodes, low_node_size);

//	printf("Before: Level %d : %d, Level %d : %d\n",
//			level+1, high_node_size,
//			level, low_node_size);

	int i=0, j=0;
	// Renumber the level of nodes for VarHigh
	for(i=0; i<high_node_size; i++) {
		node_reader* nr = initNodeReader(high_nodes[i], true);
		MEDDLY_DCASSERT(nr->getLevel()==high_level);
		MEDDLY_DCASSERT(nr->getSize()==high_size);

		for(int k=0; k<high_size; k++){
			if(getNodeLevel(nr->d(k))==level){
				// Remove the nodes corresponding to functions that
				// are independent of VarLow
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
		setNodeLevel(low_nodes[i], high_level);
	}

	// Update the variable order
	order_var[high_var] = level;
	order_var[low_var] = high_level;
	order_level[high_level] = low_var;
	order_level[level] = high_var;

	// Process the rest of nodes for VarHigh
	for(i=0; i<high_node_size; i++) {
		node_reader* high_nr = initNodeReader(high_nodes[i], true);
		node_builder& high_nb = useNodeBuilder(high_level, low_size);

		for(j=0; j<low_size; j++) {
			node_builder& low_nb = useNodeBuilder(level, high_size);
			for(int k=0; k<high_size; k++) {
				node_handle node_k=high_nr->d(k);
				node_handle node_kj=(getNodeLevel(node_k)==high_level ? getDownPtr(node_k, j) : node_k);
				low_nb.d(k)=linkNode(node_kj);
			}
			high_nb.d(j)=createReducedNode(-1, low_nb);
		}

		node_reader::recycle(high_nr);

		node_handle node=createReducedNode(-1, high_nb);
		assert(getInCount(node)==1);
		assert(getNodeLevel(node)==high_level);

		swapNodes(high_nodes[i], node);
		unlinkNode(node);
	}

	free(high_nodes);
	free(low_nodes);

//	printf("After: Level %d : %d, Level %d : %d\n",
//			level+1, unique->getNumEntries(low_var),
//			level, unique->getNumEntries(high_var));
//	printf("#Node: %d\n", getCurrentNumNodes());
}


void MEDDLY::mtmxd_forest::moveDownVariable(int high, int low)
{
	throw error(error::NOT_IMPLEMENTED);
}

void MEDDLY::mtmxd_forest::moveUpVariable(int low, int high)
{
	throw error(error::NOT_IMPLEMENTED);
}

MEDDLY::mtmxd_forest::mtmxd_iterator::~mtmxd_iterator()
{
}

bool MEDDLY::mtmxd_forest::mtmxd_iterator::start(const dd_edge &e)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmxd_forest::mtmxd_iterator::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  int k = -1;
  node_handle down = 0;
  for (;;) { 
    nzp[k]++;
    if (nzp[k] < path[k].getNNZs()) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      break;
    }
    if (k<0) {
      k = -k;
    } else {
      if (maxLevel == k) {
        level_change = k;
        return false;
      }
      k = -k-1;
    }
  }
  level_change = k;

  return first( (k>0) ? -k : -k-1, down);
}

bool MEDDLY::mtmxd_forest::mtmxd_iterator::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==down) return false;

  bool isFully = F->isFullyReduced();

  for ( ; k; k = downLevel(k) ) {
    MEDDLY_DCASSERT(down);
    int kdn = F->getNodeLevel(down);
    MEDDLY_DCASSERT(!isLevelAbove(kdn, k));

    if (isLevelAbove(k, kdn)) {
      if (k>0 || isFully) {
        F->initRedundantReader(path[k], k, down, false);
      } else {
        F->initIdentityReader(path[k], k, index[-k], down, false);
      }
    } else {
      F->initNodeReader(path[k], down, false);
    }
    nzp[k] = 0;
    index[k] = path[k].i(0);
    down = path[k].d(0);
  }
  // save the terminal value
  index[0] = down;
  return true;
}

// ******************************************************************
// *                                                                *
// *           mtmxd_forest::mtmxd_fixedrow_iter  methods           *
// *                                                                *
// ******************************************************************

MEDDLY::mtmxd_forest::
mtmxd_fixedrow_iter::mtmxd_fixedrow_iter(const expert_forest *F)
 : mt_iterator(F)
{
}

MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter::~mtmxd_fixedrow_iter()
{
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter
::start(const dd_edge &e, const int* minterm)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  for (int k=1; k<=maxLevel; k++) {
    index[k] = minterm[k];
  }
  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  node_handle down = 0;
  // Only try to advance the column, because the row is fixed.
  for (int k=-1; k>=-maxLevel; k--) { 
    for (nzp[k]++; nzp[k] < path[k].getNNZs(); nzp[k]++) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      level_change = k;
      if (first(downLevel(k), down)) return true;
    }
  } // for

  return false;
}


bool MEDDLY::mtmxd_forest::mtmxd_fixedrow_iter::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==k) {
    index[0] = down;
    return true;
  }

  // Check that this "row" node has a non-zero pointer
  // for the fixed index.
  MEDDLY_DCASSERT(k>0);
  int cdown;
  if (isLevelAbove(k, F->getNodeLevel(down))) {
    // skipped unprimed level, must be "fully" reduced
    cdown = down;
  } else {
    cdown = F->getDownPtr(down, index[k]);
  }
  if (0==cdown) return false;

  //
  // Ok, set up the "column" node below
  k = downLevel(k);
  MEDDLY_DCASSERT(k<0);

  if (isLevelAbove(k, F->getNodeLevel(cdown))) {
    // Skipped level, we can be fast about this.
    // first, recurse.
    if (!first(downLevel(k), cdown)) return false;
    // Ok, there is a valid path.
    // Set up this level.
    nzp[k] = 0;
    if (F->isFullyReduced()) {
      F->initRedundantReader(path[k], k, cdown, false);
      index[k] = 0;
    } else {
      index[k] = index[upLevel(k)];
      F->initIdentityReader(path[k], k, index[k], cdown, false);
    }
    return true;
  } 

  // Proper node here.
  // cycle through it and recurse... 

  F->initNodeReader(path[k], cdown, false);

  for (int z=0; z<path[k].getNNZs(); z++) {
    if (first(downLevel(k), path[k].d(z))) {
      nzp[k] = z;
      index[k] = path[k].i(z);
      return true;
    }
  }

  return false;
}

// ******************************************************************
// *                                                                *
// *           mtmxd_forest::mtmxd_fixedcol_iter  methods           *
// *                                                                *
// ******************************************************************

MEDDLY::mtmxd_forest::
mtmxd_fixedcol_iter::mtmxd_fixedcol_iter(const expert_forest *F)
 : mt_iterator(F)
{
}

MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter::~mtmxd_fixedcol_iter()
{
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter
::start(const dd_edge &e, const int* minterm)
{
  if (F != e.getForest()) {
    throw error(error::FOREST_MISMATCH);
  }
  
  for (int k=1; k<=maxLevel; k++) {
    index[-k] = minterm[k];
  }

  return first(maxLevel, e.getNode());
}

bool MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter::next()
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  node_handle down = 0;
  // Only try to advance the row, because the column is fixed.
  for (int k=1; k<=maxLevel; k++) { 
    for (nzp[k]++; nzp[k] < path[k].getNNZs(); nzp[k]++) {
      index[k] = path[k].i(nzp[k]);
      down = path[k].d(nzp[k]);
      MEDDLY_DCASSERT(down);
      level_change = k;
      if (first(downLevel(k), down)) return true;
    }
  } // for

  return false;
}


bool MEDDLY::mtmxd_forest::mtmxd_fixedcol_iter::first(int k, node_handle down)
{
  MEDDLY_DCASSERT(F);
  MEDDLY_DCASSERT(F->isForRelations());
  MEDDLY_DCASSERT(index);
  MEDDLY_DCASSERT(nzp);
  MEDDLY_DCASSERT(path);

  if (0==k) {
    index[0] = down;
    return true;
  }

  if (k<0) {
    // See if this "column node" has a path
    // at the specified index.
    if (isLevelAbove(k, F->getNodeLevel(down))) {
      if (!F->isFullyReduced()) {
        // Identity node here - check index
        if (index[k] != index[upLevel(k)]) return false;
      }
      return first(downLevel(k), down);
    }
    int cdown = F->getDownPtr(down, index[k]);
    if (0==cdown) return false;
    return first(downLevel(k), cdown);
  }

  // Row node.  Find an index, if any,
  // such that there is a valid path below.
  MEDDLY_DCASSERT(k>0);
  int kdn = F->getNodeLevel(down);
  if (isLevelAbove(k, kdn)) {
    // Skipped level, handle quickly
    int kpr = downLevel(k);
    if (isLevelAbove(kpr, F->getNodeLevel(kdn))) {
      // next level is also skipped.
      // See if there is a valid path below.
      if (!first(downLevel(kpr), down)) return false;
      // There's one below, set up the one at these levels.
      F->initRedundantReader(path[k], k, down, false);
      if (F->isFullyReduced()) {
        nzp[k] = 0;
        index[k] = 0;
      } else {
        nzp[k] = index[kpr];
        index[k] = index[kpr];
      }
      return true;
    }
    // next level is not skipped.
    // See if there is a valid path below.
    int cdown = F->getDownPtr(down, index[kpr]);
    if (0==cdown) return false;
    if (!first(kpr, cdown)) return false;
    F->initRedundantReader(path[k], k, down, false);
    nzp[k] = 0;
    index[k] = 0;
    return true;
  }

  // Level is not skipped.
  F->initNodeReader(path[k], down, false);
  
  for (int z=0; z<path[k].getNNZs(); z++) {
    index[k] = path[k].i(z);
    if (first(downLevel(k), path[k].d(z))) {
      nzp[k] = z;
      return true;
    }
  }
  return false;
}

