
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
			nb.computeHash();
			modifyReducedNodeInPlace(nb, high_nodes[i]);
			nb.lock=false;
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
		high_nb.computeHash();
		modifyReducedNodeInPlace(high_nb, high_nodes[i]);
		high_nb.lock=false;
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
			nb.computeHash();
			modifyReducedNodeInPlace(nb, low_nodes[i]);
			nb.lock=false;
		}
	}

	free(high_nodes);
	free(low_nodes);
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
