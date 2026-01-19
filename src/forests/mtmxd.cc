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
#include <unordered_map>

#include "mtmxd.h"
#include "../unique_table.h"

MEDDLY::mtmxd_forest
::mtmxd_forest(domain* d, range_type t, const policies &p)
 : mt_forest(d, true, t, p)
{
}

// ******************************************************************
// *                                                                *
// *              mtmxd_forest::mtmxd_iterator methods              *
// *                                                                *
// ******************************************************************

void MEDDLY::mtmxd_forest::swapAdjacentVariables(int level)
{
  if (getPolicies().isVarSwap()) {
    swapAdjacentVariablesByVarSwap(level);
  }
  else if (getPolicies().isLevelSwap()) {
    swapAdjacentVariablesByLevelSwap(level);
  }
}

void MEDDLY::mtmxd_forest::swapAdjacentVariablesByVarSwap(int level)
{
  // Swap VarHigh and VarLow
  MEDDLY_DCASSERT(level >= 1);
  MEDDLY_DCASSERT(level < getNumVariables());

  int hvar = getVarByLevel(level+1);
  int lvar = getVarByLevel(level);
  int hsize = getVariableSize(hvar);
  // int lsize = getVariableSize(lvar);

  // Renumber the level of nodes for VarHigh
  int hnum = unique->getNumEntries(hvar);
  node_handle* hnodes = new node_handle[hnum];
  unique->getItems(hvar, hnodes, hnum);
  for (int i = 0; i < hnum; i++) {
    setNodeLevel(hnodes[i], level);
  }

  // Renumber the level of nodes for VarHigh'
  int phnum = unique->getNumEntries(-hvar);
  node_handle* phnodes = new node_handle[phnum];
  unique->getItems(-hvar, phnodes, phnum);
  for (int i = 0; i < phnum; i++) {
    setNodeLevel(phnodes[i], -level);
    // Protect nodes from being destroyed
    linkNode(phnodes[i]);
  }

  // Renumber the level of nodes for VarLow
  int lnum = unique->getNumEntries(lvar);
  node_handle* lnodes = new node_handle[lnum];
  unique->getItems(lvar, lnodes, lnum);
  for (int i = 0; i < lnum; i++) {
    setNodeLevel(lnodes[i], level+1);
  }
  delete[] lnodes;

  // Renumber the level of nodes for VarLow'
  int plnum = unique->getNumEntries(-lvar);
  node_handle* plnodes = new node_handle[plnum];
  unique->getItems(-lvar, plnodes, plnum);
  for (int i = 0; i < plnum; i++) {
    setNodeLevel(plnodes[i], -(level+1));
  }
  delete[] plnodes;

  //	printf("Before: Level %d : %d, Level %d : %d, Level %d : %d, Level %d : %d\n",
  //			level+1, num_high, -(level+1), num_phigh,
  //			level, num_low, -level, num_plow);

  // Update the variable order
  std::const_pointer_cast<variable_order>(var_order)->exchange(hvar, lvar);

  std::vector<node_handle> t;
  std::unordered_map<node_handle, node_handle> dup;

  // Reconstruct nodes for VarHigh
  for (int i = 0; i < hnum; i++) {
    node_handle node = swapAdjacentVariablesOf(hnodes[i]);
    if (hnodes[i] == node) {
      // VarLow is DONT_CHANGE in the MxD
      unlinkNode(node);
    }
    else if (getNodeInCount(node) > 1) {
      MEDDLY_DCASSERT(getNodeLevel(node) == -(level+1));

      // Duplication conflict
      dup.emplace(node, hnodes[i]);
//      std::cout << "UPDATE: " << node << " -> " << hnodes[i] << std::endl;
    }
    else {
      // Newly created node
      swapNodes(hnodes[i], node);
      unlinkNode(node);
      if (getNodeLevel(hnodes[i]) == (level+1)) {
        t.push_back(hnodes[i]);
      }
    }
  }
  delete[] hnodes;

  if (!dup.empty()) {
    for (const auto& n : t) {
      MEDDLY_DCASSERT(getNodeLevel(n) == (level+1));

      // unpacked_node* nr = newUnpacked(n, FULL_ONLY);
      unpacked_node* nr = unpacked_node::newFromNode(this, n, FULL_ONLY);
      bool update = false;
      for (int i = 0; i < hsize; i++) {
        if (dup.find(nr->down(i)) != dup.end()) {
          update = true;
          break;
        }
      }

      if (update) {
        unpacked_node* nb =
            unpacked_node::newWritable(this, level + 1, hsize, FULL_ONLY);
        for (int i = 0; i < hsize; i++) {
          auto search = dup.find(nr->down(i));
          nb->setFull(i, linkNode(search == dup.end() ? nr->down(i) : search->second));
        }
        node_handle node = createReducedNode(-1, nb);
        MEDDLY_DCASSERT(getNodeInCount(node) == 1 && getNodeLevel(node) == level + 1);
        swapNodes(n, node);
        unlinkNode(node);
      }

      unpacked_node::Recycle(nr);
    }

    for (const auto& it : dup) {
      MEDDLY_DCASSERT(getNodeInCount(it.first) == 1);
      swapNodes(it.first, it.second);
      unlinkNode(it.first);
    }

    dup.clear();
  }

  {
    int j = 0;
    for (int i = 0; i < phnum; i++) {
      // Revoke the protection
      unlinkNode(phnodes[i]);
      if (isActiveNode(phnodes[i])) {
        phnodes[j++] = phnodes[i];
      }
    }
    phnum = j;
  }

  // Reconstruct nodes for VarHigh'
  for (int i = 0; i < phnum; i++) {
    MEDDLY_DCASSERT(isActiveNode(phnodes[i]));

    // unpacked_node* nr = newUnpacked(phnodes[i], FULL_ONLY);
    unpacked_node* nr = unpacked_node::newFromNode(this, phnodes[i], FULL_ONLY);
    bool skip = true;
    for (int j = 0; j < hsize; j++){
      if (!isLevelAbove(-level, getNodeLevel(nr->down(j)))) {
        skip = false;
        break;
      }
    }
    unpacked_node::Recycle(nr);

    if (skip) {
      // VarLow is DONT_CARE + DONT_CHANGE in the MxD
      continue;
    }

    node_handle node = swapAdjacentVariablesOf(phnodes[i]);
    MEDDLY_DCASSERT(phnodes[i] != node);

    if (getNodeInCount(node) > 1) {
      MEDDLY_DCASSERT(getNodeLevel(node) == -(level+1));

      // Duplication conflict
      dup.emplace(node, phnodes[i]);
//      std::cout << "UPDATE: " << node << " -> " << phnodes[i] << std::endl;
    }
    else {
      // Newly created node
      swapNodes(phnodes[i], node);
      unlinkNode(node);
      if (getNodeLevel(phnodes[i]) == (level+1)) {
        t.push_back(phnodes[i]);
      }
    }
  }
  delete[] phnodes;

  // XXX: Duplicate code
  if (!dup.empty()) {
    for (const auto& n : t) {
      MEDDLY_DCASSERT(getNodeLevel(n)==(level+1));

      // unpacked_node* nr = newUnpacked(n, FULL_ONLY);
      unpacked_node* nr = unpacked_node::newFromNode(this, n, FULL_ONLY);

      bool update = false;
      for (int i = 0; i < hsize; i++) {
        if(dup.find(nr->down(i)) != dup.end()){
          update=true;
          break;
        }
      }

      if (update) {
        unpacked_node* nb =
            unpacked_node::newWritable(this, level + 1, hsize, FULL_ONLY);
        for (int i = 0; i < hsize; i++) {
          auto search = dup.find(nr->down(i));
          nb->setFull(i, linkNode(search == dup.end() ? nr->down(i) : search->second));
        }
        node_handle node = createReducedNode(-1, nb);
        MEDDLY_DCASSERT(getNodeInCount(node) == 1 && getNodeLevel(node) == level + 1);
        swapNodes(n, node);
        unlinkNode(node);
      }

      unpacked_node::Recycle(nr);
    }

    for (const auto& it : dup) {
      MEDDLY_DCASSERT(getNodeInCount(it.first) == 1);
      swapNodes(it.first, it.second);
      unlinkNode(it.first);
    }

    dup.clear();
  }

  //	printf("After: Level %d : %d,  Level %d : %//
  // Complete adjacent variable swap by swapping two levels 4 times
  // Work for fully-fully reduction only
  //d, Level %d : %d, Level %d : %d\n",
  //			level+1, unique->getNumEntries(var_low),
  //			-(level+1), unique->getNumEntries(-var_low),
  //			level, unique->getNumEntries(var_high),
  //			-level, unique->getNumEntries(-var_high));
  //	printf("#Node: %d\n", getCurrentNumNodes());
}

MEDDLY::node_handle MEDDLY::mtmxd_forest::swapAdjacentVariablesOf(node_handle node)
{
  int level = ABS(getNodeLevel(node));
  int hvar = getVarByLevel(level);
  int lvar = getVarByLevel(level+1);
  int hsize = getVariableSize(hvar);
  int lsize = getVariableSize(lvar);

  // Unprimed high node builder
  unpacked_node* hnb =
      unpacked_node::newWritable(this, level + 1, lsize, FULL_ONLY);
  if (isFullyReduced() || isQuasiReduced()) {
    for (int m = 0; m < lsize; m++) {
      // Primed high node builder
      unpacked_node* phnb =
          unpacked_node::newWritable(this, -(level + 1), lsize, FULL_ONLY);
      for (int n = 0; n < lsize; n++) {
        // Unprimed low node builder
        unpacked_node* lnb =
            unpacked_node::newWritable(this, level, hsize, FULL_ONLY);
        for (int p = 0; p < hsize; p++) {
          // Primed low node builder
          unpacked_node* plnb =
              unpacked_node::newWritable(this, -level, hsize, FULL_ONLY);
          for (int q = 0; q < hsize; q++){
            node_handle node_p = (getNodeLevel(node) == level ? getDownPtr(node, p) : node);
            node_handle node_pq = (getNodeLevel(node_p) == -(level) ? getDownPtr(node_p, q) : node_p);
            node_handle node_pqm = (getNodeLevel(node_pq) == (level+1) ? getDownPtr(node_pq, m) : node_pq);
            plnb->setFull(q, linkNode(getNodeLevel(node_pqm) == -(level+1) ? getDownPtr(node_pqm, n) : node_pqm));
          }
          lnb->setFull(p, createReducedNode(p, plnb));
        }
        phnb->setFull(n, createReducedNode(-1, lnb));
      }
      hnb->setFull(m, createReducedNode(m, phnb));
    }
  }
  else if (isIdentityReduced()) {
    for (int m = 0; m < lsize; m++) {
      // Primed high node builder
      unpacked_node* phnb =
          unpacked_node::newWritable(this, -(level + 1), lsize, FULL_ONLY);
      for (int n = 0; n < lsize; n++) {
        // Unprimed low node builder
        unpacked_node* lnb =
            unpacked_node::newWritable(this, level, hsize, FULL_ONLY);
        for (int p = 0; p < hsize; p++) {
          // Primed low node builder
          unpacked_node* plnb =
              unpacked_node::newWritable(this, -level, hsize, FULL_ONLY);
          for (int q = 0; q < hsize; q++) {
            node_handle node_p = (getNodeLevel(node) == level ? getDownPtr(node, p) : node);
            if (getNodeLevel(node_p) != -(level) && q != p) {
              plnb->setFull(q, linkNode(getTransparentNode()));
            }
            else {
              node_handle node_pq = (getNodeLevel(node_p) == -(level) ? getDownPtr(node_p, q) : node_p);
              node_handle node_pqm = (getNodeLevel(node_pq) == (level+1) ? getDownPtr(node_pq, m) : node_pq);
              if (getNodeLevel(node_pqm) != -(level+1) && n != m) {
                plnb->setFull(q, linkNode(getTransparentNode()));
              }
              else {
                plnb->setFull(q, linkNode(getNodeLevel(node_pqm) == -(level+1) ? getDownPtr(node_pqm, n) : node_pqm));
              }
            }
          }
          lnb->setFull(p, createReducedNode(p, plnb));
        }
        phnb->setFull(n, createReducedNode(-1, lnb));
      }
      hnb->setFull(m, createReducedNode(m, phnb));
    }
  }
  else {
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
  }

  return createReducedNode(-1, hnb);
}

//
// Complete adjacent variable swap by swapping two levels four times
// Work for fully-fully and quasi-quasi reductions only
//
void MEDDLY::mtmxd_forest::swapAdjacentVariablesByLevelSwap(int level)
{
  // Swap VarHigh and VarLow
  MEDDLY_DCASSERT(level >= 1);
  MEDDLY_DCASSERT(level < getNumVariables());

  if(!isFullyReduced() && !isQuasiReduced()){
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }

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
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

//  MEDDLY_DCASSERT(ABS(level) >= 1);
//  MEDDLY_DCASSERT(ABS(level) <= getNumVariables());
//
//  int hlevel = (level<0 ? -level : (-level-1));
//  int hvar = getVarByLevel(hlevel);
//  int lvar = getVarByLevel(level);
//  int hsize = getVariableSize(ABS(hvar));
//  int lsize = getVariableSize(ABS(lvar));
//
//  int hnum = unique->getNumEntries(hvar);
//  node_handle* hnodes = new node_handle[hnum];
//  unique->getItems(hvar, hnodes, hnum);
//
//  int lnum = unique->getNumEntries(lvar);
//  node_handle* lnodes = new node_handle[lnum];
//  unique->getItems(lvar, lnodes, lnum);
//
//  //	printf("Before: Level %d : %d, Level %d : %d\n",
//  //			level+1, high_node_size,
//  //			level, low_node_size);
//
//  int num = 0;
//  // Renumber the level of nodes for VarHigh
//  for (int i = 0; i < hnum; i++) {
//    unpacked_node* nr = unpacked_node::useUnpackedNode();
//    nr->initFromNode(this, hnodes[i], true);
//
//    MEDDLY_DCASSERT(nr->getLevel() == hlevel);
//    MEDDLY_DCASSERT(nr->getSize() == hsize);
//
//    for (int j = 0; j < hsize; j++) {
//      if (getNodeLevel(nr->down(j)) == level) {
//        // Remove the nodes corresponding to functions that
//        // are independent of the variable to be moved up
//        hnodes[num++] = hnodes[i];
//        break;
//      }
//    }
//    unpacked_node::Recycle(nr);
//
//    setNodeLevel(hnodes[i], level);
//  }
//  hnum = num;
//
//  // Renumber the level of nodes for the variable to be moved up
//  for(int i = 0; i < lnum; i++) {
//    setNodeLevel(lnodes[i], hlevel);
//  }
//  delete[] lnodes;
//
//  // Update the variable order
//  order_var[hvar] = level;
//  order_var[lvar] = hlevel;
//  order_level[hlevel] = lvar;
//  order_level[level] = hvar;
//
////  std::const_pointer_cast<variable_order>(var_order)->exchange(hvar, lvar);
//
//  // Process the rest of nodes for the variable to be moved down
//  for (int i = 0; i < hnum; i++) {
//    unpacked_node* high_nr = unpacked_node::useUnpackedNode();
//    high_nr->initFromNode(this, hnodes[i], true);
//    unpacked_node* high_nb = unpacked_node::newFull(this, hlevel, lsize);
//    for (int j = 0; j < lsize; j++) {
//      unpacked_node* low_nb = unpacked_node::newFull(this, level, hsize);
//      for (int k = 0; k < hsize; k++) {
//        node_handle node_k = high_nr->down(k);
//        node_handle node_kj = (getNodeLevel(node_k) == hlevel ? getDownPtr(node_k, j) : node_k);
//        low_nb->d_ref(k) = linkNode(node_kj);
//      }
//      high_nb->d_ref(j) = createReducedNode(-1, low_nb);
//    }
//
//    unpacked_node::Recycle(high_nr);
//
//    node_handle node = createReducedNode(-1, high_nb);
//    MEDDLY_DCASSERT(getNodeInCount(node) == 1);
//    MEDDLY_DCASSERT(getNodeLevel(node) == hlevel);
//
//    swapNodes(hnodes[i], node);
//    unlinkNode(node);
//  }
//
//  delete[] hnodes;

  //	printf("After: Level %d : %d, Level %d : %d\n",
  //			level+1, unique->getNumEntries(low_var),
  //			level, unique->getNumEntries(high_var));
  //	printf("#Node: %d\n", getCurrentNumNodes());
}

void MEDDLY::mtmxd_forest::dynamicReorderVariables(int top, int bottom)
{
  MEDDLY_DCASSERT(top > bottom);
  MEDDLY_DCASSERT(top <= getNumVariables());
  MEDDLY_DCASSERT(bottom >= 1);

  removeAllComputeTableEntries();

  std::vector<int> vars;
  vars.reserve(top - bottom + 1);
  for (int i = bottom; i <= top; i++) {
    vars.push_back(getVarByLevel(i));
  }

  for (unsigned i = 0; i < vars.size(); i++) {
    unsigned max = i;
    unsigned max_num = unique->getNumEntries(vars[max]) + unique->getNumEntries(-vars[max]);
    for (unsigned j = i + 1; j < vars.size(); j++){
      unsigned num = unique->getNumEntries(vars[j]) + unique->getNumEntries(-vars[j]);
      if (num > max_num) {
        max = j;
        max_num = num;
      }
    }

    int temp = vars[max];
    vars[max] = vars[i];
    vars[i] = temp;

    sifting(vars[i], top, bottom);
  }
}

void MEDDLY::mtmxd_forest::sifting(int var, int top, int bottom)
{
  int level = getLevelByVar(var);

  MEDDLY_DCASSERT(level <= top && level >= bottom);

#ifdef DEVELOPMENT_CODE
  int num = getCurrentNumNodes();
#endif
  if(level <= (top + bottom) / 2) {
    // Move to the bottom
    while(level > bottom) {
      swapAdjacentVariables(level - 1);
      level--;
    }

    int change = 0;
    int min_level = bottom;

    MEDDLY_DCASSERT(level == bottom);
    // Move to the top
    while(level < top) {
      int high_var = getVarByLevel(level + 1);
      size_t old_sum = unique->getNumEntries(var) + unique->getNumEntries(-var)
          + unique->getNumEntries(high_var) + unique->getNumEntries(-high_var);
      swapAdjacentVariables(level);
      size_t new_sum = unique->getNumEntries(var) + unique->getNumEntries(-var)
          + unique->getNumEntries(high_var) + unique->getNumEntries(-high_var);
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
      swapAdjacentVariables(level);
      level++;
    }

    int change = 0;
    int min_level = top;
    int min = change;

    MEDDLY_DCASSERT(level == top);
    // Move to the bottom
    while(level > bottom) {
      int low_var = getVarByLevel(level - 1);
      size_t old_sum = unique->getNumEntries(var) + unique->getNumEntries(-var)
          + unique->getNumEntries(low_var) + unique->getNumEntries(-low_var);
      swapAdjacentVariables(level - 1);
      size_t new_sum = unique->getNumEntries(var) + unique->getNumEntries(-var)
          + unique->getNumEntries(low_var) + unique->getNumEntries(-low_var);
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

#ifdef DEVELOPMENT_CODE
  MEDDLY_DCASSERT(getCurrentNumNodes() <= num);
#endif
}

void MEDDLY::mtmxd_forest::moveDownVariable(int high, int low)
{
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::mtmxd_forest::moveUpVariable(int low, int high)
{
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

