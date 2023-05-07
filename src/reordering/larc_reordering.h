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

#ifndef MEDDLY_LARC_REORDERING_H
#define MEDDLY_LARC_REORDERING_H

#include "../heap.h"
#include "../unique_table.h"

namespace MEDDLY{

//
// Lowest Average Reference Count
//
class larc_reordering : public reordering_base
{
protected:
  double calculate_average_ref_count(expert_forest* forest, int level)
  {
    int lvar = forest->getVarByLevel(level);
    int lnum = get_unique_table(forest)->getNumEntries(lvar);
    if (lnum == 0) {
      return 0;
    }
    else {
      node_handle* lnodes = new node_handle[lnum];
      get_unique_table(forest)->getItems(lvar, lnodes, lnum);

      int edges = 0;
      for (int i = 0; i < lnum; i++) {
        edges += getInCount(forest, lnodes[i]);
      }

      delete[] lnodes;
      return (double) edges / lnum;
    }
  }

public:
  virtual void reorderVariables(expert_forest* forest, const int* level2var)
  {
    int size = forest->getDomain()->getNumVariables();

    // Transpose order
    int* var2level = new int[size];
    var2level[0] = 0;
    for (int i = 1; i <= size; i++) {
      var2level[level2var[i]] = i;
    }

    IndexedHeap<long, less<double>> heap(size);
    for (int i = 1; i < size; i++) {
      if (var2level[forest->getVarByLevel(i)] > var2level[forest->getVarByLevel(i + 1)]) {
        double weight = calculate_average_ref_count(forest, i);
        heap.push(i, weight);
      }
    }

    while (!heap.empty()) {
      int level = heap.top_key();
      forest->swapAdjacentVariables(level);
      heap.pop();

      if (level < size-1
          && var2level[forest->getVarByLevel(level + 1)] > var2level[forest->getVarByLevel(level + 2)]) {
        double weight = calculate_average_ref_count(forest, level + 1);
        heap.push(level + 1, weight);
      }
      if (level > 1
          && var2level[forest->getVarByLevel(level - 1)] > var2level[forest->getVarByLevel(level)]) {
        double weight = calculate_average_ref_count(forest, level - 1);
        heap.push(level - 1, weight);
      }
    }

    delete[] var2level;
  }
};

}

#endif
