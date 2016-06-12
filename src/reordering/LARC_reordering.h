#ifndef LARC_REORDERING_H
#define LARC_REORDERING_H

#include "../heap.h"

namespace MEDDLY{

//
// Lowest Average Reference Count
//
class LARC_reordering : public reordering_base
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
  virtual void reorderVariables(expert_forest* forest, const int* order)
  {
    int size = forest->getDomain()->getNumVariables();
    IndexedHeap<long, less<double> > heap(size);

    for (int i = 1; i < size; i++) {
      if (order[forest->getVarByLevel(i)] > order[forest->getVarByLevel(i+1)]) {
        double weight = calculate_average_ref_count(forest, i);
        heap.push(i, weight);
      }
    }

    while(!heap.empty()) {
      int level = heap.top_key();
      forest->swapAdjacentVariables(level);
      heap.pop();

      if (level < size-1
          && (order[forest->getVarByLevel(level+1)] > order[forest->getVarByLevel(level+2)])) {
        double weight = calculate_average_ref_count(forest, level + 1);
        heap.push(level + 1, weight);
      }
      if (level > 1
          && (order[forest->getVarByLevel(level-1)] > order[forest->getVarByLevel(level)])) {
        double weight = calculate_average_ref_count(forest, level - 1);
        heap.push(level - 1, weight);
      }
    }
  }
};

}

#endif
