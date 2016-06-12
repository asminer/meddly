#ifndef LOWEST_COST_REORDERING_H
#define LOWEST_COST_REORDERING_H

#include "../heap.h"

namespace MEDDLY{

class lowest_cost_reordering : public reordering_base
{
protected:
  long calculate_swap_cost(expert_forest* forest, int level)
  {
    int lvar = forest->getVarByLevel(level);
    int hvar = forest->getVarByLevel(level+1);
    return (long) get_unique_table(forest)->getNumEntries(hvar) * forest->getVariableSize(hvar);
  }

public:
  virtual void reorderVariables(expert_forest* forest, const int* order)
  {
    int size = forest->getDomain()->getNumVariables();
    IndexedHeap<long, less<long>> heap(size);

    for(int i=1; i<size; i++) {
      if(order[forest->getVarByLevel(i)] > order[forest->getVarByLevel(i+1)]) {
        long weight = calculate_swap_cost(forest, i);
        heap.push(i, weight);
      }
    }

    while(!heap.empty()) {
      int level = heap.top_key();
      forest->swapAdjacentVariables(level);
      heap.pop();

      if (level < size - 1
          && order[forest->getVarByLevel(level + 1)] > order[forest->getVarByLevel(level + 2)]) {
        long weight = calculate_swap_cost(forest, level + 1);
        heap.push(level+1, weight);
      }
      if (level > 1
          && order[forest->getVarByLevel(level - 1)] > order[forest->getVarByLevel(level)]) {
        long weight = calculate_swap_cost(forest, level - 1);
        heap.push(level-1, weight);
      }
    }
  }
};

}

#endif
