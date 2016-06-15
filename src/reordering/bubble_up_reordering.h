// $Id$

#ifndef BUBBLE_UP_REORDERING_H
#define BUBBLE_UP_REORDERING_H

namespace MEDDLY{

class bubble_up_reordering : public reordering_base
{
public:
  virtual void reorderVariables(expert_forest* forest, const int * order)
  {
    int size = forest->getDomain()->getNumVariables();

    // Construct the mapping from level to variable
    int* level_to_var = new int[size + 1];
    for(int i = 1; i < size+1; i++) {
      level_to_var[order[i]] = i;
    }

    for(int i = size; i > 1; i--) {
      int level = forest->getLevelByVar(level_to_var[i]);
      while (level < i) {
        forest->swapAdjacentVariables(level);
        level++;
      }
    }

    delete[] level_to_var;
  }
};

}

#endif
