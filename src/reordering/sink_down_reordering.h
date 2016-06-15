// $Id$

#ifndef SINK_DOWN_REORDERING_H
#define SINK_DOWN_REORDERING_H

namespace MEDDLY{

class sink_down_reordering : public reordering_base
{
public:
  virtual void reorderVariables(expert_forest* forest, const int * order)
  {
    int size = forest->getDomain()->getNumVariables();

    // Construct the mapping from level to variable
    int* level_to_var = new int[size + 1];
    for (int i = 1; i < size+1; i++) {
      level_to_var[order[i]] = i;
    }

    for (int i = 1; i < size; i++) {
      int level = forest->getLevelByVar(level_to_var[i]);
      while (level > i) {
        forest->swapAdjacentVariables(level-1);
        level--;
      }
    }

    delete[] level_to_var;
  }
};

}

#endif
