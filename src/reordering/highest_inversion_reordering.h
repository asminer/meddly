#ifndef HIGH_INVERSION_REORDERING_H
#define HIGH_INVERSION_REORDERING_H

namespace MEDDLY{

class highest_inversion_reordering : public reordering_base
{
public:
  virtual void reorderVariables(expert_forest* forest, const int* order)
  {
    int size = forest->getDomain()->getNumVariables();

    // The variables above ordered_level are ordered
    int ordered_level = size-1;
    int level = size-1;
    while (ordered_level > 0) {
      level = ordered_level;
      while(level<size && (order[forest->getVarByLevel(level)] > order[forest->getVarByLevel(level+1)])) {
        forest->swapAdjacentVariables(level);
        level++;
      }
      ordered_level--;
    }
  }
};

}

#endif
