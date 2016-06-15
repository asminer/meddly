// $Id$

#ifndef LOWEST_INVERSION_REORDERING_H
#define LOWEST_INVERSION_REORDERING_H

namespace MEDDLY{

class lowest_inversion_reordering : public reordering_base
{
public:
  virtual void reorderVariables(expert_forest* forest, const int* order)
  {
    int size = forest->getDomain()->getNumVariables();

    // The variables below ordered_level are ordered
    int ordered_level = 1;
    int level = 1;
    while (ordered_level < size) {
      level = ordered_level;
      while(level>0 && (order[forest->getVarByLevel(level)] > order[forest->getVarByLevel(level+1)])) {
        forest->swapAdjacentVariables(level);
        level--;
      }
      ordered_level++;
    }
  }
};

}

#endif
