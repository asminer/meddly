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

#ifndef MEDDLY_RANDOM_REORDERING_H
#define MEDDLY_RANDOM_REORDERING_H

#include <time.h>
#include <vector>

namespace MEDDLY{

class random_reordering : public reordering_base
{
public:
  virtual void reorderVariables(forest* forest, const int* level2var)
  {
    int size = forest->getDomain()->getNumVariables();

    // Transpose order
    int* var2level = new int[size];
    var2level[0] = 0;
    for (int i = 1; i <= size; i++) {
      var2level[level2var[i]] = i;
    }

    std::vector<bool> inversions(size+1, false);
    std::vector<int> levels;
    for (int i = 1; i < size; i++) {
      if (var2level[forest->getVarByLevel(i)] > var2level[forest->getVarByLevel(i + 1)]) {
        inversions[i] = true;
        levels.push_back(i);
      }
    }

    srand(time(nullptr));
    int seed = rand();
    srand(seed);

    while (!levels.empty()) {
      int index = rand() % levels.size();
      int level = levels[index];
      ASSERT(__FILE__, __LINE__, inversions[level]);
      forest->swapAdjacentVariables(level);

      inversions[level] = false;
      if (level > 1) {
        if (!inversions[level-1]
            && var2level[forest->getVarByLevel(level - 1)] > var2level[forest->getVarByLevel(level)]){
          // New inversion at lower level
          inversions[level-1] = true;
          levels.push_back(level - 1);
        }
      }
      if (level < size - 1) {
        if (!inversions[level + 1]
            && var2level[forest->getVarByLevel(level + 1)] > var2level[forest->getVarByLevel(level + 2)]) {
          // New inversion at upper level
          inversions[level+1] = true;
          levels.push_back(level + 1);
        }
      }
      levels[index] = levels.back();
      levels.pop_back();
    }

    delete[] var2level;
  }
};

}

#endif
