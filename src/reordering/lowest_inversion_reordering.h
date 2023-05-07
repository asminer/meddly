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

#ifndef MEDDLY_LOWEST_INVERSION_REORDERING_H
#define MEDDLY_LOWEST_INVERSION_REORDERING_H

namespace MEDDLY{

class lowest_inversion_reordering : public reordering_base
{
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

    // The variables below ordered_level are ordered
    int ordered_level = 1;
    int level = 1;
    while (ordered_level < size) {
      level = ordered_level;
      while (level > 0 && (var2level[forest->getVarByLevel(level)] > var2level[forest->getVarByLevel(level + 1)])) {
        forest->swapAdjacentVariables(level);
        level--;
      }
      ordered_level++;
    }

    delete[] var2level;
  }
};

}

#endif
