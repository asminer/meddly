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

#ifndef MEDDLY_SINK_DOWN_REORDERING_H
#define MEDDLY_SINK_DOWN_REORDERING_H

namespace MEDDLY{

class sink_down_reordering : public reordering_base
{
public:
  virtual void reorderVariables(expert_forest* forest, const int * level2var)
  {
    int size = forest->getDomain()->getNumVariables();

    for (int i = 1; i < size; i++) {
      int level = forest->getLevelByVar(level2var[i]);
      while (level > i) {
        forest->swapAdjacentVariables(level - 1);
        level--;
      }
    }
  }
};

}

#endif
