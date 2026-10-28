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

#ifndef MEDDLY_MTMXD_H
#define MEDDLY_MTMXD_H

#include "mt.h"

namespace MEDDLY {
  class mtmxd_forest;
};

class MEDDLY::mtmxd_forest : public mt_forest {
  public:
    mtmxd_forest(domain* d, range_type t, const policies &p);

    virtual void swapAdjacentVariables(int level);
    virtual void moveDownVariable(int high, int low);
    virtual void moveUpVariable(int low, int high);

    virtual void dynamicReorderVariables(int top, int bottom);
    void sifting(int var, int top, int bottom);

    void swapAdjacentVariablesByVarSwap(int level);
    /** Return the root node after swapping the adjacent variables
        in the MxD with the given root node.
    */
    node_handle swapAdjacentVariablesOf(node_handle node);

    void swapAdjacentVariablesByLevelSwap(int level);
    void swapAdjacentLevels(int level);
};

#endif

