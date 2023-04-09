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

#include "defines.h"
#include "varorder.h"

// ******************************************************************
// *                                                                *
// *                     variable_order methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::variable_order::variable_order(const int* order, int size) {
  MEDDLY_DCASSERT(order[0] == 0);

  level2var.assign(size + 1, 0);
  var2level.assign(size + 1, 0);
  for (int i = 1; i < size + 1; i++) {
    level2var[i] = order[i];
    var2level[order[i]] = i;
  }
}

MEDDLY::variable_order::variable_order(const variable_order& order) {
  MEDDLY_DCASSERT(order.getVarByLevel(0) == 0);

  level2var.assign(order.level2var.begin(), order.level2var.end());
  var2level.assign(order.var2level.begin(), order.var2level.end());
}

// Exchange two variables
// The two variables don't have to be adjacent
void MEDDLY::variable_order::exchange(int var1, int var2) {
  MEDDLY_DCASSERT(var1 > 0 && var2 > 0);

  level2var[var2level[var1]] = var2;
  level2var[var2level[var2]] = var1;

  int temp = var2level[var1];
  var2level[var1] = var2level[var2];
  var2level[var2] = temp;
}

bool MEDDLY::variable_order::is_compatible_with(const int* order) const {
  MEDDLY_DCASSERT(order[0] == 0);
  for (unsigned int i = 1; i < level2var.size(); i++) {
    if (level2var[i] != order[i]) {
      return false;
    }
  }
  return true;
}

bool MEDDLY::variable_order::is_compatible_with(const variable_order& order) const {
  if (this == &order) {
    return true;
  }
  if (level2var.size() != order.level2var.size()) {
    return false;
  }
  for (unsigned int i = 0; i < level2var.size(); i++) {
    if (level2var[i] != order.getVarByLevel(i)) {
      return false;
    }
  }
  return true;
}


