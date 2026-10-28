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

#ifndef MEDDLY_REORDERING_BASE_H
#define MEDDLY_REORDERING_BASE_H

namespace MEDDLY{
    class forest;

class reordering_base
{
public:
  virtual ~reordering_base() {}

protected:
  const unique_table* get_unique_table(forest* forest) const;
  int getInCount(forest* forest, node_handle p) const;

public:
  virtual void reorderVariables(forest* forest, const int* level2var) = 0;
};

inline const unique_table* reordering_base::get_unique_table(forest* forest) const
{
  return forest->getUT();
}

inline int reordering_base::getInCount(forest* forest, node_handle p) const
{
  return forest->getNodeInCount(p);
}

}

#endif
