
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


#ifndef ARRAY_GRID_H
#define ARRAY_GRID_H

namespace MEDDLY {
  class array_grid_style;
};

/**
    Factory for a memory manager based on an array plus grid.

    All details, including the actual memory manager constructed,
    are hidden in implementation file :^)

*/

class MEDDLY::array_grid_style : public memory_manager_style {
  public:
    array_grid_style(const char* n);
    virtual ~array_grid_style();

    virtual memory_manager* initManager(unsigned char granularity,
      unsigned char minsize, memstats &stats) const;
};

#endif

