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

#ifndef MEDDLY_SIMPLE_H
#define MEDDLY_SIMPLE_H

#include "../node_storage.h"

namespace MEDDLY {
  class simple_separated_style;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                  simple_separated_style class                  *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Simple storage mechanism.
    Memory management is completely separated out.
*/

class MEDDLY::simple_separated_style : public node_storage_style {
  public:
    simple_separated_style(const char* n);
    virtual ~simple_separated_style();
    virtual node_storage* createForForest(forest* f,
        const memory_manager_style* mst, memstats &ms) const;
};


#endif

