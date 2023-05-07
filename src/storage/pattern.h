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

#ifndef MEDDLY_PATTERN_H
#define MEDDLY_PATTERN_H

#include "../node_storage.h"

namespace MEDDLY {
  class pattern_storage_style;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                  pattern_storage_style class                  *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Pattern storage mechanism.
 Currently incompatible with extensible nodes.
 Memory management is completely separated out.
 */

class MEDDLY::pattern_storage_style : public node_storage_style {
public:
  pattern_storage_style(const char* n);
  virtual ~pattern_storage_style();
  virtual node_storage* createForForest(expert_forest* f,
                                        const memory_manager_style* mst) const;
};


#endif



