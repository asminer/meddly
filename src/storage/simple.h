
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



// TODO: Testing

#ifndef SIMPLE_H
#define SIMPLE_H

#include "../defines.h"
#include "../hash_stream.h"
#include "holeman.h"

namespace MEDDLY {
  class simple_separated_style;

// for now...

  class simple_grid_style;
  class simple_array_style;
  class simple_heap_style;
  class simple_none_style;
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
    virtual node_storage* createForForest(expert_forest* f,
        const memory_manager_style* mst) const;
};



//
//
//
// OLD styles below here
//
//
//


// ******************************************************************
// *                                                                *
// *                                                                *
// *                    simple_grid_style  class                    *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Simple storage using the original grid mechanism for holes.
    Should be equivalent to old_node_storage; use to check for overheads.
*/

class MEDDLY::simple_grid_style : public node_storage_style {
  public:
    simple_grid_style(const char* n);
    virtual ~simple_grid_style();
    virtual node_storage* createForForest(expert_forest* f,
        const memory_manager_style* mst) const;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                    simple_array_style class                    *
// *                                                                *
// *                                                                *
// ******************************************************************

/// Simple storage using a new array of lists mechanism for holes.
class MEDDLY::simple_array_style : public node_storage_style {
  public:
    simple_array_style(const char* n);
    virtual ~simple_array_style();
    virtual node_storage* createForForest(expert_forest* f,
        const memory_manager_style* mst) const;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                    simple_heap_style  class                    *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Simple storage using a clever grid of heaps mechanism for holes.
    This allows us to use "earliest holes first".
*/

class MEDDLY::simple_heap_style : public node_storage_style {
  public:
    simple_heap_style(const char* n);
    virtual ~simple_heap_style();
    virtual node_storage* createForForest(expert_forest* f,
        const memory_manager_style* mst) const;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                    simple_none_style  class                    *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Simple storage using no hole management whatsoever.
*/
class MEDDLY::simple_none_style : public node_storage_style {
  public:
    simple_none_style(const char* n);
    virtual ~simple_none_style();
    virtual node_storage* createForForest(expert_forest* f,
        const memory_manager_style* mst) const;
};


#endif

