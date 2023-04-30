
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


/*! \file meddly_expert.h

    Low-level MDD library interface.

    This interface is for "expert" users who want to define new
    operations, or for library developers to define the built-in operations.
    Casual users probably only need the interface provided by "meddly.h".

    The first part of the interface describes the expert interface and the
    second part contains implementations of virtual functions in the interface.

    IMPORTANT: meddly.h must be included before including this file.
    TODO: Operations are not thread-safe.
*/

#ifndef MEDDLY_EXPERT_H
#define MEDDLY_EXPERT_H

#include <string.h>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <map>
#include <set>

#include "initializer.h"

// #define DEBUG_MARK_SWEEP
// #define DEBUG_BUILDLIST

// #define KEEP_LL_COMPUTES

namespace MEDDLY {

  // classes defined here

  class expert_variable;
  class expert_domain;

  // wrapper for temporary nodes
  class unpacked_node;  // replacement for node_reader, node_builder

  class relation_node;

  /*

    class op_initializer;

    Generalized to class initializer_list.
  */

  class initializer_list;

  /*
    class cleanup_procedure;

    Subsumed by class initializer_list.
  */

  // Memory managers, for node storage and compute tables
  class memory_manager_style;
  class memory_manager;

  // Node header storage
  class node_headers;

  // Actual node storage
  class node_storage_style;
  class node_storage;

  class expert_forest;

  class opname;
  class unary_opname;
  class binary_opname;
  class specialized_opname;
  class numerical_opname;
  class satpregen_opname;
  class satotf_opname;
  class satimpl_opname;
  class sathyb_opname;
  class constrained_opname;

  class ct_initializer;
  class compute_table_style;
  class compute_table;
  class ct_entry_type;
  class ct_entry_result;

  class operation;
  class unary_operation;
  class binary_operation;
  class specialized_operation;


  // classes defined elsewhere
  class base_table;
  class unique_table;
  class impl_unique_table;

  class reordering_base;

  /// Safely destroy the given numerical operation.
  void destroyOperation(specialized_operation* &op);

  // ******************************************************************
  // *                                                                *
  // *                  library management functions                  *
  // *                                                                *
  // ******************************************************************

  /*
    /// Builds an initializer for MEDDLY's builtin operations.
    /// Use defaultInitializerList() instead
    op_initializer* makeBuiltinInitializer();

  */

}; // namespace MEDDLY


#endif

