
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

  class global_rebuilder;

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


// ******************************************************************
// *                                                                *
// *                    global_rebuilder  class                     *
// *                                                                *
// ******************************************************************

/** Rebuild the dd_edge from the source forest in the target forest.
    The source and target forests may have different variable orders.
    While rebuilding, extra nodes may be created in the source forest
    because of the restrict operation.
*/

class MEDDLY::global_rebuilder {
private:
  struct RestrictKey {
    node_handle p;
    int var;
    int idx;

    bool operator==(const RestrictKey &other) const {
      return (p == other.p && var == other.var && idx == other.idx);
    }
  };

  struct RestrictKeyHasher {
    size_t operator()(const RestrictKey &key) const;
  };

  struct TransformKey {
    int sig;
//    int var;

    bool operator==(const TransformKey &other) const {
      return sig == other.sig;
//      return (sig == other.sig && var == other.var);
    }
  };

  struct TransformEntry {
    // Partial assignment
    std::vector<int> pa;
    node_handle p;

    bool operator==(const TransformEntry &other) const {
      return p == other.p;
    }
  };

  struct TransformKeyHasher {
    size_t operator()(const TransformKey &key) const;
  };

  class SignatureGenerator {
  protected:
    global_rebuilder &_gr;

  public:
    SignatureGenerator(global_rebuilder& gr);
    virtual ~SignatureGenerator();
    virtual void precompute() = 0;
    virtual int signature(node_handle p) = 0;
  };

  class TopDownSignatureGenerator: public SignatureGenerator {
  public:
    TopDownSignatureGenerator(global_rebuilder& gr);
    void precompute() override;
    int signature(node_handle p) override;
  };

  class BottomUpSignatureGenerator: public SignatureGenerator {
  private:
    std::unordered_map<node_handle, int> _cache_sig;
    std::unordered_map<node_handle, int> _cache_rec_sig;

    int rec_signature(node_handle p);

  public:
    BottomUpSignatureGenerator(global_rebuilder& gr);
    void precompute() override;
    int signature(node_handle p) override;
  };

  std::unordered_map<RestrictKey, node_handle, RestrictKeyHasher> _computed_restrict;
  std::unordered_multimap<TransformKey, TransformEntry, TransformKeyHasher> _computed_transform;
  SignatureGenerator* _sg;

  expert_forest* _source;
  expert_forest* _target;
  node_handle _root;
  int _hit;
  int _total;

  node_handle transform(node_handle p, int target_level, std::vector<int>& pa);
  node_handle restrict(node_handle p, std::vector<int>& pa);

  bool restrict_exist(node_handle p, const std::vector<int>& pa,
      unsigned start, node_handle& result);
  int signature(node_handle p) const;

  // Return the top variable in the sub-order of the target variable order
  // starting from 0 to target_level
  // such that the given decision diagram depends on it.
  int check_dependency(node_handle p, int target_level) const;

public:
  friend class SignatureGenerator;

  global_rebuilder(expert_forest* source, expert_forest* target);
  ~global_rebuilder();

  dd_edge rebuild(const dd_edge& e);
  void clearCache();
  double hitRate() const;
};

#endif

