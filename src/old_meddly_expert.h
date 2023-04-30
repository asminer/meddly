
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

  // ******************************************************************
  // *                                                                *
  // *                      Operation management                      *
  // *                                                                *
  // ******************************************************************

  /// Remove an existing operation from the operation cache.
  void removeOperationFromCache(operation* );

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest
        @param  res     Result forest
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    expert_forest* arg, expert_forest* res);

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest from this dd_edge
        @param  res     Result forest from this dd_edge
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    const dd_edge& arg, const dd_edge& res);

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest
        @param  res     Result type
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    expert_forest* arg, opnd_type result);

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest from this dd_edge
        @param  res     Result type
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opname* code,
    const dd_edge& arg, opnd_type result);


  /** Find, or build if necessary, a binary operation.
        @param  code    Operation we want
        @param  arg1    Argument 1 forest
        @param  arg2    Argument 2 forest
        @param  res     Result forest
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  binary_operation* getOperation(const binary_opname* code,
    expert_forest* arg1, expert_forest* arg2, expert_forest* res);

  /** Find, or build if necessary, a binary operation.
        @param  code    Operation we want
        @param  arg1    Argument 1 forest taken from this dd_edge
        @param  arg2    Argument 2 forest taken from this dd_edge
        @param  res     Result forest taken from this dd_edge
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  binary_operation* getOperation(const binary_opname* code,
    const dd_edge& arg1, const dd_edge& arg2, const dd_edge& res);

  /** Safely destroy the given unary operation.
      It should be unnecessary to call this directly.
  */
  void destroyOperation(unary_operation* &op);

  /** Safely destroy the given binary operation.
      It should be unnecessary to call this directly.
  */
  void destroyOperation(binary_operation* &op);

  /// Safely destroy the given numerical operation.
  void destroyOperation(specialized_operation* &op);

  /// Should not be called directly.
  void destroyOpInternal(operation* op);

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
// *                        operation  class                        *
// *                                                                *
// ******************************************************************

// #include "compute_table.h"

/** Generic operation.
    Operations are tied to specific forests.
    Necessary for compute table entries.
*/
class MEDDLY::operation {
    const opname* theOpName;
    unsigned oplist_index;
    bool is_marked_for_deletion;

    // declared and initialized in meddly.cc
    static compute_table* Monolithic_CT;
    // declared and initialized in meddly.cc
    static operation** op_list;
    // declared and initialized in meddly.cc
    static unsigned* op_holes;
    // declared and initialized in meddly.cc
    static unsigned list_size;
    // declared and initialized in meddly.cc
    static unsigned list_alloc;
    // declared and initialized in meddly.cc
    static unsigned free_list;

    // should ONLY be called during library cleanup.
    static void destroyAllOps();

    /**
      Starting slot for entry_types, assigned
      by compute_table.
    */
    unsigned first_etid;

  protected:
    /// Compute table to use (for entry type 0), if any.
    compute_table* CT0;
    /// Array of compute tables, one per entry type.
    compute_table** CT;
    /** Array of entry types.
        Owned by the compute_table class; we have
        these pointers for convenience.
    */
    ct_entry_type** etype;
    /** Array of entry results.
        Use these during computation.
        We only ever need one result per entry type.
    */
    ct_entry_result* CTresult;

    /**
      Number of entry_types needed by this operation.
      This gives the dimension of arrays CT and etype.
    */
    unsigned num_etids;

    /// Struct for CT searches.
    // compute_table::entry_key* CTsrch;
    // for cache of operations.
    operation* next;

    virtual ~operation();
    void markForDeletion();
    void registerInForest(forest* f);
    void unregisterInForest(forest* f);
    void allocEntryForests(int nf);
    void addEntryForest(int index, expert_forest* f);
    void allocEntryObjects(int no);
    void addEntryObject(int index);

    virtual bool checkForestCompatibility() const = 0;

    void registerEntryType(unsigned slot, ct_entry_type* et);
    void buildCTs();

    friend class forest;
    friend void MEDDLY::destroyOpInternal(operation* op);
    friend void MEDDLY::cleanup();

    friend class ct_initializer;

  public:
    /**
        Constructor.
          @param  n         Operation "name"
          @param  et_slots  Number of different compute table entry types
                            used by this operation.
                            Derived class constructors must register
                            exactly this many entry types.
    */
    operation(const opname* n, unsigned et_slots);

    bool isMarkedForDeletion() const;
    void setNext(operation* n);
    operation* getNext();

    static bool usesMonolithicComputeTable();
    static void removeStalesFromMonolithic();
    static void removeAllFromMonolithic();

    /// Remove stale compute table entries for this operation.
    void removeStaleComputeTableEntries();

    /// Remove all compute table entries for this operation.
    void removeAllComputeTableEntries();

    // for compute tables.

    unsigned getIndex() const;
    static operation* getOpWithIndex(unsigned i);
    static unsigned getOpListSize();

    void setFirstETid(unsigned slot);
    unsigned getFirstETid() const;
    unsigned getNumETids() const;

    // for debugging:

    static void showMonolithicComputeTable(output &, int verbLevel);
    static void showAllComputeTables(output &, int verbLevel);
    static void countAllNodeEntries(const expert_forest* f, size_t* counts);

    void showComputeTable(output &, int verbLevel) const;
    void countCTEntries(const expert_forest* f, size_t* counts) const;

    // handy
    const char* getName() const;
    const opname* getOpName() const;
};

// ******************************************************************
// *                                                                *
// *                     unary_operation  class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a unary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::unary_operation : public operation {
  protected:
    expert_forest* argF;
    expert_forest* resF;
    opnd_type resultType;

    virtual ~unary_operation();

    virtual bool checkForestCompatibility() const;

  public:
    unary_operation(const unary_opname* code, unsigned et_slots,
      expert_forest* arg, expert_forest* res);

    unary_operation(const unary_opname* code, unsigned et_slots,
      expert_forest* arg, opnd_type res);

    bool matches(const expert_forest* arg, const expert_forest* res)
      const;

    bool matches(const expert_forest* arg, opnd_type res) const;

    // high-level front-ends

    /**
      Checks forest comatability and then calls computeDDEdge().
    */
    void compute(const dd_edge &arg, dd_edge &res);
    void computeTemp(const dd_edge &arg, dd_edge &res);

    virtual void compute(const dd_edge &arg, long &res);
    virtual void compute(const dd_edge &arg, double &res);
    virtual void compute(const dd_edge &arg, ct_object &c);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
};

// ******************************************************************
// *                                                                *
// *                     binary_operation class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a binary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::binary_operation : public operation {
  protected:
    bool can_commute;
    expert_forest* arg1F;
    expert_forest* arg2F;
    expert_forest* resF;
    opnd_type resultType;

    virtual ~binary_operation();
    void operationCommutes();

    // Check if the variables orders of relevant forests are compatible
    virtual bool checkForestCompatibility() const;

  public:
    binary_operation(const binary_opname* code, unsigned et_slots,
      expert_forest* arg1, expert_forest* arg2, expert_forest* res);

    bool matches(const expert_forest* arg1, const expert_forest* arg2,
      const expert_forest* res) const;

    // high-level front-end

    /**
      Checks forest comatability and then calls computeDDEdge().
    */
    void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);
    void computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

    virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res, bool userFlag)
      = 0;

    // low-level front ends.  TBD - REMOVE THESE BECAUSE THEY BREAK MARK & SWEEP

#ifdef KEEP_LL_COMPUTES

    /// Low-level compute on nodes a and b, return result.
    virtual node_handle compute(node_handle a, node_handle b);
    /// Low-level compute at level k on nodes a and b, return result.
    virtual node_handle compute(int k, node_handle a, node_handle b);

    /// Low-level compute on EV edges (av, ap) and (bv, bp), return result.
    virtual void compute(int av, node_handle ap, int bv, node_handle bp,
      int &cv, node_handle &cp);
    virtual void compute(long av, node_handle ap, long bv, node_handle bp,
      long &cv, node_handle &cp);
    virtual void compute(long av, node_handle ap, node_handle bp,
      long &cv, node_handle &cp);

    /// Low-level compute on EV edges (av, ap) and (bv, bp), return result.
    virtual void compute(float av, node_handle ap, float bv, node_handle bp,
      float &cv, node_handle &cp);

#endif

};

// ******************************************************************
// *                                                                *
// *                  specialized_operation  class                  *
// *                                                                *
// ******************************************************************

/** Mechanism to apply specialized operations.
*/
class MEDDLY::specialized_operation : public operation {
  public:
    specialized_operation(const specialized_opname* code, unsigned et_slots);
  protected:
    virtual ~specialized_operation();
  public:

    /** For unary (like) operations.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &arg, dd_edge &res);

    /** For unary (like) operations with boolean results.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &arg, bool &res);

    /** For binary (like) operations.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

    /** For numerical operations.
        compute y += some function of x, depending on the operation.
        Default behavior is to throw an exception.
    */
    virtual void compute(double* y, const double* x);

    /** For tenary (like) operations.
        Note that there could be other "built in" operands.
        Default behavior is to throw an exception.
    */
    virtual void compute(const dd_edge &ar1, const dd_edge &ar2, const dd_edge &ar3, dd_edge &res);
};

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

