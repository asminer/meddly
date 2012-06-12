
// $Id$

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

#include <map>
#include <vector>
#include <string.h>	// For memcpy


// Flags for development version only. Significant reduction in performance.
#ifdef DEVELOPMENT_CODE
#define RANGE_CHECK_ON
#define DCASSERTS_ON
//#define DEBUG_PRINTS_ON
#ifdef DEBUG_PRINTS_ON
#define DEBUG_MDD_H
#define TRACK_DELETIONS
#define TRACK_CACHECOUNT
#endif
#endif

// Use this for assertions that will fail only when your
// code is wrong.  Handy for debugging.
#ifdef DCASSERTS_ON
#define MEDDLY_DCASSERT(X) assert(X)
#else
#define MEDDLY_DCASSERT(X)
#endif

// Use this for range checking assertions that should succeed.
#ifdef RANGE_CHECK_ON
#define MEDDLY_CHECK_RANGE(MIN, VALUE, MAX) { assert(VALUE < MAX); assert(VALUE >= MIN); }
#else
#define MEDDLY_CHECK_RANGE(MIN, VALUE, MAX)
#endif
 
namespace MEDDLY {

  // Functions for reinterpreting an int to a float and vice-versa
  float   toFloat (int a);
  int     toInt   (float a);
  float*  toFloat (int* a);

  // classes defined here
  struct settings;
  class expert_variable;
  class expert_domain;
  class expert_forest;
  class temp_dd_edge;

  class opname;
  class unary_opname;
  class binary_opname;
  class numerical_opname;

  class compute_table;

  class operation;
  class unary_operation;
  class binary_operation;
  class numerical_operation;

  class op_initializer;

  // classes defined elsewhere
  class base_table;

  // ******************************************************************
  // *                   Named numerical operations                   *
  // ******************************************************************

  /** Computes y = y + xA.
      x and y are vectors, stored explicitly, and A is a matrix.
      x_ind and y_ind specify how minterms are mapped to indexes
      for vectors x and y, respectively.
  */
  extern const numerical_opname* VECT_MATR_MULT;

  /** Computes y = y + Ax.
      x and y are vectors, stored explicitly, and A is a matrix.
      x_ind and y_ind specify how minterms are mapped to indexes
      for vectors x and y, respectively.
  */
  extern const numerical_opname* MATR_VECT_MULT;

  // ******************************************************************
  // *                      Operation management                      *
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
  void destroyOperation(numerical_operation* &op);

  /// Should not be called directly.
  void destroyOpInternal(operation* op);
  
  /// Builds an initializer for MEDDLY's builtin operations.
  op_initializer* makeBuiltinInitializer();

}; // namespace MEDDLY


// ******************************************************************
// *                         ct_object class                        *
// ******************************************************************

/** Generic objects in compute tables.
    Used for things other than dd_edges and simple types.
    Defined in ops.cc
*/
class MEDDLY::ct_object {
  public:
    ct_object();
    virtual ~ct_object();
    virtual opnd_type getType() = 0;
};

// ******************************************************************
// *                         settings  class                        *
// ******************************************************************

/** "Global" settings for MEDDLY.
    These settings DO NOT CHANGE once the library is initialized.
  
    The compute cache by default uses a hash table with chaining (i.e. a
    list of entries at each hash location). This can be changed to a
    hash-table without chaining. The maximum size of the hash table can
    also be fixed (to limit the amount of memory it uses).
*/
struct MEDDLY::settings {
  public:
    struct computeTableSettings {
      public:
        enum styleOption {
          /// One huge hash table that uses chaining.
          MonolithicChainedHash,
          /// One huge hash table that does not use chaining.
          MonolithicUnchainedHash,
          /// A hash table (with chaining) for each operation.
          OperationChainedHash,
          /// A hash table (no chaining) for each operation.
          OperationUnchainedHash,
          /// A STL "map" for each operation.
          OperationMap
        };
        enum staleRemovalOption {
          /// Whenever we see a stale entry, remove it.
          Aggressive,
          /// Only remove stales when we need to expand the table
          Moderate,
          /// Only remove stales during Garbage Collection.
          Lazy
        };
      public:
        /// The type of compute table(s) that should be used.
        styleOption style;
        /// Maximum compute table size.
        unsigned maxSize;
        /// How aggressively should we try to eliminate stale entries.
        staleRemovalOption staleRemoval;
      public:
        /// Constructor, to set defaults.
        computeTableSettings() {
          style = MonolithicUnchainedHash;
          staleRemoval = Moderate;
          maxSize = 16777216;
        }
    };
  public:
    /// Settings for the compute table(s)
    computeTableSettings computeTable;
    /// Initializer for operations
    op_initializer* operationBuilder;
  public:
    /// Constructor, to set defaults.
    settings() : computeTable() {
      operationBuilder = makeBuiltinInitializer();
    }
    /// super handly
    inline bool usesMonolithicComputeTable() {
      return (
       computeTableSettings::MonolithicChainedHash == computeTable.style ||
       computeTableSettings::MonolithicUnchainedHash == computeTable.style 
      );
    }
};
  
// ******************************************************************
// *                     expert_variable  class                     *
// ******************************************************************

class MEDDLY::expert_variable : public variable {
  public:
    expert_variable(int b, char* n);
  private:
    virtual ~expert_variable();
  public:
    /// Update our list of domains: add \a d.
    void addToList(domain* d);
    /// Update our list of domains: remove \a d.
    void removeFromList(const domain* d);

    /** Enlarge the possible values for a variable.
      This could modify all nodes in all forests, depending on the
      choice of reduction rule.
      @param  prime   If prime is true, enlarge the bound for
                      the primed variable only, otherwise both
                      the primed and unprimed are enlarged.
      @param  b       New bound, if less than the current bound
                      an error is thrown.
    */
    void enlargeBound(bool prime, int b);

    /** Shrink the possible values for a variable.
      This could modify all nodes in all forests, depending on the
      choice of reduction rule.
      @param  b       New bound, if more than the current bound
                      an error is thrown.
      @param  force   If \a b is too small, and information will be lost,
                      proceed anyway if \a force is true, otherwise
                      return an error code.
    */
    void shrinkBound(int b, bool force);
  private:
    domain** domlist;
    int dl_alloc;
    int dl_used;
};


// ******************************************************************
// *                      expert_domain  class                      *
// ******************************************************************

class MEDDLY::expert_domain : public domain {
  public:
    expert_domain(variable**, int);

  protected:
    ~expert_domain();

  public:
    /** Create all variables at once, from the top down.
      Requires the domain to be "empty" (containing no variables or forests).
      @param  bounds  Current variable bounds.
                      bounds[0] gives the bound for the top-most variable,
                      and bounds[N-1] gives the bound for the bottom-most
                      variable.
      @param  N       Number of variables.
    */
    void createVariablesTopDown(const int* bounds, int N);

    /** Insert a new variable.
          @param  lev   Level to insert above; use 0 for a 
                        new bottom-most level.
          @param  v     Variable to insert.
    */
    void insertVariableAboveLevel(int lev, variable* v);

    /** Remove a variable at the specified level.
        An error is thrown if the variable size is not 1.
        Use shrinkVariableBound() to make the bound 1.
        All forests are modified as appropriate.
          @param  lev   Level number.
    */
    void removeVariableAtLevel(int lev);

    /** Find the level of a given variable.
          @param  v   Variable to search for.
          @return 0, if the variable was not found;
                  i, with getVar(i) == v, otherwise.
    */
    int findLevelOfVariable(const variable* v) const;

    inline expert_variable* getExpertVar(int lev) const {
      return (expert_variable*) vars[lev];
    }
    inline const expert_variable* readExpertVar(int lev) const {
      return (expert_variable*) vars[lev];
    }

    

    /** Add a new variable with bound 1.
      Can be used when the domain already has forests, in which case
      all forests are modified as appropriate.
      @param  below   Placement information: the new variable will appear
                      immediately above the level \a below.
    */
    void createVariable(int below);

    /** Add a new variable with bound 1.
      Deprecated as of version 0.5; use one paramater version instead.
    */
#ifdef _MSC_VER
    __declspec(deprecated)
#endif
#ifdef __GNUC__
    __attribute__ ((deprecated))
#endif
    inline void createVariable(int below, int &vh) {
      createVariable(below);
      vh = below+1;
    }

    /** Destroy a variable with bound 1.
        Deprecated as of version 0.5; use removeVariableAtLevel instead.
    */
#ifdef _MSC_VER
    __declspec(deprecated)
#endif
#ifdef __GNUC__
    __attribute__ ((deprecated))
#endif
    inline void destroyVariable(int vh) { removeVariableAtLevel(vh); }

    /** Get the position of the variable in this domain's variable-order.
      \a TERMINALS are considered to be at height 0.
      be at height 0.
      Deprecated as of version 0.5: level height always equals level handle.
      @param  vh      Any variable handle.
      @return         The variable at this height. 0 for \a TERMINALS.
                      If \a vh is not a valid level handle, return -1.
    */
#ifdef _MSC_VER
    __declspec(deprecated)
#endif
#ifdef __GNUC__
    __attribute__ ((deprecated))
#endif
    inline int getVariableHeight(int vh) const { return vh; }

    /** Get the variable with height \a ht in this domain's variable-order.
      \a TERMINALS are considered to be at height 0.
      Deprecated as of version 0.5: level height always equals level handle.
      @param  ht      Height of the variable.
      @return         The variable with this height. If the height is not in
                      [0, height of top variable], returns -1.
    */
#ifdef _MSC_VER
    __declspec(deprecated)
#endif
#ifdef __GNUC__
    __attribute__ ((deprecated))
#endif
    inline int getVariableWithHeight(int ht) const { return ht; }

    /** Swap the locations of variables in forests.
      I.e., changes the variable ordering of all forests with this domain.
      @param  lev1    Level of first variable.
      @param  lev2    Level of second variable.
    */
    void swapOrderOfVariables(int lev1, int lev2);

    /** Find the actual bound of a variable.
      @param  vh      Variable handle.
      @return         The smallest shrinkable bound before information loss
                      for variable \a vh. If \a vh is invalid, or TERMINALS,
                      returns 0.
    */
    int findVariableBound(int vh) const;

    /** Enlarge the possible values for a variable.
      This could modify all nodes in all forests, depending on the
      choice of reduction rule.
      @param  lev     Variable handle.
      @param  prime   If prime is true, enlarge the bound for
                      the primed variable only, otherwise both
                      the primed and unprimed are enlarged.
      @param  b       New bound, if less than the current bound
                      an error code is returned.
    */
    inline void enlargeVariableBound(int vh, bool prime, int b) {
      getExpertVar(vh)->enlargeBound(prime, b);
    }

    /** Shrink the possible values for a variable.
      This could modify all nodes in all forests, depending on the
      choice of reduction rule.
      @param  lev     Variable handle.
      @param  b       New bound, if more than the current bound
                      an error code is returned.
      @param  force   If \a b is too small, and information will be lost,
                      proceed anyway if \a force is true, otherwise
                      return an error code.
    */
    void shrinkVariableBound(int vh, int b, bool force) {
      getExpertVar(vh)->shrinkBound(b, force);
    }

    virtual void createVariablesBottomUp(const int* bounds, int N);
    virtual forest* createForest(bool rel, forest::range_type t,
      forest::edge_labeling ev);
    virtual void showInfo(FILE* strm);

    /// Free the slot that the forest is using.
    void unlinkForest(expert_forest* f, int slot);

    // --------------------------------------------------------------------

  private:
    expert_forest** forests;
    int szForests;

    /// Find a free slot for a new forest.
    int findEmptyForestSlot();

    /// Mark this domain for deletion
    void markForDeletion();

    friend void MEDDLY::destroyDomain(domain* &d);
};


// ******************************************************************
// *                      expert_forest  class                      *
// ******************************************************************

class MEDDLY::expert_forest : public forest
{
  public:
    /** Constructor.
      @param  dslot   slot used to store the forest, in the domain
      @param  d       domain to which this forest belongs to.
      @param  rel     does this forest represent a relation.
      @param  t       the range of the functions represented in this forest.
      @param  ev      edge annotation.
      @param  r       reduction rule.
      @param  s       storage rule.
      @param  ndp     node deletion policy.
    */
    expert_forest(int dslot, domain *d, bool rel, range_type t, 
      edge_labeling ev, reduction_rule r, node_storage s, 
      node_deletion_policy ndp);

  protected:
    /// Destructor.
    virtual ~expert_forest();  

    /// Phase one destruction.
    void markForDeletion();

    friend void MEDDLY::destroyForest(MEDDLY::forest* &f);
    friend class expert_domain;

  public:
    /// Returns a non-modifiable pointer to this forest's domain.
    const domain* getDomain() const;

    inline const expert_domain* getExpertDomain() const {
      return (expert_domain*) getDomain();
    }
    inline expert_domain* useExpertDomain() {
      return (expert_domain*) getDomain();
    }

    /// Returns a modifiable pointer to this forest's domain.
    domain* useDomain();

    /// Does this forest represent functions that are relations?
    bool isForRelations() const;

    /// Returns the type of range represented by functions in this forest.
    forest::range_type getRangeType() const;

    /// Returns the edge annotation mechanism used by this forest.
    forest::edge_labeling getEdgeLabeling() const;

    /// Returns the reduction rule used by this forest.
    forest::reduction_rule getReductionRule() const;

    /// Returns the storage mechanism used by this forest.
    forest::node_storage getNodeStorage() const;

    /// Returns the node deletion policy used by this forest.
    forest::node_deletion_policy getNodeDeletion() const;

    /// Sets the reduction rule for this forest.
    void setReductionRule(forest::reduction_rule r);

    /// Sets the storage mechanism for this forest.
    void setNodeStorage(forest::node_storage ns);

    /// Sets the node deletion policy for this forest.
    void setNodeDeletion(forest::node_deletion_policy np);

    /// Remove any stale compute table entries associated with this forest.
    void removeStaleComputeTableEntries();

    /// Remove all compute table entries associated with this forest.
    void removeAllComputeTableEntries();

    /// Create a temporary node -- a node that can be modified by the user.
    /// If \a clear is true, downpointers are initialized to 0.
    virtual int createTempNode(int lh, int size, bool clear = true) = 0;

    /// Create a temporary node with the maximum size allowed for this level.
    /// If \a clear is true, downpointers are initialized to 0.
    int createTempNodeMaxSize(int lh, bool clear = true);

    /// Create a temporary node with the given downpointers. Note that
    /// downPointers[i] corresponds to the downpointer at index i.
    /// IMPORTANT: The incounts for the downpointers are not incremented.
    /// The returned value is the handle for the temporary node.
    virtual int createTempNode(int lh, std::vector<int>& downPointers) = 0;

    /// Same as createTempNode(int, vector<int>) except this is for EV+MDDs.
    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<int>& edgeValues) = 0;

    /// Same as createTempNode(int, vector<int>) except this is for EV*MDDs.
    virtual int createTempNode(int lh, std::vector<int>& downPointers,
        std::vector<float>& edgeValues) = 0;

    /// Increase the size of the temporary node.
    /// The maximum size is dictated by domain to which this forest belongs to.
    virtual void resizeNode(int node, int size) = 0;

    /// Build a copy of the given node.
    /// The new node's size will be equal to max(sizeof(a), size).
    virtual int makeACopy(int node, int size = 0) = 0;

    /// The maximum size (number of indices) a node at this level can have
    int getLevelSize(int lh) const;

    /// Apply reduction rule to the temporary node and finalize it. Once
    /// a node is reduced, its contents cannot be modified.
    virtual int reduceNode(int node) = 0;

    /// Reduce and finalize an node with an incoming edge value
    virtual void normalizeAndReduceNode(int& node, int& ev) = 0;
    virtual void normalizeAndReduceNode(int& node, float& ev) = 0;

    /// A is a temporary node, and B is a reduced node.
    /// Accumulate B into A, i.e. A += B.
    /// A still remains a temporary node.
    /// B is not modified.
    /// throws error::INVALID_OPERATION if A or B are inactive nodes.
    virtual void accumulate(int& A, int B) = 0;

    /// A is a temporary node, and B is a minterm.
    /// Accumulate B into A, i.e. A += B.
    /// A still remains a temporary node.
    /// B is not modified.
    /// Assert violation will occur if A is not an active node,
    /// or if B is 0.
    /// Returns true if a new element was added to MDD, false otherwise.
    virtual bool accumulate(int& A, int* B);

    /// Accumuluate a minterm into a MXD.
    /// A is a temporary node.
    /// vlist and vplist constitute the unprimed and primed levels
    /// in the minterm.
    virtual bool accumulate(int& A, int* vlist, int* vplist);

    /// Has the node been reduced
    /// Terminal nodes are also considered to be reduced nodes.
    virtual bool isReducedNode(int node) const = 0;

    /// Is this a full or truncated-full node?
    /// A full node of size 6
    /// Example: n = [0, 3, 234, -1, 0, 344223]
    /// getFullNodeSize(n) returns 6
    /// getFullNodeDownPtr(n, 2) returns 234
    bool isFullNode(int node) const;

    /// Is this a sparse node?
    /// A sparse node of size 6
    /// Example: n = [(1:3), (2:234), (3:-1), (5:344223)]
    /// getSparseNodeSize(n) returns 4
    /// getSparseNodeIndex(n, 1) returns 2
    /// getSparseNodeDownPtr(n, 1) returns 234
    bool isSparseNode(int node) const;

    /// Is this a terminal node?
    /// Note: terminal nodes have a handle <= 0.
    /// (i) represented node handle i.
    /// In MDDs/EVMDDs: (0) = false, (-1) = true.
    /// In MTMDDs: (-i) = i; i.e. a terminal value of i is represented
    ///            internally using a node handle of -i.
    bool isTerminalNode(int node) const;

    /// Get the integer value represented by this terminal node.
    bool getBoolean(int terminalNode) const;

    /// Get the integer value represented by this terminal node.
    int getInteger(int terminalNode) const;

    /// Get the real (float) value represented by this terminal node.
    float getReal(int terminalNode) const;

    /// Get the terminal node representing this boolean value.
    int getTerminalNode(bool booleanValue) const;

    /// Get the terminal node representing this integer value.
    int getTerminalNode(int integerValue) const;

    /// Get the terminal node representing this real (float) value.
    int getTerminalNode(float realValue) const;

    /// Get the node's level handle
    int getNodeLevel(int node) const;

    /// Get the node's height
    int getNodeHeight(int node) const;

    /// Get the number of entries in a full node
    int getFullNodeSize(int node) const;

    /// Get the number of entries in a sparse nodes
    int getSparseNodeSize(int node) const;

    /// Sets the specified index of this node to point to down.
    /// The old downpointer at this index is unlinked (via unlinkNode()).
    /// This index now points to the new downpointer. The reference count
    /// to the new downpointer is incremented (via linkNode()).
    /// Note 1: Only non-reduced nodes can be modified.
    /// Note 2: Non-reduced nodes are always Full nodes regardless of the
    ///         forest's node storage policy.
    void setDownPtr(int node, int index, int down);

    /// Sets the specified index of this node to point to down.
    /// Same as setDownPtr(), except that the previous downpointer is
    /// discarded without unlinking. This is useful if you are using
    /// nodes with uninitialized downpointers created by calling
    /// createTempNode(level, size, false) (i.e. clear flag is false).
    void setDownPtrWoUnlink(int node, int index, int down);

    /// Sets the specified index's edge value -- only for Edge-valued MDDs
    /// Note 1: Only non-reduced nodes can be modified.
    /// Note 2: Non-reduced nodes are always Full nodes regardless of the
    ///         forest's node storage policy.
    void setEdgeValue(int node, int index, int edge);
    void setEdgeValue(int node, int index, float edge);

    /// Get the node pointed to at the given index -- only for Full nodes
    int getFullNodeDownPtr(int node, int index) const;

    /// Get the node pointed to at the given index -- works for Full and
    /// Sparse nodes. For Full nodes, this is the same as calling
    /// getFullNodeDownPtr(node, index). For Sparse nodes, this searches
    /// the sparse node's indexes for the given index and if found, it
    /// returns the downpointer associated with the index. Note that this
    /// is not the same as calling getSparseNodeDownPtr(int node, int index).
    int getDownPtr(int node, int index) const;

    /// Get the nodes pointed to by this node (works with Full or Sparse nodes.
    /// The vector downPointers is increased in size if necessary (but
    /// never reduced in size).
    /// Returns false is operation is not possible. Possible reasons are:
    /// node does not exist; node is a terminal; or node has not been reduced.
    bool getDownPtrs(int node, std::vector<int>& downPointers) const;

    /// Similar to getDownPtrs() except for EV+MDDs.
    /// If successful (return value true), the vectors hold the
    /// downpointers and edge-values.
    virtual bool getDownPtrsAndEdgeValues(int node,
        std::vector<int>& downPointers, std::vector<int>& edgeValues) const
        = 0;

    /// Similar to getDownPtrs() except for EV*MDDs.
    /// If successful (return value true), the vectors hold the
    /// downpointers and edge-values.
    virtual bool getDownPtrsAndEdgeValues(int node,
        std::vector<int>& downPointers, std::vector<float>& edgeValues) const
        = 0;

    /// Get the edge value for the given index -- only for Full nodes
    void getFullNodeEdgeValue(int node, int index, int& ev) const;
    void getFullNodeEdgeValue(int node, int index, float& ev) const;

    /// Get the real index at the given index -- only for Sparse nodes
    int getSparseNodeIndex(int node, int index) const;

    /// Get the node pointed to at the given index -- only for Sparse nodes
    int getSparseNodeDownPtr(int node, int index) const;

    /// Get the edge value at the given index -- only for Sparse nodes
    void getSparseNodeEdgeValue(int node, int index, int& ev) const;
    void getSparseNodeEdgeValue(int node, int index, float& ev) const;

    bool getDownPtrs(int node, const int*& dptrs) const;
    bool getSparseNodeIndexes(int node, const int*& indexes) const;
    bool getEdgeValues(int node, const float*& edgeValues) const;
    bool getEdgeValues(int node, const int*& edgeValues) const;

    bool getDownPtrs(int node, int*& dptrs);
    bool getEdgeValues(int node, float*& edgeValues);
    bool getEdgeValues(int node, int*& edgeValues);

    /// Increase the link count to this node. Call this when another node is
    /// made to point (down-pointer) to this node.
    ///   @return node, for convenience.
    int linkNode(int node);

    /// Decrease the link count to this node. If link count reduces to 0, this
    /// node may get marked for deletion. Call this when another node releases
    /// its connection to this node.
    void unlinkNode(int node);

    /// Increase the cache count for this node. Call this whenever this node
    /// is added to a cache. 
    ///   @return node, for convenience.
    int cacheNode(int node);

    /// Decrease the cache count for this node. Call this whenever this node
    /// is removed from a cache.
    void uncacheNode(int node);

    /// Returns the in-count for a node. This indicates the number of MDD
    /// nodes that link to this node. A node is never deleted when its
    /// in-count is more than zero.
    /// Note that a reference to the in-count is returned. Therefore, the
    /// in-count of the node can be modified by modifying the reference.
    int& getInCount(int node) const;

    /// Returns the cache-count for a node. This indicates the number of
    /// compute cache entries that link to this node.
    /// Note that a reference to the cache-count is returned. Therefore, the
    /// cache-count of the node can be modified by modifying the reference.
    int& getCacheCount(int node) const;

    /// A node can be discarded once it goes stale. Whether a node is
    /// considered stale depends on the forest's deletion policy.
    /// Optimistic deletion: A node is said to be stale only when both the
    ///   in-count and cache-count are zero.
    /// Pessimistic deletion: A node is said to be stale when the in-count
    ///  is zero regardless of the cache-count.
    bool isStale(int node) const;

    /// Returns the node count for this node. The node count is the number
    /// of unique nodes in the decision diagram represented by this node.
    unsigned getNodeCount(int node) const;

    /// Returns the edge count for this node. The edge count is the number
    /// of unique edges in the decision diagram represented by this node.
    unsigned getEdgeCount(int node, bool countZeroes) const;

    /// Display the contents of node
    virtual void showNode(FILE* s, int node, int verbose = 0) const = 0;
    virtual void showNodeGraph(FILE* s, int node) const = 0;

    /// Is this forest an MDD?
    bool isMdd() const;

    /// Is this forest an Multi-Terminal MDD?
    bool isMtMdd() const;

    /// Is this forest an Matrix Diagram?
    bool isMxd() const;

    /// Is this forest an Multi-terminal Matrix Diagram?
    bool isMtMxd() const;

    /// Is this forest an EV+ MDD? In EV+, edge-valued are summed
    /// as we move down a path in the MDD.
    bool isEvplusMdd() const;

    /// Is this forest an EV* MDD? In EV*, edge-valued are multiplied
    /// as we move down a path in the MDD.
    bool isEvtimesMdd() const;

    /// Does this node represent an Index Set?
    /// Note: Only applicable to EV+MDDs.
    bool isIndexSet(int node) const;

    /// Returns an Index Set's cardinality.
    /// Note: Only applicable to Index Sets (which are
    /// implemented using EV+MDDs).
    /// Returns:
    /// -1     Don't know yet
    /// 0      Invalid node
    /// > 0    Cardinality of the Index Set
    int getIndexSetCardinality(int node) const;

    /// Sets an Index Set's cardinality.
    /// Note: Only applicable to Index Sets (which are
    /// implemented using EV+MDDs).
    /// Throws INVALID_OPERATION if the operation is not valid for this forest.
    void setIndexSetCardinality(int node, int c);

    /// Is this node an active node?
    /// An active node is a node that has not yet been recycled (i.e. dead).
    /// Temporary nodes, Reduced nodes and Terminal nodes are considered
    /// to be active nodes.
    bool isActiveNode(int node) const;

    /// If the node is at a primed level, it returns (2 * node_height - 1).
    /// Otherwise, it returns (2 * node_height).
    int getMappedNodeHeight(int node) const;

    /// Register an operation with this forest.
    void registerOperation(const operation* op);

    /// Unregister an operation.
    void unregisterOperation(const operation* op);

  protected:
    // for debugging:
    void showComputeTable(FILE* s, int verbLevel) const;

    void unregisterDDEdges();

    int getInternalNodeSize(int node) const;
    int* getNodeAddress(int node) const;
    int* getAddress(int k, int offset) const;
    int getNodeOffset(int node) const;

    int mapLevel(int level) const;
    int unmapLevel(int level) const;

    bool isZombieNode(int node) const;

    bool isPessimistic() const;

    // the following virtual functions are implemented in node_manager
    virtual bool isValidNodeIndex(int node) const = 0;
    virtual bool isValidLevel(int level) const = 0;
    virtual void reclaimOrphanNode(int node) = 0;     // for linkNode()
    virtual void handleNewOrphanNode(int node) = 0;   // for unlinkNode()
    virtual void deleteOrphanNode(int node) = 0;      // for uncacheNode()
    virtual void freeZombieNode(int node) = 0;        // for uncacheNode()
    // for isStale()
    virtual bool discardTemporaryNodesFromComputeCache() const = 0;

    /// Address of each node.
    typedef struct {
      /**
        Node level
        If the node is active, this indicates node level.
        */
      int level;
      /** 
        Offset to node's data in corresponding level's data array.
        If the node is active, this is the offset (>0) in the data array.
        If the node is deleted, this is -next deleted node
        (part of the unused address list).
        */
      int offset;
      /**
        Cache count
        The number of cache entries that refer to this node (excl. unique
        table). If this node is a zombie, cache_count is negative.
        */
      int cache_count;
    } mdd_node_data;

    /// address info for nodes
    mdd_node_data *address;

    /// Level data for each level in a MDD
    typedef struct {
      /// level height
      int height;
      /// data array
      int* data;
      /// Size of data array.
      int size;
      /// Last used data slot.  Also total number of ints "allocated"
      int last;
      /// Mark for compaction
      bool compactLevel;
      /// Node representing a variable at this level pointing to terminals
      /// based on index
      int levelNode;

      // Holes grid info

      /// Pointer to top of holes grid
      int holes_top;
      /// Pointer to bottom of holes grid
      int holes_bottom;
      /// Total ints in holes
      int hole_slots;

      // performance stats

      /// Largest traversed height of holes grid
      int max_hole_chain;
      /// Number of zombie nodes
      int zombie_nodes;
      /// Number of temporary nodes -- nodes that have not been reduced
      int temp_nodes;
      /// Total number of compactions
      int num_compactions;
    } mdd_level_data;

    /// Level data. Each level maintains its own data array and hole grid.
    mdd_level_data *level;

    domain *d;
    bool isRelation;
    forest::range_type rangeType;
    forest::edge_labeling edgeLabel;
    forest::reduction_rule reductionRule;
    forest::node_storage nodeStorage;
    forest::node_deletion_policy nodeDeletionPolicy;

  private:
    int d_slot;
    bool is_marked_for_deletion;

    friend class dd_edge;
    void registerEdge(dd_edge& e);
    void unregisterEdge(dd_edge& e);


    // structure to store references to registered dd_edges.
    typedef struct {
      // If nextHole >= 0, this is a hole, in which case nextHole also
      // refers to the array index of the next hole in edge[]
      int nextHole;
      // registered dd_edge
      dd_edge* edge;
    } edge_data;

    // Array of registered dd_edges
    edge_data *edge;

    // Size of edge[]
    unsigned sz;
    // Array index: 1 + (last used slot in edge[])
    unsigned firstFree;
    // Array index: most recently created hole in edge[]
    int firstHole;

    // Used to count registered operations
    int* opCount;
    int szOpCount;

#if 0
  // TODO: 

  public:

    /** Set the specified reduction rule for one variable.
      @param  r   Desired reduction rule.
      @param  vh  Variable handle where the rule will be applied. 
      @return     An appropriate error code.
     */
    virtual error setReductionRuleForVariable(reduction_rule r, int vh) = 0;

    /** Get the current reduction rule for the specified variable.
      @param  vh  Variable handle.
      @return     The reduction rule for the variable.
      TODO: dealing with errors?
     */
    virtual reduction_rule getReductionRuleForVariable(int vh) const = 0;

    /** Set the node storage meachanism for one variable.
      @param  ns  Desired node storage mechanism.
      @param  vh  Variable handle where the meachanism will be applied. 
      @return     An appropriate error code.
     */
    virtual error setNodeStorageForVariable(node_storage ns, int vh) = 0;

    /** Get the current node storage meachanism for one variable.
      @param  vh  Variable handle.  
      @return     The reduction rule for the variable.
      TODO: dealing with errors?
     */
    virtual node_storage getNodeStorageForVariable(int vh) const = 0;

    /** Get the current memory allocation for a given variable.
      I.e., the total amount of memory allocated for nodes
      for the given variable, including memory space currently
      unused (e.g., holes from deleted nodes).
      @param  vh  Variable handle.
      @return     The memory usage in bytes, or -1 on error.
     */
    virtual int getMemoryAllocationForVariable(int vh) const = 0;

    /** Get the current memory usage for a given variable.
      Like getMemoryAllocationForVariable(), except we
      report only the memory actually in use.
      @param  vh  Variable handle.
      @return     The memory usage in bytes, or -1 on error.
     */
    virtual int getMemoryUsedForVariable(int vh) const = 0;

    /** Compact the memory for a given variable.
      Can be expensive.
      @param  vh  Variable handle.
      @return     An appropriate error code.
     */
    virtual error compactMemoryForVariable(int vh) = 0;

    /** Set a compaction threshold for a given variable.
      The threshold is set as a percentage * 100.
      Whenever the \b unused memory for this variable exceeds the threshold
      (as a percentage), a compaction operation is performed automatically.
      @param  vh  Variable handle.
      @param  p   Percentage * 100, i.e., use 50 for 50\%.
      @return     An appropriate error code.
     */    
    virtual error setCompactionThresholdForVariable(int vh, int p) = 0;

    /** Get the compaction threshold currently set for a given variable.
      @param  vh  Variable handle.
      @return     The threshold, or -1 on error.
     */    
    virtual int getCompactionThresholdForVariable(int vh) const = 0;

    /** Set the node deletion policy for a given variable.
      Determines how aggressively node memory should be reclaimed when
      nodes become disconnected.
      @param  vh  Variable handle.
      @param  np  Policy to use.
      @return     An appropriate error code.
      TODO: need to specify what happens when the policy is changed.
     */
    virtual error 
      setNodeDeletionForVariable(int vh, node_deletion_policy np) = 0;

    /** Get the node deletion policy for a given variable.
      @param  vh  Variable handle.
      @return     The current policy.
      TODO: errors?
     */
    virtual node_deletion_policy getNodeDeletionForVariable(int vh) const = 0;

#endif

};


/** A bare-bones class for the construction of temporary dd_edges.

    The motivation behind offering this class is to speed up construction
    of DDs, especially in the explicit addition of edges to the DD.

    WARNING: THE USER IS COMPLETELY IN CHARGE OF MEMORY ALLOCATION.

    Intended usage:
    (1) Create an empty temp_dd_edge.
    (2) Assign it a forest.
    (3) Assign it to a level handle (-ve for primed levels).
    (4) If this is a terminal node:
        (a) Set levelHandle to 0.
        (b) If MDD or EVMDD, the terminal is assumed to be TRUE.
        (c) If MTMDD, set terminal value using either iValue or rValue,
            depending on the range of the MTMDD.
    (5) If this is a non-terminal node:
        (a) Set levelHandle != 0.
        (b) Set size > 0, and allocate memory for downpointer[].
        (c) If EV+MDD, also allocate memory for iEdgeValue[].
        (d) If EV*MDD, also allocate memory for rEdgeValue[].
    (6) All downpointers point to either 0 (null) or other temp_dd_edges.
    (7) The downpointers may point to any temp_dd_edge as long as they
        belong to the same forest and maintain the variable ordering of the
        associated domain.
    (8) Note that the downpointers DO NOT have to point to distinct
        temp_dd_edges (i.e. sharing is allowed).
    (9) The user is in charge of creating, keeping track of and
        deleting temp_dd_edges. Meddly is not monitoring any of this.
    
    (10) When the temp_dd_edge is ready, use convertToDDEdge() to convert
         to a dd_edge based on the forest's reduction rules.
         NOTE: The temp_dd_edge is left unchanged and the user is
         responsible for deleting the temp_dd_edges after this conversion.
*/
class MEDDLY::temp_dd_edge {
  public:
    temp_dd_edge();
    // Deallocates downpointers[], iEdgeValues[] and rEdgeValues[].
    // It DOES NOT DEALLOCATE any temp_dd_edges referred to by the
    // downpointers.
    ~temp_dd_edge();

    // Adds an element to the temporary dd_edge.
    // void add(const int* vlist);
    void add(const int* vlist, const int* vplist);

    // Converts the temporary dd_edge to a reduced dd_edge.
    // Returns false if conversion fails (please refer to the intended
    // usage rules listed above).
    // 
    // NOTE: The temp_dd_edge is left unchanged and the user is
    // responsible for deleting the temp_dd_edges after this conversion.
    //
    bool convertToDDEdge(dd_edge& e) const;

    int             levelHandle;
    expert_forest*  forestHandle;
    long            iValue;
    double          rValue;
    int             size;
    temp_dd_edge**  downpointers;
    long*           iEdgeValues;
    double*         rEdgeValues;

  private:
    bool reduce(int& result) const;
    bool reduce(std::map<temp_dd_edge*, int>& ct, int zero, int& result) const;
};



// ******************************************************************
// *                          opname class                          *
// ******************************************************************

/// Class for names of operations.
class MEDDLY::opname {
    const char* name;
    int index;
    static int next_index;

    friend void MEDDLY::initialize(const settings &);
    friend void MEDDLY::cleanup();
  public:
    opname(const char* n);
    virtual ~opname();

    inline int getIndex() const         { return index; }
    inline const char* getName() const  { return name; }
};

// ******************************************************************
// *                       unary_opname class                       *
// ******************************************************************

/// Unary operation names.
class MEDDLY::unary_opname : public opname {
  public:
    unary_opname(const char* n);
    virtual ~unary_opname();

    virtual unary_operation* 
      buildOperation(expert_forest* arg, expert_forest* res) const;

    virtual unary_operation* 
      buildOperation(expert_forest* arg, opnd_type res) const;
};

// ******************************************************************
// *                      binary_opname  class                      *
// ******************************************************************

/// Binary operation names.
class MEDDLY::binary_opname : public opname {
  public:
    binary_opname(const char* n);
    virtual ~binary_opname();

    virtual binary_operation* buildOperation(expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res) const = 0;
};


// ******************************************************************
// *                     numerical_opname class                     *
// ******************************************************************

/// Numerical operation names.
class MEDDLY::numerical_opname : public opname {
  public:
    numerical_opname(const char* n);
    virtual ~numerical_opname();

    /** Note - unlike the more general binary and unary ops,
        a numerical operation might be crafted to the specific
        arguments, for speed.
        The idea is that these operations will be called several
        times (say, within a linear solver) with the same dd_edges.
    */
    virtual numerical_operation* buildOperation(const dd_edge &x_ind,
      const dd_edge &A, const dd_edge &y_ind) const = 0;
};


// ******************************************************************
// *                      compute_table  class                      *
// ******************************************************************

/** Interface for compute tables.
    Anyone implementing an operation (see below) will
    probably want to use this.
*/
class MEDDLY::compute_table {
    public:
      /// The maximum size of the hash table.
      unsigned maxSize;
      /// Do we try to eliminate stales during a "find" operation
      bool checkStalesOnFind;
      /// Do we try to eliminate stales during a "resize" operation
      bool checkStalesOnResize;

      struct stats {
        unsigned numEntries;
        unsigned hits;
        unsigned pings;
        static const int searchHistogramSize = 256;
        long searchHistogram[searchHistogramSize];
        long numLargeSearches;
        int maxSearchLength;
      };

      class search_key {
          friend class MEDDLY::base_table;
          int hashLength;
          int* data;
          bool killData;
          int* key_data;
          const operation* op;
          /// used only for range checking during "development".
          int keyLength;  
        public:
          search_key();
          ~search_key();
          inline int& key(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<keyLength);
#endif
            return key_data[i]; 
          }
          inline const int* rawData() const { return data; }
          inline int dataLength() const { return hashLength; }
          inline const operation* getOp() const { return op; }
      };

      class temp_entry {
          friend class MEDDLY::base_table;
          int handle;
          int hashLength;
          int* entry;
          int* key_entry;
          int* res_entry;
          // The remaining entries are used only in development code
          int keyLength;
          int resLength;
        public:
          inline int& key(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<keyLength);
#endif
            return key_entry[i]; 
          }
          inline int& result(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<resLength);
#endif
            return res_entry[i]; 
          }
          inline void copyResult(int i, void* data, size_t bytes) {
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i+bytes<=resLength*sizeof(int));
#endif
            memcpy(res_entry+i, data, bytes);
          }
          // The following are used by the compute table.
          inline const int* readEntry(int off) const { return entry+off; }
          inline int readHandle() const { return handle; }
          inline int readLength() const { return hashLength; }
          inline int& data(int i) {
            return entry[i];
          }
      };

    public:
      /// Constructor
      compute_table(const settings::computeTableSettings &s);

      /** Destructor. 
          Does NOT properly discard all table entries;
          use \a removeAll() for this.
      */
      virtual ~compute_table();

      /// Is this a per-operation compute table?
      virtual bool isOperationTable() const = 0;

      /// Initialize a search key for a given operation.
      virtual void initializeSearchKey(search_key &key, operation* op) = 0;

      /** Find an entry in the compute table based on the key provided.
          @param  key   Key to search for.
          @return       0, if not found;
                        otherwise, an integer array of size 
                        op->getCacheEntryLength()
      */
      virtual const int* find(const search_key &key) = 0;

      /** Start a new compute table entry.
          The operation should "fill in" the values for the entry,
          then call \a addEntry().
      */
      virtual temp_entry& startNewEntry(operation* op) = 0;

      /** Add the "current" new entry to the compute table.
          The entry may be specified by filling in the values 
          for the struct returned by \a startNewEntry().
      */
      virtual void addEntry() = 0;

      /** Remove all stale entries.
          Scans the table for entries that are no longer valid (i.e. they are
          stale, according to operation::isEntryStale) and removes them. This
          can be a time-consuming process (proportional to the number of cached
          entries).
      */
      virtual void removeStales() = 0;

      /** Removes all entries.
      */
      virtual void removeAll() = 0;

      /// Get performance stats for the table.
      inline const stats& getStats() {
        return perf;
      }

      /// For debugging.
      virtual void show(FILE *s, int verbLevel = 0) = 0;

    protected:
      stats perf;
      temp_entry currEntry;
};

// ******************************************************************
// *                        operation  class                        *
// ******************************************************************

/** Generic operation.
    Operations are tied to specific forests.
    Necessary for compute table entries.
*/
class MEDDLY::operation {
    const opname* theOpName;
    bool is_marked_for_deletion;

    friend void MEDDLY::initialize(const settings &);
    friend void MEDDLY::cleanup();

    // declared and initialized in meddly.cc
    static compute_table* Monolithic_CT;
    // declared and initialized in meddly.cc
    static operation** op_list;
    // declared and initialized in meddly.cc
    static int* op_holes;
    // declared and initialized in meddly.cc
    static int list_size;
    // declared and initialized in meddly.cc
    static int list_alloc;
    // declared and initialized in meddly.cc
    static int free_list;

    int oplist_index;

    int key_length; 
    int ans_length; 

  protected:
    /// Compute table to use, if any.
    compute_table* CT;
    /// Struct for CT searches.
    compute_table::search_key CTsrch;
    // for cache of operations.
    operation* next;
    // must stale compute table hits be discarded.
    // if the result forest is using pessimistic deletion, then true.
    // otherwise, false.  MUST BE SET BY DERIVED CLASSES.
    bool discardStaleHits;
  public:
    /// New constructor.
    /// @param  n   Operation "name"
    /// @param  kl  Key length of compute table entries.
    ///             Use 0 if this operation does not use the compute table.
    /// @param  al  Answer length of compute table entries.
    ///             Use 0 if this operation does not use the compute table.
    operation(const opname* n, int kl, int al);
    //operation(const opname* n, bool uses_CT);

  protected:
    virtual ~operation();

    inline void setAnswerForest(const expert_forest* f) {
      discardStaleHits = f 
        ?   f->getNodeDeletion() == forest::PESSIMISTIC_DELETION
        :   false;  // shouldn't be possible, so we'll do what's fastest.
    }

    void markForDeletion();

    friend class expert_forest;
    friend void MEDDLY::destroyOpInternal(operation* op);

  public:
    inline bool isMarkedForDeletion() const { return is_marked_for_deletion; }

    // should ONLY be called during library cleanup.
    static void destroyOpList();

    inline void setNext(operation* n) { next = n; }
    inline operation* getNext()       { return next; }

    inline static bool useMonolithicComputeTable() { return Monolithic_CT; }
    static void removeStalesFromMonolithic();

    /// Remove stale compute table entries for this operation.
    void removeStaleComputeTableEntries();

    /// Remove all compute table entries for this operation.
    void removeAllComputeTableEntries();

    // for compute tables.

    inline int getIndex() const { return oplist_index; }
    static inline operation* getOpWithIndex(int i) { return op_list[i]; }
    static inline int getOpListSize() { return list_size; }

    // for debugging:

    static void showMonolithicComputeTable(FILE*, int verbLevel);
    static void showAllComputeTables(FILE*, int verbLevel);
    void showComputeTable(FILE*, int verbLevel) const;

    // handy
    inline const char* getName() const { return theOpName->getName(); }
    inline const opname* getOpName() const { return theOpName; }

    /// Number of ints that make up the key (usually the operands).
    inline int getKeyLength() const { 
      return key_length; 
    }

    /// Number of ints that make up the answer (usually the results).
    inline int getAnsLength() const { 
      return ans_length; 
    }

    /// Number of ints that make up the entire record (key + answer)
    inline int getCacheEntryLength() const { 
      return key_length + ans_length; 
    }

    /// Checks if the cache entry (in entryData[]) is stale.
    inline bool isEntryStale(const int* data) {
      return (is_marked_for_deletion || isStaleEntry(data));
    }

  protected:
    virtual bool isStaleEntry(const int* entry) = 0;

  public:
    /// Removes the cache entry (in entryData[]) by informing the
    /// applicable forests that the nodes in this entry are being removed
    /// from the cache
    virtual void discardEntry(const int* entryData) = 0;

    /// Prints a string representation of this cache entry on strm (stream).
    virtual void showEntry(FILE* strm, const int *entryData) const = 0;

    inline bool shouldStaleCacheHitsBeDiscarded() const {
      return discardStaleHits;
    }

  protected:
    void allocEntryForests(int nf);
    void addEntryForest(int index, expert_forest* f);
    void allocEntryObjects(int no);
    void addEntryObject(int index);
};

// ******************************************************************
// *                     unary_operation  class                     *
// ******************************************************************

/** Mechanism to apply a unary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::unary_operation : public operation {
  protected:
    expert_forest* argF;
    expert_forest* resF;
    opnd_type resultType;
  public:
    unary_operation(const unary_opname* code, int kl, int al,
      expert_forest* arg, expert_forest* res);

    unary_operation(const unary_opname* code, int kl, int al,
      expert_forest* arg, opnd_type res);

  protected:
    virtual ~unary_operation();

  public:
    inline bool matches(const expert_forest* arg, const expert_forest* res) 
      const { 
        return (arg == argF && res == resF); 
      }

    inline bool matches(const expert_forest* arg, opnd_type res) const { 
      return (arg == argF && resultType == res); 
    }

    // high-level front-ends
    virtual void compute(const dd_edge &arg, dd_edge &res);
    virtual void compute(const dd_edge &arg, long &res);
    virtual void compute(const dd_edge &arg, double &res);
    virtual void compute(const dd_edge &arg, ct_object &c);

    // TBD: low-level front-ends?
    // e.g.,
    // virtual int computeDD(int k, int p);
    // virtual void computeEvDD(int k, int v, int p, int &w, int &q);
};

// ******************************************************************
// *                     binary_operation class                     *
// ******************************************************************

/** Mechanism to apply a binary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::binary_operation : public operation {
  protected:
    expert_forest* arg1F;
    expert_forest* arg2F;
    expert_forest* resF;
    opnd_type resultType;
  public:
    binary_operation(const binary_opname* code, int kl, int al,
      expert_forest* arg1, expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~binary_operation();

  public:
    inline bool matches(const expert_forest* arg1, const expert_forest* arg2, 
      const expert_forest* res) const { 
        return (arg1 == arg1F && arg2 == arg2F && res == resF); 
      }

    // high-level front-end
    virtual void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
      = 0;

    // low-level front ends

    /// Low-level compute on nodes a and b, return result.
    virtual int compute(int a, int b);
    /// Low-level compute at level k on nodes a and b, return result.
    virtual int compute(int k, int a, int b);
};

// ******************************************************************
// *                   numerical_operation  class                   *
// ******************************************************************

/** Mechanism to apply numerical operations to specific edges.
*/
class MEDDLY::numerical_operation : public operation {
  public:
    numerical_operation(const numerical_opname* code);
  protected:
    virtual ~numerical_operation();
  public:
    /// compute y += some function of x, depending on the operation.
    virtual void compute(double* y, const double* x);
};

// ******************************************************************
// *                      op_initializer class                      *
// ******************************************************************

/** Preferred mechanism for users to initialize their own operations.
    Derive a class from this one, provide the \a execute method.
*/
class MEDDLY::op_initializer {
  op_initializer* before;
public:
  /// Constructor.
  ///   @param  bef   initializer(s) to execute before this one.
  op_initializer(op_initializer* bef);

  virtual ~op_initializer();

  void initChain(const settings &s);
  void cleanupChain();

protected:

  virtual void init(const settings &s) = 0;
  virtual void cleanup() = 0;
};


// ****************************************************************************
// *                                                                          *
// *                          Implementation details                          *
// *                                                                          *
// *            Everything below here  can (and should) be ignored            *
// *                                                                          *
// ****************************************************************************

inline
MEDDLY::unary_operation* 
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg, 
  const dd_edge& res)
{
  return getOperation(code, 
    (expert_forest*) arg.getForest(), 
    (expert_forest*) res.getForest()
  );
}

inline
MEDDLY::unary_operation* 
MEDDLY::getOperation(const unary_opname* code, const dd_edge& arg, 
  opnd_type res)
{
  return getOperation(code, (expert_forest*) arg.getForest(), res);
}

inline
MEDDLY::binary_operation* 
MEDDLY::getOperation(const binary_opname* code, const dd_edge& arg1, 
  const dd_edge& arg2, const dd_edge& res)
{
  return getOperation(code, 
    (expert_forest*) arg1.getForest(), 
    (expert_forest*) arg2.getForest(), 
    (expert_forest*) res.getForest()
  );
}

// ****************************************************************************

inline
int MEDDLY::expert_forest::getInternalNodeSize(int p) const
{
  MEDDLY_DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return *(getNodeAddress(p) + 2);
}


inline
int* MEDDLY::expert_forest::getAddress(int k, int offset) const
{
  MEDDLY_DCASSERT(level != 0 && level[mapLevel(k)].data != 0);
  return  (level[mapLevel(k)].data + offset);
}


inline
int* MEDDLY::expert_forest::getNodeAddress(int p) const
{
  MEDDLY_DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return getAddress(getNodeLevel(p), getNodeOffset(p));
}


inline
int MEDDLY::expert_forest::getNodeOffset(int p) const
{
  MEDDLY_DCASSERT(isValidNodeIndex(p));
  return  (address[p].offset);
}

// map the level value (which could be primed) to the correct index in
// the level[]
// map logical level to physical level (index of level in the level vector)
inline
int MEDDLY::expert_forest::mapLevel(int k) const
{
  return (k >= 0)? 2 * k: ((-2 * k) - 1);
}

// map physical level to logical level (level in terms of prime, unprime)

inline
int MEDDLY::expert_forest::unmapLevel(int k) const
{
  return (k%2 == 0)? k/2: -(k + 1)/2;
}


inline
int MEDDLY::expert_forest::getMappedNodeHeight(int p) const
{
  return
    getNodeLevel(p) >= 0
    ? 2 * getNodeHeight(p)        // 2n
    : 2 * getNodeHeight(p) - 1;   // 2n-1
}


inline
bool MEDDLY::expert_forest::isActiveNode(int p) const
{
  return (isValidNodeIndex(p) && (isTerminalNode(p) || getNodeOffset(p) > 0));
}


inline
bool MEDDLY::expert_forest::isZombieNode(int p) const
{
  MEDDLY_DCASSERT(isValidNodeIndex(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  return (getCacheCount(p) < 0);
}


inline
int& MEDDLY::expert_forest::getInCount(int p) const
{
  MEDDLY_DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return *(getNodeAddress(p));
}


inline
int& MEDDLY::expert_forest::getCacheCount(int p) const
{
  MEDDLY_DCASSERT(isValidNodeIndex(p));
  MEDDLY_DCASSERT(!isTerminalNode(p));
  return address[p].cache_count;
}

inline
bool MEDDLY::expert_forest::isTerminalNode(int p) const
{
  return (p < 1);
}

inline
bool MEDDLY::expert_forest::getBoolean(int terminalNode) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN ||
      getEdgeLabeling() != forest::MULTI_TERMINAL);
  MEDDLY_DCASSERT(terminalNode == 0 || terminalNode == -1);
  return (terminalNode == 0)? false: true;
}

inline
int MEDDLY::expert_forest::getTerminalNode(bool booleanValue) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN ||
      getEdgeLabeling() != forest::MULTI_TERMINAL);
  return booleanValue? -1: 0;
}

inline
int MEDDLY::expert_forest::getInteger(int terminalNode) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);
  MEDDLY_DCASSERT(isTerminalNode(terminalNode) && terminalNode <= 0);
  // set 32nd bit based on 31st bit.  << gets rid of MSB; >> sign extends.
  return terminalNode << 1 >> 1;
}

inline
int MEDDLY::expert_forest::getTerminalNode(int integerValue) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::INTEGER);
  // value has to fit within 31 bits (incl. sign)
  // int(0xc0000000) == -1073741824
  // int(0x3fffffff) == +1073741823
  // MEDDLY_DCASSERT(-1073741824 <= integerValue && integerValue <= 1073741823);
  MEDDLY_DCASSERT(-1073741824 <= integerValue && integerValue <= 1073741823);
  return integerValue == 0? 0: integerValue | 0x80000000;
}

#define INLINED_REALS 1 // Inlining getReal and getTerminalNode(float)

#ifdef INLINED_REALS
// Warning: Optimizing (-O2) with inlined functions with gcc vers 4.1.3
//          produced errors in output. Assembly output showed incorrect
//          sequence of operations that caused errors. This does not occur
//          with lower optimization (-O0 and -O1). This also does not occur
//          if these two operations (getReal(int), getTerminalNode(float))
//          are made non-inlined. To overcome this issue (while keeping
//          -O2 optimization and inlined functions), memcpy is used instead
//          of typecasting.

inline
float MEDDLY::expert_forest::getReal(int term) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::REAL);
#if 0
  // Warning: does not work when compiled with -O2 with gcc 4.1.3
  return *(float*)((void*)&(term <<= 1));
#else
  if (term == 0) return 0.0;
  term <<= 1;
#if 1
  float ret;
  memcpy(&ret, &term, sizeof(float));
  return ret;
#else
  return toFloat(term);
#endif
#endif
}

inline
int MEDDLY::expert_forest::getTerminalNode(float a) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::REAL);
#if 0
  // Warning: does not work when compiled with -O2 with gcc 4.1.3
  return (a == 0.0)? 0: (*((int*)((void*)&a)) >> 1) | 0x80000000;
#else
  if (a == 0.0) return 0;
#if 1
  int node;
  memcpy(&node, &a, sizeof(int));
  // printf("%x\n", node);
  return (node >> 1) | 0x80000000;
#else
  return (toInt(a) >> 1) | 0x80000000;
#endif
#endif
}

#endif

inline
int MEDDLY::expert_forest::createTempNodeMaxSize(int lh, bool clear)
{
  return createTempNode(lh, getLevelSize(lh), clear);
}

inline
int MEDDLY::expert_forest::getLevelSize(int lh) const {
  MEDDLY_DCASSERT(isValidLevel(lh));
  MEDDLY_DCASSERT(lh == 0 || level[mapLevel(lh)].data != NULL);
  if (lh < 0) {
    return static_cast<expert_domain*>(d)->getVariableBound(-lh, true);
  } else {
    return static_cast<expert_domain*>(d)->getVariableBound(lh, false);
  }
}

inline
int MEDDLY::expert_forest::getNodeLevel(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  MEDDLY_DCASSERT(isActiveNode(p) || isZombieNode(p));
  return (isTerminalNode(p)? 0: address[p].level);
}


inline
int MEDDLY::expert_forest::getNodeHeight(int p) const
{
  MEDDLY_DCASSERT(isActiveNode(p));
  return level[mapLevel(getNodeLevel(p))].height;
}


inline
bool MEDDLY::expert_forest::isFullNode(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d, size: %d\n", __func__, p, getInternalNodeSize(p));
#endif
  MEDDLY_DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return (getInternalNodeSize(p) > 0);
}


inline
bool MEDDLY::expert_forest::isSparseNode(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d, size: %d\n", __func__, p, getInternalNodeSize(p));
#endif
  MEDDLY_DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return (getInternalNodeSize(p) < 0);
}


inline
int MEDDLY::expert_forest::getFullNodeSize(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  MEDDLY_DCASSERT(isFullNode(p));
  return getInternalNodeSize(p);
}


inline
int MEDDLY::expert_forest::getSparseNodeSize(int p) const
{
#ifdef DEBUG_MDD_H
printf("%s: p: %d\n", __func__, p);
#endif
  MEDDLY_DCASSERT(isSparseNode(p));
  return -getInternalNodeSize(p);
}


inline
void MEDDLY::expert_forest::setDownPtr(int p, int i, int value)
{
  MEDDLY_DCASSERT(!isReducedNode(p));
  MEDDLY_DCASSERT(isFullNode(p));
  MEDDLY_DCASSERT(isActiveNode(value));
  MEDDLY_CHECK_RANGE(0, i, getFullNodeSize(p));
  int temp = *(getNodeAddress(p) + 3 + i);
  // linkNode to new node
  linkNode(value);
  *(getNodeAddress(p) + 3 + i) = value;
  // unlinkNode old node
  unlinkNode(temp);
}


inline
void MEDDLY::expert_forest::setDownPtrWoUnlink(int p, int i, int value)
{
  MEDDLY_DCASSERT(!isReducedNode(p));
  MEDDLY_DCASSERT(isFullNode(p));
  MEDDLY_DCASSERT(isActiveNode(value));
  MEDDLY_CHECK_RANGE(0, i, getFullNodeSize(p));
  // linkNode to new node
  linkNode(value);
  *(getNodeAddress(p) + 3 + i) = value;
}


inline
void MEDDLY::expert_forest::setEdgeValue(int node, int index, int ev)
{
  MEDDLY_DCASSERT(!isReducedNode(node));
  MEDDLY_DCASSERT(isFullNode(node));
  MEDDLY_CHECK_RANGE(0, index, getFullNodeSize(node));
  *(getNodeAddress(node) + 3 + getFullNodeSize(node) + index) = ev;
}


inline
void MEDDLY::expert_forest::setEdgeValue(int node, int index, float ev)
{
  MEDDLY_DCASSERT(!isReducedNode(node));
  MEDDLY_DCASSERT(isFullNode(node));
  MEDDLY_CHECK_RANGE(0, index, getFullNodeSize(node));
  *(getNodeAddress(node) + 3 + getFullNodeSize(node) + index) = toInt(ev);
}


inline
bool MEDDLY::expert_forest::getDownPtrs(int node, int*& dptrs)
{
  MEDDLY_DCASSERT(isActiveNode(node));
  if (isTerminalNode(node) || isReducedNode(node)) return false;
  MEDDLY_DCASSERT(isFullNode(node));
  dptrs = getNodeAddress(node) + 3;
  return true;
}


inline
bool MEDDLY::expert_forest::getEdgeValues(int node, int*& evs)
{
  MEDDLY_DCASSERT(isActiveNode(node));
  if (isTerminalNode(node) || isReducedNode(node)) return false;
  MEDDLY_DCASSERT(isFullNode(node));
  evs = getNodeAddress(node) + 3 + getFullNodeSize(node);
  return true;
}


inline
bool MEDDLY::expert_forest::getEdgeValues(int node, float*& evs)
{
  MEDDLY_DCASSERT(isActiveNode(node));
  if (isTerminalNode(node) || isReducedNode(node)) return false;
  MEDDLY_DCASSERT(isFullNode(node));
  evs = toFloat(getNodeAddress(node) + 3 + getFullNodeSize(node));
  return true;
}


inline
bool MEDDLY::expert_forest::getDownPtrs(int node, const int*& dptrs) const
{
  if (isTerminalNode(node)) return false;
  dptrs = isFullNode(node)
    ? getNodeAddress(node) + 3
    : getNodeAddress(node) + 3 + getSparseNodeSize(node);
  return true;
}


inline
bool MEDDLY::expert_forest
::getSparseNodeIndexes(int node, const int*& indexes) const
{
  MEDDLY_DCASSERT(isSparseNode(node));
  if (!isSparseNode(node)) return false;
  indexes = getNodeAddress(node) + 3;
  return true;
}


inline
bool MEDDLY::expert_forest::getEdgeValues(int node, const int*& evs) const
{
  MEDDLY_DCASSERT(isReducedNode(node));
  if (isTerminalNode(node)) return false;
  evs = isFullNode(node)
    ? getNodeAddress(node) + 3 + getFullNodeSize(node)
    : getNodeAddress(node) + 3 + (getSparseNodeSize(node) * 2);
  return true;
}


inline
bool MEDDLY::expert_forest::getEdgeValues(int node, const float*& evs) const
{
  MEDDLY_DCASSERT(isReducedNode(node));
  if (isTerminalNode(node)) return false;
  evs = toFloat(isFullNode(node)
    ? getNodeAddress(node) + 3 + getFullNodeSize(node)
    : getNodeAddress(node) + 3 + (getSparseNodeSize(node) * 2));
  return true;
}


inline
int MEDDLY::expert_forest::getFullNodeDownPtr(int p, int i) const
{
  MEDDLY_DCASSERT(isFullNode(p));
  MEDDLY_CHECK_RANGE(0, i, getFullNodeSize(p));
  return getNodeAddress(p)[3 + i];
}


inline
int MEDDLY::expert_forest::getSparseNodeDownPtr(int p, int i) const
{
  MEDDLY_DCASSERT(isSparseNode(p));
  MEDDLY_CHECK_RANGE(0, i, getSparseNodeSize(p));
  return *(getNodeAddress(p) + 3 + getSparseNodeSize(p) + i);
}


inline
int MEDDLY::expert_forest::getSparseNodeIndex(int p, int i) const
{
  MEDDLY_DCASSERT(isSparseNode(p));
  MEDDLY_CHECK_RANGE(0, i, getSparseNodeSize(p));
  return *(getNodeAddress(p) + 3 + i);
}


inline
void MEDDLY::expert_forest
::getFullNodeEdgeValue(int node, int index, int& ev) const
{
  MEDDLY_DCASSERT(isFullNode(node));
  MEDDLY_CHECK_RANGE(0, index, getFullNodeSize(node));
  ev = *(getNodeAddress(node) + 3 + getFullNodeSize(node) + index);
}


inline
void MEDDLY::expert_forest
::getSparseNodeEdgeValue(int node, int index, int& ev) const
{
  MEDDLY_DCASSERT(isSparseNode(node));
  MEDDLY_CHECK_RANGE(0, index, getSparseNodeSize(node));
  ev = *(getNodeAddress(node) + 3 + (getSparseNodeSize(node) * 2) + index);
}


inline
void MEDDLY::expert_forest
::getFullNodeEdgeValue(int node, int index, float& ev) const
{
  MEDDLY_DCASSERT(isFullNode(node));
  MEDDLY_CHECK_RANGE(0, index, getFullNodeSize(node));
  ev = toFloat(*(getNodeAddress(node) + 3 + getFullNodeSize(node) + index));
}


inline
void MEDDLY::expert_forest
::getSparseNodeEdgeValue(int node, int index, float& ev) const
{
  MEDDLY_DCASSERT(isSparseNode(node));
  MEDDLY_CHECK_RANGE(0, index, getSparseNodeSize(node));
  ev = toFloat(*(getNodeAddress(node) + 3 +
        (getSparseNodeSize(node) * 2) + index));
}


inline
int MEDDLY::expert_forest::linkNode(int p)
{ 
  MEDDLY_DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return p;
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));

  // increase incount
  ++getInCount(p);

  if (getInCount(p) == 1) {
    reclaimOrphanNode(p);
  }

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t+Node %d count now %d\n", p, getInCount(p));
  fflush(stdout);
#endif
  return p;
}  



inline
void MEDDLY::expert_forest::unlinkNode(int p)
{
  MEDDLY_DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return;
  MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
  MEDDLY_DCASSERT(getInCount(p) > 0);
  
  // decrement incoming count
  --getInCount(p);

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t-Node %d count now %d\n", p, getInCount(p));
  fflush(stdout);
#endif

  if (getInCount(p) == 0) { handleNewOrphanNode(p); }
}



inline
int MEDDLY::expert_forest::cacheNode(int p)
{
  MEDDLY_DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return p;
  MEDDLY_DCASSERT(isReducedNode(p));
  getCacheCount(p)++;
#ifdef TRACK_CACHECOUNT
  fprintf(stdout, "\t+Node %d is in %d caches\n", p, getCacheCount(p));
  fflush(stdout);
#endif
  return p;
}



inline
void MEDDLY::expert_forest::uncacheNode(int p)
{
  if (isTerminalNode(p)) return;
  MEDDLY_DCASSERT(isActiveNode(p) ||
      (!isActiveNode(p) && isPessimistic() && isZombieNode(p)));

  if (isPessimistic() && isZombieNode(p)) {
    MEDDLY_DCASSERT(getCacheCount(p) < 0);
    getCacheCount(p)++;                           // special case; stored -ve
    if (getCacheCount(p) == 0) {
      freeZombieNode(p);
    }
    return;
  }

  MEDDLY_DCASSERT(getCacheCount(p) > 0);
  getCacheCount(p)--;
#ifdef TRACK_CACHECOUNT
  fprintf(stdout, "\t-Node %d is in %d caches\n", p, getCacheCount(p));
  fflush(stdout);
#endif

  if (getCacheCount(p) == 0 && getInCount(p) == 0) {
    deleteOrphanNode(p);
  }
}


inline
const MEDDLY::domain* MEDDLY::expert_forest::getDomain() const
{
  return d;
}


inline
MEDDLY::domain* MEDDLY::expert_forest::useDomain()
{
  return d;
}


inline
bool MEDDLY::expert_forest::isForRelations() const
{
  return isRelation;
}


inline
MEDDLY::forest::range_type MEDDLY::expert_forest::getRangeType() const
{
  return rangeType;
}


inline
MEDDLY::forest::edge_labeling MEDDLY::expert_forest::getEdgeLabeling() const
{
  return edgeLabel;
}


inline
MEDDLY::forest::reduction_rule MEDDLY::expert_forest::getReductionRule() const
{
  return reductionRule;
}


inline
MEDDLY::forest::node_storage MEDDLY::expert_forest::getNodeStorage() const
{
  return nodeStorage;
}


inline
MEDDLY::forest::node_deletion_policy 
MEDDLY::expert_forest::getNodeDeletion() const
{
  return nodeDeletionPolicy;
}


inline
void MEDDLY::expert_forest::setNodeDeletion(node_deletion_policy np)
{
  if (np == forest::NEVER_DELETE)
    throw error(error::NOT_IMPLEMENTED);
  if (getCurrentNumNodes() > 0)
    throw error(error::INVALID_OPERATION);
  nodeDeletionPolicy = np;
}

inline
void MEDDLY::expert_forest::setNodeStorage(node_storage ns)
{
  if (nodeStorage == ns) return;
  if (getCurrentNumNodes() > 0)
    throw error(error::INVALID_OPERATION);
  nodeStorage = ns;
}

inline
void MEDDLY::expert_forest::setReductionRule(reduction_rule r)
{
  if (reductionRule == r) return;
  if (getCurrentNumNodes() > 0)
    throw error(error::INVALID_OPERATION);
  if (!isForRelations() && r == forest::IDENTITY_REDUCED) {
    // cannot have IDENTITY reduced for non-relation forests.
    throw error(error::INVALID_OPERATION);
  }
  reductionRule = r;
}


inline
bool MEDDLY::expert_forest::isPessimistic() const {
  return nodeDeletionPolicy == forest::PESSIMISTIC_DELETION;
}


inline
bool MEDDLY::expert_forest::isMdd() const {
  return !isForRelations() &&
         getRangeType() == forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool MEDDLY::expert_forest::isMtMdd() const {
  return !isForRelations() &&
         // same as == INTEGER || == REAL
         getRangeType() != forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool MEDDLY::expert_forest::isMxd() const {
  return isForRelations() &&
         getRangeType() == forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool MEDDLY::expert_forest::isMtMxd() const {
  return isForRelations() &&
         // same as == INTEGER || == REAL
         getRangeType() != forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool MEDDLY::expert_forest::isEvplusMdd() const {
  return !isForRelations() &&
         getEdgeLabeling() == forest::EVPLUS;
}


inline
bool MEDDLY::expert_forest::isEvtimesMdd() const {
  return !isForRelations() &&
         getEdgeLabeling() == forest::EVTIMES;
}


inline
bool MEDDLY::expert_forest::isIndexSet(int node) const {
  return getIndexSetCardinality(node) > 0;
}


inline
int MEDDLY::expert_forest::getIndexSetCardinality(int node) const {
  MEDDLY_DCASSERT(isEvplusMdd());
  MEDDLY_DCASSERT(isActiveNode(node));
  // Cardinality is stored just after the downpointers and edge-values
  return
    isTerminalNode(node)
    ? node == 0? 0: 1
    : *(
        getNodeAddress(node) +
        3 +
        (isFullNode(node)
         ? 2 * getFullNodeSize(node)
         : 3 * getSparseNodeSize(node))
       );
}


inline
void MEDDLY::expert_forest::setIndexSetCardinality(int node, int c) {
  MEDDLY_DCASSERT(isEvplusMdd());
  MEDDLY_DCASSERT(isActiveNode(node));
  if (isEvplusMdd() && isActiveNode(node) && !isTerminalNode(node)) {
    *(
        getNodeAddress(node) +
        3 +
        (isFullNode(node)
         ? 2 * getFullNodeSize(node)
         : 3 * getSparseNodeSize(node))
     ) = c;
     return;
  }
  throw error(error::INVALID_OPERATION);
}


inline
bool MEDDLY::expert_forest::accumulate(int& A, int* vlist, int* vplist) {
  return false;
}

inline
bool MEDDLY::expert_forest::accumulate(int& A, int* B) {
  return false;
}


inline
float MEDDLY::toFloat(int a) {
  union { int i; float f; } n = {a};
  return n.f;
}


inline
int MEDDLY::toInt(float a) {
  union { float f; int i; } n = {a};
  return n.i;
}


inline
float* MEDDLY::toFloat(int* a) {
  union { int* i; float* f; } n = {a};
  return n.f;
}


#endif
