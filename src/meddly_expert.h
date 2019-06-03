
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

#define OLD_NODE_HEADERS
// #define COMPACTED_HEADERS

namespace MEDDLY {

  // classes defined here

  class expert_variable;
  class expert_domain;

  // wrapper for temporary nodes
  class unpacked_node;  // replacement for node_reader, node_builder

  // EXPERIMENTAL - matrix wrappers for unprimed, primed pairs of nodes
  // class unpacked_matrix;
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
  class constrained_opname;

  class ct_initializer;
  class compute_table_style;
  class compute_table;

  class operation;
  class unary_operation;
  class binary_operation;
  class specialized_operation;

  class global_rebuilder;

  // classes defined elsewhere
  class base_table;
  class unique_table;

  class reordering_base;

  // ******************************************************************
  // *                                                                *
  // *                   Named numerical operations                   *
  // *                                                                *
  // ******************************************************************

  /** Computes y = y + xA.
      x and y are vectors, stored explicitly, and A is a matrix.
      x_ind and y_ind specify how minterms are mapped to indexes
      for vectors x and y, respectively.
  */
  extern const numerical_opname* EXPLVECT_MATR_MULT;
  // extern const numerical_opname* VECT_MATR_MULT; // renamed!

  /** Computes y = y + Ax.
      x and y are vectors, stored explicitly, and A is a matrix.
      x_ind and y_ind specify how minterms are mapped to indexes
      for vectors x and y, respectively.
  */
  extern const numerical_opname* MATR_EXPLVECT_MULT;
  // extern const numerical_opname* MATR_VECT_MULT; // renamed!

  // ******************************************************************
  // *                                                                *
  // *                  Named saturation operations                   *
  // *                                                                *
  // ******************************************************************

  /** Forward reachability using saturation.
      Transition relation is already known.
  */
  extern const satpregen_opname* SATURATION_FORWARD;

  /** Backward reachability using saturation.
      Transition relation is already known.
  */
  extern const satpregen_opname* SATURATION_BACKWARD;

  /** Forward reachability using saturation.
      Transition relation is not completely known,
      will be built along with reachability set.
  */
  extern const satotf_opname* SATURATION_OTF_FORWARD;

  /** Forward reachability using saturation.
      Transition relation is specified implicitly.
  */
  extern const satimpl_opname* SATURATION_IMPL_FORWARD;

  /** Minimum-witness operations.
  */
  extern const constrained_opname* CONSTRAINED_BACKWARD_BFS;
  extern const constrained_opname* CONSTRAINED_FORWARD_DFS;
  extern const constrained_opname* CONSTRAINED_BACKWARD_DFS;
  extern const constrained_opname* TRANSITIVE_CLOSURE_DFS;

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

  /**
    Build list of initializers for Meddly.
    Custom-built initialization lists will usually include this list.
      @param    prev    Initializers to execute before the default list;
                        can be null.

      @return   List of initializers.
  */
  initializer_list* defaultInitializerList(initializer_list* prev);


  /** Initialize the library with custom settings.
      Should be called before using any other functions.
        @param  L   List of initializers.  Will execute the "setup()"
                    methods in order now, and the "cleanup()" methods
                    in reverse order on library cleanup.
  */
  void initialize(initializer_list* L);


}; // namespace MEDDLY



// ******************************************************************
// *                                                                *
// *                     expert_variable  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_variable : public variable {
  public:
    expert_variable(int b, char* n);

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
                      If bound<=0, the variable is marked as extensible,
                      with initial bound as abs(bound).
                      Note: an extensible variable has a range [1 .. +infinity].
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

    virtual ~expert_variable();
};


// ******************************************************************
// *                                                                *
// *                      expert_domain  class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_domain : public domain {
  public:
    expert_domain(variable**, int);

    virtual void createVariablesBottomUp(const int* bounds, int N);

    /** Create all variables at once, from the top down.
      Requires the domain to be "empty" (containing no variables or forests).
      @param  bounds  Current variable bounds.
                      bounds[0] gives the bound for the top-most variable,
                      and bounds[N-1] gives the bound for the bottom-most
                      variable.
                      If bound<=0, the variable is marked as extensible,
                      with initial bound as abs(bound).
                      Note: an extensible variable has a range [1 .. +infinity].
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

    expert_variable* getExpertVar(int lev) const;
    const expert_variable* readExpertVar(int lev) const;

    /** Add a new variable with bound 1.
      Can be used when the domain already has forests, in which case
      all forests are modified as appropriate.
      @param  below   Placement information: the new variable will appear
                      immediately above the level \a below.
    */
    void createVariable(int below);


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
      @param  vh     Variable handle.
      @param  prime   If prime is true, enlarge the bound for
                      the primed variable only, otherwise both
                      the primed and unprimed are enlarged.
      @param  b       New bound, if less than the current bound
                      an error code is returned.
                      If bound<=0, the variable is marked as extensible,
                      with initial bound as abs(bound).
                      Note: an extensible variable has a range [1 .. +infinity].
    */
    void enlargeVariableBound(int vh, bool prime, int b);

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
    void shrinkVariableBound(int vh, int b, bool force);

    virtual void write(output &s) const;
    virtual void read(input &s);

  protected:
    ~expert_domain();
};


// ******************************************************************
// *                                                                *
// *                      relation_node  class                      *
// *                                                                *
// ******************************************************************

/** Pieces of an implicit relation.
 
 Each piece (this class) is a function of a single variable.
 The function specifies the "next state" for the given
 current state, but we are only allowed to depend on one
 state variable.
 
 This is an abstract base class.  The function is specified
 by deriving a class from this one, and specifying method
 nextOf().
 TBD - Need a bogus value for nextOf()...
 
 Additionally, you must specify method equals(), which is used
 to detect when two nodes actually represent the same function.
 
 */
typedef int rel_node_handle;
class MEDDLY::relation_node {
public:
  /** Constructor.
   @param  signature   Hash for this node, such that
   two equal nodes must have the same
   signature.
   @param  level          Level affected.
   @param  down           Handle to a relation node below us. 
   */
  relation_node(unsigned long signature, int level, rel_node_handle down);
  virtual ~relation_node();
  
  // the following should be inlined in meddly_expert.hh
  
  /** A signature for this function.
   This helps class implicit_relation to detect duplicate
   functions (using a hash table where the signature
   is taken as the hash value).
   */
  unsigned long getSignature() const;
  
  /** The state variable affected by this part of the relation.
   */
  int getLevel() const;
  
  /** Pointer to the (ID of the) next piece of the relation.
   */
  rel_node_handle getDown() const;
  
  /** The unique ID for this piece.
   */
  rel_node_handle getID() const;
  
  /** Set the unique ID for this piece.
   */
  void setID(rel_node_handle ID);
  
  /** The token_update array for this piece.
   */
  long* getTokenUpdate() const;
  
  /** Set the token_update array for this piece.
   */
  void setTokenUpdate(long* token_update);
  
  /** The size of token_update array for this piece.
   */
  long getPieceSize() const;
  
  /** Set the size of token_update array for this piece.
   */
  void setPieceSize(long pS);
  
  /** Expand the tokenUpdate array as the variable increases
   */
  void expandTokenUpdate(long i);
  
  /** Set the tokenUpdate array at location i to val
   */
  void setTokenUpdateAtIndex(long i,long val);
  
  // the following must be provided in derived classes.
  
  /** If the variable at this level has value i,
   what should the new value be?
   */
  virtual long nextOf(long i);
  
  /** Determine if this node is equal to another one.
   */
  virtual bool equals(const relation_node* n) const;
  
private:
  unsigned long signature;
  int level;
  rel_node_handle down;
  rel_node_handle ID;
  long* token_update;
  long piece_size;        
  
  // used by the hash table in implicit_relation
  relation_node* hash_chain;
  
  // friend class implicit_relation;
};  // class relation_node

// ******************************************************************
// *                                                                *
// *                      unpacked_node  class                      *
// *                                                                *
// ******************************************************************

/** Class for reading nodes.
    Ideally - used anywhere we want to read node data.
    Backend implementation may change :^)
    Implemented in node_wrappers.cc.
    Readers may be "full" or "sparse",
    regardless of how the actual node is stored.

    Currently, a node is:
      array of node_handles, for the downward pointers
      array of integers, for indexes (sparse only)
      "compact" chunk of memory for edge values
*/
class MEDDLY::unpacked_node { 
  public:
    /**
        Options for filling an unpacked node from an existing one.
    */
    enum storage_style {
      /// Unpacked node should be stored as truncated full
      FULL_NODE,
      /// Unpacked node should be stored sparsely
      SPARSE_NODE,
      /// Unpacked node should be stored same as packed node
      AS_STORED
    };

  public:
    /** Constructor.
     The class must be "filled" by a forest before
     it can be used, however.
     */
    unpacked_node();

    /// Destructor.
    ~unpacked_node();

    /// Free memory, but don't delete.
    void clear();

  public:
  /* Initialization methods, primarily for reading */

    void initFromNode(const expert_forest *f, node_handle node, bool full);
    void initFromNode(const expert_forest *f, node_handle node, storage_style st2);

    void initRedundant(const expert_forest *f, int k, node_handle node, bool full);
    void initRedundant(const expert_forest *f, int k, int ev, node_handle node, bool full);
    void initRedundant(const expert_forest *f, int k, long ev, node_handle node, bool full);
    void initRedundant(const expert_forest *f, int k, float ev, node_handle node, bool full);

    void initIdentity(const expert_forest *f, int k, int i, node_handle node, bool full);
    void initIdentity(const expert_forest *f, int k, int i, int ev, node_handle node, bool full);
    void initIdentity(const expert_forest *f, int k, int i, long ev, node_handle node, bool full);
    void initIdentity(const expert_forest *f, int k, int i, float ev, node_handle node, bool full);

  /* Create blank node, primarily for writing */
    
    void initFull(const expert_forest *f, int level, int tsz); 
    void initSparse(const expert_forest *f, int level, int nnz);

  public:
  /* For convenience: get recycled instance and initialize */

    static unpacked_node* newFromNode(const expert_forest *f, node_handle node, bool full);
    static unpacked_node* newFromNode(const expert_forest *f, node_handle node, storage_style st2);

    static unpacked_node* newRedundant(const expert_forest *f, int k, node_handle node, bool full);
    static unpacked_node* newRedundant(const expert_forest *f, int k, long ev, node_handle node, bool full);
    static unpacked_node* newRedundant(const expert_forest *f, int k, float ev, node_handle node, bool full);

    static unpacked_node* newIdentity(const expert_forest *f, int k, int i, node_handle node, bool full);
    static unpacked_node* newIdentity(const expert_forest *f, int k, int i, long ev, node_handle node, bool full);
    static unpacked_node* newIdentity(const expert_forest *f, int k, int i, float ev, node_handle node, bool full);

    static unpacked_node* newFull(const expert_forest *f, int level, int tsz);
    static unpacked_node* newSparse(const expert_forest *f, int level, int nnz);

  public:
  /* Display / access methods */

    /** Write a node in human-readable format.

        @param  s       Output stream.
        @param  details Should we show "details" or not.
    */
    void show(output &s, bool details) const;

    /** Write a node in machine-readable format.

        @param  s       Output stream.
        @param  map     Translation to use on node handles.
                        Allows us to renumber nodes as we write them.
    */
    void write(output &s, const node_handle* map) const;

    /// Get a pointer to the unhashed header data.
    const void* UHptr() const;

    /// Modify a pointer to the unhashed header data
    void* UHdata();

    /// Get the number of bytes of unhashed header data.
    int UHbytes() const;

    /// Get a pointer to the hashed header data.
    const void* HHptr() const;

    /// Modify a pointer to the hashed header data.
    void* HHdata();

    /// Get the number of bytes of hashed header data.
    int HHbytes() const;

    /** Get a downward pointer.
          @param  n   Which pointer.
          @return     If this is a full reader,
                      return pointer with index n.
                      If this is a sparse reader,
                      return the nth non-zero pointer.
    */
    node_handle d(int n) const;

    /** Reference to a downward pointer.
          @param  n   Which pointer.
          @return     If this is a full reader,
                      modify pointer with index n.
                      If this is a sparse reader,
                      modify the nth non-zero pointer.
    */
    node_handle& d_ref(int n);

    /** Get the index of the nth non-zero pointer.
        Use only for sparse readers.
     */
    int i(int n) const;

    /** Modify the index of the nth non-zero pointer.
        Use only for sparse readers.
    */
    int& i_ref(int n);

    /// Get a pointer to an edge
    const void* eptr(int i) const;

    /// Modify pointer to an edge
    void* eptr_write(int i);

    /// Get the edge value, as an integer.
    void getEdge(int i, long& ev) const;

    /// Get the edge value, as a float.
    void getEdge(int i, float& ev) const;

    /// Set the edge value, as an integer.
    void setEdge(int i, long ev);

    /// Set the edge value, as a float.
    void setEdge(int i, float ev);

    /// Get the edge value, as an integer.
    long ei(int i) const;

    /// Get the edge value, as a float.
    float ef(int i) const;

    // -------------------------------------------------------------------------
    // Methods to access the extensible portion of the node
    //
    // Note: Only nodes that belong to an extensible level can be extensible.
    // 
    // Extensible node  : (Regular part, Extensible part).
    // Regular Part     : same as non-extensible nodes.
    // Extensible Part  : a single edge, i.e. a tuple <index, node, edge-value>,
    //                    for all indices in [extensible_index, +infinity].
    // 
    // Example:
    // The node
    //    [0:n0:ev0, 1:n1:ev1, 2:n2:ev2, 3:n2:ev2, ..., +infinity:n2:ev2]
    // is represented as the following extensible node:
    //    ([0:n0:ev0, 1:n1:ev1], Extensible: [2:n2:ev2])
    // with,
    //    size of node           : 2
    //    extensible index       : 2
    //    extensible node_handle : n2
    //    extensible edge-value  : ev2
    // -------------------------------------------------------------------------

    /// Does this reader store an extensible edge?
    /// Note: in an extensible node, all edges starting from
    ///       index to +infinity refer to the same edge, i.e. (node, edge-value).
    bool isExtensible() const;

    void markAsExtensible();
    void markAsNotExtensible();

    int ext_d() const;
    int ext_i() const;
    int ext_ei() const;
    float ext_ef() const;

    // -------------------- End of extensible portion --------------------------

    /// Get the level number of this node.
    int getLevel() const;

    /// Set the level number of this node.
    void setLevel(int k);

    /// Get the size of this node (full readers only).
    int getSize() const;

    /// Get the number of nonzeroes of this node (sparse readers only).
    int getNNZs() const;

    /// Is this a sparse reader?
    bool isSparse() const;

    /// Is this a full reader?
    bool isFull() const;

    /// Does this node have edge values?
    bool hasEdges() const;

    /// Number of bytes per edge
    int edgeBytes() const;

    // For debugging unique table.
    unsigned hash() const;

    void setHash(unsigned H);

    void computeHash();

    /// Removes redundant trailing edges.
    /// If the unpacked node is sparse, it assumes its indices to be in ascending order.
    void trim();

    /// If the unpacked node is sparse, it is sorted so that the
    /// indices are in ascending order.
    void sort();

    /// Checks if the node is has no trailing redundant edges
    bool isTrim() const;

    // checks if the node indices are in ascending order
    bool isSorted() const;


  public:

    /// Change the size of a node
    void resize(int ns);

    /// Shrink the size of a (truncated) full node
    void shrinkFull(int ns);

    /// Shrink the size of a sparse node
    void shrinkSparse(int ns);

    /// Called within expert_forest to allocate space.
    ///   @param  p     Parent.
    ///   @param  k     Level number.
    ///   @param  ns    Size of node.
    ///   @param  full  If true, we'll be filling a full reader.
    ///                 Otherwise it is a sparse one.
    void bind_to_forest(const expert_forest* p, int k, int ns, bool full);

    /// Called by node_storage when building an unpacked
    /// node based on how it's stored.
    void bind_as_full(bool full);

  public:
    // Centralized recycling
    static unpacked_node* useUnpackedNode();
    static void recycle(unpacked_node* r);
    static void freeRecycled();


  private:
    const expert_forest* parent;
    static unpacked_node* freeList;
    unpacked_node* next; // for recycled list
    /*
      TBD - extra info that is not hashed
    */
    void* extra_unhashed;
    int ext_uh_alloc;
    int ext_uh_size;
    /*
     Extra info that is hashed
     */
    void* extra_hashed;
    int ext_h_alloc;
    int ext_h_size;
    /*
     Down pointers, indexes, edge values.
     */
    node_handle* down;
    int* index;
    void* edge;
    bool is_extensible;
    int alloc;
    int ealloc;
    int size;
    int nnzs;
    int level;
    unsigned h;
    char edge_bytes; // number of bytes for an edge value.
    bool is_full;
#ifdef DEVELOPMENT_CODE
    bool has_hash;
#endif
};




// ******************************************************************
// *                                                                *
// *                     unpacked_matrix  class                     *
// *                                                                *
// ******************************************************************

/** Class for unpacked matrices: unprimed, primed "2 levels" of nodes.
    Implemented in node_wrappers.cc.

    TBD : inline what we should
*/

#if 0

class MEDDLY::unpacked_matrix {

  public:
    /**
      Different internal storage types.
      Affects the allowed read/write access methods.
    */
    enum storage_type {
      FULL_FULL     = 0,  // Full storage; acts like a 2-d array
      FULL_SPARSE   = 1,  // Compressed row storage; each row is sparse
      SPARSE_SPARSE = 2   // Coordinate list storage; list of (i, j, down)
    };

  public:
    unpacked_matrix();
    ~unpacked_matrix();

    /**
        An instance must be initialized, using init(), before use.
        Ties the instance to a particular level of a forest.
        Can be re-initialized with a different level or forest.
          @param  k   Level.
          @param  p   Forest.
    */
    void init(int k, const expert_forest* p);
    
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      Basic matrix information
        
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /// How the matrix can be accessed.
    storage_type howStored() const;

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      Methods for filling the unpacked matrix from nodes,
      primarily used by forests.
        
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /// Clear node and use FULL_FULL storage
    void initFullFull(int rows, int cols);

    /// Clear node and use FULL_SPARSE storage
    void initFullSparse(int rows, int maxnzs);

    /// Clear node and use SPARSE_SPARSE storage
    void initSparseSparse(int maxnzs);

    /** Add a non-zero element to the matrix.
        Must be called in lexicographical order of rows, columns.
        I.e., enumerate the rows in increasing order, and for each row,
        enumerate the columns in increasing order.
    */
    void appendElement(int i, int j, node_handle d);

    // tbd - edge version of appendElement.

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      FULL_FULL storage access:  methods for reading and modifying
      matrix entries.
        
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    node_handle down(int i, int j) const;

    void set_down(int i, int j, node_handle d);


    // tbd - edges?


    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      FULL_SPARSE and SPARSE_SPARSE storage access:  
      methods for reading matrix entries.
        
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /** 
      For FULL_SPARSE: Get nonzero index for first element in row i.
      Equivalently, get one plus the nonzero index for last element
      in row i-1.
        @param  i   Row index, between 0 and number of rows (inclusive).
        @return     A nonzero index, can be used in methods col_index()
                    and down().
    */
    int row_start(int i) const; 

    /// For SPARSE_SPARSE: Get row for nonzero number z.
    int row_index(int z) const;

    /// Get column for nonzero number z.
    int col_index(int z) const;

    /// Get downward pointer for nonzero number z.
    node_handle down(int z) const;

    // tbd - edges?

  private:
    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      Implementation details.  Nothing to see here.
        
    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


};  // end of unpacked_matrix class

#endif

// ******************************************************************
// *                                                                *
// *                     initializer_list class                     *
// *                                                                *
// ******************************************************************

/** Mechanism for initializing and/or cleaning up library structures.
    Any user additions to the library should utilize this class.
    Derive a class from this one, provide the \a setup and \a cleanup
    methods.
    Implementation in meddly.cc
*/
class MEDDLY::initializer_list {
  public:
    /**
        Constructor.
        Takes the initializer(s) to run before this one.
        Cleanup runs in the reverse order.
    */
    initializer_list(initializer_list* previous);
    virtual ~initializer_list();

    /**
        Run all setup methods for the list of initializers,
        "previous first".
    */
    void setupAll();

    /**
        Run all cleanup methods for the list of initializers,
        "previous last".
    */
    void cleanupAll();

  protected:
    virtual void setup() = 0;
    virtual void cleanup() = 0;

  private:
    initializer_list* previous;
};

// ******************************************************************
// *                                                                *
// *                   memory_manager_style class                   *
// *                                                                *
// ******************************************************************

/** Abstract base class for memory manager factories.
    
    This allows us to use specialized implementations of
    memory managers (say, using templates) based on the granularity.

    Implementation is in memory_managers/base_manager.cc
*/
class MEDDLY::memory_manager_style {
    const char* name; 
  public:
    memory_manager_style(const char* n);
    virtual ~memory_manager_style();

    /**
        Build a new memory manager.
        
          @param  granularity   Unit of storage, in bytes.
                                Must be greater than 0.
                                All sizes specified to the memory manager,
                                for allocating and freeing chunks, are in
                                terms of the granularity.  For example, to
                                manage a collection of arrays of integers,
                                set the granularity to be sizeof(int) and
                                use sizes equal to the number of integers.
                                For behavior exactly the same as malloc,
                                use a granularity of 1.


          @param  minsize       The smallest size chunk that will ever be
                                requested.  Must be greater than 0.
                                This is specified here in case
                                that information can help the memory manager.

          @param  stats         Structure to use for updating memory stats.

          @return   A pointer to a new instance of a memory manager, or 0
                    if some error occurred, for example if the requested
                    granularity cannot be supported by this type of memory
                    manager.
    */
    virtual memory_manager* initManager(unsigned char granularity, 
      unsigned char minsize, memstats& stats) const = 0;


    /**
        Human readable name.
        Used for debugging and reporting.
    */
    const char* getName() const;
};

// ******************************************************************
// *                                                                *
// *                      memory_manager class                      *
// *                                                                *
// ******************************************************************

/** 
    Interface for memory managers.

    Implementation is in memory_managers/base_manager.cc
*/
class MEDDLY::memory_manager {

  public:
    memory_manager(const char* sn, memstats& stats);
    virtual ~memory_manager();

    /**
        Is this memory manager unable to free everything on its own?

          @return   True, if destroying the memory manager DOES NOT
                    automatically recycle all non-freed requested chunks
                    (because the memory manager does not track them).

                    False, if destroying the memory manager DOES
                    automatically recycle all non-freed requested chunks.
    */
    virtual bool mustRecycleManually() const = 0;

    /**
        Does the memory manager require that the most significant bit
        of the first slot is cleared?

        (If true, the memory manager sets this bit for recycled chunks.)

          @return   True, if the first slot of a chunk must hold a value 
                    such that the most significant bit is cleared.
                    If this cannot be guaranteed by what is stored there,
                    then you must allocate an extra slot and not use the
                    first one.

                    False, if there are no restrictions on what may be
                    stored in the first slot of a chunk.
    */
    virtual bool firstSlotMustClearMSB() const = 0;

    /**
        Does the memory manager require that the most significant bit
        of the last slot is cleared?

        (If true, the memory manager sets this bit for recycled chunks.)

          @return   True, if the last slot of a chunk must hold a value 
                    such that the most significant bit is cleared.
                    If this cannot be guaranteed by what is stored there,
                    then you must allocate an extra slot and not use the
                    last one.

                    False, if there are no restrictions on what may be
                    stored in the last slot of a chunk.
    */
    virtual bool lastSlotMustClearMSB() const = 0;


    /**
        Request a chunk of memory.

          @param  numSlots    Number of slots.  
                              INPUT: number of requested slots.
                              OUTPUT: number of slots in the given chunk,
                              might be larger than the number requested.
                              Will not be smaller than the number requested,
                              unless a failure occurred, in which case 
                              it will be set to zero.

          @return   A non-zero handle for a new chunk of memory, containing
                    numSlots slots (and requiring numSlots * granularity
                    bytes), on success.
                    Zero, on failure.

    */
    virtual node_address requestChunk(size_t &numSlots) = 0;
    
    /**
        Recycle a chunk of memory.
        
          @param  h           Handle of the chunk, as returned by 
                              method requestChunk().

          @param  numSlots    Total number of slots in the chunk.
    */
    virtual void recycleChunk(node_address h, size_t numSlots) = 0;

    /**
        Convert a handle to an actual pointer we can use.

        A default, fast, inlined implementation that will work for
        most memory managers is implemented here, based on
          address = base + m * h
        where m*h is the number of bytes to shift base by.
        This requires derived classes to maintain the pointer "base"
        and multiplier "m" by calling protected methods
            void setChunkBase(void* base)
            void setChunkMultiplier(unsigned int m)
        Note that if m is zero (its default value), this method
        will fall back to slowChunkAddress().

          @param  h     Non-null handle of the chunk, as returned by requestChunk().

          @return       If h is 0, or an invalid handle, then the result is 
                        undefined.  Otherwise, we return a pointer to the 
                        chunk given by handle h.
                        This pointer is guaranteed to be fixed until the next
                        call to requestChunk() or recycleChunk(); after that,
                        the pointer for a handle could change.
    */
    void* getChunkAddress(node_address h) const;


  protected:
    /**
        See getChunkAddress.
        This method (default behavior is assert(false))
        should be overridden in derived classes when
        the mapping from node address h to pointer
        does not follow the formula used by getChunkAddress().

        Convert a handle to an actual pointer we can use.

          @param  h     Handle of the chunk, as returned by requestChunk(), or 0.

          @return       If h is 0, then we return 0.  Otherwise, we return
                        a pointer to the chunk given by handle h.
                        This pointer is guaranteed to be fixed until the next
                        call to requestChunk() or recycleChunk(); after that,
                        the pointer for a handle could change.
    */
    virtual void* slowChunkAddress(node_address h) const;


  public:
    /**
        Check if a handle is valid (non-null or otherwise).
        Since this is not always possible,
        this method is conservative.

        @return   false if h is null or definitely invalid;
                  true otherwise.
    */
    virtual bool isValidHandle(node_address h) const = 0;

    /** Show various statistics.
          @param  s         Output stream to write to
          @param  pad       Padding string, written at the start of
                            each output line.
          @param  human     If false, just display raw byte counts. 
                            If true, use units (e.g., Mbytes, Kbytes).
          @param  details   If false, just display basic statistics.
                            If true, display details.
    */
    virtual void reportStats(output &s, const char* pad, 
      bool human, bool details) const = 0;


    /** Display manager-specific internals.
        For debugging.
          @param  s       Output stream to use
    */
    virtual void dumpInternal(output &s) const = 0;


  
    /** Get first address.
        Used to cycle through all addresses, if we can.
        If we cannot, then always return 0.
    */
    virtual node_address getFirstAddress() const = 0;

    /** Is a given address in use?
        Best effort answer only.
        If unsure, return false.
    */
    virtual bool isAddressInUse(node_address addr) const = 0;

    /** Get the next address.
        Used to cycle through all addresses, if we can.
        If we cannot, then always return 0.
          @param  addr    An address given by any of the methods that
                          "cycle through addresses", and one where
                          isAddressInUse() returned false.

          @return   Next address to check,  Presumably this address
                    is unused and is a hole of some kind.
    */
    virtual node_address getNextAddress(node_address addr) const = 0;

    /** Show information about an unused address.
        If we do not track unused addresses, then do nothing. 
          @param  addr    An address given by any of the methods that
                          "cycle through addresses", and one where
                          isAddressInUse() returned false.
    */
    virtual void dumpInternalUnused(output &s, node_address addr) const = 0;

  protected:
    void incMemUsed(size_t b);
    void decMemUsed(size_t b);
    void incMemAlloc(size_t b);
    void decMemAlloc(size_t b);

    void zeroMemUsed();
    void zeroMemAlloc();

    /**
        Set base pointer used for fast getChunkAddress().
    */
    void setChunkBase(void* p);

    /*
        Set multiplier used for fast getChunkAddress().
        If zero, we call getSlowChunkAddress().
    */
    void setChunkMultiplier(unsigned int m);

  public:
    /**
        Return the name of the style that created us.
    */
    const char* getStyleName() const;

  private:
    /// Name of the style that invoked us
    const char* style_name;
    memstats &my_mem;

    /// Base pointer for getChunkAddress
    char* chunk_base;

    /// Handle multiplier for getChunkAddress; if zero must call virtual function
    unsigned int chunk_multiplier;
};


// ******************************************************************
// *                                                                *
// *                       node_headers class                       *
// *                                                                *
// ******************************************************************

/**
    Node header information, for a collection of nodes.

    This is used within a forest to store meta-information
    about each node.  The actual node itself (i.e., the downward
    pointers and any edge values) is stored in a separate data
    structure (the node_storage class) using various encoding schemes.

    Conceptually, this class is simply an array of structs,
    plus a mechanism for recycling unused elements.  
    Node handle 0 is reserved for special use, so
    "real" nodes must have a non-zero handle.

    For each node, its header contains:
      node_address  address     The "address" of the node in the node_storage class.
      integer       level       Level number of the node.  0 indicates that the
                                node has been deleted, but we cannot recycle the
                                node handle yet (probably because it might be
                                contained in a compute table somewhere).
      natural       cache_count Optional (i.e., can be turned on/off for all nodes).
                                Number of compute table references to this handle.
      natural       incoming    Optional (i.e., can be turned on/off for all nodes).
                                Number of incoming edges to the node referred to
                                by this handle.

    In practice, we may use a different data structure for speed
    or (more likely) to save space.

    
    For now, we are using an array of structs, but this may change.


    Inlined methods are found in meddly_expert.hh.
    Non-inlined methods are found in node_headers.cc.
*/
class MEDDLY::node_headers {
  public:
    node_headers(expert_forest &P);
    ~node_headers();


    /** Show various memory stats.
          @param  s       Output stream to write to
          @param  pad     Padding string, written at the start of
                          each output line.
          @param  flags   Which stats to display, as "flags";
                          use bitwise or to combine values.
    */
    void reportStats(output &s, const char* pad, unsigned flags) const;

    /**
        Indicate that we don't need to track cache counts.
        The default is that we will track cache counts.
        For efficiency, this should be called soon after 
        construction (otherwise we may allocate space for nothing).
    */
    void turnOffCacheCounts();

    /**
        Indicate that we don't need to track incoming counts.
        The default is that we will track incoming counts.
        For efficiency, this should be called soon after 
        construction (otherwise we may allocate space for nothing).
    */
    void turnOffIncomingCounts();

    /**
        Set node recycling to pessimistic or optimistic.
        Pessimistic: disconnected nodes are recycled immediately.
        Optimistic:  disconnected nodes are recycled only when
                      they no longer appear in any caches.
    */
    void setPessimistic(bool pess);
    
  public: // node handle management

    /**
        Get an unused node handle.
        This is either a recycled one or
        the next one in the available pool
        (which will be expanded if necessary).
    */
    node_handle getFreeNodeHandle();

    /**
        Recycle a used node handle.
        The recycled handle can eventually be
        reused when returned by a call to
        getFreeNodeHandle().
    */
    void recycleNodeHandle(node_handle p);

    /**
        Get largest handle of active nodes.
    */
    node_handle lastUsedHandle() const;

    /**
        Swap two nodes.
        Used for in-place forest reordering.
          @param  p               First node
          @param  q               Second node
          @param  swap_incounts   If true, swap everything;
                                  if false, don't swap incoming counts.
    */
    void swapNodes(node_handle p, node_handle q, bool swap_incounts);

  public: // node status
    bool isActive(node_handle p) const;
    /// Is this a zombie node (dead but not able to be deleted yet)
    bool isZombie(node_handle p) const;
    /// Is this a deleted node
    bool isDeleted(node_handle p) const;
    /// Deactivated: 0 level  
    bool isDeactivated(node_handle p) const;


  public: // address stuff

    /// Get the address for node p.
    node_address getNodeAddress(node_handle p) const;

    /// Set the address for node p to a.
    void setNodeAddress(node_handle p, node_address a);

    /** Change the address for node p.
        Same as setNodeAddress, except we can verify the old address
        as a sanity check.
    */
    void moveNodeAddress(node_handle p, node_address old_addr, node_address new_addr);

  public: // level stuff

    /// Get the level for node p.
    int getNodeLevel(node_handle p) const;

    /// Set the level for node p to k.
    void setNodeLevel(node_handle p, int k);

  public: // cache count stuff

    /// Are we tracking cache counts
    bool trackingCacheCounts() const;

    /// Get the cache count for node p.
    unsigned long getNodeCacheCount(node_handle p) const;

    /// Increment the cache count for node p and return p.
    void cacheNode(node_handle p);

    /// Decrement the cache count for node p.
    void uncacheNode(node_handle p);
    
  public: // incoming count stuff

    /// Are we tracking incoming counts
    bool trackingIncomingCounts() const;

    /// Get the incoming count for node p.
    unsigned long getIncomingCount(node_handle p) const;

    /// Increment the incoming count for node p and return p.
    node_handle linkNode(node_handle p);

    /// Decrement the incoming count for node p.
    void unlinkNode(node_handle p);
    
  public: // implicit stuff
    
    /// Get whether node p is implicit.
    int getNodeImplicitFlag(node_handle p) const;
    
    /// Set as true if node p is implicit.
    void setNodeImplicitFlag(node_handle p,bool flag);


  public: // for debugging

    void dumpInternal(output &s) const;


  private:  // for debugging help
    void validateFreeLists() const;

#ifndef OLD_NODE_HEADERS
  public: // interface for node header size changes
    void changeHeaderSize(unsigned oldbits, unsigned newbits);
#endif

  private:  // helper methods

    /// Increase the number of node handles.
    void expandHandleList();

    /// Decrease the number of node handles.
    void shrinkHandleList();

    void deactivate(node_handle p);

#ifdef OLD_NODE_HEADERS

    node_handle getNextOf(node_handle p) const;
    void setNextOf(node_handle p, node_handle n);

  private:
    
    // Nothing to see here.  Move along.
    // Anything below here is subject to change without notice.


    struct node_header {
          /** Offset to node's data in the corresponding node storage structure.
              If the node is active, this is the offset (>0) in the data array.
              If the node is deleted, this is the next deleted node
              (part of the unused address list).
          */
          node_address offset;

          /** Node level
              If the node is active, this indicates node level.
          */
          int level;

          /** Cache count
              The number of cache entries that refer to this node (excl. unique
              table). 
          */
          unsigned int cache_count;

          /** Incoming count
              The number of incoming edges to this node.
          */
          unsigned int incoming_count;

          /// Is node marked.  This is pretty horrible, but is only temporary
          // bool marked;
       
          /** Implicit node
              This indicates if node is implicit
           */
          bool is_implicit = false;

    };

    /// address info for nodes
    node_header *address;
    /// Size of address/next array.
    node_handle a_size;
    /// Last used address.
    node_handle a_last;
    /// Number of recycled addresses
    node_handle a_freed;
    /// Pointer to unused address lists, based on size
    node_handle a_unused[8];  // number of bytes per handle
    /// Lowest non-empty address list
    char a_lowest_index;
    /// Next time we shink the address list.
    node_handle a_next_shrink;

    /// Do we track cache counts
    bool usesCacheCounts;
    /// Do we track incoming counts
    bool usesIncomingCounts;

    /// Are we using the pessimistic strategy?
    bool pessimistic;

    /// Parent forest, needed for recycling
    expert_forest &parent;

    static const int a_min_size = 1024;

#else

    size_t getNextOf(size_t p) const;
    void setNextOf(size_t p, size_t n);

  private:

    class level_array {
        node_headers &parent;
#ifdef COMPACTED_HEADERS
        char* data8;
        short* data16;
#endif
        int* data32;
        size_t size;
        unsigned char bytes;
      public:
        level_array(node_headers &p, int max_level);
        ~level_array();

        void expand(size_t ns);
        void shrink(size_t ns);

        int get(size_t i) const;
        void set(size_t i, int v);
        void swap(size_t i, size_t j);

        void show(output &s, size_t first, size_t last, int width) const;
        size_t entry_bits() const;
    };

    class counter_array {
        node_headers &parent;
#ifdef COMPACTED_HEADERS
        unsigned char* data8;
        unsigned short* data16;
#endif
        unsigned int* data32;
        size_t size;
#ifdef COMPACTED_HEADERS
        size_t counts_09bit;  // number of counts requiring at least 9 bits
        size_t counts_17bit;  // number of counts requiring at least 17 bits
#endif
        unsigned char bytes;
      public:
        counter_array(node_headers &p);
        ~counter_array();

        void expand(size_t ns);
        void shrink(size_t ns);

        unsigned int get(size_t i) const;
        void swap(size_t i, size_t j);
        void increment(size_t i);
        void decrement(size_t i);

        bool isZeroBeforeIncrement(size_t i);
        bool isPositiveAfterDecrement(size_t i);

        void show(output &s, size_t first, size_t last, int width) const;
        size_t entry_bits() const;

#ifdef COMPACTED_HEADERS
        // Expand from 8-bit to 16-bit entries because of element i
        void expand8to16(size_t i); 
        // Expand from 16-bit to 32-bit entries because of element i
        void expand16to32(size_t i); 

        void shrink16to8(size_t ns); 

        void shrink32to16(size_t ns); 
        void shrink32to8(size_t ns); 
#endif
    };

    class address_array {
        node_headers &parent;
#ifdef COMPACTED_HEADERS
        unsigned int* data32;
#endif
        unsigned long* data64;
        size_t size;
#ifdef COMPACTED_HEADERS
        size_t num_large_elements;
#endif
        unsigned char bytes;
      public:
        address_array(node_headers &p);
        ~address_array();

        void expand(size_t ns);
        void shrink(size_t ns);

        unsigned long get(size_t i) const;
        void set(size_t i, unsigned long v);
        void swap(size_t i, size_t j);

        void show(output &s, size_t first, size_t last, int width) const;
        size_t entry_bits() const;

#ifdef COMPACTED_HEADERS
        void expand32to64();
        void shrink64to32(size_t ns);
#endif
    };

    class bitvector {
        node_headers &parent;
        bool* data;
        size_t size;
      public:
        bitvector(node_headers &p);
        ~bitvector();

        void expand(size_t ns);
        void shrink(size_t ns);

        bool get(size_t i) const;
        void set(size_t i, bool v);
        void swap(size_t i, size_t j);
        size_t entry_bits() const;
    };

  private:
    address_array* addresses;
    level_array* levels;
    counter_array* cache_counts;
    counter_array* incoming_counts;
    bitvector* implicit_bits;

    /// Last used address.
    size_t a_last;

    /// Allocated sizes of arrays.
    size_t a_size;

    /// Next time we shink the address list.
    size_t a_next_shrink;

    /// Number of addresses in free lists
    size_t a_freed;

    /// Pointer to unused address lists, based on size
    size_t a_unused[8];  // number of bytes per handle

    /// Current size of node "header"
    unsigned h_bits;

    /// Lowest non-empty address list
    char a_lowest_index;

    /// Are we using the pessimistic strategy?
    bool pessimistic;

    /// Parent forest, needed for recycling
    expert_forest &parent;

#endif

};


// ******************************************************************
// *                                                                *
// *                    node_storage_style class                    *
// *                                                                *
// ******************************************************************

/** Abstract base class for node storage factories.

    The base class is implemented in node_wrappers.cc;
    various backends are implemented in directory storage/.

*/
class MEDDLY::node_storage_style {
    const char* name;
  public:
    node_storage_style(const char* n);
    virtual ~node_storage_style();


    /** Build a new node storage mechanism, bound to the given forest.

          @param  f   Forest to bind to
          @param  mm  Memory manager style to use

          @return     A pointer to a node storage class,
                      initialized for forest f.
    */
    virtual node_storage* createForForest(expert_forest* f, 
        const memory_manager_style* mmst) const = 0;

    const char* getName() const;
};

// ******************************************************************
// *                                                                *
// *                       node_storage class                       *
// *                                                                *
// ******************************************************************

/** Abstract base class for node storage.

    The base class is implemented in node_wrappers.cc;
    various backends are implemented in directory storage/,
    or you may implement your own scheme :^)

    Nodes are represented by an address, which for a valid
    node, will be greater than 0.  The first valid address
    must be 1.

    Whatever scheme is used to store nodes internally,
    it must be possible to set pointers \a count and \a next
    such that count[addr] and next[addr] give the incoming
    count and next pointer for the given node addr.
    Derived classes are responsible for setting up
    and maintaining these pointers.
*/
class MEDDLY::node_storage {
  public:
    node_storage(const char* sn, expert_forest* f);
    virtual ~node_storage();

    /** Go through and collect any garbage.

          @param  shrink    If true, we will shrink data structures
                            as we free space; otherwise, we won't.
    */
    virtual void collectGarbage(bool shrink) = 0;

    /** Show various stats.
          @param  s       Output stream to write to
          @param  pad     Padding string, written at the start of
                          each output line.
          @param  flags   Controls what is displayed.
    */
    virtual void reportStats(output &s, const char* pad, unsigned flags) const = 0;

    /** Dump the internal storage details.
        Primarily used for debugging.

          @param  s       Output stream to use
          @param  flags   What to show.
                            0x01  Show active memory
                            0x02  Show memory "holes"

        TBD - remove this
    */
    void dumpInternal(output &s, unsigned flags) const;

    /** Allocate space for, and store, a node.
        I.e., create a new node that is a copy of the given one.
        The node might be "compressed" in various ways to reduce
        storage requirements.  (Indeed, that is the whole point
        of the node_storage class.)
            @param  p     Node handle number, in case it is used
            @param  nb    Node data is copied from here.
            @param  opt   Ways we can store the node.
            @return       The "address" of the new node.
    */
    virtual node_address makeNode(node_handle p, const unpacked_node &nb,
                                  node_storage_flags opt) = 0;

    /** Destroy a node.
        Unlink the downward pointers, and recycle the memory
        used by the node.
            @param  addr    Address of the node.
    */
    virtual void unlinkDownAndRecycle(node_address addr) = 0;


    // various ways to read a node

    /** Check for duplicates.
          @param  addr    Node address in this structure
          @param  nr      Node to compare against

          @return true    iff the nodes are duplicates
    */
    virtual bool areDuplicates(node_address addr, const unpacked_node &nr) 
        const = 0;

    /**
        Copy the node at the specified address, into an unpacked node.
        Useful for reading an entire node.
          @param  un      Result will be stored here.  Will be resized if needed.
          @param  addr    Node address in this structure.
    */
    virtual void fillUnpacked(unpacked_node &un, node_address addr, 
        unpacked_node::storage_style st2) const = 0;

    /** Compute the hash value for a node.
        Should give the same answer as filling a unpacked_node
        and computing the hash on the unpacked_node.

          @param  levl  Level of the node of interest
          @param  addr  Address of the node of interest
    */
    virtual unsigned hashNode(int level, node_address addr) const = 0;

    /** Determine if this is an extensible node.
          @param  addr    Node Address
          @return         True if the node stores an extensible edge,
                          False otherwise.
    */
    virtual bool isExtensible(node_address addr) const = 0;

    /** Determine if this is a singleton node.
        Used for identity reductions.
          @param  addr    Address of the node we care about
          @param  down    Output:
                          The singleton downward pointer, or undefined.

          @return   If the node has only one non-zero downward pointer,
                    then return the index for that pointer.
                    Otherwise, return a negative value.
    */
    virtual int getSingletonIndex(node_address addr, node_handle &down) 
        const = 0;


    /** Get the specified downward pointer for a node.
        Fast if we just want one.
          @param  addr    Address of the node we care about
          @param  index   Index of downward pointer
          @return         Desired pointer
          @throw          INVALID_VARIABLE, if index is negative.
    */
    virtual node_handle getDownPtr(node_address addr, int index) const = 0;

    /** Get the specified outgoing edge for a node.
        Fast if we just want one.

          @param  addr    Address of the node we care about
          @param  ind     Index of the pointer we want.
          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    virtual void getDownPtr(node_address addr, int ind, int& ev,
          node_handle& dn) const = 0;
    virtual void getDownPtr(node_address addr, int ind, long& ev,
              node_handle& dn) const = 0;

    /** Get the specified outgoing edge for a node.
        Fast if we just want one.

          @param  addr    Address of the node we care about
          @param  ind     Index of the pointer we want.
          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    virtual void getDownPtr(node_address addr, int ind, float& ev,
          node_handle& dn) const = 0;


    /** Read the unhashed header portion of a node.

          @param  addr    Address of the node we care about
    */
    virtual const void* getUnhashedHeaderOf(node_address addr) const = 0;

    /** Read the hashed header portion of a node.

          @param  addr    Address of the node we care about
    */
    virtual const void* getHashedHeaderOf(node_address addr) const = 0;


    /**
        Get next pointer for node at this address.
        Used by unique table in case of chaining.
          @param  addr    Address of node
          @return         Handle of next node, or 0.
    */
    virtual node_handle getNextOf(node_address addr) const = 0;


    /**
        Set next pointer for node at this address.
        Used by unique table in case of chaining.
          @param  addr    Address of node
          @param  n       Non-negative.  Handle of
                          next node in a unique table chain,
                          or 0 if none.
    */
    virtual void setNextOf(node_address addr, node_handle n) = 0;


    /**
        Return the name of the style that created us.
    */
    const char* getStyleName() const;

  protected:
    /// TBD - remove this
    /// Dump information not related to individual nodes.
    virtual void dumpInternalInfo(output &s) const = 0;

    /** Get the first interesting address.
        If this cannot be determined, return 0.
        TBD - remove this
    */
    virtual node_address firstNodeAddress() const = 0;

    /** Dump the node/hole information at the given address.
          @param  s       Output stream to use
          @param  addr    Address
          @param  flags   What chunks should be displayed

          @return   Next interesting address, if we can determine this; otherwise 0.

        TBD - make this public
    */
    virtual node_address dumpInternalNode(output &s, node_address addr,
        unsigned flags) const = 0;

    /// TBD - remove this
    /// Dump final info (after node info)
    virtual void dumpInternalTail(output &s) const = 0;

    // Hooks from other classes, so we don't need to make
    // all the derived classes "friends".

    void moveNodeOffset(node_handle node, node_address old_addr, 
        node_address new_addr);

    //
    // Methods for derived classes to deal with
    // members owned by the base class
    //

    const expert_forest* getParent() const;
    expert_forest* getParent();

  private:
    /// Name of the style that invoked us
    const char* style_name;

    /// Parent forest.
    expert_forest* parent;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      expert_forest  class                      *
// *                                                                *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_forest: public forest
{
    // flags for reporting; DO NOT rely on specific values
  public:
    /// Should memory be reported in a human readable format
    static const unsigned HUMAN_READABLE_MEMORY;
    /// Basic forest stats
    static const unsigned BASIC_STATS;
    static const unsigned EXTRA_STATS;
    /// Specific forest stats, dependent on forest type
    static const unsigned FOREST_STATS;
    /// Stats specific to the node storage mechanism.
    static const unsigned STORAGE_STATS;
    /// Detailed stats for the node storage mechanism.
    static const unsigned STORAGE_DETAILED;
    /// Stats specific to the unique table.
    static const unsigned UNIQUE_TABLE_STATS;
    /// Stats specific to the unique table.
    static const unsigned UNIQUE_TABLE_DETAILED;
    /// Stats specific to the hole manager.
    static const unsigned HOLE_MANAGER_STATS;
    /// Stats specific to the hole manager.
    static const unsigned HOLE_MANAGER_DETAILED;

    // ************************************************************
    // *                                                          *
    // *               Constants for  showing nodes               *
    // *                                                          *
    // ************************************************************

    /// Display deleted nodes
    static const unsigned int SHOW_DELETED;
    /// Display zombie nodes
    static const unsigned int SHOW_ZOMBIE;
    /// Display node details
    static const unsigned int SHOW_DETAILS;
    /// Display the node index
    static const unsigned int SHOW_INDEX;
    /// Display terminal nodes
    static const unsigned int SHOW_TERMINALS;

    // ************************************************************
    // *                                                          *
    // *    Preferred way to encode and decode terminal values    *
    // *    (classes so we can use them in template functions)    *
    // *                                                          *
    // ************************************************************

    /** Encoding for booleans into (terminal) node handles */
    class bool_Tencoder {
        // -1 true
        //  0 false
      public:
        static node_handle value2handle(bool v);
        static bool handle2value(node_handle h);
        static void show(output &s, node_handle h);
        static void write(output &s, node_handle h);
        static node_handle read(input &s);
    };
    /** Encoding for integers into (terminal) node handles */
    class int_Tencoder {
      public:
        static node_handle value2handle(int v);
        static int handle2value(node_handle h);
        static void show(output &s, node_handle h);
        static void write(output &s, node_handle h);
        static node_handle read(input &s);
    };
    /** Encoding for floats into (terminal) node handles */
    class float_Tencoder {
        union intfloat {
            float real;
            int integer;
        };
      public:
        static node_handle value2handle(float v);
        static float handle2value(node_handle h);
        static void show(output &s, node_handle h);
        static void write(output &s, node_handle h);
        static node_handle read(input &s);
    };
    // preferred way to encode and decode edge values
    // (classes so we can use them in template functions)
    template<typename T>
    class EVencoder {
      public:
        static size_t edgeBytes();
        static void writeValue(void* ptr, T val);
        static void readValue(const void* ptr, T &val);
        static void show(output &s, const void* ptr);
        static void write(output &s, const void* ptr);
        static void read(input &s, void* ptr);
    };


    friend class reordering_base;

    /** Constructor.
      @param  dslot   slot used to store the forest, in the domain
      @param  d       domain to which this forest belongs to.
      @param  rel     does this forest represent a relation.
      @param  t       the range of the functions represented in this forest.
      @param  ev      edge annotation.
      @param  p       Polcies for reduction, storage, deletion.
      @param  level_reduction_rule       Rules for reduction on different levels.
    */
    expert_forest(int dslot, domain *d, bool rel, range_type t,
                  edge_labeling ev, const policies &p, int* level_reduction_rule);

  // ------------------------------------------------------------
  // inlined helpers.

    /**
        Convenience function.
        Based on the forest type, convert the desired value
        into a terminal node handle.
          @param  v   Value to encode
          @return     Handle for terminal node
    */
    template <typename T>
    node_handle handleForValue(T v) const;

    /**
        Convenience function.
        Based on the forest type, convert the terminal node handle
        into its encoded value.
          @param  n   Node handle
          @param  v   Output: encoded value
    */
    template <typename T>
    void getValueFromHandle(node_handle n, T& v) const;

    /**
        Convenience function.
        Based on the forest type, convert the terminal node handle
        into its encoded boolean value.
          @param  n   Node handle
    */
    bool getBooleanFromHandle(node_handle n) const;

    /**
        Convenience function.
        Based on the forest type, convert the terminal node handle
        into its encoded integer value.
          @param  n   Node handle
    */
    int getIntegerFromHandle(node_handle n) const;

    /**
        Convenience function.
        Based on the forest type, convert the terminal node handle
        into its encoded real (float) value.
          @param  n   Node handle
    */
    float getRealFromHandle(node_handle n) const;

    memstats& changeMemStats();
    /// Number of bytes for an edge value.
    char edgeBytes() const;
    /// Are edge values included when computing the hash.
    bool areEdgeValuesHashed() const;
    /// Extra bytes per node, not hashed.
    char unhashedHeaderBytes() const;
    /// Extra bytes per node, hashed.
    char hashedHeaderBytes() const;

    const expert_domain* getExpertDomain() const;
    expert_domain* useExpertDomain();

  // --------------------------------------------------
  // Node address information
  // --------------------------------------------------
  protected:
    node_address getNodeAddress(node_handle p) const;
    void setNodeAddress(node_handle p, node_address a);

  // --------------------------------------------------
  // Node level information
  // --------------------------------------------------
  public:
    /**
        Negative values are used for primed levels or variables.
    */
    int getVarByLevel(int level) const {
      return level > 0 ? var_order->getVarByLevel(level) : -var_order->getVarByLevel(-level);
    }
    int getLevelByVar(int var) const {
      return var > 0 ? var_order->getLevelByVar(var) : -var_order->getLevelByVar(-var);
    }

    int getNodeLevel(node_handle p) const;
    bool isPrimedNode(node_handle p) const;
    bool isUnprimedNode(node_handle p) const;
    int getNumVariables() const;
    // returns 0 or -K 
    int getMinLevelIndex() const;
    bool isValidLevel(int k) const;
    bool isExtensibleLevel(int k) const;

    /// The maximum size (number of indices) a node at this level can have
    int getLevelSize(int lh) const;
    // The maximum size (number of indices) a variable can have.
    int getVariableSize(int var) const;

  protected:
    void setNodeLevel(node_handle p, int level);

  // --------------------------------------------------
  // Managing incoming edge counts
  // --------------------------------------------------
  public:
    /// Returns true if we are tracking incoming counts
    bool trackingInCounts() const;

    /// Returns the in-count for a node.
    unsigned long getNodeInCount(node_handle p) const;

    /** Increase the link count to this node. Call this when another node is
        made to point to this node.
          @return p, for convenience.
    */
    node_handle linkNode(node_handle p);
    void getDownPtr(node_handle p, int index, long& ev, node_handle& dn) const;

    /** Decrease the link count to this node. If link count reduces to 0, this
        node may get marked for deletion. Call this when another node releases
        its connection to this node.
    */
    void unlinkNode(node_handle p);

  // --------------------------------------------------
  // Managing cache counts
  // --------------------------------------------------
  public:
    /// Returns true if we are tracking incoming counts
    bool trackingCacheCounts() const;

    /** Increase the cache count for this node. Call this whenever this node
        is added to a cache.
          @param  p     Node we care about.
          @return p, for convenience.
    */
    void cacheNode(node_handle p);

    /** Increase the cache count for this node. Call this whenever this node
        is added to a cache.
          @param  p     Node we care about.
    */
    void uncacheNode(node_handle p);

  // --------------------------------------------------
  // Node status
  // --------------------------------------------------
  public:
    bool isActiveNode(node_handle p) const;
    bool isZombieNode(node_handle p) const;
    bool isDeletedNode(node_handle p) const;
    static bool isTerminalNode(node_handle p);
    /// Sanity check: is this a valid nonterminal node index.
    bool isValidNonterminalIndex(node_handle p) const;
    /// Sanity check: is this a valid node index.
    bool isValidNodeIndex(node_handle p) const;
    node_handle getLastNode() const;

    // --------------------------------------------------
    // Extensible Node Information:
    // --------------------------------------------------
    /// Is this an extensible node
    bool isExtensible(node_handle p) const;
    /// If \a p is an extensible node, this returns an index
    /// that can be used by getDownPtr(...) to access the
    /// extensible portion of this node.
    int getExtensibleIndex(node_handle p) const;

  public:

    /// Get the cardinality of an Index Set.
    int getIndexSetCardinality(node_handle node) const;

    // --------------------------------------------------
    // Used by the unique table
    // --------------------------------------------------
    node_handle getNext(node_handle p) const;
    bool isImplicit(node_handle p) const;
    void setNext(node_handle p, node_handle n);
    unsigned hash(node_handle p) const;


    /// A node can be discarded once it goes stale. Whether a node is
    /// considered stale depends on the forest's deletion policy.
    /// Optimistic deletion: A node is said to be stale only when both the
    ///   in-count and cache-count are zero.
    /// Pessimistic deletion: A node is said to be stale when the in-count
    ///  is zero regardless of the cache-count.
    MEDDLY::forest::node_status getNodeStatus(node_handle node) const;

  // ------------------------------------------------------------
  // non-virtual, handy methods for debugging or logging.

    /**
        Display all nodes in the forest.
          @param  s       File stream to write to
          @param  flags   Switches to control output;
                          see constants "SHOW_DETAILED", etc.
    */
    void dump(output &s, unsigned int flags) const;
    void dumpInternal(output &s) const;
    void dumpUniqueTable(output &s) const;
    void validateIncounts(bool exact);
    void validateCacheCounts() const;
    void countNodesByLevel(long* active) const;



  // ------------------------------------------------------------
  // non-virtual, handy methods.

    /** Build a list of nodes in the subgraph below the given node.
        This for example is used to determine which nodes must
        be printed to display a subgraph.
        Terminal nodes are NOT included.

          @param  roots   Array of root nodes in the forest.
                          Each root node will be included in the list,
                          except for terminal nodes.

          @param  N       Dimension of \a roots array.

          @param  sort    If true, the list will be in increasing order.
                          Otherwise, the list will be in some convenient order
                          (currently, it is the order that nodes
                          are discovered).

          @return   A malloc'd array of non-terminal nodes, terminated by 0.
                    Or, a null pointer, if the list is empty.
    */
    node_handle*
    markNodesInSubgraph(const node_handle* roots, int N, bool sort) const;

    /** Count and return the number of non-terminal nodes
        in the subgraph below the given node.
    */
    long getNodeCount(node_handle node) const;

    /** Count and return the number of non-terminal nodes
        in the subgraph below the given nodes.
    */
    long getNodeCount(const node_handle* roots, int N) const;

    /** Count and return the number of edges
        in the subgraph below the given node.
    */
    long getEdgeCount(node_handle node, bool countZeroes) const;

    /** Display the contents of a single node.
          @param  s       File stream to write to.
          @param  node    Node to display.
          @param  flags   Switches to control output;
                          see constants "SHOW_DETAILED", etc.

          @return true, iff we displayed anything
    */
    bool showNode(output &s, node_handle node, unsigned int flags = 0) const;

    /// Show all the nodes in the subgraph below the given nodes.
    void showNodeGraph(output &s, const node_handle* node, int n) const;

    /// Write all the nodes in the subgraph below the given nodes
    /// in a graphical format specified by the extension.
    void writeNodeGraphPicture(const char* filename, const char *extension,
        const node_handle* node, int n) const;


    /** Show various stats for this forest.
          @param  s       Output stream to write to
          @param  pad     Padding string, written at the start of
                          each output line.
          @param  flags   Which stats to display, as "flags";
                          use bitwise or to combine values.
                          For example, BASIC_STATS | FOREST_STATS.
    */
    void reportStats(output &s, const char* pad, unsigned flags) const;


    /// Compute a hash for a node.
    unsigned hashNode(node_handle p) const;

    /** Check and find the index of a single downward pointer.

          @param  node    Node we care about
          @param  down    Output:
                          The singleton downward pointer, or undefined.

          @return   If the node has only one non-zero downward pointer,
                    then return the index for that pointer.
                    Otherwise, return a negative value.
    */
    int getSingletonIndex(node_handle p, node_handle &down) const;

    /** Check and get a single downward pointer.

          @param  node    Node we care about
          @param  index   Index we're trying to match

          @return   If the only non-zero downward pointer for
                    this node happens at \a index, then return the pointer.
                    Otherwise, return 0.
    */
    node_handle getSingletonDown(node_handle node, int index) const;

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        unpacked_node should be used instead.

          @param  p       Node to look at
          @param  index   Index of the pointer we want.

          @return         The downward pointer at that index.
    */
    node_handle getDownPtr(node_handle p, int index) const;

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        unpacked_node should be used instead.

          @param  p       Node to look at
          @param  index   Index of the pointer we want.

          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    void getDownPtr(node_handle p, int index, int& ev, node_handle& dn) const;

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        unpacked_node should be used instead.

          @param  p       Node to look at
          @param  index   Index of the pointer we want.

          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    void getDownPtr(node_handle p, int index, float& ev, node_handle& dn) const;

    /**
        Get the transparent node.
        This is the default node value for "skipped" edges in sparse nodes.
    */
    node_handle getTransparentNode() const;

  // ------------------------------------------------------------
  // Copy a node into an unpacked node

  void fillUnpacked(unpacked_node &un, node_handle node, unpacked_node::storage_style st2) const;

  /**   Return a forest node equal to the one given.
        The node is constructed as necessary.
        This version should be used only for
        multi terminal forests.
        The unpacked node is recycled.
          @param  in    Incoming pointer index;
                        used for identity reductions.
          @param  un    Temporary node; will be recycled.

          @return       A node handle equivalent
                        to \a un, taking into account
                        the forest reduction rules
                        and if a duplicate node exists.
    
  */
  node_handle createReducedNode(int in, unpacked_node *un);
  
  /** Return a forest node for implicit node equal to the one given.
   The implicit node is already constructed inside satimpl_opname.
   This version should be used only for
   multi terminal forests.
   @param  un    Implicit relation node
   
   @return       A node handle equivalent
   to \a un.
   
   */
  node_handle createRelationNode(MEDDLY::relation_node *un);

  /** Return a forest node equal to the one given.
      The node is constructed as necessary.
      This version should be used only for
      edge valuded forests.
      The node builder nb is recycled.
        @param  in    Incoming pointer index;
                      used for identity reductions.
        @param  nb    Constructed node.
        @param  ev    Output: edge value
        @param  node  Output: node handle.
                      On exit, the edge value and the node
                      handle together are equivalent to nb;
                      taking into account the forest reduction rules
                      and if a duplicate node exists.
  */
  template <class T>
  void createReducedNode(int in, unpacked_node* nb, T& ev, node_handle& node);


    /** Swap the content of nodes.
        Do not update their parents and inCount.
    */
    void swapNodes(node_handle p, node_handle q);

    /*
     * Modify a node in place.
     * Does not check if the modified node is duplicate or redundant.
     * The level of the node may change.
     * Keep the reference number and the cache count of the node.
     */
    node_handle modifyReducedNodeInPlace(unpacked_node* un, node_handle p);

  // ------------------------------------------------------------
  // virtual in the base class, but implemented here.
  // See meddly.h for descriptions of these methods.

    virtual void writeEdges(output &s, const dd_edge* E, int n) const;
    virtual void readEdges(input &s, dd_edge* E, int n);
    virtual void garbageCollect();
    virtual void compactMemory();
    virtual void showInfo(output &strm, int verbosity);

  // ------------------------------------------------------------
  // abstract virtual, must be overridden.
  //

    /**
        Is the given edge transparent?
        If so it may be "skipped" in a sparse node.
        Should not be called for multi-terminal; it is better to simply
        compare with the transparent node as given by getTransparentNode().
          @param  ep    Node part of the edge to check
          @param  ev    Value part of the edge to check
    */
    virtual bool isTransparentEdge(node_handle ep, const void* ev) const = 0;

    /**
        Get the transparent edge value.
        This is the default edge value for "skipped" edges in sparse nodes.
        Copy the transparent edge value into the parameters.
          @param  ep      Put node part of the edge here.
          @param  ev      Put value part of the edge here.
    */
    virtual void getTransparentEdge(node_handle &ep, void* ptr) const = 0;

    /** Are the given edge values "duplicates".
        I.e., when determining if two nodes are duplicates,
        do the given edge values qualify as duplicate values.
          @param  eva     Pointer to the first edge value
          @param  evb     Pointer to the second edge value

          @return     true, iff the edge values are "equal".
          @throws     A TYPE_MISMATCH error if the forest
                      does not store edge values.
    */
    virtual bool areEdgeValuesEqual(const void* eva, const void* evb) const = 0;

    virtual void normalize(unpacked_node &nb, long& ev) const;

    /** Discover duplicate nodes.
        Right now, used for sanity checks only.
          @param  node    Handle to a node.
          @param  nr      Some other node.

          @return   true, iff the nodes are duplicates.
    */
    bool areDuplicates(node_handle node, const unpacked_node &nr) const;

    /** Is this a redundant node that can be eliminated?
        Must be implemented in derived forests
        to deal with the default edge value.
          @param  nb    Node we're trying to build.

          @return   True, if nr is a redundant node
                          AND it should be eliminated.
    */
    virtual bool isRedundant(const unpacked_node &nb) const = 0;

    /** Is the specified edge an identity edge, that can be eliminated?
        Must be implemented in derived forests
        to deal with the default edge value.

          @param  nr    Node we're trying to build.
                        We know there is a single non-zero downward pointer.

          @param  i     Candidate edge (or edge index for sparse nodes).

          @return True, if nr[i] is an identity edge, and the
                        identity node should be eliminated.
    */
    virtual bool isIdentityEdge(const unpacked_node &nb, int i) const = 0;


    /**
        Build an iterator.
        Used by class enumerator.
    */
    virtual enumerator::iterator* makeFullIter() const = 0;

    /**
        Build an iterator with a fixed row.
        Default behavior - throw an "INVALID_FOREST" error.
    */
    virtual enumerator::iterator* makeFixedRowIter() const;

    /**
        Build an iterator with a fixed column.
        Default behavior - throw an "INVALID_FOREST" error.
    */
    virtual enumerator::iterator* makeFixedColumnIter() const;

    /*
     * Reorganize the variables in a certain order.
     */
    void reorderVariables(const int* level2var);

    void getVariableOrder(int* level2var) const;

    std::shared_ptr<const variable_order> variableOrder() const;

    /*
     * Swap the variables at level and level+1.
     * This method should only be called by expert_domain.
     */
    virtual void swapAdjacentVariables(int level) = 0;

    /*
     * Move the variable at level high down to level low.
     * The variables from level low to level high-1 will be moved one level up.
     */
    virtual void moveDownVariable(int high, int low) = 0;

    /*
     * Move the variable at level low up to level high.
     * The variables from level low+1 to level high will be moved one level down.
     */
    virtual void moveUpVariable(int low, int high) = 0;

    virtual void dynamicReorderVariables(int top, int bottom) {
    	throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    }

    /** Show a terminal node.
          @param  s       Stream to write to.
          @param  tnode   Handle to a terminal node.
    */
    virtual void showTerminal(output &s, node_handle tnode) const;

    /** Write a terminal node in machine-readable format.
          @param  s       Stream to write to.
          @param  tnode   Handle to a terminal node.
    */
    virtual void writeTerminal(output &s, node_handle tnode) const;

    /** Read a terminal node in machine-readable format.
          @param  s       Stream to read from.
          @return         Handle to a terminal node.
          @throws         An invalid file exception if the stream does not
                          contain a valid terminal node.
    */
    virtual node_handle readTerminal(input &s);

    /** Show an edge value.
          @param  s       Stream to write to.
          @param  edge    Pointer to edge value chunk
    */
    virtual void showEdgeValue(output &s, const void* edge) const;

    /** Write an edge value in machine-readable format.
          @param  s       Stream to write to.
          @param  edge    Pointer to edge value chunk
    */
    virtual void writeEdgeValue(output &s, const void* edge) const;

    /** Read an edge value in machine-readable format.
          @param  s       Stream to read from.
          @param  edge    Pointer to edge value chunk
    */
    virtual void readEdgeValue(input &s, void* edge);

    /** Show the hashed header values.
          @param  s       Stream to write to.
          @param  hh      Pointer to hashed header data.
    */
    virtual void showHashedHeader(output &s, const void* hh) const;

    /** Write the hashed header in machine-readable format.
          @param  s       Stream to write to.
          @param  hh      Pointer to hashed header data.
    */
    virtual void writeHashedHeader(output &s, const void* hh) const;

    /** Read the hashed header in machine-readable format.
          @param  s       Stream to write to.
          @param  nb      Node we're building.
    */
    virtual void readHashedHeader(input &s, unpacked_node &nb) const;

    /** Show the unhashed header values.
          @param  s       Stream to write to.
          @param  uh      Array of all unhashed header values.
    */
    virtual void showUnhashedHeader(output &s, const void* uh) const;

    /** Write the unhashed header in machine-readable format.
          @param  s       Stream to write to.
          @param  hh      Pointer to unhashed header data.
    */
    virtual void writeUnhashedHeader(output &s, const void* uh) const;

    /** Read the unhashed header in machine-readable format.
          @param  s       Stream to write to.
          @param  nb      Node we're building.
    */
    virtual void readUnhashedHeader(input &s, unpacked_node &nb) const;



  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                       protected  methods                       |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:

    /// Destructor.
    virtual ~expert_forest();

    /** Initialize data.
        Should be called in the child class constructors.
        Allows us to use class properties to initialize the data.
    */
    void initializeForest();


  // ------------------------------------------------------------
  // inlined setters for derived classes to use.

    void setEdgeSize(char ebytes, bool hashed);
    void setUnhashedSize(char ubytes);
    void setHashedSize(char hbytes);

  // ------------------------------------------------------------
  // virtual, with default implementation.
  // Should be overridden in appropriate derived classes.

    /// Character sequence used when writing forests to files.
    virtual const char* codeChars() const;


    /** Normalize a node.
        Used only for "edge valued" DDs with range type: integer.
        Different forest types will have different normalization rules,
        so the default behavior given here (throw an error) will need
        to be overridden by all edge-valued forests.

          @param  nb    Array of downward pointers and edge values;
                        may be modified.
          @param  ev    The incoming edge value, may be modified
                        as appropriate to normalize the node.
    */
    virtual void normalize(unpacked_node &nb, int& ev) const;

    /** Normalize a node.
        Used only for "edge valued" DDs with range type: real.
        Different forest types will have different normalization rules,
        so the default behavior given here (throw an error) will need
        to be overridden by all edge-valued forests.

          @param  nb    Array of downward pointers and edge values;
                        may be modified.
          @param  ev    The incoming edge value, may be modified
                        as appropriate to normalize the node.
    */
    virtual void normalize(unpacked_node &nb, float& ev) const;

    /** Show forest-specific stats.
          @param  s     Output stream to use
          @param  pad   String to display at the beginning of each line.
    */
    virtual void reportForestStats(output &s, const char* pad) const;



  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                        private  methods                        |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:

  // ------------------------------------------------------------
  // inlined helpers for this class

    bool isTimeToGc() const;

    /** Change the location of a node.
        Used by node_storage during compaction.
        Should not be called by anything else.
          @param  node        Node we're moving
          @param  old_addr    Current address of node, for sanity check
          @param  new_addr    Where we're moving the node to
    */
    void moveNodeOffset(node_handle node, node_address old_addr, node_address new_addr);
    friend class MEDDLY::node_storage;
    friend class MEDDLY::node_headers;  // calls deleteNode().

    friend class MEDDLY::global_rebuilder;

  // ------------------------------------------------------------
  // helpers for this class

    /**
        Disconnects all downward pointers from p,
        and removes p from the unique table.
    */
    void deleteNode(node_handle p);

    /** Apply reduction rule to the temporary node and finalize it. 
        Once a node is reduced, its contents cannot be modified.
          @param  in    Incoming index, used only for identity reduction;
                        Or -1.
          @param  un    Unpacked node.
          @return       Handle to a node that encodes the same thing.
    */
    node_handle createReducedHelper(int in, unpacked_node &nb);
  
    /** Create implicit node in the forest. Just add a handle and it points to original location
     @param  un    Relation node.
     @return       Handle to a node that encodes the same thing.
     */
    node_handle createImplicitNode(MEDDLY::relation_node &nb);


    /** Apply reduction rule to the temporary extensible node and finalize it. 
        Once a node is reduced, its contents cannot be modified.
          @param  in    Incoming index, used only for identity reduction;
                        Or -1.
          @param  un    Unpacked extensible node. Must be sorted by indices if sparse.
          @return       Handle to a node that encodes the same thing.
    */
    node_handle createReducedExtensibleNodeHelper(int in, unpacked_node &nb);

    // Sanity check; used in development code.
    void validateDownPointers(const unpacked_node &nb) const;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                              Data                              |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  protected:
    /// uniqueness table, still used by derived classes.
    unique_table* unique;

    /// Should a terminal node be considered a stale entry in the compute table.
    /// per-forest policy, derived classes may change as appropriate.
    MEDDLY::forest::node_status terminalNodesStatus;

    /// Class that stores nodes.
    node_storage* nodeMan;

    /// Transparent value
    node_handle transparent;

    std::shared_ptr<const variable_order> var_order;

  private:
    // Garbage collection in progress
    bool performing_gc;

    // memory for validating incounts
    node_handle* in_validate;
    int  in_val_size;
    // depth of delete/zombie stack; validate when 0
    int delete_depth;

    /// Node header information
    node_headers nodeHeaders;


    /// Number of bytes for an edge
    char edge_bytes;
    /// Are edge values hashed?
    bool hash_edge_values;
    /// Number of bytes of unhashed header
    char unhashed_bytes;
    /// Number of bytes of hashed header
    char hashed_bytes;

    class nodecounter;
};
// end of expert_forest class.


// ******************************************************************
// *                                                                *
// *                expert_forest::nodecounter class                *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_forest::nodecounter: public edge_visitor {
    expert_forest* parent;
    int* counts;
  public:
    nodecounter(expert_forest*p, int* c);
    virtual ~nodecounter();
    virtual void visit(dd_edge &e);
};


// ******************************************************************
// *                                                                *
// *                          opname class                          *
// *                                                                *
// ******************************************************************

/// Class for names of operations.
class MEDDLY::opname {
    const char* name;
    int index;
    static int next_index;

    friend void MEDDLY::initialize(initializer_list *);
    friend void MEDDLY::cleanup();
  public:
    opname(const char* n);
    virtual ~opname();

    int getIndex() const;
    const char* getName() const;
};

// ******************************************************************
// *                                                                *
// *                       unary_opname class                       *
// *                                                                *
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
// *                                                                *
// *                      binary_opname  class                      *
// *                                                                *
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
// *                                                                *
// *                    specialized_opname class                    *
// *                                                                *
// ******************************************************************

/// Specialized operation names.
class MEDDLY::specialized_opname : public opname {
  public:
    /**
      Abstract base class for arguments to buildOperation().
      Derived operation names must provide derived classes for arguments.
    */
    class arguments {
      public:
        arguments();
        virtual ~arguments();

        /**
            Specify if arguments should be destroyed or not.
            If yes (the default), the operation will destroy
            the arguments once they are no longer needed.
        */
        void setAutoDestroy(bool destroy);
        bool autoDestroy() const;

      private:
        bool destroyWhenDone;
    };

    specialized_opname(const char* n);
    virtual ~specialized_opname();

    /** Note - unlike the more general binary and unary ops,
        a specialized operation might be crafted to the specific
        arguments (passed as an abstract class).
        Examples are:
          - operations that will be called several times with
            several of the arguments unchanged, and we want
            to do some preprocessing on the unchanged arguments.

          - operations with very bizarre and/or user-defined
            parameters (like on-the-fly saturation).

        @param  a   Arguments.  Will be destroyed when we are finished,
                    if autoDestroy() is set for the arguments.
    */
    virtual specialized_operation* buildOperation(arguments* a) const = 0;
};

// ******************************************************************
// *                                                                *
// *                     numerical_opname class                     *
// *                                                                *
// ******************************************************************

/// Numerical operation names.
class MEDDLY::numerical_opname : public specialized_opname {
  public:
    class numerical_args : public specialized_opname::arguments {
      public:
        const dd_edge &x_ind;
        const dd_edge &A;
        const dd_edge &y_ind;

        numerical_args(const dd_edge &xi, const dd_edge &a, const dd_edge &yi);
        virtual ~numerical_args();
    };

    numerical_opname(const char* n);
    virtual ~numerical_opname();
    virtual specialized_operation* buildOperation(arguments* a) const = 0;

    /// For convenience, and backward compatability :^)
    specialized_operation* buildOperation(const dd_edge &x_ind,
      const dd_edge &A, const dd_edge &y_ind) const;
};


// ******************************************************************
// *                                                                *
// *                     satpregen_opname class                     *
// *                                                                *
// ******************************************************************

/// Saturation, with already generated transition relations, operation names.
class MEDDLY::satpregen_opname : public specialized_opname {
  public:
    satpregen_opname(const char* n);
    virtual ~satpregen_opname();

    /// Arguments should have type "pregen_relation".
    virtual specialized_operation* buildOperation(arguments* a) const = 0;

    /** Class for a partitioned transition relation, already known
        The relation can be partitioned "by events" or "by levels".
        In the case of "by events", we can have more than one relation
        per level; otherwise, there is at most one relation per level.
    */
    class pregen_relation : public specialized_opname::arguments {
      public:
        /** Constructor, by events
              @param  inmdd       MDD forest containing initial states
              @param  mxd         MxD forest containing relations
              @param  outmdd      MDD forest containing result
              @param  num_events  Number of events; specifies the maximum
                                  number of calls to addToRelation().
        */
        pregen_relation(forest* inmdd, forest* mxd, forest* outmdd,
          int num_events);

        /** Constructor, by levels
              @param  inmdd       MDD forest containing initial states
              @param  mxd         MxD forest containing relations
              @param  outmdd      MDD forest containing result
        */
        pregen_relation(forest* inmdd, forest* mxd, forest* outmdd);

        virtual ~pregen_relation();
        void addToRelation(const dd_edge &r);

        // Options for controlling the amount of processing performed by
        // \a finalize(splittingOption).
        enum splittingOption {
          // None.
          None,
          // Transitions from level K that do not effect level K,
          // are moved to a lower level.
          SplitOnly,
          // SplitOnly + duplicated transitions between adjacent levels
          // are removed from the higher level.
          SplitSubtract,
          // SplitOnly + all duplicate transitions are removed.
          SplitSubtractAll,
          // Same as SplitSubtractAll, but using an algorithm that
          // first combines all transitions before splitting it up per level.
          MonolithicSplit
        };

        /** To be called after all events have been added to
            the transition relation.
            This method modifies the decision diagrams stored at different
            levels, to reduce duplicated transitions.
              @param  split       This parameter only applies to "by levels",
                                  and it controls the amount of processing
                                  that is performed.
                                  Please refer to splittingOption for details.
        */
        void finalize(splittingOption split = SplitSubtract);

        bool isFinalized() const;

        forest* getInForest() const;
        forest* getRelForest() const;
        forest* getOutForest() const;

        // the following methods assume the relation has been finalized.
        node_handle* arrayForLevel(int k) const;
        int lengthForLevel(int k) const;

      private:
        // helper for constructors
        void setForests(forest* inf, forest* mxd, forest* outf);
        // helper for finalize,
        // find intersection of diagonals of events[k],
        // subtracts the intersection of events[k] and adds it to events[k-1].
        void splitMxd(splittingOption split);
        // helper for finalize
        // adds all event[k]; sets all event[k] to 0;
        // sets events[level(sum)] = sum
        void unionLevels();

        forest* insetF;
        expert_forest* mxdF;
        forest* outsetF;
        int K;
        // array of sub-relations
        node_handle* events;
        // next pointers, unless we're finalized
        int* next;
        // size of events array
        int num_events;
        // last used element of events array
        int last_event;

        // If null, then we are "by levels".  Otherwise, we are "by events",
        // and before we're finalized, level_index[k] points to a linked-list
        // of sub-relations that affect level k.
        // after we're finalized, the events array is sorted, so
        // level_index[k] is the index of the first event affecting level k.
        // Dimension is number of variables + 1.
        int* level_index;
    };

};



// ******************************************************************
// *                                                                *
// *                      satotf_opname  class                      *
// *                                                                *
// ******************************************************************

/// Saturation, transition relations built on the fly, operation names.
class MEDDLY::satotf_opname : public specialized_opname {
  public:
    satotf_opname(const char* n);
    virtual ~satotf_opname();

    /// Arguments should have type "otf_relation", below
    virtual specialized_operation* buildOperation(arguments* a) const = 0;

    class otf_relation;

    // ============================================================

    /**
        User must derive a subclass from this.
        Part of an enabling or updating function.
        It knows what variables it depends on, and how to build itself
        (provided by the user).
    */
    class subevent {
      public:
        /// Constructor, specify variables that this function depends on,
        /// and if it is a firing or enabling event.
        subevent(forest* f, int* v, int nv, bool firing);
        virtual ~subevent();

        /// Get the forest to which this function belongs to.
        expert_forest* getForest();

        /// Get number of variables this function depends on.
        int getNumVars() const;

        /// Get array of variables this function depends on.
        const int* getVars() const;

        /// Get the DD encoding of this function
        const dd_edge& getRoot() const;

        /// Get the "top" variable for this function
        int getTop() const;

        /// Is this a firing subevent?
        bool isFiring() const;

        /// Is this an enabling subevent
        bool isEnabling() const;

        /**
          Rebuild the function to include the
          local state "index" for the variable "v".
          Updates root with the updated function.
          User MUST provide this method.
        */
        virtual void confirm(otf_relation &rel, int v, int index) = 0;

        /// If num_minterms > 0,
        ///   Add all minterms to the root
        ///   Delete all minterms.
        void buildRoot();

        /// Debugging info
        void showInfo(output& out) const;

        long mintermMemoryUsage() const;
        void clearMinterms();

      protected:
        bool addMinterm(const int* from, const int* to);
        bool usesExtensibleVariables() const;

        int* vars;
        int num_vars;
        dd_edge root;
        int top;
        expert_forest* f;
        int** unpminterms;
        int** pminterms;
        int num_minterms;
        int size_minterms;
        bool is_firing;
        bool uses_extensible_variables;

    };  // end of class subevent

    // ============================================================

    /**
        An "event".
        Produces part of the transition relation, from its sub-functions.

        TBD - do we need to split the enabling and updating sub-functions,
        or will one giant list work fine?
    */
    class event {
        // TBD - put a list of events that have priority over this one

        // TBD - for priority - when is this event enabled?
      public:
        event(subevent** se, int nse);
        virtual ~event();

        /// Get the forest to which the subevents belong to
        inline expert_forest* getForest() { return f; } 

        /// Get number of subevents
        inline int getNumOfSubevents() const { return num_subevents; }

        /// Get array of subevents
        inline subevent** getSubevents() const { return subevents; }

        /// Get the "top" variable for this event
        inline int getTop() const { return top; }

        /// Get the number of variables that are effected by this event
        inline int getNumVars() const { return num_vars; }

        /// Get a (sorted) array of variables that are effected by this event
        inline const int* getVars() const { return vars; }

        inline const dd_edge& getRoot() const { return root; }

        inline bool isDisabled() const { return is_disabled; }

        inline bool needsRebuilding() const { return needs_rebuilding; }

        inline void markForRebuilding() { needs_rebuilding = true; }

        /**
            If this event has been marked for rebuilding:
              Build this event as a conjunction of its sub-events.

            @return               true, if the event needed rebuilding and
                                        the rebuilding modified the root.
                                  false, otherwise.
        */
        virtual bool rebuild();

        /// Enlarges the "from" variable to be the same size as the "to" variable
        void enlargeVariables();

        /// Debugging info
        void showInfo(output& out) const;

        long mintermMemoryUsage() const;

      protected:
        void buildEventMask();

      private:
        subevent** subevents;
        int num_subevents;
        int top;
        int num_vars;
        int* vars;
        dd_edge root;
        bool needs_rebuilding;
        expert_forest* f;

        bool is_disabled;
        int num_firing_vars;
        int* firing_vars;
        dd_edge event_mask;
        int* event_mask_from_minterm;
        int* event_mask_to_minterm;

    };  // end of class event

    // ============================================================

    /**
        Overall relation.
        This includes all events, and keeping track of which local
        variables are confirmed.

        TBD.
    */
    class otf_relation : public specialized_opname::arguments {
      public:
        /** Constructor.
              @param  inmdd       MDD forest containing initial states
              @param  mxd         MxD forest containing relations
              @param  outmdd      MDD forest containing result
              @param  E           List of events
              @param  nE          Number of events
        */
        otf_relation(forest* inmdd, forest* mxd, forest* outmdd,
          event** E, int ne);

        virtual ~otf_relation();

        /// Returns the MDD forest that stores the initial set of states
        expert_forest* getInForest() const;

        /// Returns the MXD forest that stores the events
        expert_forest* getRelForest() const;

        /// Returns the MDD forest that stores the resultant set of states
        expert_forest* getOutForest() const;

        /// Returns true if the local state is already confirmed.
        bool isConfirmed(int level, int index) const;

        /// Returns an array of local states for this level, such that
        /// result[i] == isConfirmed(level, i).
        const bool* getLocalStates(int level);

        /// Returns the number of confirmed states at this level
        int getNumConfirmed(int level) const;

        /// Confirms all local states enabled in the given MDD
        void confirm(const dd_edge& set);

        /** Confirm a variable's previously unconfirmed state.
            Any event that is dependent on this variable is marked
            as "stale" --- so that it is rebuilt before use.

            @param  level       variable's level
            @param  index       the state of the variable being confirmed.
            @return             false: if state was previously confirmed.
                                true: if state was previously unconfirmed.
         */
        bool confirm(int level, int index);

        /** Get the number of events at whose "top" is this level.

            @param  level       level for the events.
            @return             number of events whose "top" is this level.
         */
        int getNumOfEvents(int level) const;

        /** Gets an event from the set of events whose "top" is this level.

            @param  level       level for the events.
            @param  i           index of the event.
            @return             if 0 <= i < getNumOfEvents(level),
                                the ith event at this level;
                                otherwise, 0.
         */
        node_handle getEvent(int level, int i);

        /** Rebuild an event.

            @param  i           index of the event.
            @return             true, if event was updated.
          */
        bool rebuildEvent(int level, int i);

        /** Build a Monolithic Next State Function that is equivalent to
            the union of all events while restricting the size of each
            variable to that of the largest confirmed index.

            @return             union of bounded OTF events.
        */
        node_handle getBoundedMonolithicNSF();

        /** Bound all extensible variables
            using the maximum confirmed local state as the bound.
        */
        void bindExtensibleVariables();

        /** Get the number of arcs in the OTF relation
            restricted by the confirmed local states.

            Only works with non-extensible variables (call
            bindExtensibleVariables() prior to calling this).

            @param    count_duplicates  if false, counts arcs that are common
                                          to mutiple transitions as one.
            @return                     the number of arcs in the OTF relation.
        */
        double getArcCount(const dd_edge& mask, bool count_duplicates);

        /// For Debugging
        void showInfo(output &strm) const;

        long mintermMemoryUsage() const;

        void clearMinterms();

      protected:
        void enlargeConfirmedArrays(int level, int sz);
        node_handle getBoundedMxd(node_handle mxd, const int* bounds_array, int sz,
            std::unordered_map<node_handle, node_handle>& cache);

      private:
        expert_forest* insetF;
        expert_forest* mxdF;
        expert_forest* outsetF;
        int num_levels;

        // All events that begin at level i,
        // are listed in events_by_top_level[i].
        // An event will appear in only one list
        // (as there is only one top level per event).
        // num_events_by_top_level[i] gives the size of events_by_top_level[i]
        event*** events_by_top_level;
        int *num_events_by_top_level;

        // All events that depend on a level i,
        // are listed in events_by_level[i]
        // Therefore, an event that depends on n levels,
        // will appear in n lists
        // num_events_by_level[i] gives the size of events_by_level[i]
        event*** events_by_level;
        int *num_events_by_level;

        // All subevents that depend on a level i,
        // are listed in subevents_by_level[i]
        // Therefore, an subevent that depends on n levels,
        // will appear in n lists
        // num_subevents_by_level[i] gives the size of subevents_by_level[i]
        subevent*** subevents_by_level;
        int *num_subevents_by_level;

        // List of confirmed local states at each level
        bool** confirmed;
        int* size_confirmed;
        int* num_confirmed;

    };  // end of class otf_relation

};  // end of class satotf_opname

// ******************************************************************
// *                                                                *
// *                      satimpl_opname class                      *
// *                                                                *
// ******************************************************************

/// Saturation, transition relations stored implcitly, operation names.
class MEDDLY::satimpl_opname: public specialized_opname {
  public:

    satimpl_opname(const char* n);
    virtual ~satimpl_opname();

    /// Arguments should have type "implicit_relation", below
    virtual specialized_operation* buildOperation(arguments* a) const;

  public:

    /** An implicit relation, as a DAG of relation_nodes.

        The relation is partitioned by "events", where each event
        is the conjunction of local functions, and each local function
        is specified as a single relation_node.  The relation_nodes
        are chained together with at most one relation_node per state
        variable, and any skipped variables are taken to be unchanged
        by the event.

        If the bottom portion (suffix) of two events are identical,
        then they are merged.  This is done by "registering" nodes
        which assigns a unique ID to each node, not unlike an MDD forest.

        Note: node handles 0 and 1 are reserved.
        0 means null node.
        1 means special bottom-level "terminal" node
        (in case we need to distinguish 0 and 1).
    */
    class implicit_relation : public specialized_opname::arguments {
      public:
        /** Constructor.

            @param  inmdd       MDD forest containing initial states
            @param  outmdd      MDD forest containing result

            Not 100% sure we need these...
        */
        implicit_relation(forest* inmdd, forest* relmxd, forest* outmdd);
        virtual ~implicit_relation();
      
        /// Returns the Relation forest that stores the mix of relation nodes and mxd nodes
        expert_forest* getMixRelForest() const;


        /// Returns the MDD forest that stores the initial set of states
        expert_forest* getInForest() const;


        /// Returns the MDD forest that stores the resultant set of states
        expert_forest* getOutForest() const;

        /// Returns true iff a state in \a constraint is reachable
        /// from the states in \a initial_states
        /// Note: \a constraint can be an XDD
        bool isReachable(const dd_edge& initial_states, const dd_edge& constraint);

        /** Register a relation node.

            If we have seen an equivalent node before, then
            return its handle and destroy n.
            Otherwise, add n to the unique table, assign it a unique
            identifier, and return that identifier.

            @param  is_event_top    If true, this is also the top
                                    node of some event; register it
                                    in the list of events.

            @param  n               The relation node to register.

            @return Unique identifier to use to refer to n.
        */
        rel_node_handle registerNode(bool is_event_top, relation_node* n);

        /** Check if the relation node is unique
            @param n  The relation node.
            @return   If unique, 0
                      Else, existing node handle
        */
        rel_node_handle isUniqueNode(relation_node* n);


        /** Indicate that there will be no more registered nodes.
            Allows us to preprocess the events or cleanup or convert
            to a more useful representation for saturation.
        */

        //void finalizeNodes();

        /** Get the relation node associated with the given handle.

            Should be a fast, inlined implementation.
        */
        relation_node* nodeExists(rel_node_handle n);

        /** Get the relation node associated with the given handle.

            Should be a fast, inlined implementation.
        */
        bool isReserved(rel_node_handle n);

      private:
        expert_forest* insetF;
        expert_forest* outsetF;
        expert_forest* mixRelF;
        
        int num_levels;

      private:

        /// Last used ID of \a relation node.
        long last_in_node_array;

      private:
        // TBD - add a data structure for the "uniqueness table"
        // of relation_nodes, so if we register a node that
        // is already present in a node_array, we can detect it.

        std::unordered_map<rel_node_handle, relation_node*> impl_unique;

      private:
        // TBD - add a data structure for list of events with top level k,
        // for all possible k.
        // Possibly this data structure is built by method
        // finalizeNodes().

        rel_node_handle** event_list;
        long* event_list_alloc; // allocated space
        long* event_added; //how many events added so far

        long* confirm_states; //total no. of confirmed states of a level
        bool** confirmed; // stores whether a particular local state is confirmed
        long* confirmed_array_size; // stores size of confirmed array


      public:

        /// Get total number of events upto given level
        long getTotalEvent(int level);

        /// Resizes the Event List
        void resizeEventArray(int level);

        /// Returns the number of events that have this level as top
        long lengthForLevel(int level) const;

        /// Returns the array of events that have this level as top
        rel_node_handle* arrayForLevel(int level) const;

        /// Returns the number of confirmed states at a level
        long getConfirmedStates(int level) const;

        /// Confirms the local states at a level
        void setConfirmedStates(int level, int i);

        /// Confirms the local states in the given MDD
        void setConfirmedStates(const dd_edge &set);


        /// Checks if i is confirmed
        bool isConfirmedState(int level, int i);

        /// Expand confirm array
        void resizeConfirmedArray(int level, int index);

        /** Bound all extensible variables
            using the maximum confirmed local state as the bound.
        */
        void bindExtensibleVariables();

      public:
        /// Prints the implicit relation
        void show();

        /// Build mxd forest
        MEDDLY::node_handle buildMxdForest();

        /// Build each event_mxd
        dd_edge buildEventMxd(rel_node_handle event_top, forest *mxd);

        /// Get relation forest  
        expert_forest* getRelForest() const;


      private:
        expert_forest* mxdF;
      
     public:
      
      /*
       Group the list of events at a given level by same next-of values 
       @param level    level at which saturation is called
       @param i        number of tokens before event is fired
       @param R        array of events at level top
       
       Return the map.
       */
      std::unordered_map<long,std::vector<rel_node_handle>> getListOfNexts(int level, long i, relation_node **R);
      
      /*
       Returns whether there exist a possibility of doing union 
       @param level    level at which saturation is called
       @param i        number of tokens before event is fired
       @param R        array of events at level top
       
       Return bool.
       */
      bool isUnionPossible(int level, long i, relation_node **R);
      
    };  // class implicit_relation



};
class MEDDLY::constrained_opname : public specialized_opname {
public:
	constrained_opname(const char* n);

	class constrained_args : public specialized_opname::arguments {
	public:
	  forest* consForest;
	  forest* inForest;
	  forest* relForest;
	  forest* outForest;

	  constrained_args(forest* consF, forest* inF, forest* relF, forest* outF)
	    : consForest(consF), inForest(inF), relForest(relF), outForest(outF)
	  {
	  }
	};
};

// ******************************************************************
// *                                                                *
// *                         ct_object class                        *
// *                                                                *
// ******************************************************************

/** Generic objects in compute tables.
    Used for things other than dd_edges and simple types.
    Defined in compute_table.cc
*/
class MEDDLY::ct_object {
  public:
    ct_object();
    virtual ~ct_object();
    virtual opnd_type getType() = 0;
};


// ******************************************************************
// *                                                                *
// *                      ct_initializer  class                     *
// *                                                                *
// ******************************************************************

/** Interface for initializing Meddly's compute table(s).
    Implemented in compute_table.cc.
    Note - this is a singleton class but this is not enforced.
  
    This is exposed here because it allows us to avoid a
    "chicken and egg" problem:  to initialize the library, we want to 
    set the compute table style, but we cannot guarantee that those
    pointers are set up before we initialize the library.
    So, settings for compute tables should be made as follows.

    (1) call defaultInitializerList(), to build an instance of this class,
        and save the result.  That will set up the default settings.

    (2) change settings using static members

    (3) initialize Meddly using the saved initializer list.

*/
class MEDDLY::ct_initializer : public initializer_list {
  public:
    enum staleRemovalOption {
      /// Whenever we see a stale entry, remove it.
      Aggressive,
      /// Only remove stales when we need to expand the table
      Moderate,
      /// Only remove stales during Garbage Collection.
      Lazy
    };

    enum builtinCTstyle {
      /// One huge hash table that uses chaining.
      MonolithicChainedHash,

      /// One huge hash table that does not use chaining.
      MonolithicUnchainedHash,

      /// A hash table (with chaining) for each operation.
      OperationChainedHash,

      /// A hash table (no chaining) for each operation.
      OperationUnchainedHash,
    };

    enum compressionOption {
      /// No compression at all
      None,
      /// Compression based on item type
      TypeBased
      // TBD - others
    };

    struct settings {
      public:
        /// Memory manager to use for compute table entries
        const memory_manager_style *MMS;
        /// Maximum compute table size
        size_t maxSize;
        /// Stale removal policy
        staleRemovalOption staleRemoval;
        /// Compression policy
        compressionOption compression;
      public:
        settings() {
          MMS = 0;
        }
    };

  public:
    ct_initializer(initializer_list* previous);
    virtual ~ct_initializer();

  protected:
    virtual void setup();
    virtual void cleanup();
    static void setMemoryManager(const memory_manager_style*);

  // use these to change defaults, before library initialization
  public:
    static void setStaleRemoval(staleRemovalOption sro);
    static void setMaxSize(unsigned ms);
    static void setBuiltinStyle(builtinCTstyle cts);
    static void setUserStyle(const compute_table_style*);
    static void setCompression(compressionOption co);

    // for convenience
    static compute_table* createForOp(operation* op, unsigned slot);

  private:
    static settings the_settings;
    static const compute_table_style* ct_factory;
    static compute_table_style* builtin_ct_factory;
};

// ******************************************************************
// *                                                                *
// *                    compute_table_style class                   *
// *                                                                *
// ******************************************************************

/** Interface for building compute tables.
*/
class MEDDLY::compute_table_style {
  public:
    compute_table_style();
    virtual ~compute_table_style();

    /** Build a new, monolithic table.
        Monolithic means that the table stores entries for several
        (ideally, all) operations.

        Default throws an error.
    */
    virtual compute_table* create(const ct_initializer::settings &s)
      const;


    /**
        Build a new table for a single operation.
        Default throws an error.
    */
    virtual compute_table* create(const ct_initializer::settings &s,
      operation* op, unsigned slot) const;


    /**
        Does this style build monolithic CTs?
    */
    virtual bool usesMonolithic() const = 0;
};

// ******************************************************************
// *                                                                *
// *                      compute_table  class                      *
// *                                                                *
// ******************************************************************

/** Interface for compute tables.
    Anyone implementing an operation (see below) will
    probably want to use this.
    Implementation is in compute_table.cc.
*/
class MEDDLY::compute_table {
    protected:
      /// The maximum size of the hash table.
      unsigned maxSize;
      /// Do we try to eliminate stales during a "find" operation
      bool checkStalesOnFind;
      /// Do we try to eliminate stales during a "resize" operation
      bool checkStalesOnResize;

    public:

      //
      // ******************************************************************
      //

      struct stats {
        unsigned long numEntries;
        unsigned long hits;
        unsigned long pings;
        static const unsigned searchHistogramSize = 256;
        unsigned long searchHistogram[searchHistogramSize];
        unsigned long numLargeSearches;
        unsigned maxSearchLength;
        unsigned long completedScans;
        unsigned long resizeScans;
      };

      //
      // ******************************************************************
      //

      enum typeID {
        ERROR = 0,
        NODE = 1,
        INTEGER = 2,
        LONG = 3,
        FLOAT = 4,
        DOUBLE = 5,
        GENERIC = 6
      };

      //
      // ******************************************************************
      //

      union entry_item {
        int I;
        unsigned int U;
        long L;
        unsigned long UL;
        node_handle N;
        float F;
        double D;
        ct_object* G;
      };

      //
      // ******************************************************************
      //

      /**
        Type information about entries.
        Usually there is one type of entry for each operation,
        but there could be more than one type.

        These are built by operations and then registered
        with the compute table.
      */
      class entry_type {
        public:
          /**
            Constructor.
              @param  name    Name of the entry type; used only for displaying
                              CT entries (usually while debugging).

              @param  pattern Pattern for an entry.  The following characters
                              are supported in this string:
                                'N': node (in a forest)
                                'I': int 
                                'L': long 
                                'F': float
                                'D': double
                                'G': pointer to a ct_object
                                ':': separates key portion from result portion;
                                     must appear exactly once
                                '.': for repeating entries; can appear at most once.
                                     Everything between '.' and ':' can repeat
                                     zero or more times.

              @throws INVALID_ARGUMENT if pattern is illegal
          */
          entry_type(const char* name, const char* pattern);
          ~entry_type();

          unsigned getID() const;

          /**
            Set the forest for 'N' items in the pattern.
              @param  i   Slot.  Character i in the pattern must be 'N'.
              @param  f   Forest.
          */
          void setForestForSlot(unsigned i, expert_forest* f);

          /**
            Results might be overwritten.
            Indicate that in these entries, the result portion of the
            entry might be updated.
            The CT will make storage decisions based on this.
          */
          void mightUpdateResults();
          
          /**
              Is the result portion updatable?
          */
          bool isResultUpdatable() const;

          //
          // The remaining interface is for use by the compute table.
          // All these should be inlined for speed (see meddly_expert.hh)
          //

          const char* getName() const;

          /**
              Does this entry type allow repetitions in the key?
              I.e., was there a '.' in the pattern?
          */
          bool isRepeating() const;

          /**
              Get the number of items in the key.
                @param  reps  Number of repetitions.
                              If this is not a repeating type,
                              then this is ignored.

                @return Total number of slots in the key.
          */
          unsigned getKeySize(unsigned reps) const;

          /**
              Get the number of bytes in the key.
                @param  reps  Number of repetitions.
                              If this is not a repeating type,
                              then this is ignored.

                @return Total number of bytes required for the key.
          */
          unsigned getKeyBytes(unsigned reps) const;

          /**
              Get the type for item i in the key.
              Automatically handles repetitions.
                @param  i   Slot number, between 0 and getKeySize().   

                @param  t   On output, the type for item i.
                @param  f   If t is 'N', the forest for item i.
                            Otherwise, null.
          */
          void getKeyType(unsigned i, typeID &t, expert_forest* &f) const;

          /**
              Get the type for item i in the key.
              Automatically handles repetitions.
                @param  i   Slot number, between 0 and getKeySize().   
          */
          typeID getKeyType(unsigned i) const;

          /**
              Get the forest for item i in the key.
              Automatically handles repetitions.
                @param  i   Slot number, between 0 and getKeySize().
                @return     Forest for that slot, or 0 if the type
                            is not 'N'.
          */
          expert_forest* getKeyForest(unsigned i) const;

          /**
              Get the number of items in the result
          */
          unsigned getResultSize() const;

          /**
              Get the number of bytes in the result
          */
          unsigned getResultBytes() const;

          /**
              Get the type for item i in the result.
                @param  i   Slot number, between 0 and getResultSize().

                @param  t   On output, the type for item i.
                @param  f   If t is 'N', the forest for item i.
                            Otherwise, null.
          */
          void getResultType(unsigned i, typeID &t, expert_forest* &f) const;

          /**
              Get the type for item i in the result.
                @param  i   Slot number, between 0 and getResultSize().
          */
          typeID getResultType(unsigned i) const;

          /**
              Get the forest for item i in the result.
                @param  i   Slot number, between 0 and getResultSize().
                @return     Forest for that slot, or 0 if the type
                            is not 'N'.
          */
          expert_forest* getResultForest(unsigned i) const;

          /// Mark for deletion
          void markForDeletion();
 
          /// Unmark for deletion
          void unmarkForDeletion();
 
          /// Should we remove all CT entries of this type?
          bool isMarkedForDeletion() const;
        private:
          /// Unique ID, set by compute table
          unsigned etID;

          const char* name;

          /// Starting portion of key pattern.
          typeID* ks_type;
          /// Forests in starting portion of key.
          expert_forest** ks_forest;
          /// Length of ks_type and ks_forest arrays.
          unsigned len_ks_type;
          /// Total bytes in the starting portion of the key.
          unsigned ks_bytes;

          /// Repeating portion of key pattern (or null for no repeats).
          typeID* kr_type;
          /// Forests in repeating portion of key (or null).
          expert_forest** kr_forest;
          /// Length of kr_type and kr_forest arrays (zero if no repeats).
          unsigned len_kr_type;
          /// Total bytes in the repeating portion of the key.
          unsigned kr_bytes;

          /// Result pattern
          typeID* r_type;
          /// Forests in result
          expert_forest** r_forest;
          /// Length of r_type and r_forest arrays.
          unsigned len_r_type;
          /// Total bytes in the result.
          unsigned r_bytes;

          bool updatable_result;

          bool is_marked_for_deletion;

          friend class compute_table;
      };

      //
      // ******************************************************************
      //

      /** 
        The key portion of an entry.
        Internally, in the compute table, we may store
        entries differently.  This class is used to build
        keys for searching and to construct CT entries.
      */
      class entry_key {
        public:
          entry_key();
          ~entry_key();

        protected:
          /// Start using for this operation
          void setup(const compute_table::entry_type* et, unsigned repeats);

        public:
          const compute_table::entry_type* getET() const;

          // interface, for operations.  All inlined in meddly_expert.hh
          void writeN(node_handle nh);
          void writeI(int i);
          void writeL(long i);
          void writeF(float f);
          // For templates
          inline void write_ev(long i)  { writeL(i); }
          inline void write_ev(float f) { writeF(f); }

        public: 
          // interface, for compute_table.  All inlined in meddly_expert.hh
          const entry_item* rawData() const;
          unsigned dataLength() const;
          unsigned numRepeats() const;

          const void* readTempData() const;
          unsigned numTempBytes() const;
          void* allocTempData(unsigned bytes);
          /// Increase cache counters for nodes in this portion of the entry.
          void cacheNodes() const;
          unsigned getHash() const;

        protected:
          // protected interface, for compute_table.  All inlined in meddly_expert.hh
          void setHash(unsigned h);

        private:
          typeID theSlotType() const;

        private:
          const compute_table::entry_type* etype;
          entry_item* data;
          void* temp_data; 
          unsigned temp_bytes;
          unsigned temp_alloc;
          unsigned num_repeats;
          unsigned hash_value;
          unsigned data_alloc;

          unsigned currslot;
          unsigned total_slots;
#ifdef DEVELOPMENT_CODE
          bool has_hash;
#endif
        protected:
          /// Used for linked-list of recycled search keys in compute_table
          entry_key* next;

        friend class compute_table;
      };

      //
      // ******************************************************************
      //

      /** 
        The result portion of an entry.
        Internally, in the compute table, we may store
        entries differently.  This class is used to return 
        results from searches and to construct CT entries.
      */
      class entry_result {
        public:
          entry_result();
          ~entry_result();

        public:
          // For delayed construction
          void initialize(const entry_type* et);

          // interface, for operations (reading).
          node_handle readN();
          int readI();
          float readF();
          long readL();
          double readD();
          ct_object* readG();
          // for templates
          void read_ev(long &l)   { l = readL(); }
          void read_ev(float &f)  { f = readF(); }

          // interface, for operations (building). 
          void reset();
          void writeN(node_handle nh);
          void writeI(int i);
          void writeF(float f);
          void writeL(long L);
          void writeD(double D);
          void writeG(ct_object* G);

          // interface, for compute tables.
          void setValid();
          void setValid(const entry_item* d);
          void setInvalid();
          operator bool() const;
          /// Increase cache counters for nodes in this portion of the entry.
          void cacheNodes() const;

          const entry_item* rawData() const;
          unsigned dataLength() const;


        private:
          const entry_type* etype;
          entry_item* build;
          const entry_item* data;
          bool is_valid;
          unsigned currslot;
      };

      //
      // ******************************************************************
      //

      // convenience methods, for grabbing edge values
      static void readEV(const node_handle* p, int &ev);
      static void readEV(const node_handle* p, long &ev);
      static void readEV(const node_handle* p, float &ev);

      /// Constructor
      compute_table(const ct_initializer::settings &s);

      /** Destructor.
          Does NOT properly discard all table entries;
          use \a removeAll() for this.
      */
      virtual ~compute_table();

      /**
          Start using an entry_key for the given operation.
      */
      static entry_key* useEntryKey(const entry_type* et, unsigned repeats);

      /**
          Done using an entry_key.
      */
      static void recycle(entry_key* k);


      /// Is this a per-operation compute table?
      virtual bool isOperationTable() const = 0;

      /** Find an entry in the compute table based on the key provided.
          @param  key   Key to search for.
          @param  res   Where to store the result, if any.
      */
      virtual void find(entry_key* key, entry_result &res) = 0;

      /**
          Add an entry (key plus result) to the compute table.
            @param  key   Key portion of the entry.  Will be recycled.
            @param  res   Result portion of the entry.
      */
      virtual void addEntry(entry_key* key, const entry_result &res) = 0;

      /**
          Update an existing entry in the compute table.
            @param  key   Key portion of the entry.  Will be recycled.
            @param  res   Updated result portion of the entry.
      */
      virtual void updateEntry(entry_key* key, const entry_result &res) = 0;

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
      const stats& getStats();

      /// For debugging.
      virtual void show(output &s, int verbLevel = 0) = 0;

      /** Also for debugging.
          Examine all entries, and for each pointer to forest f node p,
          increment counts[p].
      */
      virtual void countNodeEntries(const expert_forest* f, size_t* counts) const = 0;


      static void initialize();
      static void destroy();

    protected:
      /** Register an operation.
          Sets aside a number of entry_type slots for the operation.
      */
      static void registerOp(operation* op, unsigned num_ids);

      /// Register an entry_type.
      static void registerEntryType(unsigned etid, entry_type* et);

      /** Unregister an operation.
          Frees the entry_type slots for the operation.
      */
      static void unregisterOp(operation* op, unsigned num_ids);

    public:
      /// Find entry_type for operation and slot number.
      static const entry_type* getEntryType(operation* op, unsigned slot);

      /// Find entry type for given entryID
      static const entry_type* getEntryType(unsigned etID);

    protected:
      void setHash(entry_key *k, unsigned h);

    protected:
      stats perf;

    private:
      static entry_type** entryInfo;
      static unsigned entryInfoAlloc;
      static unsigned entryInfoSize;

    private:
      static entry_key* free_keys;

    friend class operation;
};

// ******************************************************************
// *                                                                *
// *                        operation  class                        *
// *                                                                *
// ******************************************************************

/** Generic operation.
    Operations are tied to specific forests.
    Necessary for compute table entries.
*/
class MEDDLY::operation {
    const opname* theOpName;
    int oplist_index;
    bool is_marked_for_deletion;

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
    compute_table::entry_type** etype;
    /** Array of entry results.
        Use these during computation.
        We only ever need one result per entry type.
    */
    compute_table::entry_result* CTresult;

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

    void registerEntryType(unsigned slot, compute_table::entry_type* et);
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

    int getIndex() const;
    static operation* getOpWithIndex(int i);
    static int getOpListSize();

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
    virtual void compute(const dd_edge &arg, dd_edge &res);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res);
    virtual void compute(const dd_edge &arg, long &res);
    virtual void compute(const dd_edge &arg, double &res);
    virtual void compute(const dd_edge &arg, ct_object &c);

    // TBD: low-level front-ends?
    // e.g.,
    // virtual int computeDD(int k, int p);
    // virtual void computeEvDD(int k, int v, int p, int &w, int &q);
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
    virtual void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);
    virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
      = 0;

    // low-level front ends

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

  bool restrict_exist(node_handle p, const std::vector<int>& pa, int start,
      node_handle& result);
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

#include "meddly_expert.hh"
#endif

