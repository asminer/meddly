
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

// Flags for development version only. Significant reduction in performance.
#ifdef DEVELOPMENT_CODE
#define RANGE_CHECK_ON
#define DCASSERTS_ON
#endif

// #define DEBUG_MDD_H
// #define TRACK_DELETIONS
// #define TRACK_CACHECOUNT

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

//
// Design decision: should we remember the hashes for a reduced node?
// Tests so far indicate: doesn't save much, if any, time; and has
// a noticeable increase in memory.
//
// #define SAVE_HASHES

// #define NODE_STORAGE_PER_LEVEL
 
// #define ONE_NODE_FREELIST

namespace MEDDLY {

  // Functions for reinterpreting an int to a float and vice-versa
  float   toFloat (int a);
  int     toInt   (float a);

  // classes defined here
  struct settings;
  class expert_variable;
  class expert_domain;

  // wrappers for nodes
  class node_reader;
  class node_builder;

  // Actual node storage
  struct node_header;
  class node_storage;

  class expert_forest;

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
  class unique_table;

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
  extern const numerical_opname* VECT_MATR_MULT;

  /** Computes y = y + Ax.
      x and y are vectors, stored explicitly, and A is a matrix.
      x_ind and y_ind specify how minterms are mapped to indexes
      for vectors x and y, respectively.
  */
  extern const numerical_opname* MATR_VECT_MULT;

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
  void destroyOperation(numerical_operation* &op);

  /// Should not be called directly.
  void destroyOpInternal(operation* op);
  
  /// Builds an initializer for MEDDLY's builtin operations.
  op_initializer* makeBuiltinInitializer();

}; // namespace MEDDLY


// ******************************************************************
// *                                                                *
// *                         ct_object class                        *
// *                                                                *
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
// *                                                                *
// *                         settings  class                        *
// *                                                                *
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
    /// Default forest policies for "sets"
    forest::policies mddDefaults;
    /// Default forest policies for "relations"
    forest::policies mxdDefaults;
  public:
    /// Constructor, to set defaults.
    settings() : computeTable(), mddDefaults(0), mxdDefaults(1) {
      operationBuilder = makeBuiltinInitializer();
    }
    /// Copy constructor.
    settings(const settings &s) : mddDefaults(0), mxdDefaults(1) { 
      init(s); 
    }
    /// Destructor
    ~settings() { clear(); }
    inline void operator=(const settings &s) {
      if (&s != this) {
        clear();
        init(s);
      }
    }
    /// super handly
    inline bool usesMonolithicComputeTable() {
      return (
       computeTableSettings::MonolithicChainedHash == computeTable.style ||
       computeTableSettings::MonolithicUnchainedHash == computeTable.style 
      );
    }
  private:
    void init(const settings &s);
    void clear();
};
  
// ******************************************************************
// *                                                                *
// *                     expert_variable  class                     *
// *                                                                *
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
// *                                                                *
// *                      expert_domain  class                      *
// *                                                                *
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
};


// ******************************************************************
// *                                                                *
// *                       node_reader  class                       *
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
class MEDDLY::node_reader {
    public:
        /** Constructor.
            The class must be "filled" by a forest befure
            it can be used, however.
        */  
        node_reader();

        /// Destructor.
        ~node_reader();

        /// Free memory, but don't delete.
        void clear();

        /** Get a downward pointer.
              @param  n   Which pointer.
              @return     If this is a full reader, 
                          return pointer with index n.
                          If this is a sparse reader,
                          return the nth non-zero pointer.
        */
        inline node_handle d(int n) const {
          MEDDLY_DCASSERT(down);
          MEDDLY_CHECK_RANGE(0, n, (is_full ? size : nnzs));
          return down[n];
        }

        /** Get the index of the nth non-zero pointer.
            Use only for sparse readers.
        */
        inline int i(int n) const {
          MEDDLY_DCASSERT(index);
          MEDDLY_DCASSERT(!is_full);
          MEDDLY_CHECK_RANGE(0, n, nnzs);
          return index[n];
        }
        
        /// Get the edge value, as an integer.
        inline int ei(int n) const {
          MEDDLY_DCASSERT(edge);
          MEDDLY_CHECK_RANGE(0, n, (is_full ? size : nnzs));
          MEDDLY_DCASSERT(sizeof(int) == edge_bytes);
          return ((int*)edge)[n];
        }

        /// Get the edge value, as a float.
        inline float ef(int n) const {
          MEDDLY_DCASSERT(edge);
          MEDDLY_CHECK_RANGE(0, n, (is_full ? size : nnzs));
          MEDDLY_DCASSERT(sizeof(float) == edge_bytes);
          return ((float*)edge)[n];
        }

        /// Get a pointer to an edge
        inline const void* eptr(int i) const {
          MEDDLY_DCASSERT(edge);
          MEDDLY_CHECK_RANGE(0, i, size);
          return ((char*)edge) + i * edge_bytes;
        }

        /// Get the raw edge array.
        inline const void* rawEdges() const {
          return edge;
        }

        /// Get the level number of this node.
        inline int getLevel() const {
          return level;
        }

        /// Get the size of this node (full readers only).
        inline int getSize() const {
          MEDDLY_DCASSERT(is_full);
          return size;
        }
        /// Get the number of nonzeroes of this node (sparse readers only).
        inline int getNNZs() const {
          MEDDLY_DCASSERT(!is_full);
          return nnzs;
        }
        /// Is this a sparse reader?
        inline bool isSparse() const { 
          return !is_full;
        }
        /// Is this a full reader?
        inline bool isFull() const {
          return is_full;
        }
        /// Does this node have edge values?
        inline bool hasEdges() const {
          return edge_bytes;
        }
        /// Number of bytes per edge
        inline int edgeBytes() const {
          return edge_bytes;
        }

        // For debugging unique table.
        inline unsigned hash() const {
#ifdef DEVELOPMENT_CODE
          MEDDLY_DCASSERT(has_hash);
#endif
          return h;
        }
        inline void setHash(unsigned H) {
#ifdef DEVELOPMENT_CODE
          MEDDLY_DCASSERT(!has_hash);
          has_hash = true;
#endif
          h = H;
        };

        void computeHash(bool hashEdgeValues);

    // Centralized recycling
    public:
        inline static node_reader* useReader() {
          node_reader* nr;
          if (freeList) {
            nr = freeList;
            freeList = nr->next;
          } else {
            nr = new node_reader;
          }
#ifdef DEVELOPMENT_CODE
          nr->has_hash = false;
#endif
          return nr;
        }
        inline static void recycle(node_reader* r) {
          if (r) {
            r->next = freeList;
            freeList = r;
          }
        }
    private:
        inline static void freeRecycled() {
          while (freeList) {
            node_reader* n = freeList->next;
            delete freeList;
            freeList = n;
          }
        }

        static node_reader* freeList;
    private:
        node_reader* next;    // for recycled list
        /* 
          Extra info that is hashed
        */
        void* extra_hashed;
        int ext_alloc;
        int ext_size;
        /*
          Down pointers, indexes, edge values.
        */
        node_handle* down;
        int* index;
        void* edge;
        int alloc;
        int ealloc;
        int size;
        int nnzs;
        int level;
        unsigned h;
        char edge_bytes;  // number of bytes for an edge value.
        bool is_full;
#ifdef DEVELOPMENT_CODE
        bool has_hash;
#endif

        /// Called within expert_forest to allocate space.
        ///   @param  p     Parent.
        ///   @param  k     Level number.
        ///   @param  ns    Size of node.
        ///   @param  eb    Bytes for each edge.
        ///   @param  full  If true, we'll be filling a full reader.
        ///                 Otherwise it is a sparse one.
        void resize(int k, int ns, char eb, bool full);

        /// Allocate space for extra hashed info
        void resize_header(int extra_bytes);

        friend class expert_forest; 
        friend class node_storage; 
        friend void cleanup();

};


// ******************************************************************
// *                                                                *
// *                       node_builder class                       *
// *                                                                *
// ******************************************************************

/** Class for building nodes.
    Effectively, a reserved chunk of memory for storing down pointers
    and edge values.
    Implemented in node_wrappers.cc.
*/
class MEDDLY::node_builder {
        const expert_forest* parent;
        void* extra_hashed;
        void* extra_unhashed;
        node_handle* down;
        int* indexes;
        void* edge;
        int level;  // level of this node.
        int size;
        int alloc;
        unsigned h;
        char edge_bytes;
#ifdef DEVELOPMENT_CODE
        bool has_hash;
#endif
        bool is_sparse;
    private:
        // These cheat; let's see if any compiler complains.
        inline void* raw_hh() const {
          MEDDLY_DCASSERT(extra_hashed);
          return extra_hashed;
        }
        inline void* raw_uh() const {
          MEDDLY_DCASSERT(extra_unhashed);
          return extra_unhashed;
        }
        inline node_handle& raw_d(int i) const {
          MEDDLY_DCASSERT(down);
          MEDDLY_CHECK_RANGE(0, i, size);
          return down[i];
        }
        inline int& raw_i(int i) const {
          MEDDLY_DCASSERT(indexes);
          MEDDLY_DCASSERT(is_sparse);
          MEDDLY_CHECK_RANGE(0, i, size);
          return indexes[i];
        }
        inline int& raw_ei(int i) const {
          MEDDLY_DCASSERT(edge);
          MEDDLY_CHECK_RANGE(0, i, size);
          return ((int*)edge)[i];
        }
        inline float& raw_ef(int i) const {
          MEDDLY_DCASSERT(edge);
          MEDDLY_CHECK_RANGE(0, i, size);
          return ((float*)edge)[i];
        }
    public:
        bool lock;
    public:
        node_builder();
        ~node_builder();
        void init(int k, const expert_forest* p);
        inline bool hasEdges() const {
          return edge_bytes > 0;
        }
        inline void resize(int s) {
          is_sparse = false;
#ifdef DEVELOPMENT_CODE
          has_hash = false;
#endif
          size = s;
          if (size > alloc) enlarge();
        }
        inline void resparse(int s) {
          is_sparse = true;
#ifdef DEVELOPMENT_CODE
          has_hash = false;
#endif
          size = s;
          if (size > alloc) enlarge();
        }
        // used when we don't know exactly the sparse size
        inline void shrinkSparse(int ns) {
          MEDDLY_DCASSERT(is_sparse);
          MEDDLY_CHECK_RANGE(0, ns, size+1);
          size = ns;
        }
        inline bool isSparse() const {
          return is_sparse;
        }
        inline bool isFull() const {
          return !is_sparse;
        }
        inline int rawSize() const {
          return size;
        }
        inline int getSize() const { 
          MEDDLY_DCASSERT(!is_sparse);
          return size;
        }
        inline int getNNZs() const {
          MEDDLY_DCASSERT(is_sparse);
          return size;
        }
        inline int getLevel() const { 
          return level;
        }
        void setUH(const void* x);
        void getUH(void* x) const;
        void setHH(const void* x);
        void getHH(void* x) const;
        inline node_handle& d(int i)        { return  raw_d(i); }
        inline node_handle  d(int i) const  { return  raw_d(i); }
        inline int&  i(int i)               { return  raw_i(i); }
        inline int   i(int i) const         { return  raw_i(i); }
        inline int& ei(int i)               { return raw_ei(i); }
        inline int  ei(int i) const         { return raw_ei(i); }
        inline float& ef(int i)             { return raw_ef(i); }
        inline float  ef(int i) const       { return raw_ef(i); }
        // edge setting - for templates
        inline void setEdge(int i, int ev) {
          ei(i) = ev;
        }
        inline void setEdge(int i, float ev) {
          ef(i) = ev;
        }
        // raw edge info
        inline const void* eptr(int i) const {
          MEDDLY_DCASSERT(edge);
          MEDDLY_CHECK_RANGE(0, i, size);
          return ((char*)edge) + i * edge_bytes;
        }
    public: // for unique table
        inline unsigned hash() const {
#ifdef DEVELOPMENT_CODE
          MEDDLY_DCASSERT(has_hash);
#endif
          return h;
        }
        void computeHash();

    protected:
        void enlarge();

}; // end of node_builder class




// ******************************************************************
// *                                                                *
// *                       node_header  class                       *
// *                                                                *
// ******************************************************************

/** Header information for each node.
    The node handles (integers) need to keep some
    information about the node they point to,
    both for addressing and for bookkeeping purposes.
    This struct holds that information.
*/
struct MEDDLY::node_header {
    /** Offset to node's data in the corresponding node storage structure.
        If the node is active, this is the offset (>0) in the data array.
        If the node is deleted, this is -next deleted node
        (part of the unused address list).
    */
    node_address offset;

    /** Node level
        If the node is active, this indicates node level.
    */
    int level;

    /** Cache count
        The number of cache entries that refer to this node (excl. unique
        table). If this node is a zombie, cache_count is negative.
    */
    int cache_count;

#ifdef SAVE_HASHES
    /// Remember the hash for speed.
    unsigned hash;
#endif

    // Handy functions, in case our internal storage changes.

    inline bool isActive() const  { return offset > 0; }
    inline bool isZombie() const  { return cache_count < 0; }

    inline bool isDeleted() const { return 0 == level; }
    inline void setDeleted()      { level = 0; }
    inline void setNotDeleted()   { level = 1; }
      
    inline int getNextDeleted() const { return -offset; }
    inline void setNextDeleted(int n) { offset = -n; }

    inline void makeZombie() { 
      MEDDLY_DCASSERT(cache_count > 0);
      cache_count *= -1; 
      offset = 0;
    }
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
    node_storage();
    virtual ~node_storage();

    /** Build a new mechanism, bound to the given forest.

        Any instance serves also as a factory to build new instances.
        However, an instance is not initialized until it is bound
        to a given forest.

          @param  f   Forest to bind to

          @return     A pointer to a node storage class,
                      initialized for forest f.
                      Most likely, we are returning the same
                      instance as ourself, but this is not required.
    */
    virtual node_storage* createForForest(expert_forest* f) const = 0;

    /** Should be called by createForForest in derived classes.
        Initializes the base class for a given forest.

          @param  f   Forest to bind to
    */
    void initForForest(expert_forest* f);

    /** Go through and collect any garbage.
        
          @param  shrink    If true, we will shrink data structures
                            as we free space; otherwise, we won't.
    */
    virtual void collectGarbage(bool shrink) = 0;

    /** Show various stats.
          @param  s       Output stream to write to
          @param  pad     Padding string, written at the start of 
                          each output line.
          @param  flags   Controls what is displated.
    */
    virtual void reportStats(FILE* s, const char* pad, unsigned flags) const = 0;

    /** Write a node in human-readable format.

        Ideally, the output format for each node is the same
        regardless of how it is stored.
          
        @param  s       Output stream.
        @param  addr    Address of the node we care about.
        @param  verb    Verbose output or not.
    */
    virtual void showNode(FILE* s, node_address addr, bool verb) const = 0;


    /** Dump the internal storage details.
        Primarily used for debugging.

          @param  s       Output stream to use
          @param  flags   What to show.               
                            0x01  Show active memory
                            0x02  Show memory "holes"
    */
    void dumpInternal(FILE* s, unsigned flags) const;

 

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
    virtual node_address makeNode(node_handle p, const node_builder &nb, 
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
          @param  nb      Node to compare against

          @return true    iff the nodes are duplicates
    */
    virtual bool areDuplicates(node_address addr, const node_builder &nb) const = 0;

    /** Check for duplicates.
          @param  addr    Node address in this structure
          @param  nr      Node to compare against

          @return true    iff the nodes are duplicates
    */
    virtual bool areDuplicates(node_address addr, const node_reader &nr) const = 0;

    /** Fill a node reader.
        Useful if we want to read an entire node.
          @param  addr    Node address in this structure.
          @param  nr      Node reader to fill, has already
                          been "resized".  Full or sparse
                          will be decided from \a nr settings.
    */
    virtual void fillReader(node_address addr, node_reader &nr) const = 0;

    /** Compute the hash value for a node.
        This requires an actual "node", not just an address,
        because the hash contains the node's level.
        Should give the same answer as filling a node_reader
        and computing the hash on the node_reader.

          @param  p     Node of interest.
    */
    virtual unsigned hashNode(const node_header& p) const = 0;


    /** Determine if this is a singleton node.
        Used for identity reductions.
          @param  addr    Address of the node we care about
          @param  down    Output:
                          The singleton downward pointer, or undefined.

          @return   If the node has only one non-zero downward pointer,
                    then return the index for that pointer.
                    Otherwise, return a negative value.
    */
    virtual int getSingletonIndex(node_address addr, node_handle &down) const = 0;


    /** Get the specified downward pointer for a node.
        Fast if we just want one.
          @param  addr    Address of the node we care about
          @param  index   Index of downward pointer
          @return         Desired pointer
    */
    virtual node_handle getDownPtr(node_address addr, int index) const = 0;

    /** Get the specified outgoing edge for a node.
        Fast if we just want one.

          @param  addr    Address of the node we care about
          @param  ind     Index of the pointer we want.
          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    virtual void getDownPtr(node_address addr, int ind, int& ev, node_handle& dn) const = 0;

    /** Get the specified outgoing edge for a node.
        Fast if we just want one.

          @param  addr    Address of the node we care about
          @param  ind     Index of the pointer we want.
          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    virtual void getDownPtr(node_address addr, int ind, float& ev, node_handle& dn) const = 0;


    /** Read the unhashed header portion of a node.

          @param  addr    Address of the node we care about
    */
    virtual const void* getUnhashedHeaderOf(node_address addr) const = 0;

    /** Read the hashed header portion of a node.

          @param  addr    Address of the node we care about
    */
    virtual const void* getHashedHeaderOf(node_address addr) const = 0;

  // Inlines, for speed
  public:
    // --------------------------------------------------
    // incoming count data
    // --------------------------------------------------

    /// Get the number of incoming pointers to a node.
    inline node_handle getCountOf(node_address addr) const {
      MEDDLY_DCASSERT(counts);
      MEDDLY_DCASSERT(addr>0);
      return counts[addr];
    }

    /// Set the number of incoming pointers to a node.
    inline void setCountOf(node_address addr, node_handle c) {
      MEDDLY_DCASSERT(counts);
      MEDDLY_DCASSERT(addr>0);
      counts[addr] = c;
    }

    /// Increment (and return) the number of incoming pointers to a node.
    inline node_handle incCountOf(node_address addr) {
      MEDDLY_DCASSERT(counts);
      MEDDLY_DCASSERT(addr>0);
      return ++counts[addr];
    };

    /// Decrement (and return) the number of incoming pointers to a node.
    inline node_handle decCountOf(node_address addr) {
      MEDDLY_DCASSERT(counts);
      MEDDLY_DCASSERT(addr>0);
      return --counts[addr];
    }

    // --------------------------------------------------
    // next pointer data
    // --------------------------------------------------

    inline node_handle getNextOf(node_address addr) const {
      MEDDLY_DCASSERT(nexts);
      MEDDLY_DCASSERT(addr>0);
      return nexts[addr];
    }

    inline void setNextOf(node_address addr, node_handle n) {
      MEDDLY_DCASSERT(nexts);
      MEDDLY_DCASSERT(addr>0);
      nexts[addr] = n;
    }
    
  protected:
    /// Dump information not related to individual nodes.
    virtual void dumpInternalInfo(FILE*) const = 0;

    /** Dump the node/hole information at the given address.
          @param  s       Output stream to use
          @param  addr    Address
          @param  flags   What chunks should be displayed

          @return   Next interesting address.
    */
    virtual node_address 
    dumpInternalNode(FILE* s, node_address addr, unsigned flags) const = 0;

    /// Dump final info (after node info)
    virtual void dumpInternalTail(FILE*) const = 0;


  // Hooks from other classes, so we don't need to make
  // all the derived classes "friends".
  protected:
    // inlined, below
    void moveNodeOffset(node_handle node, node_address old_addr, node_address new_addr); 

    static inline void resize_header(node_reader& nr, int extra_slots) {
      nr.resize_header(extra_slots);
    }
    static inline void* extra_hashed(node_reader& nr) {
      return nr.extra_hashed;
    }
    static inline node_handle* down_of(node_reader& nr) {
      return nr.down;
    }
    static inline int* index_of(node_reader& nr) {
      return nr.index;
    }
    static inline char* edge_of(node_reader& nr) {
      return (char*) nr.edge;
    }
    static inline int& nnzs_of(node_reader& nr) {
      return nr.nnzs;
    }

  //
  // Methods for derived classes to deal with
  // members owned by the base class
  //
  protected:
    inline const expert_forest* getParent() const {
      MEDDLY_DCASSERT(parent);
      return parent;
    }
    inline expert_forest* getParent() {
      MEDDLY_DCASSERT(parent);
      return parent;
    }

    inline void incMemUsed(long delta) {
      if (stats) stats->incMemUsed(delta);
    }
    inline void decMemUsed(long delta) {
      if (stats) stats->decMemUsed(delta);
    }
    inline void incMemAlloc(long delta) {
      if (stats) stats->incMemAlloc(delta);
    }
    inline void decMemAlloc(long delta) {
      if (stats) stats->decMemAlloc(delta);
    }
    inline void incCompactions() {
      if (stats) stats->num_compactions++;
    }
    inline void updateCountArray(node_handle* cptr) {
      counts = cptr;
    }
    inline void updateNextArray(node_handle* nptr) {
      nexts = nptr;
    }

  //
  // Hooks for hole managers
  //
  protected:
    // Change the data array
    virtual void updateData(node_handle* data) = 0;

    // How small can a node be?
    virtual int smallestNode() const = 0;

    friend class holeman;

  private:
    /// Parent forest.
    expert_forest* parent;

    /// Memory stats
    forest::statset* stats;

    /// Count array, so that counts[addr] gives the count for node at addr.
    node_handle* counts;

    /// Next array, so that nexts[addr] gives the next value for node at addr.
    node_handle* nexts;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                      expert_forest  class                      *
// *                                                                *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_forest : public forest
{
  // flags for reporting; DO NOT rely on specific values
  public:
      /// Should memory be reported in a human readable format
      static const unsigned HUMAN_READABLE_MEMORY   = 0x0001;
      /// Basic forest stats
      static const unsigned BASIC_STATS             = 0x0002;
      /// Extra forest stats
      static const unsigned EXTRA_STATS             = 0x0004;
      /// Specific forest stats, dependent on forest type
      static const unsigned FOREST_STATS            = 0x0008;
      /// Stats specific to the node storage mechanism.
      static const unsigned STORAGE_STATS           = 0x0010;
      /// Detailed stats for the node storage mechanism.
      static const unsigned STORAGE_DETAILED        = 0x0020;
      /// Stats specific to the unique table.
      static const unsigned UNIQUE_TABLE_STATS      = 0x0040;
      /// Stats specific to the unique table.
      static const unsigned UNIQUE_TABLE_DETAILED   = 0x0080;
      /// Stats specific to the hole manager.
      static const unsigned HOLE_MANAGER_STATS      = 0x0100;
      /// Stats specific to the hole manager.
      static const unsigned HOLE_MANAGER_DETAILED   = 0x0200;

  public:
    /** Constructor.
      @param  dslot   slot used to store the forest, in the domain
      @param  d       domain to which this forest belongs to.
      @param  rel     does this forest represent a relation.
      @param  t       the range of the functions represented in this forest.
      @param  ev      edge annotation.
      @param  p       Polcies for reduction, storage, deletion.
    */
    expert_forest(int dslot, domain *d, bool rel, range_type t, 
      edge_labeling ev, const policies &p);

  protected:
    /// Destructor.
    virtual ~expert_forest();  

    /** Initialize data.
        Should be called in the child class constructors.
        Allows us to use class properties to initialize the data.
    */
    void initializeForest();
    
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                         public methods                         |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ------------------------------------------------------------
  // inlined helpers.
  public:
    inline statset& changeStats() {
      return stats;
    }
    /// Number of bytes for an edge value.
    inline char edgeBytes() const {
      return edge_bytes;
    }
    /// Are edge values included when computing the hash.
    inline bool areEdgeValuesHashed() const {
      return hash_edge_values;
    }
    /// Extra bytes per node, not hashed.
    inline char unhashedHeaderBytes() const {
      return unhashed_bytes;
    }
    /// Extra bytes per node, hashed.
    inline char hashedHeaderBytes() const {
      return hashed_bytes;
    }


// TBD: MOVE THESE TO MT CLASS?

    /// Get the integer value represented by this terminal node.
    static inline bool getBoolean(node_handle terminalNode) {
      MEDDLY_DCASSERT(terminalNode <= 0);
      return terminalNode;
    }

    /// Get the terminal node representing this boolean value.
    static inline node_handle getTerminalNode(bool booleanValue) {
      return booleanValue ? -1 : 0;
    }

    /// Get the integer value represented by this terminal node.
    static inline int getInteger(node_handle terminalNode) {
      MEDDLY_DCASSERT(terminalNode <= 0);
      // set 32nd bit based on 31st bit.  << gets rid of MSB; >> sign extends.
      return terminalNode << 1 >> 1;
    }

    /// Get the terminal node representing this integer value.
    static inline node_handle getTerminalNode(int integerValue) {
      // value has to fit within 31 bits (incl. sign)
      // int(0xc0000000) == -1073741824
      // int(0x3fffffff) == +1073741823
      // MEDDLY_DCASSERT(-1073741824 <= integerValue && integerValue <= 1073741823);
      MEDDLY_DCASSERT(-1073741824 <= integerValue && integerValue <= 1073741823);
      return integerValue == 0? 0: integerValue | 0x80000000;
    }

    /// Get the real (float) value represented by this terminal node.
    static inline float getReal(int term) {
      MEDDLY_DCASSERT(term <= 0);
      MEDDLY_DCASSERT(sizeof(float) <= sizeof(int));
      if (0 == term) return 0.0;
      term <<= 1;
      float ret;
      memcpy(&ret, &term, sizeof(float));
      return ret;
    }

    /// Get the terminal node representing this real (float) value.
    static inline node_handle getTerminalNode(float a) {
      if (0.0 == a) return 0;
      int node;
      memcpy(&node, &a, sizeof(node_handle));
      return (node >> 1) | 0x80000000;
    }
// TO HERE



    inline const expert_domain* getExpertDomain() const {
      return (expert_domain*) getDomain();
    }
    inline expert_domain* useExpertDomain() {
      return (expert_domain*) useDomain();
    }
    /// Ignores prime/unprime.
    inline int getNumVariables() const {
      return getDomain()->getNumVariables();
    }
    inline int getMinLevelIndex() const {
      return isForRelations() ? -getNumVariables() : 0;
    }
    inline bool isValidLevel(int k) const {
      return (k>=getMinLevelIndex()) && (k<=getNumVariables());
    }
    /// The maximum size (number of indices) a node at this level can have
    inline int getLevelSize(int lh) const {
      MEDDLY_DCASSERT(isValidLevel(lh));
      if (lh < 0) {
        return getDomain()->getVariableBound(-lh, true);
      } else {
        return getDomain()->getVariableBound(lh, false);
      }
    }
    /// used by node_storage, for displaying.
#ifdef NODE_STORAGE_PER_LEVEL
    inline long getLevelNumber(const node_storage* L) const {
      for (long i=L; i>=-L; i--) {
        if (levels[i] == L) return i;
      }
      return 0;
    }
#endif
    /// Is this a terminal node?
    inline static bool isTerminalNode(node_handle p) {
      return (p < 1);
    }
    /// Sanity check: is this a valid nonterminal node index.
    inline bool isValidNonterminalIndex(node_handle node) const {
      return (node>0) && (node <= a_last);
    }
    /// Sanity check: is this a valid node index.
    inline bool isValidNodeIndex(node_handle node) const {
      return node <= a_last;
    }
    inline node_handle getLastNode() const {
      return a_last;
    }
    /// Get data for a given nonterminal node.
    inline const node_header& getNode(node_handle p) const {
      MEDDLY_DCASSERT(address);
      MEDDLY_CHECK_RANGE(1, p, 1+a_last);
      return address[p];
    }
    /** Get the node's level as an integer.
        Negative values are used for primed levels.
    */
    inline int getNodeLevel(node_handle p) const {
      if (isTerminalNode(p)) return 0;
      MEDDLY_DCASSERT(address);
      MEDDLY_CHECK_RANGE(1, p, 1+a_last);
      return address[p].level;
    }
    inline bool isPrimedNode(node_handle p) const {
      return getNodeLevel(p) < 0;
    }
    inline bool isUnprimedNode(node_handle p) const {
      return getNodeLevel(p) > 0;
    }

    /** Get the node's height.
        For convenience.  Height is the total
        number of levels until terminal nodes.
    */
    inline int getNodeHeight(node_handle p) const {
      if (isForRelations()) {
        int k = getNodeLevel(p);
        return (k<0) ? (-2*k-1) : 2*k;
      } else {
        return getNodeLevel(p);
      }
    }

    /** Does this node represent an Index Set?
        Note: Only applicable to EV+MDDs.
    */
    inline bool isIndexSet(node_handle node) const {
      if (!isEVPlus()) return false;
      if (isTerminalNode(node)) return node;
      // yes iff the unhashed extra header is non-zero.
      const node_header& nd = getNode(node);
#ifdef NODE_STORAGE_PER_LEVEL
      const node_storage* nodeMan = levels[nd.level];
#endif
      const int* uhh = (const int*) nodeMan->getUnhashedHeaderOf(nd.offset);
      return *uhh > 0;
    }

    /// Get the cardinality of an Index Set.
    inline int getIndexSetCardinality(node_handle node) const {
      MEDDLY_DCASSERT(isEVPlus());
      if (isTerminalNode(node)) return (node != 0) ? 1 : 0;
      // yes iff the unhashed extra header is non-zero.
      const node_header& nd = getNode(node);
#ifdef NODE_STORAGE_PER_LEVEL
      const node_storage* nodeMan = levels[nd.level];
#endif
      const int* uhh = (const int*) nodeMan->getUnhashedHeaderOf(nd.offset);
      MEDDLY_DCASSERT(*uhh > 0);
      return *uhh;
    }

    // --------------------------------------------------
    // Used by the unique table
    // --------------------------------------------------
    inline node_handle getNext(node_handle p) const {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
#ifdef NODE_STORAGE_PER_LEVEL
      MEDDLY_DCASSERT(levels);
      MEDDLY_DCASSERT(address[p].level);
      return levels[address[p].level]->getNextOf(address[p].offset);
#else
      return nodeMan->getNextOf(address[p].offset);
#endif
    }
    inline void setNext(node_handle p, node_handle n) {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
#ifdef NODE_STORAGE_PER_LEVEL
      MEDDLY_DCASSERT(levels);
      MEDDLY_DCASSERT(address[p].level);
      levels[address[p].level]->setNextOf(address[p].offset, n);
#else
      nodeMan->setNextOf(address[p].offset, n);
#endif
    }
    inline unsigned hash(node_handle p) const {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
#ifdef SAVE_HASHES
      return address[p].hash;
#else
      return hashNode(p);
#endif
    }

    // --------------------------------------------------
    // Managing reference counts
    // --------------------------------------------------

    /// Returns the in-count for a node.
    inline long readInCount(node_handle p) const {
#ifdef NODE_STORAGE_PER_LEVEL
      const node_header& node = getNode(p);
      return levels[node.level]->getCountOf(node.offset);
#else
      return nodeMan->getCountOf(getNode(p).offset);
#endif
    }

    /** Increase the link count to this node. Call this when another node is
        made to point to this node.
          @return p, for convenience.
    */
    inline long linkNode(node_handle p) {
        MEDDLY_DCASSERT(isActiveNode(p));
        if (isTerminalNode(p)) return p;
        MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));

        long count = incInCount(p);
        if (1 == count) {
          // Reclaim an orphan node
          stats.reclaimed_nodes++;
          stats.orphan_nodes--;
        }
#ifdef TRACK_DELETIONS
        fprintf(stdout, "\t+Node %d count now %ld\n", p, count);
        fflush(stdout);
#endif
        return p;
    }

    /** Decrease the link count to this node. If link count reduces to 0, this
        node may get marked for deletion. Call this when another node releases
        its connection to this node.
    */
    inline void unlinkNode(node_handle p) {
        MEDDLY_DCASSERT(isActiveNode(p));
        if (isTerminalNode(p)) return;
        MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
        MEDDLY_DCASSERT(getInCount(p) > 0);

        long count = decInCount(p);
  
#ifdef TRACK_DELETIONS
        fprintf(stdout, "\t-Node %d count now %ld\n", p, count);
        fflush(stdout);
#endif
        if (count) return;

        handleNewOrphanNode(p);
    }

    // --------------------------------------------------
    // Managing cache entries
    // --------------------------------------------------

    /** Increase the cache count for this node. Call this whenever this node
        is added to a cache. 
          @param  p     Node we care about.
          @return p, for convenience.
    */
    inline node_handle cacheNode(node_handle p) {
      MEDDLY_DCASSERT(isActiveNode(p));
      if (isTerminalNode(p)) return p;
      cacheCount(p)++;
#ifdef TRACK_CACHECOUNT
      fprintf(stdout, "\t+Node %d is in %d caches\n", p, getCacheCount(p));
      fflush(stdout);
#endif
      return p;
    }

    /** Increase the cache count for this node. Call this whenever this node
        is added to a cache. 
          @param  p     Node we care about.
    */
    inline void uncacheNode(node_handle p) {
      if (isTerminalNode(p)) return;
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
      MEDDLY_DCASSERT(isActiveNode(p) ||
          (!isActiveNode(p) && isPessimistic() && isZombieNode(p)));
      int& cc = cacheCount(p);
      if (isPessimistic() && isZombieNode(p)) {
        // special case: we store the negative of the count.
        MEDDLY_DCASSERT(cc < 0);
        cc++;
        if (0 == cc) {
          stats.zombie_nodes--;
          recycleNodeHandle(p);
        }
        return;
      }
      // we store the actual count.
      MEDDLY_DCASSERT(cc > 0);
      cc--;
#ifdef TRACK_CACHECOUNT
      fprintf(stdout, "\t-Node %d is in %d caches\n", p, cc);
      fflush(stdout);
#endif

      if (cc == 0 && readInCount(p) == 0) {
        MEDDLY_DCASSERT(!isPessimistic());
        stats.orphan_nodes--;
        deleteNode(p);
      }
    }


    /// A node can be discarded once it goes stale. Whether a node is
    /// considered stale depends on the forest's deletion policy.
    /// Optimistic deletion: A node is said to be stale only when both the
    ///   in-count and cache-count are zero.
    /// Pessimistic deletion: A node is said to be stale when the in-count
    ///  is zero regardless of the cache-count.
    inline bool isStale(node_handle node) const {
      return
        isMarkedForDeletion() || (
          isTerminalNode(node)
          ? terminalNodesAreStale
          : isPessimistic()
            ? isZombieNode(node)
            : (readInCount(node) == 0)
        );
    }


  // ------------------------------------------------------------
  // non-virtual, handy methods for debugging.
  public:
    void dump(FILE *s) const; 
    void dumpInternal(FILE *s) const; 
    void dumpUniqueTable(FILE *s) const;
    void validateIncounts(bool exact);


  // ------------------------------------------------------------
  // non-virtual, handy methods.
  public:
    /** Build a list of nodes in the subgraph below the given node.
        This for example is used to determine which nodes must
        be printed to display a subgraph.
        Terminal nodes are NOT included.

          @param  root    Root node in the forest.
                          Will be included in the list, unless it
                          is a terminal node.

          @param  sort    If true, the list will be in increasing order.
                          Otherwise, the list will be in some convenient order
                          (currently, it is the order that nodes 
                          are discovered).

          @return   A malloc'd array of non-terminal nodes, terminated by 0.
                    Or, a null pointer, if the list is empty.
    */
    node_handle* markNodesInSubgraph(node_handle root, bool sort) const;

    /** Count and return the number of non-terminal nodes
        in the subgraph below the given node.
    */
    long getNodeCount(node_handle node) const;

    /** Count and return the number of edges
        in the subgraph below the given node.
    */
    long getEdgeCount(node_handle node, bool countZeroes) const;

    /// Display the contents of a single node.
    void showNode(FILE* s, node_handle node, int verbose = 0) const;

    /// Show all the nodes in the subgraph below the given node.
    void showNodeGraph(FILE* s, node_handle node) const;


    /** Show various stats for this forest.
          @param  s       Output stream to write to
          @param  pad     Padding string, written at the start of 
                          each output line.
          @param  flags   Which stats to display, as "flags";
                          use bitwise or to combine values.
                          For example, BASIC_STATS | FOREST_STATS.
    */
    void reportStats(FILE* s, const char* pad, unsigned flags) const;


    /** Show stats about memory usage for this forest.
          @param  s     Output stream to write to
          @param  pad   Padding string, written at the start of 
                        each output line.
          @param  verb  Level of detail, between 0 (least detailed)
                        and 9 (most detailed).
    */
    // void reportMemoryUsage(FILE * s, const char* pad, int verb) const;

    /// Compute a hash for a node.
    inline unsigned hashNode(node_handle p) const {
#ifdef NODE_STORAGE_PER_LEVEL
        const node_header& node = getNode(p);
        MEDDLY_DCASSERT(node.level);
        return levels[node.level]->hashNode(node);
#else
        return nodeMan->hashNode(getNode(p));
#endif
    }


    /** Check and find the index of a single downward pointer.

          @param  node    Node we care about
          @param  down    Output:
                          The singleton downward pointer, or undefined.

          @return   If the node has only one non-zero downward pointer,
                    then return the index for that pointer.
                    Otherwise, return a negative value.
    */
    inline int getSingletonIndex(node_handle p, node_handle &down) const {
#ifdef NODE_STORAGE_PER_LEVEL
        const node_header& node = getNode(p);
        MEDDLY_DCASSERT(node.level);
        return levels[node.level]->getSingletonIndex(node.offset, down);
#else
        return nodeMan->getSingletonIndex(getNode(p).offset, down);
#endif

    }

    /** Check and get a single downward pointer.

          @param  node    Node we care about
          @param  index   Index we're trying to match

          @return   If the only non-zero downward pointer for
                    this node happens at \a index, then return the pointer.
                    Otherwise, return 0.
    */
    inline node_handle getSingletonDown(node_handle node, int index) const {
      node_handle down;
      if (getSingletonIndex(node, down) == index) return down;
      return 0;
    }

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        node_reader should be used instead.

          @param  p       Node to look at
          @param  index   Index of the pointer we want.

          @return         The downward pointer at that index.
    */
    inline node_handle getDownPtr(node_handle p, int index) const {
#ifdef NODE_STORAGE_PER_LEVEL
        const node_header& node = getNode(p);
        MEDDLY_DCASSERT(node.level);
        return levels[node.level]->getDownPtr(node.offset, index);
#else
        return nodeMan->getDownPtr(getNode(p).offset, index);
#endif
    }

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        node_reader should be used instead.

          @param  p       Node to look at
          @param  index   Index of the pointer we want.

          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    inline void getDownPtr(node_handle p, int index, int& ev, node_handle& dn) const {
#ifdef NODE_STORAGE_PER_LEVEL
        const node_header& node = getNode(p);
        MEDDLY_DCASSERT(node.level);
        levels[node.level]->getDownPtr(node.offset, index, ev, dn);
#else
        nodeMan->getDownPtr(getNode(p).offset, index, ev, dn);
#endif
    }

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        node_reader should be used instead.

          @param  p       Node to look at
          @param  index   Index of the pointer we want.

          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    inline void getDownPtr(node_handle p, int index, float& ev, node_handle& dn) const {
#ifdef NODE_STORAGE_PER_LEVEL
        const node_header& node = getNode(p);
        MEDDLY_DCASSERT(node.level);
        levels[node.level]->getDownPtr(node.offset, index, ev, dn);
#else
        nodeMan->getDownPtr(getNode(p).offset, index, ev, dn);
#endif
    }


  // ------------------------------------------------------------
  // Preferred mechanism for reading nodes
  public:
    /** Initialize a node reader.
          @param  nr      Node reader to fill.
          @param  node    The node to use.
          @param  full    true:   Use a full reader.
                          false:  Use a sparse reader.
    */
    inline void initNodeReader(node_reader &nr, node_handle node, bool full) const {
      const node_header &n = getNode(node);
#ifdef NODE_STORAGE_PER_LEVEL
      node_storage& nodeMan = levels[n.level];
#endif
      nr.resize(n.level, getLevelSize(n.level), 
        edgeBytes(), full);
      nodeMan->fillReader(n.offset, nr);
    }

    /// Allocate and initialize a node reader.
    inline node_reader* initNodeReader(node_handle node, bool full) const {
      node_reader* nr = node_reader::useReader();
      MEDDLY_DCASSERT(nr);
      initNodeReader(*nr, node, full);
      return nr;
    }

    /** Initialize a redundant node reader.
        Use for multi-terminal forests.
        For convenience.
          @param  nr      Node reader to fill.
          @param  k       Level that was skipped.
          @param  node    Downward pointer to use.
          @param  full    Use a full reader or sparse.
    */
    void initRedundantReader(node_reader &nr, int k, node_handle node, bool full)
    const;

    /// Allocate and initialize a redundant node reader.
    inline node_reader* initRedundantReader(int k, node_handle node, bool full) const {
      node_reader* nr = node_reader::useReader();
      MEDDLY_DCASSERT(nr);
      initRedundantReader(*nr, k, node, full);
      return nr;
    }

    /** Initialize a redundant node reader.
        Use for edge-valued forests, whose edge values
        require a single integer slot.
        For convenience.
          @param  nr      Node reader to fill.
          @param  k       Level that was skipped.
          @param  ev      Edge value to use.
          @param  node    Downward pointer to use.
          @param  full    Use a full reader or sparse.
    */
    void initRedundantReader(node_reader &nr, int k, int ev, node_handle node, 
      bool full) const;

    /// Allocate and initialize a redundant node reader.
    inline node_reader* initRedundantReader(int k, int ev, node_handle nd, bool full) 
    const {
      node_reader* nr = node_reader::useReader();
      MEDDLY_DCASSERT(nr);
      initRedundantReader(*nr, k, ev, nd, full);
      return nr;
    }

    /** Initialize an identity node reader.
        Use for multi-terminal forests.
        For convenience.
          @param  nr      Node reader to fill.
          @param  k       Level that was skipped.
          @param  i       Index of identity reduction
          @param  n       Downward pointer to use.
          @param  f       Use a full reader or sparse.
    */
    void initIdentityReader(node_reader &nr, int k, int i, node_handle n, bool f) 
    const;

    /** Initialize an identity node reader.
        Use for edge-valued forests, whose edge values
        require a single integer slot.
        For convenience.
          @param  nr      Node reader to fill.
          @param  k       Level that was skipped.
          @param  i       Index of identity reduction
          @param  ev      Edge value.
          @param  n       Downward pointer to use.
          @param  f       Use a full reader or sparse.
    */
    void initIdentityReader(node_reader &nr, int k, int i, int ev, node_handle n, 
      bool f) const;

    /// Allocate and initialize an identity node reader.
    inline node_reader* initIdentityReader(int k, int i, node_handle node, 
      bool full) const 
    {
      node_reader* nr = node_reader::useReader();
      MEDDLY_DCASSERT(nr);
      initIdentityReader(*nr, k, i, node, full);
      return nr;
    }

  // ------------------------------------------------------------
  // Preferred mechanism for building nodes
  public:
    inline node_builder& useNodeBuilder(int level, int tsz) {
      MEDDLY_DCASSERT(isValidLevel(level));
      MEDDLY_DCASSERT(!builders[level].lock);
      builders[level].resize(tsz);
      builders[level].lock = true;
      return builders[level];
    }
    inline node_builder& useSparseBuilder(int level, int nnz) {
      MEDDLY_DCASSERT(isValidLevel(level));
      MEDDLY_DCASSERT(!builders[level].lock);
      builders[level].resparse(nnz);
      builders[level].lock = true;
      return builders[level];
    }
    inline void doneNodeBuilder(node_builder& nb) {
      MEDDLY_DCASSERT(nb.lock);
      nb.lock = false;
    }

    /** Return a forest node equal to the one given.
        The node is constructed as necessary.
        This version should be used only for
        multi terminal forests.
          @param  in    Incoming pointer index;
                        used for identity reductions.
          @param  nb    Constructed node.

          @return       A node handle equivalent
                        to nb, taking into account
                        the forest reduction rules
                        and if a duplicate node exists.
    */
    inline node_handle createReducedNode(int in, node_builder& nb) {
      nb.computeHash();
      node_handle q = createReducedHelper(in, nb);
      MEDDLY_DCASSERT(nb.lock);
      nb.lock = false;
#ifdef TRACK_DELETIONS
      printf("Created node %d\n", q);
#endif
      return q;
    }

    /** Return a forest node equal to the one given.
        The node is constructed as necessary.
        This version should be used only for
        edge valuded forests.
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
    inline void createReducedNode(int in, node_builder& nb, T& ev, node_handle& node) {
      normalize(nb, ev);
      nb.computeHash();
      node = createReducedHelper(in, nb);
      MEDDLY_DCASSERT(nb.lock);
      nb.lock = false;
#ifdef TRACK_DELETIONS
      printf("Created node %d\n", node);
#endif
    }

    
  // ------------------------------------------------------------
  // virtual in the base class, but implemented here.
  // See meddly.h for descriptions of these methods.
  public:
    virtual void garbageCollect();
    virtual void compactMemory();
    virtual void showInfo(FILE* strm, int verbosity);

  // ------------------------------------------------------------
  // abstract virtual, must be overridden.
  // 
  public:

    /** Are the given edge values "duplicates".
        I.e., when determining if two nodes are duplicates,
        do the given edge values qualify as duplicate values.
          @param  eva     Pointer to the first edge value
          @param  evb     Pointer to the second edge value

          @return     true, iff the edge values are "equal".
          @throws     A TYPE_MISMATCH error if the forest 
                      does not store edge values.
    */
    virtual bool areEdgeValuesEqual(const void* eva, const void* evb) const;

    /** Discover duplicate nodes.
        Required for the unique table. 
          @param  node    Handle to a node.
          @param  nb      Node we've just built.
          
          @return   true, iff the nodes are duplicates.
    */
    inline bool areDuplicates(long node, const node_builder &nb) const {
      MEDDLY_DCASSERT(node>0);
      MEDDLY_DCASSERT(address);
      MEDDLY_CHECK_RANGE(1, node, 1+a_last);
      if (address[node].level != nb.getLevel()) return false;
#ifdef NODE_STORAGE_PER_LEVEL
      return levels[node.level]->areDuplicates(address[node].offset, nb);
#else
      return nodeMan->areDuplicates(address[node].offset, nb);
#endif
    };

    /** Discover duplicate nodes.
        Right now, used for sanity checks only.
          @param  node    Handle to a node.
          @param  nr      Some other node.
          
          @return   true, iff the nodes are duplicates.
    */
    inline bool areDuplicates(long node, const node_reader &nr) const {
      MEDDLY_DCASSERT(node>0);
      MEDDLY_DCASSERT(address);
      MEDDLY_CHECK_RANGE(1, node, 1+a_last);
      if (address[node].level != nr.getLevel()) return false;
#ifdef NODE_STORAGE_PER_LEVEL
      return levels[node.level]->areDuplicates(address[node].offset, nr);
#else
      return nodeMan->areDuplicates(address[node].offset, nr);
#endif
    }

    /** Is this a redundant node that can be eliminated?
        Must be implemented in derived forests
        to deal with the default edge value.
          @param  nr    Node we're trying to build.

          @return   True, if nr is a redundant node 
                          AND it should be eliminated.
    */
    virtual bool isRedundant(const node_builder &nb) const = 0;

    /** Is the specified edge an identity edge, that can be eliminated?
        Must be implemented in derived forests
        to deal with the default edge value.

          @param  nr    Node we're trying to build.
                        We know there is a single non-zero downward pointer.

          @param  i     Candidate edge (or edge index for sparse nodes).

          @return True, if nr[i] is an identity edge, and the
                        identity node should be eliminated.
    */
    virtual bool isIdentityEdge(const node_builder &nb, int i) const = 0;
        


  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                       protected  methods                       |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ------------------------------------------------------------
  // inlined setters for derived classes to use.
  protected:
    inline void setEdgeSize(char ebytes, bool hashed) {
      MEDDLY_DCASSERT(0==edge_bytes);
      edge_bytes = ebytes;
      hash_edge_values = hashed;
    };
    inline void setUnhashedSize(char ubytes) {
      MEDDLY_DCASSERT(0==unhashed_bytes);
      unhashed_bytes = ubytes;
    }
    inline void setHashedSize(char hbytes) {
      MEDDLY_DCASSERT(0==hashed_bytes);
      hashed_bytes = hbytes;
    }

  // ------------------------------------------------------------
  // inlined helpers.
  protected:
    inline bool isZombieNode(long p) const {
      MEDDLY_DCASSERT(isValidNodeIndex(p));
      MEDDLY_DCASSERT(!isTerminalNode(p));
      return (getCacheCount(p) < 0);
    }

    inline bool isActiveNode(long p) const {
      return 
      (isValidNodeIndex(p) && (isTerminalNode(p) || getNode(p).offset > 0));
    }

    inline bool isDeletedNode(long p) const {
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
      return !(isActiveNode(p) || isZombieNode(p));
    }

  // ------------------------------------------------------------
  // virtual, with default implementation.
  // Should be overridden in appropriate derived classes.
  protected:
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
    virtual void normalize(node_builder &nb, int& ev) const;

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
    virtual void normalize(node_builder &nb, float& ev) const;

    /** Show forest-specific stats.
          @param  s     Output stream to use
          @param  pad   String to display at the beginning of each line.
    */
    virtual void reportForestStats(FILE* s, const char* pad) const;

  public:
    /** Show a terminal node.
          @param  s       Stream to write to.
          @param  tnode   Handle to a terminal node.
    */
    virtual void showTerminal(FILE* s, node_handle tnode) const;

    /** Show an edge value.
          @param  s       Stream to write to.
          @param  edge    Array of all edge values.
          @param  i       Index of the one we want to display.
    */
    virtual void showEdgeValue(FILE* s, const void* edge, int i) const;

    /** Show the hashed header values.
          @param  s       Stream to write to.
          @param  hh      Array of all hashed header values.
    */
    virtual void showHashedHeader(FILE* s, const void* hh) const;

    /** Show the unhashed header values.
          @param  s       Stream to write to.
          @param  uh      Array of all unhashed header values.
    */
    virtual void showUnhashedHeader(FILE* s, const void* uh) const;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                        private  methods                        |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ------------------------------------------------------------
  // inlined helpers for this class
  private:
    inline bool isTimeToGc() const {
      return isPessimistic() 
        ? (stats.zombie_nodes > deflt.zombieTrigger)
        : (stats.orphan_nodes > deflt.orphanTrigger);
    }
    /// Returns the in-count for a node.
    inline long getInCount(node_handle p) {
#ifdef NODE_STORAGE_PER_LEVEL
      const node_header& node = getNode(p);
      return levels[node.level]->getCountOf(node.offset);
#else
      return nodeMan->getCountOf(getNode(p).offset);
#endif
    }
    /// Increment and return the in-count for a node
    inline long incInCount(node_handle p) {
#ifdef NODE_STORAGE_PER_LEVEL
      const node_header& node = getNode(p);
      return levels[node.level]->incCountOf(node.offset);
#else
      return nodeMan->incCountOf(getNode(p).offset);
#endif
    }
    /// Decrement and return the in-count for a node
    inline long decInCount(node_handle p) {
#ifdef NODE_STORAGE_PER_LEVEL
      const node_header& node = getNode(p);
      return levels[node.level]->decCountOf(node.offset);
#else
      return nodeMan->decCountOf(getNode(p).offset);
#endif
    }

    /// Returns the (modifiable) cache-count for a node
    inline int& cacheCount(node_handle p) {
      MEDDLY_DCASSERT(isValidNodeIndex(p));
      MEDDLY_DCASSERT(!isTerminalNode(p));
      return address[p].cache_count;
    }

    /// Returns the cache-count for a node
    inline int getCacheCount(node_handle p) const {
      MEDDLY_DCASSERT(isValidNodeIndex(p));
      MEDDLY_DCASSERT(!isTerminalNode(p));
      return address[p].cache_count;
    }


    /** Change the location of a node.
        Used by node_storage during compaction.
        Should not be called by anything else.
          @param  node        Node we're moving
          @param  old_addr    Current address of node, for sanity check
          @param  new_addr    Where we're moving the node to
    */
    inline void moveNodeOffset(node_handle node, node_address old_addr, node_address new_addr) 
    {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(old_addr == address[node].offset);
      address[node].offset = new_addr;
    }
    friend class MEDDLY::node_storage;


  // ------------------------------------------------------------
  // helpers for this class
  private:
    void handleNewOrphanNode(node_handle node);   
    void deleteNode(node_handle p);
    void zombifyNode(node_handle p);

    /// Determine a node handle that we can use.
    node_handle getFreeNodeHandle();

    /// Release a node handle back to the free pool.
    void recycleNodeHandle(node_handle p);

    /** Apply reduction rule to the temporary node and finalize it. 
        Once a node is reduced, its contents cannot be modified.
          @param  in    Incoming index, used only for identity reduction;
                        Or -1.
          @param  nb    Array of downward pointers.
          @return       Handle to a node that encodes the same thing.
    */
    node_handle createReducedHelper(int in, const node_builder &nb);

    // Sanity check; used in development code.
    void validateDownPointers(const node_builder &nb) const;

    /// Increase the number of node handles.
    void expandHandleList();

    /// Decrease the number of node handles.
    void shrinkHandleList();

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
    bool terminalNodesAreStale;

#ifdef NODE_STORAGE_PER_LEVEL
    /// Level data. Each level maintains its own data array and hole grid.
    /// The array is shifted, so we can use level[k] with negative k.
    node_storage** levels;
#else
    node_storage* nodeMan;
#endif

  private:
    // Keep a node_builder for each level.
    node_builder* raw_builders;
    node_builder* builders;

    // Garbage collection in progress
    bool performing_gc;

    // memory for validating incounts
    node_handle* in_validate;
    int  in_val_size;

#ifdef NODE_STORAGE_PER_LEVEL
    // Raw level data
    node_storage** raw_levels;
#endif


    /// address info for nodes
    node_header *address;
    /// Size of address/next array.
    node_handle a_size;
    /// Last used address.
    node_handle a_last;
#ifdef ONE_NODE_FREELIST
    /// Pointer to unused address list.
    node_handle a_unused;
#else
    /// Pointer to unused address lists, based on size
    node_handle a_unused[8];  // number of bytes per handle
    /// Lowest non-empty address list
    char a_lowest_index;
#endif
    /// Next time we shink the address list.
    node_handle a_next_shrink;


    /// Number of bytes for an edge
    char edge_bytes;
    /// Are edge values hashed?
    bool hash_edge_values;
    /// Number of bytes of unhashed header
    char unhashed_bytes;
    /// Number of bytes of hashed header
    char hashed_bytes;

  private:
    class nodecounter;
};
// end of expert_forest class.


// ******************************************************************
// *                                                                *
// *                expert_forest::nodecounter class                *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_forest::nodecounter : public edge_visitor {
    expert_forest* parent;
    int* counts;
  public:
    nodecounter(expert_forest*p, int* c);
    virtual ~nodecounter();
    virtual void visit(dd_edge &e);
};





// ******************************************************************
// *                                                                *
// *                  inlined node_storage methods                  *
// *                                                                *
// ******************************************************************


inline void MEDDLY::node_storage
::moveNodeOffset(node_handle node, node_address old_addr, node_address new_addr)
{
  MEDDLY_DCASSERT(parent);
  parent->moveNodeOffset(node, old_addr, new_addr);
}







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

    friend void MEDDLY::initialize(const settings &);
    friend void MEDDLY::cleanup();
  public:
    opname(const char* n);
    virtual ~opname();

    inline int getIndex() const         { return index; }
    inline const char* getName() const  { return name; }
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
// *                     numerical_opname class                     *
// *                                                                *
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
// *                                                                *
// *                      compute_table  class                      *
// *                                                                *
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
          node_handle* data;
          bool killData;
          node_handle* key_data;
          const operation* op;
          /// used only for range checking during "development".
          int keyLength;  
        public:
          search_key();
          ~search_key();
          inline node_handle& key(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<keyLength);
#endif
            return key_data[i]; 
          }
          inline const node_handle* rawData() const { return data; }
          inline int dataLength() const { return hashLength; }
          inline const operation* getOp() const { return op; }
      };

      class temp_entry {
          friend class MEDDLY::base_table;
          int handle;
          int hashLength;
          node_handle* entry;
          node_handle* key_entry;
          node_handle* res_entry;
          // The remaining entries are used only in development code
          int keyLength;
          int resLength;
        public:
          inline node_handle& key(int i) { 
#ifdef DEVELOPMENT_CODE
            assert(i>=0);
            assert(i<keyLength);
#endif
            return key_entry[i]; 
          }
          inline node_handle& result(int i) { 
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
          inline const node_handle* readEntry(int off) const { return entry+off; }
          inline int readHandle() const { return handle; }
          inline int readLength() const { return hashLength; }
          inline node_handle& data(int i) {
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
      virtual const node_handle* find(const search_key &key) = 0;

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
        ?   f->getNodeDeletion() == forest::policies::PESSIMISTIC_DELETION
        :   false;  // shouldn't be possible, so we'll do what's fastest.
    }

    void markForDeletion();

    friend class forest;
    friend void MEDDLY::initialize(const settings &);
    friend void MEDDLY::destroyOpInternal(operation* op);
    friend void MEDDLY::cleanup();

    inline void registerInForest(forest* f) { 
      if (f) f->registerOperation(this); 
    }

    inline void unregisterInForest(forest* f) { 
      if (f) f->unregisterOperation(this); 
    }

  private:
    // should ONLY be called during library cleanup.
    static void destroyAllOps();

  public:
    inline bool isMarkedForDeletion() const { return is_marked_for_deletion; }

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
    inline bool isEntryStale(const node_handle* data) {
      return (is_marked_for_deletion || isStaleEntry(data));
    }

  protected:
    virtual bool isStaleEntry(const node_handle* entry) = 0;

  public:
    /// Removes the cache entry (in entryData[]) by informing the
    /// applicable forests that the nodes in this entry are being removed
    /// from the cache
    virtual void discardEntry(const node_handle* entryData) = 0;

    /// Prints a string representation of this cache entry on strm (stream).
    virtual void showEntry(FILE* strm, const node_handle *entryData) const = 0;

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
    virtual node_handle compute(node_handle a, node_handle b);
    /// Low-level compute at level k on nodes a and b, return result.
    virtual node_handle compute(int k, node_handle a, node_handle b);

  protected:
    inline void operationCommutes() {
      can_commute = (arg1F == arg2F);
    }
};

// ******************************************************************
// *                                                                *
// *                   numerical_operation  class                   *
// *                                                                *
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
// *                                                                *
// *                      op_initializer class                      *
// *                                                                *
// ******************************************************************

/** Preferred mechanism for users to initialize their own operations.
    Derive a class from this one, provide the \a execute method.
    Implementation in ops.cc
*/
class MEDDLY::op_initializer {
  unsigned refcount;
  op_initializer* before;
public:
  /// Constructor.
  ///   @param  bef   initializer(s) to execute before this one.
  op_initializer(op_initializer* bef);

  void initChain(const settings &s);
  void cleanupChain();

  inline static void recycle(op_initializer *I) {
    if (0==I) return;
    MEDDLY_DCASSERT(I->refcount);
    I->refcount--;
    if (0==I->refcount) delete I;
  }

  inline static op_initializer* copy(op_initializer *I) { 
    if (I) I->refcount++;
    return I;
  }

protected:
  virtual ~op_initializer();

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

inline void MEDDLY::node_builder::setUH(const void* x) 
{
  memcpy(raw_uh(), x, parent->unhashedHeaderBytes());
}

inline void MEDDLY::node_builder::getUH(void* x) const 
{
  memcpy(x, raw_uh(), parent->unhashedHeaderBytes());
}

inline void MEDDLY::node_builder::setHH(const void* x) 
{
  memcpy(raw_hh(), x, parent->hashedHeaderBytes());
}

inline void MEDDLY::node_builder::getHH(void* x) const 
{
  memcpy(x, raw_hh(), parent->hashedHeaderBytes());
}

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
float MEDDLY::toFloat(int a) {
  union { int i; float f; } n = {a};
  return n.f;
}


inline
int MEDDLY::toInt(float a) {
  union { float f; int i; } n = {a};
  return n.i;
}

#endif
