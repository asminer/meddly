
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
// #define DEBUG_MDD_H
// #define TRACK_DELETIONS
// #define TRACK_CACHECOUNT
#endif

// #define TRACK_DELETIONS

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
//
// #define SAVE_HASHES
 
namespace MEDDLY {

  // Functions for reinterpreting an int to a float and vice-versa
  float   toFloat (int a);
  int     toInt   (float a);
  // float*  toFloat (int* a);

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
        inline int d(int n) const {
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

        // For debugging
        // void dump(FILE*) const;

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
        int* down;
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
        // void resize(const expert_forest* p, int k, int ns, bool full);
        void resize(int k, int ns, char eb, bool full);

        friend class expert_forest; 
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
*/
class MEDDLY::node_builder {
        const expert_forest* parent;
        int* extra_hashed;
        int* extra_unhashed;
        int* down;
        int* indexes;
        void* edge;
        int level;  // level of this node.
        int hhsize; // extra_hashed size.
        int uhsize; // extra_unhashed size.
        int size;
        int alloc;
        // TBD: need "extra" header stuff
        unsigned h;
        char edge_bytes;
#ifdef DEVELOPMENT_CODE
        bool has_hash;
#endif
        bool is_sparse;
        bool hashEdgeValues;
    private:
        // These cheat; let's see if any compiler complains.
        inline int& raw_hh(int i) const {
          MEDDLY_CHECK_RANGE(0, i, hhsize);
          MEDDLY_DCASSERT(extra_hashed);
          return extra_hashed[i];
        }
        inline int& raw_uh(int i) const {
          MEDDLY_CHECK_RANGE(0, i, uhsize);
          MEDDLY_DCASSERT(extra_unhashed);
          return extra_unhashed[i];
        }
        inline int& raw_d(int i) const {
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
        inline int& hh(int i)         { return raw_hh(i); }
        inline int  hh(int i) const   { return raw_hh(i); }
        inline int& uh(int i)         { return raw_uh(i); }
        inline int  uh(int i) const   { return raw_uh(i); }
        inline int&  d(int i)         { return  raw_d(i); }
        inline int   d(int i) const   { return  raw_d(i); }
        inline int&  i(int i)         { return  raw_i(i); }
        inline int   i(int i) const   { return  raw_i(i); }
        inline int& ei(int i)         { return raw_ei(i); }
        inline int  ei(int i) const   { return raw_ei(i); }
        inline float& ef(int i)       { return raw_ef(i); }
        inline float  ef(int i) const { return raw_ef(i); }
        // edge setting - for templates
        inline void setEdge(int i, int ev) {
          ei(i) = ev;
        }
        inline void setEdge(int i, float ev) {
          ef(i) = ev;
        }
    public: // for unique table
        inline unsigned hash() const {
#ifdef DEVELOPMENT_CODE
          MEDDLY_DCASSERT(has_hash);
#endif
          return h;
        }
        // bool equals(int p) const;
        void computeHash();

    public: // helpers for copying
        void copyIntoFull(int* down, int N) const;
        void copyIntoFull(int* down, void* edge, int N) const;
        void copyIntoSparse(int* down, int* indexes, int Z) const;
        void copyIntoSparse(int* down, int* indexes, void* edge, int N) const;

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
    /** Node level
        If the node is active, this indicates node level.
    */
    int level;

    /** Offset to node's data in the corresponding node storage structure.
        If the node is active, this is the offset (>0) in the data array.
        If the node is deleted, this is -next deleted node
        (part of the unused address list).
    */
    int offset;

    /** Cache count
        The number of cache entries that refer to this node (excl. unique
        table). If this node is a zombie, cache_count is negative.
    */
    int cache_count;

#ifdef SAVE_HASH
    /// Remember the hash for speed.
    unsigned hash;
#endif

    // Handy functions, in case our internal storage changes.

    inline bool isActive() const  { return offset > 0; }
    inline bool isZombie() const  { return cache_count < 0; }

    inline bool isDeleted() const { return 0 == level; }
    inline void setDeleted()      { level = 0; }
      
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

/** Actual node storage in a forest, by level.
    Implemented in forest.cc
    Details of node storage are left to the derived forests.
    However, every active node is stored in the following format.
      

      common   {  slot[0] : incoming count, >= 0.
      header --{  slot[1] : next pointer in unique table or special value.
               {  slot[2] : size.  >=0 for full storage, <0 for sparse.

      unhashed    {       : slots used for any extra information
      header    --{       : as needed on a forest-by-forest basis.
      (optional)  {       : Info does NOT affect node uniqueness.

      hashed      {       : slots used for any extra information
      header    --{       : as needed on a forest-by-forest basis.
      (optional)  {       : Info DOES affect node uniqueness.

                  {       : Downward pointers.
                  {       : If full storage, there are size pointers
      down -------{       : and entry i gives downward pointer i.
                  {       : If sparse storage, there are -size pointers
                  {       : and entry i gives a pointer but the index
                  {       : corresponds to index[i], below.
          
                  {       : Index entries.
                  {       : Unused for full storage.
      index ------{       : If sparse storage, entry i gives the
      (sparse)    {       : index for outgoing edge i, and there are
                  {       : -size entries.

                  {       : Edge values.
      edge        {       : If full storage, there are size * edgeSize
      values -----{       : slots; otherwise there are -size * edgeSize
                  {       : slots.  Derived forests are responsible
                  {       : for packing information into these slots.

                { -padlen : Any node padding to allow for future
                {         : expansion, or for memory management purposes
                {         : (e.g., memory hold is larger than requested).
      padding --{         : padlen is number of padded slots.  If the
                {         : first entry after the node proper is negative,
                {         : then it specifies the number of padding slots;
                {         : otherwise, there is no padding.

      tail    --{ slot[L] : The last slot gives the forest node number,
                            guaranteed to be non-negative.


    When nodes are deleted, the memory slots are marked as a "hole",
    using the following format.

          slot[0] : -numslots, the number of slots in the hole
            .
            .
            .
          slot[L] : -numslots, with L = numslots-1
        
    The first few slots of the hole are used for a hole management
    data structure, whose details are important only for implementation
    of this class.
    Note that a hole is guaranteed to be at least 5 slots long
    (assuming a node of size 0, with no extra header info, is impossible).
*/
class MEDDLY::node_storage {
      /// Growth parameter
      static const int add_size = 1024;

      /// Special values
      static const int non_index_hole = -2;
      static const int temp_node_value = -5;

      // header indexes
      static const int count_index = 0;
      static const int next_index = 1;    
      static const int size_index = 2;

      // Counts for extra slots
      static const int commonHeaderLength = 3;
      static const int commonTailLength = 1;
      static const int commonExtra = commonHeaderLength + commonTailLength;

      // Parent forest.
      expert_forest* parent;

      /// data array
      int* data;
      /// shifted data array, for "node start"
      int* data_down;
      /// Size of data array.
      int size;
      /// Last used data slot.  Also total number of ints "allocated"
      int last;

  // Holes grid info
  private:
      /// Pointer to top of holes grid
      int holes_top;
      /// Pointer to bottom of holes grid
      int holes_bottom;
      /// Total ints in holes
      int hole_slots;

  // performance stats
  public:
      /// Largest traversed height of holes grid
      int max_hole_chain;
      /// Number of zombie nodes
      int zombie_nodes;
      /// Total number of compactions
      // int num_compactions;

  // header sizes; vary by forest.
  public:
      /// Size of each outgoing edge's value (can be 0).
      char edgeSize;
      /// Size of extra unhashed data (typically 0).
      char unhashedHeader;
      /// Size of extra hashed data (typically 0).
      char hashedHeader;

      /// Mark for compaction
      bool compactLevel;

  // --------------------------------------------------------
  // |  Public interface.
  public:
      node_storage();
      ~node_storage();

      /// Bind this to a particular level in a particular forest.
      void init(expert_forest* p, int k);

      /** Allocate enough slots to store a node with given size.
          Also, stores the node size in the node.
            @param  sz      negative for sparse storage, otherwise full.
            @param  tail    Node id
            @param  clear   Should the node be zeroed.
            @return         Offset in the data array.
      */
      int allocNode(int sz, int tail, bool clear);

      /// Unlink the downward pointers.
      void unlinkDown(int addr);

      /// Compact this level.  (Rearrange, to remove all holes.)
      void compact(node_header* address);

      /// For debugging.
      void dumpInternal(FILE* s) const;

      /// For performance stats.
      void addToChainCounts(std::map<int, int> &chainLengths) const;

  // --------------------------------------------------------
  // |  helpers for inlines.
  private:
      inline int& rawCountOf(int addr) const {
          MEDDLY_DCASSERT(data);
          MEDDLY_CHECK_RANGE(1, addr, last+1);
          return data[addr + count_index];
      }
      inline int& rawNextOf(int addr) const { 
          MEDDLY_DCASSERT(data);
          MEDDLY_CHECK_RANGE(1, addr, last+1);
          return data[addr + next_index];
      }
      inline int& rawSizeOf(int addr) const {
          MEDDLY_DCASSERT(data);
          MEDDLY_CHECK_RANGE(1, addr, last+1);
          return data[addr + size_index];
      }
      inline int* UH(int addr) const {
          MEDDLY_DCASSERT(data);
          return data + addr + commonHeaderLength;
      }
      inline int* HH(int addr) const {
          MEDDLY_DCASSERT(data);
          return data + addr + commonHeaderLength + unhashedHeader;
      }
      inline int* FD(int addr) const {
          MEDDLY_DCASSERT(data);
          MEDDLY_DCASSERT(isFull(addr));
          return data_down + addr;
      }
      inline int* FE(int addr) const {
          return FD(addr) + rawSizeOf(addr);
      }
      inline int* SD(int addr) const {
          MEDDLY_DCASSERT(data);
          MEDDLY_DCASSERT(isSparse(addr));
          return data_down + addr;
      }
      inline int* SI(int addr) const {
          MEDDLY_DCASSERT(isSparse(addr));
          return SD(addr) - rawSizeOf(addr);  // sparse size is negative
      }
      inline int* SE(int addr) const {
          MEDDLY_DCASSERT(isSparse(addr));
          return SD(addr) - 2*rawSizeOf(addr);  // sparse size is negative
      }

  // --------------------------------------------------------
  // |  inlines.
  public:
      // TBD:
      inline int trailerSize() const {
          return unhashedHeader + hashedHeader;
      }

      /// Recycle a node stored at the given offset.
      inline void recycleNode(int off) {
          makeHole(off, activeNodeActualSlots(off));
      }

      inline bool needsCompaction(int compactPercent) const {
          MEDDLY_DCASSERT(parent);
          if (hole_slots <= 100)  return false;
          if (hole_slots > 10000) return true;
          return hole_slots * 100 > last * compactPercent;
      }

      /// How many slots would be required for a node with given size.
      ///   @param  sz  negative for sparse storage, otherwise full.
      inline int slotsForNode(int sz) const {
          int edges = (sz<0) ? (2+edgeSize) * (-sz) : (1+edgeSize) * sz;
          return commonExtra + unhashedHeader + hashedHeader + edges;
      }

      inline int getHoleSlots() const { return hole_slots; }

      inline int& countOf(int addr)       { return rawCountOf(addr); }
      inline int  countOf(int addr) const { return rawCountOf(addr); }

      inline int& nextOf(int addr)       { return rawNextOf(addr); }
      inline int  nextOf(int addr) const { return rawNextOf(addr); }

      inline bool isReduced(int addr) const {
          MEDDLY_DCASSERT(data);
          MEDDLY_CHECK_RANGE(1, addr, last+1);
          return data[addr + next_index] >= 0;
      }
      inline bool isTempNode(int addr) const {
          return temp_node_value == nextOf(addr);
      }
      inline void setTempNode(int addr) { 
          nextOf(addr) = temp_node_value; 
      }

      inline int& sizeOf(int addr)              { return rawSizeOf(addr); }
      inline int  sizeOf(int addr) const        { return rawSizeOf(addr); }

      inline int fullSizeOf(int addr) const     { return rawSizeOf(addr); }
      inline int sparseSizeOf(int addr) const   { return -rawSizeOf(addr); }

      inline int isFull(int addr) const       { return sizeOf(addr) >= 0; }
      inline int isSparse(int addr) const     { return sizeOf(addr) < 0;  }

      inline       int* unhashedHeaderOf(int addr)       { return UH(addr); }
      inline const int* unhashedHeaderOf(int addr) const { return UH(addr); }

      inline       int* hashedHeaderOf(int addr)       { return HH(addr); }
      inline const int* hashedHeaderOf(int addr) const { return HH(addr); }

      inline       int* fullDownOf(int addr)       { return FD(addr); }
      inline const int* fullDownOf(int addr) const { return FD(addr); }

      inline       int* fullEdgeOf(int addr)       { return FE(addr); }
      inline const int* fullEdgeOf(int addr) const { return FE(addr); }

      inline       int* sparseDownOf(int addr)       { return SD(addr); }
      inline const int* sparseDownOf(int addr) const { return SD(addr); }

      inline       int* sparseIndexesOf(int addr)       { return SI(addr); }
      inline const int* sparseIndexesOf(int addr) const { return SI(addr); }

      inline       int* sparseEdgeOf(int addr)       { return SE(addr); }
      inline const int* sparseEdgeOf(int addr) const { return SE(addr); }

  // --------------------------------------------------------
  // |  Hole management helpers.
  private:
      /// Find actual number of slots used for this active node.
      inline int activeNodeActualSlots(int off) const {
          MEDDLY_DCASSERT(data);
          MEDDLY_CHECK_RANGE(1, off, last+1);
          MEDDLY_DCASSERT(data[off + count_index] >= 0);
          int end = off + slotsForNode(sizeOf(off))-1;
          // account for any padding
          if (data[end] < 0) {
            end -= data[end];
          }
          return end - off + 1;
      }

      // returns offset to the hole found in level
      int getHole(int slots);

      // makes a hole of size == slots, at the specified offset
      void makeHole(int p_offset, int slots);

      // add a hole to the hole grid
      void gridInsert(int p_offset);

      inline bool isHoleNonIndex(int p_offset) {
          return (data[p_offset + 1] == non_index_hole);
      }

      // remove a non-index hole from the hole grid
      void midRemove(int p_offset);

      // remove an index hole from the hole grid
      void indexRemove(int p_offset);

      // resize the data array.
      void resize(int new_slots);
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

// TBD: MOVE THESE TO MT CLASS?

    /// Get the integer value represented by this terminal node.
    static inline bool getBoolean(int terminalNode) {
      MEDDLY_DCASSERT(terminalNode <= 0);
      return terminalNode;
    }

    /// Get the terminal node representing this boolean value.
    static inline int getTerminalNode(bool booleanValue) {
      return booleanValue ? -1 : 0;
    }

    /// Get the integer value represented by this terminal node.
    static inline int getInteger(int terminalNode) {
      MEDDLY_DCASSERT(terminalNode <= 0);
      // set 32nd bit based on 31st bit.  << gets rid of MSB; >> sign extends.
      return terminalNode << 1 >> 1;
    }

    /// Get the terminal node representing this integer value.
    static inline int getTerminalNode(int integerValue) {
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
    static inline int getTerminalNode(float a) {
      if (0.0 == a) return 0;
      int node;
      memcpy(&node, &a, sizeof(int));
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
    inline long getLevelNumber(const node_storage* L) const {
      return long(L - levels);
    }
    /// Is this a terminal node?
    inline static bool isTerminalNode(int p) {
      return (p < 1);
    }
    /// Sanity check: is this a valid nonterminal node index.
    inline bool isValidNonterminalIndex(int node) const {
      return (node>0) && (node <= a_last);
    }
    /// Sanity check: is this a valid node index.
    inline bool isValidNodeIndex(int node) const {
      return node <= a_last;
    }
    inline int getLastNode() const {
      return a_last;
    }
    /// Get data for a given nonterminal node.
    inline const node_header& getNode(int p) const {
      MEDDLY_DCASSERT(address);
      MEDDLY_CHECK_RANGE(1, p, 1+a_last);
      return address[p];
    }
    /** Get the node's level as an integer.
        Negative values are used for primed levels.
    */
    inline int getNodeLevel(int p) const {
      if (isTerminalNode(p)) return 0;
      MEDDLY_DCASSERT(address);
      MEDDLY_CHECK_RANGE(1, p, 1+a_last);
      return address[p].level;
    }
    inline bool isPrimedNode(int p) const {
      return getNodeLevel(p) < 0;
    }
    inline bool isUnprimedNode(int p) const {
      return getNodeLevel(p) > 0;
    }

    /** Get the node's height.
        For convenience.  Height is the total
        number of levels until terminal nodes.
    */
    inline int getNodeHeight(int p) const {
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
    inline bool isIndexSet(int node) const {
      if (!isEVPlus()) return false;
      if (isTerminalNode(node)) return node;
      // yes iff the unhashed extra header is non-zero.
      const node_header& nd = getNode(node);
      const node_storage& ld = levels[nd.level];
      return (ld.unhashedHeaderOf(nd.offset))[0] > 0;
    }

    /// Get the cardinality of an Index Set.
    inline int getIndexSetCardinality(int node) const {
      MEDDLY_DCASSERT(isEVPlus());
      if (isTerminalNode(node)) return (node != 0) ? 1 : 0;
      // yes iff the unhashed extra header is non-zero.
      const node_header& nd = getNode(node);
      const node_storage& ld = levels[nd.level];
      MEDDLY_DCASSERT((ld.unhashedHeaderOf(nd.offset))[0] > 0);
      return (ld.unhashedHeaderOf(nd.offset))[0];
    }

    // --------------------------------------------------
    // Used by the unique table
    // --------------------------------------------------
    inline int getNext(int p) const {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(levels);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
      MEDDLY_DCASSERT(address[p].level);
      return levels[address[p].level].nextOf(address[p].offset);
    }
    inline void setNext(int p, int n) {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(levels);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
      MEDDLY_DCASSERT(address[p].level);
      levels[address[p].level].nextOf(address[p].offset) = n;
    }
    inline unsigned hash(int p) const {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
#ifdef SAVE_HASH
      return address[p].hash;
#else
      return hashNode(p);
#endif
    }

    // --------------------------------------------------
    // Managing reference counts
    // --------------------------------------------------

    /** Increase the link count to this node. Call this when another node is
        made to point to this node.
          @return p, for convenience.
    */
    inline int linkNode(int p) {
        MEDDLY_DCASSERT(isActiveNode(p));
        if (isTerminalNode(p)) return p;
        MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));

        int& count = getInCount(p);
        count++;

        if (1 == count) {
          // Reclaim an orphan node
          stats.reclaimed_nodes++;
          stats.orphan_nodes--;
        }
#ifdef TRACK_DELETIONS
        fprintf(stdout, "\t+Node %d count now %d\n", p, count);
        fflush(stdout);
#endif
        return p;
    }

    /** Decrease the link count to this node. If link count reduces to 0, this
        node may get marked for deletion. Call this when another node releases
        its connection to this node.
    */
    inline void unlinkNode(int p) {
        MEDDLY_DCASSERT(isActiveNode(p));
        if (isTerminalNode(p)) return;
        MEDDLY_DCASSERT(!isPessimistic() || !isZombieNode(p));
        MEDDLY_DCASSERT(getInCount(p) > 0);

        int& count = getInCount(p);
        count--;
  
#ifdef TRACK_DELETIONS
        fprintf(stdout, "\t-Node %d count now %d\n", p, count);
        fflush(stdout);
#endif
        if (count) return;

        handleNewOrphanNode(p);
    }

    // --------------------------------------------------
    // Managing cache entries
    // --------------------------------------------------

    /** Returns the cache-count for a node. This indicates the number of
        compute cache entries that link to this node.
        This version returns a reference to the cache-count. Therefore, the
        cache-count of the node can be modified by modifying the reference.
    */
    inline int& getCacheCount(int p) {
      MEDDLY_DCASSERT(isValidNodeIndex(p));
      MEDDLY_DCASSERT(!isTerminalNode(p));
      return address[p].cache_count;
    }

    /** Returns the cache-count for a node. This indicates the number of
        compute cache entries that link to this node.
        This version returns a copy of the count.
    */
    inline int getCacheCount(int p) const {
      MEDDLY_DCASSERT(isValidNodeIndex(p));
      MEDDLY_DCASSERT(!isTerminalNode(p));
      return address[p].cache_count;
    }

    /** Increase the cache count for this node. Call this whenever this node
        is added to a cache. 
          @param  p     Node we care about.
          @return p, for convenience.
    */
    inline int cacheNode(int p) {
      MEDDLY_DCASSERT(isActiveNode(p));
      if (isTerminalNode(p)) return p;
      getCacheCount(p)++;
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
    inline void uncacheNode(int p) {
      if (isTerminalNode(p)) return;
      MEDDLY_DCASSERT(isActiveNode(p) ||
          (!isActiveNode(p) && isPessimistic() && isZombieNode(p)));
      int& cc = getCacheCount(p);
      if (isPessimistic() && isZombieNode(p)) {
        // special case: we store the negative of the count.
        MEDDLY_DCASSERT(cc < 0);
        cc++;
        if (0 == cc) {
          freeZombieNode(p);
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
        deleteOrphanNode(p);
      }
    }


    /// A node can be discarded once it goes stale. Whether a node is
    /// considered stale depends on the forest's deletion policy.
    /// Optimistic deletion: A node is said to be stale only when both the
    ///   in-count and cache-count are zero.
    /// Pessimistic deletion: A node is said to be stale when the in-count
    ///  is zero regardless of the cache-count.
    inline bool isStale(int node) const {
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
    void dumpInternalLevel(FILE *s, int k) const; 
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
    int* markNodesInSubgraph(int root, bool sort) const;

    /** Count and return the number of non-terminal nodes
        in the subgraph below the given node.
    */
    int getNodeCount(int node) const;

    /** Count and return the number of edges
        in the subgraph below the given node.
    */
    int getEdgeCount(int node, bool countZeroes) const;

    /// Display the contents of a single node.
    void showNode(FILE* s, int node, int verbose = 0) const;

    /// Show all the nodes in the subgraph below the given node.
    void showNodeGraph(FILE* s, int node) const;

    /** Show stats about memory usage for this forest.
          @param  s     Output stream to write to
          @param  pad   Padding string, written at the start of 
                        each output line.
          @param  verb  Level of detail, between 0 (least detailed)
                        and 9 (most detailed).
    */
    void reportMemoryUsage(FILE * s, const char* pad, int verb) const;

    /// Compute a hash for a node.
    unsigned hashNode(int p) const;

    /** Check and find the index of a single downward pointer.

          @param  node    Node we care about
          @param  down    Output:
                          The singleton downward pointer, or undefined.

          @return   If the node has only one non-zero downward pointer,
                    then return the index for that pointer.
                    Otherwise, return a negative value.
    */
    int getSingletonIndex(int node, int &down) const;

    /** Check and get a single downward pointer.

          @param  node    Node we care about
          @param  index   Index we're trying to match

          @return   If the only non-zero downward pointer for
                    this node happens at \a index, then return the pointer.
                    Otherwise, return 0.
    */
    inline int getSingletonDown(int node, int index) const {
      int down;
      if (getSingletonIndex(node, down) == index) return down;
      return 0;
    }

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        node_reader should be used instead.

          @param  node    Node to look at
          @param  index   Index of the pointer we want.

          @return         The downward pointer at that index.
    */
    int getDownPtr(int node, int index) const;

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        node_reader should be used instead.

          @param  node    Node to look at
          @param  index   Index of the pointer we want.

          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    void getDownPtr(int node, int index, int& ev, int& dn) const;

    /** For a given node, get a specified downward pointer.

        This is designed to be used for one or two indexes only.
        For reading all or several downward pointers, a
        node_reader should be used instead.

          @param  node    Node to look at
          @param  index   Index of the pointer we want.

          @param  ev      Output: edge value at that index.
          @param  dn      Output: downward pointer at that index.
    */
    void getDownPtr(int node, int index, float& ev, int& dn) const;


  // ------------------------------------------------------------
  // Preferred mechanism for reading nodes
  public:
    /** Initialize a node reader.
          @param  nr      Node reader to fill.
          @param  node    The node to use.
          @param  full    true:   Use a full reader.
                          false:  Use a sparse reader.
    */
    void initNodeReader(node_reader &nr, int node, bool full) const;

    /// Allocate and initialize a node reader.
    inline node_reader* initNodeReader(int node, bool full) const {
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
    void initRedundantReader(node_reader &nr, int k, int node, bool full) const;

    /// Allocate and initialize a redundant node reader.
    inline node_reader* initRedundantReader(int k, int node, bool full) const {
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
    void initRedundantReader(node_reader &nr, int k, int ev, int node, 
      bool full) const;

    /// Allocate and initialize a redundant node reader.
    inline node_reader* initRedundantReader(int k, int ev, int nd, bool full) const {
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
    void initIdentityReader(node_reader &nr, int k, int i, int n, bool f) const;

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
    void initIdentityReader(node_reader &nr, int k, int i, int ev, int n, 
      bool f) const;

    /// Allocate and initialize an identity node reader.
    inline node_reader* initIdentityReader(int k, int i, int node, 
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
    inline int createReducedNode(int in, node_builder& nb) {
      nb.computeHash();
      int q = createReducedHelper(in, nb);
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
    inline void createReducedNode(int in, node_builder& nb, T& ev, int& node) {
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
    /// Number of integer slots for an edge at level k.
    virtual char edgeSize(int k) const = 0;
    /// Extra unhashed node slots at level k.
    virtual char unhashedHeaderSize(int k) const = 0;
    /// Extra hashed node slots at level k.
    virtual char hashedHeaderSize(int k) const = 0;
    /// Are edge values included when computing the hash, at level k.
    virtual bool areEdgeValuesHashed(int k) const = 0;

    /** Discover duplicate nodes.
        Required for the unique table. 
        Default behavior is to throw an exception.
          @param  node    Handle to a node.
          @param  nb      Node we've just built.
          
          @return   true, iff the nodes are duplicates.
    */
    virtual bool areDuplicates(int node, const node_builder &nb) const = 0;

    /** Discover duplicate nodes.
        Right now, used for sanity checks only.
        Default behavior is to throw an exception.
          @param  node    Handle to a node.
          @param  nr      Some other node.
          
          @return   true, iff the nodes are duplicates.
    */
    virtual bool areDuplicates(int node, const node_reader &nr) const = 0;

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
  // inlined helpers.
  protected:
    inline bool isZombieNode(int p) const {
      MEDDLY_DCASSERT(isValidNodeIndex(p));
      MEDDLY_DCASSERT(!isTerminalNode(p));
      return (getCacheCount(p) < 0);
    }

    inline bool isActiveNode(int p) const {
      return 
      (isValidNodeIndex(p) && (isTerminalNode(p) || getNode(p).offset > 0));
    }

    inline bool isDeletedNode(int p) const {
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
      return !(isActiveNode(p) || isZombieNode(p));
    }

  // ----------------------------------------------------------------- 
  // | 
  // |  Management of node headers.
  // | 
  protected:
    inline int getFreeNodeHandle() {
      MEDDLY_DCASSERT(address);
      stats.incActive(1);
      stats.incMemUsed(sizeof(node_header));
      while (a_unused > a_last) {
        a_unused = address[a_unused].getNextDeleted();
      }
      if (a_unused) {     // get a recycled one
        int p = a_unused;
        MEDDLY_DCASSERT(address[p].isDeleted());
        a_unused = address[p].getNextDeleted();
        return p;
      }
      a_last++;
      if (a_last >= a_size) {
        expandHandleList();
      }
      MEDDLY_DCASSERT(a_last < a_size);
      return a_last;
    }

    inline void freeActiveNode(int p) {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
      MEDDLY_DCASSERT(address[p].isActive());
      stats.decActive(1);
      recycleNodeHandle(p);
    }

    inline void freeZombieNode(int p) {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(isValidNonterminalIndex(p));
      MEDDLY_DCASSERT(address[p].level);
      stats.zombie_nodes--;
      levels[address[p].level].zombie_nodes--;
      recycleNodeHandle(p);
    }

    inline void recycleNodeHandle(int p) {
      MEDDLY_DCASSERT(address);
      MEDDLY_DCASSERT(p>0);
      MEDDLY_DCASSERT(0==address[p].cache_count);
      stats.decMemUsed(sizeof(node_header));
      address[p].setDeleted();
      address[p].setNextDeleted(a_unused);
      a_unused = p;
      if (p == a_last) {
        while (a_last && address[a_last].isDeleted()) {
          a_last--;
        }
        if (a_last < a_next_shrink) shrinkHandleList();
      }
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

    /** Show a terminal node.
          @param  s       Stream to write to.
          @param  tnode   Handle to a terminal node.
    */
    virtual void showTerminal(FILE* s, int tnode) const;

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
    virtual void showHashedHeader(FILE* s, const int* hh) const;

    /** Show the unhashed header values.
          @param  s       Stream to write to.
          @param  uh      Array of all unhashed header values.
    */
    virtual void showUnhashedHeader(FILE* s, const int* uh) const;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // |                                                                |
  // |                        private  methods                        |
  // |                                                                |
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // ------------------------------------------------------------
  // helpers for this class
  private:
    inline bool isTimeToGc() const {
      return isPessimistic() 
        ? (stats.zombie_nodes > deflt.zombieTrigger)
        : (stats.orphan_nodes > deflt.orphanTrigger);
    }

    void handleNewOrphanNode(int node);   
    void deleteOrphanNode(int node);     
    void deleteNode(int p);
    void zombifyNode(int p);

    /** Apply reduction rule to the temporary node and finalize it. 
        Once a node is reduced, its contents cannot be modified.
          @param  in    Incoming index, used only for identity reduction;
                        Or -1.
          @param  nb    Array of downward pointers.
          @return       Handle to a node that encodes the same thing.
    */
    int createReducedHelper(int in, const node_builder &nb);

    // Sanity check; used in development code.
    void validateDownPointers(const node_builder &nb) const;

  // ----------------------------------------------------------------- 
  // | 
  // |  Dealing with "node handles".
  // | 
  private:
    /// address info for nodes
    node_header *address;
    /// Size of address/next array.
    int a_size;
    /// Last used address.
    int a_last;
    /// Pointer to unused address list.
    int a_unused;
    /// Next time we shink the address list.
    int a_next_shrink;


  private:
    void expandHandleList();
    void shrinkHandleList();

  private:
    node_storage* raw_levels;

  protected:
    /// Level data. Each level maintains its own data array and hole grid.
    /// The array is shifted, so we can use level[k] with negative k.
    node_storage *levels;


  // ------------------------------------------------------------
  // |
  // |  Critical helpers.
  // |
  protected:  // TBD - make this private
    /// Returns the in-count for a node.
    inline int& getInCount(int p) {
      const node_header& node = getNode(p);
      return levels[node.level].countOf(node.offset);
    }

  public:
    /// Returns the in-count for a node.
    inline int readInCount(int p) const {
      const node_header& node = getNode(p);
      return levels[node.level].countOf(node.offset);
    }


  

  // ------------------------------------------------------------
  // |
  // |  Miscellaneous data
  // |
  protected:  
    /// uniqueness table, still used by derived classes.
    unique_table* unique;

    /// Should a terminal node be considered a stale entry in the compute table.
    /// per-forest policy, derived classes may change as appropriate.
    bool terminalNodesAreStale;

  private:
    // Keep a node_builder for each level.
    node_builder* raw_builders;
    node_builder* builders;

    // Garbage collection in progress
    bool performing_gc;

    // memory for validating incounts
    int* in_validate;
    int  in_val_size;

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
    virtual int compute(int a, int b);
    /// Low-level compute at level k on nodes a and b, return result.
    virtual int compute(int k, int a, int b);

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
