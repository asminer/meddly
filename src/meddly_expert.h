
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
#define DCASSERT(X) assert(X)
#else
#define DCASSERT(X)
#endif

#ifdef RANGE_CHECK_ON
#define CHECK_RANGE(MIN, VALUE, MAX)  {assert(VALUE<MAX);assert(VALUE>=MIN);}
#else
#define CHECK_RANGE(MIN, VALUE, MAX)
#endif

namespace MEDDLY {

  // Functions for reinterpreting an int to a float and vice-versa
  float   toFloat (int a);
  int     toInt   (float a);
  float*  toFloat (int* a);


  // Forward declarations
  class old_operation;
  class compute_cache;
  class compute_table;
  class expert_domain;
  class expert_forest;

  // ******************************************************************
  // *                     expert_variable  class                     *
  // ******************************************************************

  class expert_variable : public variable {
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

  class expert_domain : public domain {
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

      // functions inherited from class domain
      virtual int getNumForests() const;
      virtual void createVariablesBottomUp(const int* bounds, int N);
      virtual forest* createForest(bool rel, forest::range_type t,
        forest::edge_labeling ev);
      virtual void showInfo(FILE* strm);
  
      // --------------------------------------------------------------------

    private:
      int nForests;
      expert_forest** forests;
      int szForests;

      // Forests may be deleted either by calling the destructor for the
      // forest, or by destroying the domain. The domain maintains a list of
      // forest handles to delete the forests linked to it. To make sure that a
      // forest is not deleted more than once, the list of forest handles in
      // the domain is updated by a call to unlinkForest from the forest
      // destructor.
      void unlinkForest(forest *);

      friend void MEDDLY::destroyForest(forest* &f);
      friend void MEDDLY::destroyDomain(domain* &d);
  };


  // ******************************************************************
  // *                        operation  class                        *
  // ******************************************************************
  /** Generic operation.
      Necessary for compute table entries.
  */
  class operation {
      // declared and initialized in meddly.cc
      static bool& useMonolithicCT;
      // declared and initialized in meddly.cc
      static compute_table* Monolithic_CT;
      friend void MEDDLY::initialize(MEDDLY::settings);
    protected:
      int key_length; 
      int ans_length; 
      compute_table* CT;
    public:
      operation();

      // tbd: privatize this
      virtual ~operation();

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

      // are these used?

      inline int getKeyLengthInBytes() const { 
        return sizeof(int) * key_length;
      }
      inline int getAnsLengthInBytes() const {
        return sizeof(int) * ans_length;
      }
      inline int getCacheEntryLengthInBytes() const {
        return sizeof(int) * (key_length + ans_length);
      }

      /// Checks if the cache entry (in entryData[]) is stale.
      virtual bool isEntryStale(const int* entryData) = 0;

      /// Removes the cache entry (in entryData[]) by informing the
      /// applicable forests that the nodes in this entry are being removed
      /// from the cache
      virtual void discardEntry(const int* entryData) = 0;

      /// Prints a string representation of this cache entry on strm (stream).
      virtual void showEntry(FILE* strm, const int *entryData) const = 0;
  };

  // ******************************************************************
  // *                     unary_operation  class                     *
  // ******************************************************************

  /** Mechanism to apply a unary operation in a specific forest.
      Specific operations will be derived from this class.
  */
  class unary_operation : public operation {
      // not sure if this is necessary:
      const unary_opcode* opcode;
      // for cache of operations.
      unary_operation* next;
    protected:
      expert_forest* argF;
      expert_forest* resF;
      opnd_type resultType;
    public:
      unary_operation(const unary_opcode* code, 
        expert_forest* arg, expert_forest* res);

      unary_operation(const unary_opcode* code,
        expert_forest* arg, opnd_type res);

      // tbd: privatize this...
      virtual ~unary_operation();

      // handy
      const char* getName() const;

      // used by "compute manager"

      inline bool matches(const expert_forest* arg, const expert_forest* res) 
        const { 
          return (arg == argF && res == resF); 
        }

      inline bool matches(const expert_forest* arg, opnd_type res) const { 
        return (arg == argF && resultType == res); 
      }

      inline unary_operation* getNext()         { return next; }
      inline void setNext(unary_operation* n)   { next = n; }

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
  // *                      unary_builder  class                      *
  // ******************************************************************

  /** Mechanism for building unary_operations.
      This class will have many derived classes.
      Each instance of this class will know how to build a particular
      type of unary operation.
      The goal here is to simplify the logic for building operations.
  */
  class unary_builder {
      unary_builder* next;
    protected:
      const unary_opcode* opcode;
    public:
      unary_builder(unary_opcode* oc);
      virtual ~unary_builder();

      inline unary_builder* getNext() { return next; }

      /** Does this builder know how to build
          the operation for the given forests?
          Default behavior is to always return false.
      */
      virtual bool canBuild(const forest* arg, const forest* res) const;

      /** Does this builder know how to build
          the operation for the given forest and return type?
          Default behavior is to always return false.
      */
      virtual bool canBuild(const forest* arg, opnd_type res) const;

      /** Build the requested operation (it is assumed that we can.)
          Default behavior is to throw an "unknown operation" error.
      */
      virtual unary_operation* build(const forest* ar, const forest* rs) const;

      /** Build the requested operation (it is assumed that we can.)
          Default behavior is to throw an "unknown operation" error.
      */
      virtual unary_operation* build(const forest* ar, opnd_type res) const;
  };


  // ******************************************************************
  // *                       unary_opcode class                       *
  // ******************************************************************

  /// Wrapper for unary operations.
  class unary_opcode {
      const char* name;
      unary_builder* builders;
      int index;
      unary_opcode* next;
      static int next_index;
      static unary_opcode* list;

      friend void MEDDLY::initialize(settings s);
      friend void MEDDLY::cleanup();
      friend class unary_builder;

    public:
      unary_opcode(const char* n);
      ~unary_opcode();

      inline int getIndex() const { return index; }
      inline const char* getName() const { return name; }
      
      unary_operation* buildOperation(const forest* arg, const forest* res)
        const;
      unary_operation* buildOperation(const forest* arg, opnd_type res) const;
  };

  // ******************************************************************
  // *                        Unary operations                        *
  // ******************************************************************

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest
        @param  res     Result forest
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opcode* code, 
    const forest* arg, const forest* res);

  /** Find, or build if necessary, a unary operation.
        @param  code    Operation we want
        @param  arg     Argument forest
        @param  res     Result type
        @return         The matching operation, if it already exists;
                        a new operation, otherwise.
  */
  unary_operation* getOperation(const unary_opcode* code,
    const forest* arg, opnd_type result);

  // ******************************************************************
  // *                      expert_forest  class                      *
  // ******************************************************************

class expert_forest : public forest
{
  public:
    /** Constructor.
      @param  d     domain to which this forest belongs to.
      @param  rel   does this forest represent a relation.
      @param  t     the range of the functions represented in this forest.
      @param  ev    edge annotation.
      @param  r     reduction rule.
      @param  s     storage rule.
      @param  ndp   node deletion policy.
    */
    expert_forest(domain *d, bool rel, range_type t, edge_labeling ev,
      reduction_rule r, node_storage s, node_deletion_policy ndp);

  protected:
    /// Destructor.
    virtual ~expert_forest();  

  public:
    /// Returns a non-modifiable pointer to this forest's domain.
    const domain* getDomain() const;

    inline const expert_domain* getExpertDomain() const {
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
    void linkNode(int node);

    /// Decrease the link count to this node. If link count reduces to 0, this
    /// node may get marked for deletion. Call this when another node releases
    /// its connection to this node.
    void unlinkNode(int node);

    /// Increase the cache count for this node. Call this whenever this node
    /// is added to a cache.
    void cacheNode(int node);

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

#if 0
    /// Returns the cardinality of the node. Cardinality is the number of
    /// paths in the graph rooted at this node that end in a valid terminal.
    /// Paths that end at the terminal Zero (or False) are not counted.
    virtual double getCardinality(int node) const = 0;
#endif

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

  protected:

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
    friend class dd_edge;
    friend void MEDDLY::destroyDomain(domain* &d);
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
class temp_dd_edge {
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



/** Operation class.
    Abstract base class.

    TODO description.
    All concrete classes derived from operation must be Singleton classes.
    Otherwise the compute cache will treat entries from two different
    instances of the same operation (operating on the same forests) as
    entries from different operations -- which we consider to be
    inefficient behaviour.
*/


class old_operation {

  public:

    old_operation();

    virtual ~old_operation();

    /// Operation description
    virtual const char* getName() const = 0;

    /// Number of ints that make up the key (usually the operands).
    virtual int getKeyLength() const = 0;

    /// Number of ints that make up the answer (usually the results).
    virtual int getAnsLength() const = 0;

    /// Number of ints that make up the entire record (key + answer)
    virtual int getCacheEntryLength() const = 0;

    /// keyLength * sizeof(int) -- for speed store this with object
    virtual int getKeyLengthInBytes() const = 0;

    /// ansLength * sizeof(int) -- for speed store this with object
    virtual int getAnsLengthInBytes() const = 0;

    /// cacheEntryLength * sizeof(int) -- for speed store this with object
    virtual int getCacheEntryLengthInBytes() const = 0;

    /// Checks if this operation is compatible with the forests in owner.
    virtual void typeCheck(const op_info* owner) = 0;

    /// Checks if the cache entry (in entryData[]) is stale.
    virtual bool isEntryStale(const op_info* owner, const int* entryData) = 0;

    /// Removes the cache entry (in entryData[]) by informing the
    /// applicable forests that the nodes in this entry are being removed
    /// from the cache
    virtual void discardEntry(op_info* owner, const int* entryData) = 0;

    /// Prints a string representation of this cache entry on strm (stream).
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const = 0;

    /// Compute the result of this operation on the nodes provided by
    /// operands[]. operands[] contains the input as well as output dd_edges.
    /// Refer to the derived operation for information on the order of the
    /// operands in operands[].
    virtual void compute(op_info* cc, dd_edge** operands);
    
    /// Compute the result of this operation on \a a and store the result in
    /// \a b.
    virtual void compute(op_info* cc, const dd_edge& a,
      dd_edge& b);

    /// Compute the result of this operation on \a a and store the result in
    /// \a b.
    virtual void compute(op_info* cc, const dd_edge& a,
      long& b);

    /// Compute the result of this operation on \a a and store the result in
    /// \a b.
    virtual void compute(op_info* cc, const dd_edge& a,
      double& b);

    /// Compute the result of this operation on \a a and store the result in
    /// \a b.
    virtual void compute(op_info* cc, const dd_edge& a,
      ct_object &b);

    /// Compute the result of this operation on \a a and \a b and store the
    /// result in \a c.
    virtual void compute(op_info* cc, const dd_edge& a,
      const dd_edge& b, dd_edge& c);

};


/** Operation parameter class.

    Parameter information, either a forest or a value.
*/
class op_param {
  public:
    enum type {
      FOREST      = 0,
      BOOLEAN     = 1,
      INTEGER     = 2,
      REAL        = 3,
      HUGEINT     = 4,
      FLOATVECT   = 5,
      DOUBLEVECT  = 6
    };
  public:
    op_param();
    op_param(const op_param& p);

    void set(const dd_edge&);
    void set(forest*);
    void set(long);
    void set(double);
    void set(ct_object &);
    void set(const float*);
    void set(const double*);

    type getType() const;
    bool isForest() const;
    expert_forest* getForest();
    const expert_forest* readForest() const;
    void print(FILE* s) const;
    bool operator==(const op_param& a) const;
    bool operator!=(const op_param& a) const;
    bool operator<(const op_param& a) const;
    bool operator>(const op_param& a) const;
    // for convenience
    bool isForestOf(bool r, forest::range_type t, forest::edge_labeling e) const;
    bool isMxd() const;
    bool isMdd() const;
    bool isMT() const;
    bool isEVPlus() const;
    bool isEVTimes() const;
    bool isBoolForest() const;
    bool isIntForest() const;
    bool isRealForest() const;
    const domain* getDomain() const;
  private:
    expert_forest* f;
    type my_type;
}; // param class


/** Operation Info class.

    Each operation is uniquely identified by its operation pointer combined
    with the data it operates on. The data it operates on is determined
    by the forest(s) it belongs to.

    Forests, operations and compute caches are designed so that they are
    loosely-coupled. A forest does not need to know the operations that are
    based on it. But, an operation needs to know its forests as well as where
    to store its computations. A compute cache does not need to know what
    its cache entries mean but it does need someone (in this case operation)
    to decide (during garbage collection) if a cache entry should be retained
    in the cache or discarded.

    op_info stores the tuple {operation*, forest**, nForests, compute_cache*}
    to help with the above interactions.
*/
class op_info {
  public:
    op_info();
    ~op_info();
    op_info(old_operation *oper, op_param* plist, int n, compute_cache* cache);
    op_info(const op_info& a);
    op_info& operator=(const op_info& a);
    bool operator==(const op_info& a) const;
    void print(FILE* strm) const;
    bool areAllForests() const;

    old_operation* op;
    op_param* p;
    int nParams;
    compute_cache* cc;
};


/** Generic objects in compute tables.
    Used for things other than dd_edges and simple types.

    Defined in operation.cc
*/
class ct_object {
  public:
    ct_object();
    virtual ~ct_object();
    virtual op_param::type getType() = 0;
};


/** Expert Compute Manager class.

    Concrete class implements compute_manager and provides an interface
    to user-defined operations.

    TODO: At the moment, all (concrete) operations are singleton classes.
*/

class builtin_op_key;
class custom_op_key;
class expert_compute_manager : public compute_manager {
  public:
    expert_compute_manager(const settings &s);
    virtual ~expert_compute_manager();

    virtual void clearComputeTable();
    virtual void showComputeTable(FILE* strm) const;
    virtual long getNumCacheEntries() const;

    /** Removes all the stale entries in the compute table that belong to
        the specified operation handle. Note that an operation handle
        is assigned to each distinct tuple {operation code, forests[]}.

        If op is 0, ALL the stale entries in the compute table are removed.

        @param  op    Operation handle.
    */
    virtual void removeStales(op_info* op = 0);

    /** Removes all stales entries in the compute table that belong to the
        specified forest handle.

        @param  f     Forest handle.
    */
    virtual void removeStales(expert_forest* f);

    /** Removes all entries in the compute table that belong
        to the specified operation handle.

        Note that an operation handle is assigned to each distinct
        tuple {operation code, forests[]}.

        If op is 0, ALL the entries in the compute table are removed.
        In other words, it is the same as calling clearComputeTable().

        @param  op    Operation handle.
    */
    virtual void removeEntries(op_info* op = 0);

    virtual const char* getOperationName(op_code op) const;
    virtual void apply(op_code op, const dd_edge &a, dd_edge &b);
    virtual void apply(op_code op, const dd_edge &a, long &c);
    virtual void apply(op_code op, const dd_edge &a, double &c);
    virtual void apply(op_code op, const dd_edge &a, ct_object &c);

    virtual void apply(op_code op, const dd_edge &a, const dd_edge &b,
        dd_edge &c);

    virtual void vectorMatrixMultiply(double* y, const dd_edge &y_ind,
                      const double* x, const dd_edge &x_ind, const dd_edge &A);

    virtual void matrixVectorMultiply(double* y, const dd_edge &y_ind,
                      const dd_edge &A, const double* x, const dd_edge &x_ind);

    /** Same as apply(op_code, dd_edge&, dd_edge&, dd_edge&) except with
        op_info.

        When call repeatedly this reduces the overhead of calling apply()
        with op_code. Calling op_code makes the compute_manager search
        for the concrete operation handle (located within op_info). By
        calling getOpInfo(op_code,...) and then calling apply(op_info,...)
        the compute_manager needs to search for the concrete handle just
        once (when calling getOpInfo(op_code,...)).
    */
    virtual void apply(op_info* op, const dd_edge &a, const dd_edge &b,
        dd_edge &c);
    virtual void apply(op_info* op, const dd_edge &a, dd_edge &b);
    virtual void apply(op_info* op, const dd_edge &a, long &b);
    virtual void apply(op_info* op, const dd_edge &a, double &b);
    virtual void apply(op_info* op, const dd_edge &a, ct_object &b);

    /** Obtain a concrete handle to the built-in operation.
        If a built-in operation is going to be called repeatedly for the
        same set of forests, using this function with the corresponding
        apply() will improve efficiency.
        @param  op    Operation handle.
        @param  f     Forests on which the operation will be performed.
                      If you will be calling apply(UNION, a, b, c), the
                      corresponding f[] = {a's forest, b's forest, c's forest}.
        @param  N     size of f[].
        @return       A concrete-handle to the given operation.
    */
    virtual op_info* getOpInfo(op_code op, const op_param* const plist,
        int N);
    virtual op_info* getOpInfo(const old_operation* op, 
        const op_param* const p, int N);

  private:

    void addBuiltinOp(const builtin_op_key& key, const old_operation* op,
        const op_param* plist, int n); 

    // compute cache
    compute_cache* cc;

    // store for op_info entries
    // using stl map for storing op_info entries
    std::map<builtin_op_key, op_info>* builtinOpEntries;
    std::map<custom_op_key, op_info>* customOpEntries;

    const bool& useCTchaining;
    const unsigned& maxCTsize;
};


















































class builtin_op_key {
  public:
    builtin_op_key(compute_manager::op_code op, const op_param* const p, int n);
    ~builtin_op_key();
    builtin_op_key(const builtin_op_key& a);
    builtin_op_key& operator=(const builtin_op_key& a);

    bool operator<(const builtin_op_key& key) const;

    void print(FILE* strm) const;

  private:
    compute_manager::op_code opCode;
    op_param* plist;
    int nParams;
};


class custom_op_key {
  public:
    custom_op_key(const old_operation* op, const op_param* const p, int n);
    ~custom_op_key();
    custom_op_key(const custom_op_key& a);
    custom_op_key& operator=(const custom_op_key& a);

    bool operator<(const custom_op_key& key) const;

    void print(FILE* strm) const;

  private:
    old_operation* op;
    op_param* plist;
    int nParams;
};



// ****************************************************************************
// *         Implementation details -- not meant for the casual user!         *
// ****************************************************************************

// ---------------------- op_param ------------------------

inline
op_param::op_param()
{
  f = 0;
  my_type = FOREST;
}

inline
op_param::op_param(const op_param& p)
{
  f = p.f;
  my_type = p.my_type;
}

inline
void op_param::set(const dd_edge &e)
{
  f = (expert_forest*)(e.getForest());
  my_type = FOREST;
}

inline
void op_param::set(forest* _f)
{
  f = (expert_forest*)(_f);
  my_type = FOREST;
}

inline
void op_param::set(long)
{
  f = 0;
  my_type = INTEGER;
}

inline
void op_param::set(double)
{
  f = 0;
  my_type = REAL;
}

inline
void op_param::set(ct_object &x)
{
  f = 0;
  my_type = x.getType();
}

inline
void op_param::set(const float*)
{
  f = 0;
  my_type = FLOATVECT;
}

inline
void op_param::set(const double*)
{
  f = 0;
  my_type = DOUBLEVECT;
}

inline 
op_param::type op_param::getType() const 
{ 
  return my_type; 
}

inline
bool op_param::isForest() const
{
  return f;
}

inline
expert_forest* op_param::getForest() 
{ 
  return f; 
}

inline
const expert_forest* op_param::readForest() const 
{ 
  return f; 
}

inline
void op_param::print(FILE* s) const
{
  switch (my_type) {
    case INTEGER :
        fprintf(s, "INT");
        return;
    
    case REAL :
        fprintf(s, "REAL");
        return;

    case HUGEINT :
        fprintf(s, "MPZ");
        return;

    case FOREST:
        fprintf(s, "FOREST %p", f);
        return;

    case FLOATVECT:
        fprintf(s, "float[]");
        return;

    case DOUBLEVECT:
        fprintf(s, "double[]");
        return;

    default:
        fprintf(s, "Unknown param");
  }
}

inline
bool op_param::operator==(const op_param& p) const
{
  return (f == p.f) && (my_type == p.my_type);
}

inline
bool op_param::operator!=(const op_param& p) const
{
  return (f != p.f) || (my_type != p.my_type);
}

inline
bool op_param::operator<(const op_param& p) const
{
  if (my_type < p.my_type) return true;
  if (my_type > p.my_type) return false;
  return f < p.f;
}

inline
bool op_param::operator>(const op_param& p) const
{
  if (my_type > p.my_type) return true;
  if (my_type < p.my_type) return false;
  return f > p.f;
}

inline
bool op_param::isForestOf(bool r, forest::range_type t, forest::edge_labeling e) const
{
  return  (f->isForRelations() == r)  &&
          (f->getRangeType() == t)    &&
          (f->getEdgeLabeling() == e);
}

inline
bool op_param::isMxd() const
{
  return f->isForRelations();
}

inline
bool op_param::isMdd() const
{
  return !f->isForRelations();
}

inline
bool op_param::isMT() const
{
  return f->getEdgeLabeling() == forest::MULTI_TERMINAL;
}

inline
bool op_param::isEVPlus() const
{
  return f->getEdgeLabeling() == forest::EVPLUS;
}

inline
bool op_param::isEVTimes() const
{
  return f->getEdgeLabeling() == forest::EVTIMES;
}

inline
bool op_param::isBoolForest() const
{
  return f->getRangeType() == forest::BOOLEAN;
}

inline
bool op_param::isIntForest() const
{
  return f->getRangeType() == forest::INTEGER;
}

inline
bool op_param::isRealForest() const
{
  return f->getRangeType() == forest::REAL;
}

inline
const domain* op_param::getDomain() const
{
  return f ? f->getDomain() : 0;
}

// ---------------------- op_info ------------------------


inline
bool op_info::operator==(const op_info& a) const
{
  if (this == &a) return true;
  if (op == a.op && nParams == a.nParams && cc == a.cc) {
    int i = 0;
    for ( ; i < nParams && p[i] == a.p[i]; ++i)
      ;
    if (i == nParams) return true;
  }
  return false;
}


inline
void op_info::print(FILE* s) const
{
  fprintf(s, "{op = %p, nParams = %d, cc = %p, f = {", op, nParams, cc);
  if (nParams > 0) {
    p[0].print(s);
    for (int i = 1; i < nParams; ++i) {
      fprintf(s, ", ");
      p[i].print(s);
    }
  }
  fprintf(s, "}\n");
}

inline
bool op_info::areAllForests() const
{
  for (int i = 0; i < nParams; ++i) {
    if (! p[i].isForest()) return false;
  }
  return true;
}

// ---------------------- builtin_op_key ------------------------


inline
bool builtin_op_key::operator<(const builtin_op_key& key) const
{
  if (opCode < key.opCode)
    return true;
  else if (opCode > key.opCode)
    return false;
  else if (nParams < key.nParams)
    return true;
  else if (nParams > key.nParams)
    return false;
  else {
    for (int i = 0; i < nParams; ++i)
    {
      if (plist[i] < key.plist[i])
        return true;
      else if (plist[i] > key.plist[i])
        return false;
    }
  }
  return false;
}


inline
void builtin_op_key::print(FILE* s) const
{
  fprintf(s, "{opCode = %d, nParams = %d, plist = {", opCode, nParams);
  if (nParams > 0) {
    plist[0].print(s);
    for (int i = 1; i < nParams; ++i) {
      fprintf(s, ", ");
      plist[i].print(s);
    }
  }
  fprintf(s, "}\n");
}


// ---------------------- custom_op_key ------------------------


inline
bool custom_op_key::operator<(const custom_op_key& key) const
{
  if (op < key.op)
    return true;
  else if (op > key.op)
    return false;
  else if (nParams < key.nParams)
    return true;
  else if (nParams > key.nParams)
    return false;
  else {
    for (int i = 0; i < nParams; ++i)
    {
      if (plist[i] < key.plist[i])
        return true;
      else if (plist[i] > key.plist[i])
        return false;
    }
  }
  return false;
}


inline
void custom_op_key::print(FILE* s) const
{
  fprintf(s, "{op = %p, nParams = %d, plist = {", op, nParams);
  if (nParams > 0) {
    plist[0].print(s);
    for (int i = 1; i < nParams; ++i) {
      fprintf(s, ", ");
      plist[i].print(s);
    }
  }
  fprintf(s, "}\n");
}


// ************************* expert_domain ************************************


inline int expert_domain::getNumForests() const
{
  return nForests;
}


// **************************** operation *************************************


#if 0
// Operation description
inline const char* operation::getName() const
{
  return name;
}


// Number of ints that make up the key (usually the operands).
inline int operation::getKeyLength() const
{
  return keyLength;
}


// Number of ints that make up the answer (usually the results).
inline int operation::getAnsLength() const
{
  return ansLength;
}


// Number of ints that make up the entire record (key + answer)
inline int operation::getCacheEntryLength() const
{
  return cacheEntryLength;
}


// keyLength * sizeof(int)
inline int operation::getKeyLengthInBytes() const
{
  return keyLengthInBytes;
}


// ansLength * sizeof(int)
inline int operation::getAnsLengthInBytes() const
{
  return ansLengthInBytes;
}


// cacheEntryLength * sizeof(int)
inline int operation::getCacheEntryLengthInBytes() const
{
  return cacheEntryLengthInBytes;
}

#endif

// **************************** expert_forest *********************************


inline
int expert_forest::getInternalNodeSize(int p) const
{
  DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return *(getNodeAddress(p) + 2);
}


inline
int* expert_forest::getAddress(int k, int offset) const
{
  DCASSERT(level != 0 && level[mapLevel(k)].data != 0);
  return  (level[mapLevel(k)].data + offset);
}


inline
int* expert_forest::getNodeAddress(int p) const
{
  DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return getAddress(getNodeLevel(p), getNodeOffset(p));
}


inline
int expert_forest::getNodeOffset(int p) const
{
  DCASSERT(isValidNodeIndex(p));
  return  (address[p].offset);
}

// map the level value (which could be primed) to the correct index in
// the level[]
// map logical level to physical level (index of level in the level vector)
inline
int expert_forest::mapLevel(int k) const
{
  return (k >= 0)? 2 * k: ((-2 * k) - 1);
}

// map physical level to logical level (level in terms of prime, unprime)

inline
int expert_forest::unmapLevel(int k) const
{
  return (k%2 == 0)? k/2: -(k + 1)/2;
}


inline
int expert_forest::getMappedNodeHeight(int p) const
{
  return
    getNodeLevel(p) >= 0
    ? 2 * getNodeHeight(p)        // 2n
    : 2 * getNodeHeight(p) - 1;   // 2n-1
}


inline
bool expert_forest::isActiveNode(int p) const
{
  return (isValidNodeIndex(p) && (isTerminalNode(p) || getNodeOffset(p) > 0));
}


inline
bool expert_forest::isZombieNode(int p) const
{
  DCASSERT(isValidNodeIndex(p));
  DCASSERT(!isTerminalNode(p));
  return (getCacheCount(p) < 0);
}


inline
int& expert_forest::getInCount(int p) const
{
  DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return *(getNodeAddress(p));
}


inline
int& expert_forest::getCacheCount(int p) const
{
  DCASSERT(isValidNodeIndex(p));
  DCASSERT(!isTerminalNode(p));
  return address[p].cache_count;
}

inline
bool expert_forest::isTerminalNode(int p) const
{
  return (p < 1);
}

inline
bool expert_forest::getBoolean(int terminalNode) const
{
  DCASSERT(getRangeType() == forest::BOOLEAN ||
      getEdgeLabeling() != forest::MULTI_TERMINAL);
  DCASSERT(terminalNode == 0 || terminalNode == -1);
  return (terminalNode == 0)? false: true;
}

inline
int expert_forest::getTerminalNode(bool booleanValue) const
{
  DCASSERT(getRangeType() == forest::BOOLEAN ||
      getEdgeLabeling() != forest::MULTI_TERMINAL);
  return booleanValue? -1: 0;
}

inline
int expert_forest::getInteger(int terminalNode) const
{
  DCASSERT(getRangeType() == forest::INTEGER);
  DCASSERT(isTerminalNode(terminalNode) && terminalNode <= 0);
  // set 32nd bit based on 31st bit.  << gets rid of MSB; >> sign extends.
  return terminalNode << 1 >> 1;
}

inline
int expert_forest::getTerminalNode(int integerValue) const
{
  DCASSERT(getRangeType() == forest::INTEGER);
  // value has to fit within 31 bits (incl. sign)
  // int(0xc0000000) == -1073741824
  // int(0x3fffffff) == +1073741823
  // DCASSERT(-1073741824 <= integerValue && integerValue <= 1073741823);
  DCASSERT(-1073741824 <= integerValue && integerValue <= 1073741823);
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
float expert_forest::getReal(int term) const
{
  DCASSERT(getRangeType() == forest::REAL);
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
int expert_forest::getTerminalNode(float a) const
{
  DCASSERT(getRangeType() == forest::REAL);
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
int expert_forest::createTempNodeMaxSize(int lh, bool clear)
{
  return createTempNode(lh, getLevelSize(lh), clear);
}

inline
int expert_forest::getLevelSize(int lh) const {
  DCASSERT(isValidLevel(lh));
  DCASSERT(lh == 0 || level[mapLevel(lh)].data != NULL);
  if (lh < 0) {
    return static_cast<expert_domain*>(d)->getVariableBound(-lh, true);
  } else {
    return static_cast<expert_domain*>(d)->getVariableBound(lh, false);
  }
}

inline
int expert_forest::getNodeLevel(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  DCASSERT(isActiveNode(p) || isZombieNode(p));
  return (isTerminalNode(p)? 0: address[p].level);
}


inline
int expert_forest::getNodeHeight(int p) const
{
  DCASSERT(isActiveNode(p));
  return level[mapLevel(getNodeLevel(p))].height;
}


inline
bool expert_forest::isFullNode(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d, size: %d\n", __func__, p, getInternalNodeSize(p));
#endif
  DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return (getInternalNodeSize(p) > 0);
}


inline
bool expert_forest::isSparseNode(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d, size: %d\n", __func__, p, getInternalNodeSize(p));
#endif
  DCASSERT(isActiveNode(p) && !isTerminalNode(p));
  return (getInternalNodeSize(p) < 0);
}


inline
int expert_forest::getFullNodeSize(int p) const
{
#ifdef DEBUG_MDD_H
  printf("%s: p: %d\n", __func__, p);
#endif
  DCASSERT(isFullNode(p));
  return getInternalNodeSize(p);
}


inline
int expert_forest::getSparseNodeSize(int p) const
{
#ifdef DEBUG_MDD_H
printf("%s: p: %d\n", __func__, p);
#endif
  DCASSERT(isSparseNode(p));
  return -getInternalNodeSize(p);
}


inline
void expert_forest::setDownPtr(int p, int i, int value)
{
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(isActiveNode(value));
  CHECK_RANGE(0, i, getFullNodeSize(p));
  int temp = *(getNodeAddress(p) + 3 + i);
  // linkNode to new node
  linkNode(value);
  *(getNodeAddress(p) + 3 + i) = value;
  // unlinkNode old node
  unlinkNode(temp);
}


inline
void expert_forest::setDownPtrWoUnlink(int p, int i, int value)
{
  DCASSERT(!isReducedNode(p));
  DCASSERT(isFullNode(p));
  DCASSERT(isActiveNode(value));
  CHECK_RANGE(0, i, getFullNodeSize(p));
  // linkNode to new node
  linkNode(value);
  *(getNodeAddress(p) + 3 + i) = value;
}


inline
void expert_forest::setEdgeValue(int node, int index, int ev)
{
  DCASSERT(!isReducedNode(node));
  DCASSERT(isFullNode(node));
  CHECK_RANGE(0, index, getFullNodeSize(node));
  *(getNodeAddress(node) + 3 + getFullNodeSize(node) + index) = ev;
}


inline
void expert_forest::setEdgeValue(int node, int index, float ev)
{
  DCASSERT(!isReducedNode(node));
  DCASSERT(isFullNode(node));
  CHECK_RANGE(0, index, getFullNodeSize(node));
  *(getNodeAddress(node) + 3 + getFullNodeSize(node) + index) = toInt(ev);
}


inline
bool expert_forest::getDownPtrs(int node, int*& dptrs)
{
  DCASSERT(isActiveNode(node));
  if (isTerminalNode(node) || isReducedNode(node)) return false;
  DCASSERT(isFullNode(node));
  dptrs = getNodeAddress(node) + 3;
  return true;
}


inline
bool expert_forest::getEdgeValues(int node, int*& evs)
{
  DCASSERT(isActiveNode(node));
  if (isTerminalNode(node) || isReducedNode(node)) return false;
  DCASSERT(isFullNode(node));
  evs = getNodeAddress(node) + 3 + getFullNodeSize(node);
  return true;
}


inline
bool expert_forest::getEdgeValues(int node, float*& evs)
{
  DCASSERT(isActiveNode(node));
  if (isTerminalNode(node) || isReducedNode(node)) return false;
  DCASSERT(isFullNode(node));
  evs = toFloat(getNodeAddress(node) + 3 + getFullNodeSize(node));
  return true;
}


inline
bool expert_forest::getDownPtrs(int node, const int*& dptrs) const
{
  if (isTerminalNode(node)) return false;
  dptrs = isFullNode(node)
    ? getNodeAddress(node) + 3
    : getNodeAddress(node) + 3 + getSparseNodeSize(node);
  return true;
}


inline
bool expert_forest::getSparseNodeIndexes(int node, const int*& indexes) const
{
  DCASSERT(isSparseNode(node));
  if (!isSparseNode(node)) return false;
  indexes = getNodeAddress(node) + 3;
  return true;
}


inline
bool expert_forest::getEdgeValues(int node, const int*& evs) const
{
  DCASSERT(isReducedNode(node));
  if (isTerminalNode(node)) return false;
  evs = isFullNode(node)
    ? getNodeAddress(node) + 3 + getFullNodeSize(node)
    : getNodeAddress(node) + 3 + (getSparseNodeSize(node) * 2);
  return true;
}


inline
bool expert_forest::getEdgeValues(int node, const float*& evs) const
{
  DCASSERT(isReducedNode(node));
  if (isTerminalNode(node)) return false;
  evs = toFloat(isFullNode(node)
    ? getNodeAddress(node) + 3 + getFullNodeSize(node)
    : getNodeAddress(node) + 3 + (getSparseNodeSize(node) * 2));
  return true;
}


inline
int expert_forest::getFullNodeDownPtr(int p, int i) const
{
  DCASSERT(isFullNode(p));
  CHECK_RANGE(0, i, getFullNodeSize(p));
  return getNodeAddress(p)[3 + i];
}


inline
int expert_forest::getSparseNodeDownPtr(int p, int i) const
{
  DCASSERT(isSparseNode(p));
  CHECK_RANGE(0, i, getSparseNodeSize(p));
  return *(getNodeAddress(p) + 3 + getSparseNodeSize(p) + i);
}


inline
int expert_forest::getSparseNodeIndex(int p, int i) const
{
  DCASSERT(isSparseNode(p));
  CHECK_RANGE(0, i, getSparseNodeSize(p));
  return *(getNodeAddress(p) + 3 + i);
}


inline
void expert_forest::getFullNodeEdgeValue(int node, int index, int& ev) const
{
  DCASSERT(isFullNode(node));
  CHECK_RANGE(0, index, getFullNodeSize(node));
  ev = *(getNodeAddress(node) + 3 + getFullNodeSize(node) + index);
}


inline
void expert_forest::getSparseNodeEdgeValue(int node, int index, int& ev) const
{
  DCASSERT(isSparseNode(node));
  CHECK_RANGE(0, index, getSparseNodeSize(node));
  ev = *(getNodeAddress(node) + 3 + (getSparseNodeSize(node) * 2) + index);
}


inline
void expert_forest::getFullNodeEdgeValue(int node, int index, float& ev) const
{
  DCASSERT(isFullNode(node));
  CHECK_RANGE(0, index, getFullNodeSize(node));
  ev = toFloat(*(getNodeAddress(node) + 3 + getFullNodeSize(node) + index));
}


inline
void expert_forest::getSparseNodeEdgeValue(int node, int index, float& ev) const
{
  DCASSERT(isSparseNode(node));
  CHECK_RANGE(0, index, getSparseNodeSize(node));
  ev = toFloat(*(getNodeAddress(node) + 3 +
        (getSparseNodeSize(node) * 2) + index));
}


inline
void expert_forest::linkNode(int p)
{ 
  DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return;
  DCASSERT(!isPessimistic() || !isZombieNode(p));

  // increase incount
  ++getInCount(p);

  if (getInCount(p) == 1) {
    reclaimOrphanNode(p);
  }

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t+Node %d count now %d\n", p, getInCount(p));
  fflush(stdout);
#endif
}  



inline
void expert_forest::unlinkNode(int p)
{
  DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return;
  DCASSERT(!isPessimistic() || !isZombieNode(p));
  DCASSERT(getInCount(p) > 0);
  
  // decrement incoming count
  --getInCount(p);

#ifdef TRACK_DELETIONS
  fprintf(stdout, "\t-Node %d count now %d\n", p, getInCount(p));
  fflush(stdout);
#endif

  if (getInCount(p) == 0) { handleNewOrphanNode(p); }
}



inline
void expert_forest::cacheNode(int p)
{
  DCASSERT(isActiveNode(p));
  if (isTerminalNode(p)) return;
  DCASSERT(isReducedNode(p));
  getCacheCount(p)++;
#ifdef TRACK_CACHECOUNT
  fprintf(stdout, "\t+Node %d is in %d caches\n", p, getCacheCount(p));
  fflush(stdout);
#endif
}



inline
void expert_forest::uncacheNode(int p)
{
  if (isTerminalNode(p)) return;
  DCASSERT(isActiveNode(p) ||
      (!isActiveNode(p) && isPessimistic() && isZombieNode(p)));

  if (isPessimistic() && isZombieNode(p)) {
    DCASSERT(getCacheCount(p) < 0);
    getCacheCount(p)++;                           // special case; stored -ve
    if (getCacheCount(p) == 0) {
      freeZombieNode(p);
    }
    return;
  }

  DCASSERT(getCacheCount(p) > 0);
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
const domain* expert_forest::getDomain() const
{
  return d;
}


inline
domain* expert_forest::useDomain()
{
  return d;
}


inline
bool expert_forest::isForRelations() const
{
  return isRelation;
}


inline
forest::range_type expert_forest::getRangeType() const
{
  return rangeType;
}


inline
forest::edge_labeling expert_forest::getEdgeLabeling() const
{
  return edgeLabel;
}


inline
forest::reduction_rule expert_forest::getReductionRule() const
{
  return reductionRule;
}


inline
forest::node_storage expert_forest::getNodeStorage() const
{
  return nodeStorage;
}


inline
forest::node_deletion_policy expert_forest::getNodeDeletion() const
{
  return nodeDeletionPolicy;
}


inline
void expert_forest::setNodeDeletion(node_deletion_policy np)
{
  if (np == forest::NEVER_DELETE)
    throw error(error::NOT_IMPLEMENTED);
  if (getCurrentNumNodes() > 0)
    throw error(error::INVALID_OPERATION);
  nodeDeletionPolicy = np;
}

inline
void expert_forest::setNodeStorage(node_storage ns)
{
  if (nodeStorage == ns) return;
  if (getCurrentNumNodes() > 0)
    throw error(error::INVALID_OPERATION);
  nodeStorage = ns;
}

inline
void expert_forest::setReductionRule(reduction_rule r)
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
bool expert_forest::isPessimistic() const {
  return nodeDeletionPolicy == forest::PESSIMISTIC_DELETION;
}


inline
bool expert_forest::isMdd() const {
  return !isForRelations() &&
         getRangeType() == forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool expert_forest::isMtMdd() const {
  return !isForRelations() &&
         // same as == INTEGER || == REAL
         getRangeType() != forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool expert_forest::isMxd() const {
  return isForRelations() &&
         getRangeType() == forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool expert_forest::isMtMxd() const {
  return isForRelations() &&
         // same as == INTEGER || == REAL
         getRangeType() != forest::BOOLEAN &&
         getEdgeLabeling() == forest::MULTI_TERMINAL;
}


inline
bool expert_forest::isEvplusMdd() const {
  return !isForRelations() &&
         getEdgeLabeling() == forest::EVPLUS;
}


inline
bool expert_forest::isEvtimesMdd() const {
  return !isForRelations() &&
         getEdgeLabeling() == forest::EVTIMES;
}


inline
bool expert_forest::isIndexSet(int node) const {
  return getIndexSetCardinality(node) > 0;
}


inline
int expert_forest::getIndexSetCardinality(int node) const {
  DCASSERT(isEvplusMdd());
  DCASSERT(isActiveNode(node));
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
void expert_forest::setIndexSetCardinality(int node, int c) {
  DCASSERT(isEvplusMdd());
  DCASSERT(isActiveNode(node));
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
bool expert_forest::accumulate(int& A, int* vlist, int* vplist) {
  return false;
}

inline
bool expert_forest::accumulate(int& A, int* B) {
  return false;
}


inline
float toFloat(int a) {
  union { int i; float f; } n = {a};
  return n.f;
}


inline
int toInt(float a) {
  union { float f; int i; } n = {a};
  return n.i;
}


inline
float* toFloat(int* a) {
  union { int* i; float* f; } n = {a};
  return n.f;
}

inline 
const char* unary_operation::getName() const { return opcode->getName(); }

} // namespace MEDDLY

#endif
