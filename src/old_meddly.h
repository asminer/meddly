
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


/*! \file meddly.h

    MDD library interface.

    This interface is enough for "casual" users, i.e., users
    who are happy to use only the built-in operations.
    There is also an "expert" interface for users who would
    like to define new operations (in fact, the built-in
    operations use this same interface).

    General naming principles:
    Class names are of the form class_name.
    Method names are of the form verbObjectDetails.

*/

#ifndef MEDDLY_H
#define MEDDLY_H

#include <vector>
#include <memory>

namespace MEDDLY {
  /** Special value for minterms: don't care what this variable does.
      I.e., do the same thing for all possible assignments for a variable.
  */
  const int DONT_CARE  = -1;
  /** Special value for primed minterms: don't change this variable.
      Forces the primed value to equal the unprimed value for a variable.
      Undefined for unprimed variables.
  */
  const int DONT_CHANGE = -2;

  // Typedefs
  // Classes

  class memstats;
  class initializer_list;
  class input;
//  class FILE_input;
  // class istream_input;
  class output;
  // class FILE_output;
  // class ostream_output;

  class forest;
  class expert_forest;
  class unpacked_node;

  class memory_manager_style;
  class node_storage_style;

  class variable;
  class variable_order;
  class domain;
  class dd_edge;
  class enumerator;
  class ct_object;
  class unary_opname;
  class binary_opname;
  class operation;
  class unary_operation;
  class binary_operation;

  /// Argument and result types for apply operations.
  enum class opnd_type {
    FOREST      = 0,
    BOOLEAN     = 1,
    INTEGER     = 2,
    REAL        = 3,
    HUGEINT     = 4,
    FLOATVECT   = 5,
    DOUBLEVECT  = 6
  };

  // ******************************************************************
  // *                    miscellaneous  functions                    *
  // ******************************************************************

#ifdef __GMP_H__
  ct_object& get_mpz_wrapper();
  void unwrap(const ct_object &, mpz_t &value);
#endif

  // ******************************************************************
  // *                    Memory management schemes                   *
  // ******************************************************************

  extern const memory_manager_style* ORIGINAL_GRID;
  extern const memory_manager_style* ARRAY_PLUS_GRID;
  extern const memory_manager_style* MALLOC_MANAGER;
  extern const memory_manager_style* HEAP_MANAGER;
  extern const memory_manager_style* FREELISTS;   // used for compute tables

  // ******************************************************************
  // *                     Node storage mechanisms                    *
  // ******************************************************************

  /**
    New, "simple" style with memory manager removed.
    This node storage mechanism relies on the
    memory_manager_style for memory management.
  */
  extern const node_storage_style* SIMPLE_STORAGE;


   extern const node_storage_style* PATTERN_STORAGE;
   extern const node_storage_style* BEST_STORAGE;


  //
  // From here are "old" mechanisms for node storage,
  // with built-in memory managers.
  //

  /** "Classic" node storage mechanism.
      The node storage mechanism from early versions of this library.
  */
  extern const node_storage_style* CLASSIC_STORAGE;

  /** Similar to classic.
      Differences are in class design, so we can measure overhead
      (if any) of class design differences.
  */
  extern const node_storage_style* SIMPLE_GRID;

  /** Like classic, but use an array of lists for hole management.
  */
  extern const node_storage_style* SIMPLE_ARRAY;

  /** Like classic, but use heaps for hole management.
  */
  extern const node_storage_style* SIMPLE_HEAP;

  /** Classic node storage but no hole management.
  */
  extern const node_storage_style* SIMPLE_NONE;

  /** Nodes are stored in a compressed form.
      Holes are managed using the original grid structure.
  */
  extern const node_storage_style* COMPACT_GRID;


  // ******************************************************************
  // *                     Named unary operations                     *
  // ******************************************************************

  /** Create a copy of a dd_edge.
      The copy may be stored in any forest as long as it belongs to the
      same domain as the original and the transformation is valid.
      Copying is valid with the following:
      MDD to MTMDD, MTMDD to MDD, MXD to MTMXD, MTMXD to MXD.
  */
  extern const unary_opname* COPY;

  /// Unary operation.  Return the number of variable assignments
  /// so that the function evaluates to non-zero.
  extern const unary_opname* CARDINALITY;

  /// For BOOLEAN forests, flip the return values.
  extern const unary_opname* COMPLEMENT;

  /// Find the largest value returned by the function.
  extern const unary_opname* MAX_RANGE;

  /// Find the smallest value returned by the function.
  extern const unary_opname* MIN_RANGE;

  /// Convert MDD to EV+MDD index set.  A special case of COPY, really.
  extern const unary_opname* CONVERT_TO_INDEX_SET;

  /// Extract cycles (EV+MDD) from transitive closure (EV+MxD)
  extern const unary_opname* CYCLE;

  /// Randomly select one state from a set of states
  extern const unary_opname* SELECT;

  // ******************************************************************
  // *                    Named  binary operations                    *
  // ******************************************************************

  /// Set operation for forests with range_type of BOOLEAN. All operands
  /// must belong to the same forest.
  extern const binary_opname* UNION;
  /// Set operation for forests with range_type of BOOLEAN. All operands
  /// must belong to the same forest.
  extern const binary_opname* INTERSECTION;
  /// Set operation for forests with range_type of BOOLEAN. All operands
  /// must belong to the same forest.
  extern const binary_opname* DIFFERENCE;

  /// Binary operation.  Combines two functions into a single one,
  /// where the operands are MDDs and the result is an MXD.
  /// Specifically, for MDD operands f and g, produces MXD h where
  /// h(xn, x'n, ..., x1, x'1) = f(xn, ..., x1) * g(x'n, ..., x'1)
  /// Works for BOOLEAN forests.
  extern const binary_opname* CROSS;

  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* MINIMUM;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* MAXIMUM;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* PLUS;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* MINUS;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* MULTIPLY;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* DIVIDE;
  /// For forests with range_type of INTEGER. All operands must
  /// belong to the same forest.
  extern const binary_opname* MODULO;

  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* EQUAL;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* NOT_EQUAL;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* LESS_THAN;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* LESS_THAN_EQUAL;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* GREATER_THAN;
  /// For forests with range_type of INTEGER and REAL. All operands must
  /// belong to the same forest.
  extern const binary_opname* GREATER_THAN_EQUAL;

  /** Plus operation used to compute transitive closure and further
      minimum witness. The first operand must be an EV+MxD and the second
      operand must be an EV+MDD. The result is an EV+MxD.
   */
  extern const binary_opname* PRE_PLUS;
  extern const binary_opname* POST_PLUS;

  /** Image operations on a set-of-states with a transition function.
      The first operand must be the set-of-states and the second operand
      must be the transition function. The result is a set-of-states that
      must be stored in the same forest as the first operand.

      Applies to:
      PRE_IMAGE, POST_IMAGE,
      TC_POST_IMAGE,
      REACHABLE_STATES_DFS, REACHABLE_STATES_BFS,
      REVERSE_REACHABLE_DFS, REVERSE_REACHABLE_BFS.
  */
  extern const binary_opname* PRE_IMAGE;
  extern const binary_opname* POST_IMAGE;
  extern const binary_opname* TC_POST_IMAGE;
  extern const binary_opname* REACHABLE_STATES_DFS;
  extern const binary_opname* REACHABLE_STATES_BFS;
  extern const binary_opname* REVERSE_REACHABLE_DFS;
  extern const binary_opname* REVERSE_REACHABLE_BFS;

  /** Vector matrix multiply, where the first argument is vector (MDD),
      the second argument is a matrix (MXD), and the result is a vector (MDD).
  */
  extern const binary_opname* VM_MULTIPLY;

  /** Matrix vector multiply, where the first argument is a matrix (MXD),
      the second argument is a vector (MDD), and the result is a vector (MDD).
  */
  extern const binary_opname* MV_MULTIPLY;

  /** Matrix multiplication, where the first argument is a matrix (MXD),
      the second argument is a matrix (MXD), and the result is a matrix (MXD),
      such that, C[m][n] += A[m][i] * B[i][n], for all m, n and i.
  */
  extern const binary_opname* MM_MULTIPLY;

  // ******************************************************************
  // *                  library management functions                  *
  // ******************************************************************

  /** Initialize the library with default settings.
      See meddly_expert.h for functions to initialize
      the library with non-default settings.
      Should be called before using any other functions.
  */
  void initialize();

  /** Clean up the library.
      Can be called to free memory used by the library;
      after it is called, the library may be initialized again.
  */
  void cleanup();

  /** Get the information about the library.
      @param  what  Determines the type of information to obtain.
      @return A human-readable information string.
              The string depends on parameter \a what, as follows.
              0: Title string, e.g., "MDD library version 0.0.0"
              1: Copyright string
              2: License string
              3: Reference url
              4: String with library features
              5: Release date
              Anything else: null string.
  */
  const char* getLibraryInfo(int what = 0);

  // ******************************************************************
  // *                   object creation  functions                   *
  // ******************************************************************

  /** Front-end function to create a variable.
      This is required because variable is an abstract class.
        @param  bound   The initial bound for the variable.
                        If bound<=0, the variable is marked as extensible,
                        with initial bound as abs(bound).
                        Note: an extensible variable has a range [1 .. +infinity].
        @param  name    Variable name (used only in display / debugging), or 0.
        @return A new variable, or 0 on error.
  */
  variable* createVariable(int bound, char* name);

  /** Front-end function to create a domain with the given variables.
        @param  vars    List of variables, in order.
                        vars[i] gives the variable at level i.
                        Note that vars[0] should be 0.
        @param  N       Number of variables.
                        vars[N] refers to the top-most variable.

        @return A new domain.
  */
  domain* createDomain(variable** vars, int N);

  /** Front-end function to create an empty domain.
      This is required because domain is an abstract class.
  */
  domain* createDomain();

  /** Front-end function to create a domain with given variable bounds.
      Equivalent to creating an empty domain and then building the
      domain bottom up.

        @param  bounds  variable bounds.
                        bounds[i] gives the bound for the variable
                        at level i+1.
                        If bound<=0, the variable is marked as extensible,
                        with initial bound as abs(bound).
                        Note: an extensible variable has a range [1 .. +infinity].
        @param  N       Number of variables.

        @return A new domain.
  */
  domain* createDomainBottomUp(const int* bounds, int N);

/* Commented out as of version 0.10
#ifdef _MSC_VER
  __declspec(deprecated)
#endif
#ifdef __GNUC__
  __attribute__ ((deprecated))
#endif
  /// This function is deprecated as of version 0.4;
  /// use "createDomain" instead.
  domain* MEDDLY_createDomain();
*/


  // ******************************************************************
  // *                  object destruction functions                  *
  // ******************************************************************

  /** Front-end function to destroy a domain.
      For consistency.
  */
  void destroyDomain(domain* &d);

  /** Front-end function to destroy a forest.
  */
  void destroyForest(forest* &f);

  // ******************************************************************
  // *                          Unary  apply                          *
  // ******************************************************************

  /** Apply a unary operator.
      The operand and the result are not necessarily in the same forest,
      but they must belong to forests that share the same domain.
      This is useful, for instance, for copying a function to a new forest.
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Output parameter: the result, where \a c = \a op \a a.
  */
  void apply(const unary_opname* op, const dd_edge &a, dd_edge &c);

  /** Apply a unary operator.
      For operators whose result is an integer.
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Output parameter: the result, where \a c = \a op \a a.
  */
  void apply(const unary_opname* op, const dd_edge &a, long &c);

  /** Apply a unary operator.
      For operators whose result is a real.
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Output parameter: the result, where \a c = \a op \a a.
  */
  void apply(const unary_opname* op, const dd_edge &a, double &c);

  void apply(const unary_opname* op, const dd_edge &a, opnd_type cr,
    ct_object &c);

#ifdef __GMP_H__
  /** Apply a unary operator.
      For operators whose result is an arbitrary-precision integer
      (as supplied by the GNU MP library).
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Input: an initialized MP integer.
                      Output parameter: the result, where \a c = \a op \a a.
  */
  void apply(const unary_opname* op, const dd_edge &a, mpz_t &c);
#endif

  // ******************************************************************
  // *                          Binary apply                          *
  // ******************************************************************

  /** Apply a binary operator.
      \a a, \a b and \a c are not required to be in the same forest,
      but they must have the same domain. The result will be in the
      same forest as \a result. The operator decides the type of forest
      for each \a dd_edge.
      Useful, for example, for constructing comparisons
      where the resulting type is "boolean" but the operators are not,
      e.g., c = f EQUALS g.
      @param  op    Operator handle.
      @param  a     First operand.
      @param  b     Second operand.
      @param  c     Output parameter: the result,
                    where \a c = \a a \a op \a b.
  */
  void apply(const binary_opname* op, const dd_edge &a, const dd_edge &b,
    dd_edge &c);


};  // namespace MEDDLY


// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************



// ******************************************************************
// ******************************************************************
// ******************************************************************
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        enumerator class                        *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Class for enumerating values encoded by a dd-edge.
    Effectively, a single class encapsulating
    several possible iterators.
*/
class MEDDLY::enumerator {
  public:
    class iterator {
      public:
        iterator(const expert_forest* F);
        virtual ~iterator();

        bool build_error() const;

        virtual bool start(const dd_edge &e);
        virtual bool start(const dd_edge &e, const int* m);
        virtual bool next() = 0;

        /**
            Return the highest level changed during the last increment.
        */
        int levelChanged() const;

        /** Get the current variable assignments.
            For variable i, use index i for the
            unprimed variable, and index -i for the primed variable.
        */
        const int* getAssignments() const;

        /** Get primed assignments.
            It is much faster to use getAssigments()
            and look at the negative indexes;
            however, this works.
        */
        const int* getPrimedAssignments();

        /// For integer-ranged edges, get the current non-zero value.
        virtual void getValue(int& edgeValue) const;

        /// For integer-ranged edges, get the current non-zero value.
        virtual void getValue(long& edgeValue) const;

        /// For real-ranged edges, get the current non-zero value.
        virtual void getValue(float& edgeValue) const;


      protected:
        // Current parent forest.
        const expert_forest* F;
        // Path, as list of unpacked nodes
        unpacked_node*    rawpath;
        unpacked_node*    path;   // rawpath, shifted so we can use path[-k]
        // Path nnz pointers
        int*      rawnzp;
        int*      nzp;   // rawnzp, shifted so we can use nzp[-k]
        // Path indexes
        int*      rawindex;
        int*      index;  // rawindex, shifted so we can use index[-k]
        //
        int       minLevel; // 1 or -#vars, depending.
        int       maxLevel; // #vars
        //
        int       level_change;

      private:
        // Used only by getPrimedAssignments.
        int*      prindex;
    };

  public:
    /// Enumerator type.
    enum type {
      EMPTY,
      FULL,
      ROW_FIXED,
      COL_FIXED
    };

  public:
    /// Empty constructor.
    enumerator();

    /// Proper constructor.
    enumerator(type t, const forest* F);

    /** What is usually wanted constructor.
        Equivalent to enumerator(FULL, e.getForest())
        followed by start(e).
    */
    enumerator(const dd_edge &e);

    /// Destructor.
    ~enumerator();

    /// Re-initialize
    void init(type t, const forest* F);

    operator bool() const;

    void operator++();

    /**
        Start iterating through edge e.
    */
    void start(const dd_edge &e);

    /** Start iterating through edge e.
        The unprimed variables will be fixed to
        the given values.
          @param  e         Edge to iterate.
                            Must be a relation.
          @param  minterm   Array of dimension 1+vars in e.
                            minterm[k] gives the fixed variable
                            assignment for (unprimed) variable k.
    */
    void startFixedRow(const dd_edge &e, const int* minterm);

    /** Start iterating through edge e.
        The primed variables will be fixed to
        the given values.
          @param  e         Edge to iterate.
                            Must be a relation.
          @param  minterm   Array of dimension 1+vars in e.
                            minterm[k] gives the fixed variable
                            assignment for (unprimed) variable k.
    */
    void startFixedColumn(const dd_edge &e, const int* minterm);


    /** Get the current variable assignments.
        For variable i, use index i for the
        unprimed variable, and index -i for the primed variable.
    */
    const int* getAssignments() const;

    /// Get the current primed variable assignments.
    const int* getPrimedAssignments() const;

    void getValue(int &v) const;
    void getValue(long &v) const;
    void getValue(float &v) const;

    int levelChanged() const;

    type getType() const;

  private:
    iterator* I;
    bool is_valid;
    type T;
};

#endif
