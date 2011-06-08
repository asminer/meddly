
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
    Global functions have names prefixed with "MEDDLY_".

*/

#ifndef MEDDLY_H
#define MEDDLY_H

#include <cstdio>
#include <cassert>

// Forward declarations
class domain;
class dd_edge;
class op_info;
class ct_object;

#ifdef __GMP_H__
ct_object& MEDDLY_get_mpz_wrapper();
void MEDDLY_unwrap(const ct_object &, mpz_t &value);
#endif


// ******************************************************************
// *                                                                *
// *                      library  information                      *
// *                                                                *
// ******************************************************************

/** Get the information about the library.
    @param  what  Determines the type of information to obtain.
    @return A human-readable information string.
            The string depends on parameter \a what, as follows.
            0: Title string, e.g., "MDD library version 0.0.0"
            1: Copyright string
            2: License string
            3: Reference url
            4: String with library features
            Anything else: null string.
*/
const char* MEDDLY_getLibraryInfo(int what = 0);


// ******************************************************************
// *                                                                *
// *                                                                *
// *                          forest class                          *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Forest class.
    Abstract base class.
    
    A data structure for managing collections of functions
    (or sets, or vectors, or relations, or matrices,
    depending on your conceptual view) represented in a single
    decision diagram forest over a common domain.
    
    Some interesting features of forests:
    - The reduction rule can be different for different variables.
    - The storage rule can be different for different variables.

    TBD: discussion of garbage collection.

    When a forest is destroyed, all of the corresponding dd_edges
    are also destroyed, as are any compute table entries for the forest.
*/
class forest {
  public:
    /** Types of values that we can currently store in forests.
        I.e., if every node in a forest is a function,
        these are the possible ranges for a function.
    */
    enum range_type {
      /// boolean-valued functions.
      BOOLEAN,
      /// integer-valued functions.
      INTEGER,
      /// real-valued functions.
      REAL
    };

    /** Edge annotation mechanism.
    */
    enum edge_labeling {
      /// Edges unlabeled, all values stored in distinct terminals.
      MULTI_TERMINAL,
      /// Edges labeled, values summed along path.
      EVPLUS,
      /// Edges labeled, values multiplied along path.
      EVTIMES
      // TBD: there may be others in the future :^)
    };

    /** Supported node reduction rules.
        Currently, the following reduction rules are allowed:
        - Fully reduced, meaning that duplicate and redundant nodes are
          eliminated.
        - Quasi reduced, meaning that duplicate nodes are eliminated.
        - Identity reduced, for relations only, meaning that duplicate
          nodes are eliminated, as are "identity" pairs of primed, unprimed
          variables.
    */
    enum reduction_rule {
      /// Nodes are fully reduced.
      FULLY_REDUCED,
      /// Nodes are quasi-reduced.
      QUASI_REDUCED,
      /// Nodes are identity-reduced.
      IDENTITY_REDUCED
    };

    /** Supported node storage meachanisms.
    */
    enum node_storage {
      /// Truncated full storage.
      FULL_STORAGE,
      /// Sparse storage.
      SPARSE_STORAGE,
      /// For each node, either full or sparse, whichever is more compact.
      FULL_OR_SPARSE_STORAGE
    };

    /** Supported node deletion policies.
    */
    enum node_deletion_policy {
      /** Never delete nodes.
          Useful for debugging the garbage collector.
      */
      NEVER_DELETE,
      /** Nodes are marked for deletion.
          We assume that a deleted node might be used again,
          so we keep it until it no longer appears in any compute table.
      */
      OPTIMISTIC_DELETION,
      /// Nodes are removed as soon as they become disconnected.
      PESSIMISTIC_DELETION
    };

    /** Error codes for forest methods.
        "SUCCESS" is fixed to the code 0, so we can
        test for errors using:  if (err) { HandleError(); }
    */
    enum error {
      /// Returned whenever an operation is successful.
      SUCCESS=0,
      /// Returned whenever the requested operation is not yet implemented!
      NOT_IMPLEMENTED,
      /// Returned whenever an operation fails due to lack of memory.
      INSUFFICIENT_MEMORY,
      /// Returned when an operation is not supported for the given forest.
      INVALID_OPERATION,
      /// Returned if a provided variable handle is erroneous.
      INVALID_VARIABLE,
      /// Returned if a provided variable value is out of range.
      INVALID_ASSIGNMENT
    };

    /// Constructor -- this class cannot be instantiated.
    forest();

    /// Destructor.  Will request the domain to destroy the forest.
    virtual ~forest();  

    /** Produce a human-readable name for an error code.
        @param  e   Error code.
        @return     A string that describes the error, suitable for printing.
    */
    static const char* getErrorCodeName(error e);

    /// Returns a non-modifiable pointer to this forest's domain.
    virtual const domain* getDomain() const = 0;

    /// Returns a pointer to this forest's domain.
    virtual domain* useDomain() = 0;

    /// Does this forest represent relations or matrices?
    virtual bool isForRelations() const = 0;

    /// Returns the range type.
    virtual range_type getRangeType() const = 0;

    /// Returns the edge labeling mechanism.
    virtual edge_labeling getEdgeLabeling() const = 0;

    /** Set the specified reduction rule for all variables in this forest.
        @param  r   Desired reduction rule.
        @return     An appropriate error code.
    */
    virtual error setReductionRule(reduction_rule r) = 0;

    /** Set the node storage mechanism for all variables in this forest.
        The same as calling setNodeStorageForVariable(), for all variables.
        @param  ns  Desired node storage mechanism.
        @return     An appropriate error code.
    */
    virtual error setNodeStorage(node_storage ns) = 0;

    /** Set the node deletion policy for all variables in this forest.
        Determines how aggressively node memory should be reclaimed when
        nodes become disconnected.
        @param  np  Policy to use.
        @return     An appropriate error code.
        TODO: need to specify what happens when the policy is changed.
    */
    virtual error setNodeDeletion(node_deletion_policy np) = 0;

    /** Set a compaction threshold for all variables in this forest.
        The threshold is set as a percentage * 100.
        Whenever the \b unused memory for this variable exceeds the threshold
        (as a percentage), a compaction operation is performed automatically.
        @param  p   Percentage * 100, i.e., use 50 for 50\%.
        @return     An appropriate error code.
    */
    virtual error setCompactionThreshold(unsigned p) = 0;

    /** Force garbage collection.
        All disconnected nodes in this forest are discarded along with any
        compute table entries that may include them.
        @return     An apprpropriate error code.
    */
    virtual error garbageCollect() = 0;

    /** Compact the memory for all variables in this forest.
        This is not the same as garbage collection.
        @return     An apprpropriate error code.
    */
    virtual error compactMemory() = 0;

    /** Get the current number of nodes in the forest, at all levels.
        @return     The current number of nodes, not counting deleted or
                    marked for deletion nodes.
    */
    virtual long getCurrentNumNodes() const = 0;

    /** Get the current total memory used by the forest.
        This should be equal to summing getMemoryUsedForVariable()
        over all variables.
        @return     Current memory used by the forest.
    */
    virtual long getCurrentMemoryUsed() const = 0;

    /** Get the current total memory allocated by the forest.
        This should be equal to summing getMemoryAllocatedForVariable()
        over all variables.
        @return     Current total memory allocated by the forest.
    */
    virtual long getCurrentMemoryAllocated() const = 0;

    /** Get the peak number of nodes in the forest, at all levels.
        This will be at least as large as calling getNumNodes() after
        every operation and maintaining the maximum.
        @return     The peak number of nodes that existed at one time,
                    in the forest.
    */
    virtual long getPeakNumNodes() const = 0;

    /** Get the peak memory used by the forest.
        @return     Peak total memory used by the forest.
    */
    virtual long getPeakMemoryUsed() const = 0;

    /** Get the peak memory allocated by the forest.
        @return     Peak memory allocated by the forest.
    */
    virtual long getPeakMemoryAllocated() const = 0;


    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = i for 0 <= i < size(vh).

        For example, in a forest with range_type INTEGER, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
         X   0   X  0
         X   1   X  1
         X   2   X  2
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  primedLevel
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = i for 0 <= i < size(vh).
        @return       INVALID_VARIABLE: if vh is invalid.
                      INVALID_OPERATION: if dd_edge does not belong to this
                      forest; or if size(vh) != 2 for a BOOLEAN forest.
                      INVALID_ASSIGNMENT: if this forest is not a relation,
                      and primedLevel is set to true.
                      SUCCESS: dd_edge created successfully.
    */
    virtual error createEdgeForVar(int vh, bool primedLevel, dd_edge& a) = 0;


    /** Create an edge representing the subset of a Matrix Diagram.

        size(vlist) == number of variables in the forest + 1 (for terminals)
        size(vlist[i]) == size (i.e. bound) of variable i
        size(vlist) == size(vplist)
        size(vlist[i]) == size(vplist[i])
        If vlist[i][j] is true, that index is included in the mask
        If vlist[i][j] is false, that index is excluded from the mask.
        
        TODO: write a better description (an example might also help).
    */
    virtual error createSubMatrix(const bool* const* vlist,
        const bool* const* vplist, const dd_edge a, dd_edge& b) = 0;


    /** Returns element \a e at index \a i from an Index Set EV+MDD.
    
        size(e) = number of variables in the forest + 1 (for terminals).
        TODO: complete this description
    */
    virtual error getElement(const dd_edge& a, int index, int* e) = 0;


    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = terms[i] for 0 <= i < size(vh).

        For example, in a forest with range_type BOOLEAN, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
         X   0   X  terms[0]
         X   1   X  terms[1]
         X   2   X  terms[2]
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  primedLevel
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  terms Array of boolean terminal values.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = terms[i]
                      for 0 <= i < size(vh).
        @return       INVALID_VARIABLE: if vh is invalid.
                      INVALID_OPERATION: if dd_edge does not belong to this
                      forest; or if this is not a BOOLEAN forest.
                      INVALID_ASSIGNMENT: if this forest is not a relation,
                      and primedLevel is set to true.
                      SUCCESS: dd_edge created successfully.
    */
    virtual error createEdgeForVar(int vh, bool primedLevel,
       bool* terms, dd_edge& a) = 0;

    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = terms[i] for 0 <= i < size(vh).

        For example, in a forest with range_type INTEGER, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
         X   0   X  terms[0]
         X   1   X  terms[1]
         X   2   X  terms[2]
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  primedLevel
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  terms Array of boolean terminal values.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = terms[i]
                      for 0 <= i < size(vh).
        @return       INVALID_VARIABLE: if vh is invalid.
                      INVALID_OPERATION: if dd_edge does not belong to this
                      forest; or if this is not a INTEGER forest.
                      INVALID_ASSIGNMENT: if this forest is not a relation,
                      and primedLevel is set to true.
                      SUCCESS: dd_edge created successfully.
    */
    virtual error createEdgeForVar(int vh, bool primedLevel,
       int* terms, dd_edge& a) = 0;

    /** Create an edge such that
        f(v_1, ..., vh=i, ..., v_n) = terms[i] for 0 <= i < size(vh).

        For example, in a forest with range_type REAL, with 3 variables,
        all of size 3, if vh == 2. An edge is created such that
        v_1 v_2 v_3 TERM
         X   0   X  terms[0]
         X   1   X  terms[1]
         X   2   X  terms[2]
        where X represents "all possible".

        \a primedLevel is useful only with forests that store relations. Here,
        \a primedLevel is used to indicate whether the edge is to be created
        for the primed or the unprimed level.

        @param  vh    Variable handle.
        @param  primedLevel
                      true: creates node for the primed vh variable.
                      false: creates node for the unprimed vh variable.
        @param  terms Array of boolean terminal values.
        @param  a     return a handle to a node in the forest such that
                      f(v_1, ..., vh=i, ..., v_n) = terms[i]
                      for 0 <= i < size(vh).
        @return       INVALID_VARIABLE: if vh is invalid.
                      INVALID_OPERATION: if dd_edge does not belong to this
                      forest; or if this is not a REAL forest.
                      INVALID_ASSIGNMENT: if this forest is not a relation,
                      and primedLevel is set to true.
                      SUCCESS: dd_edge created successfully.
    */
    virtual error createEdgeForVar(int vh, bool primedLevel,
       float* terms, dd_edge& a) = 0;

    /** Create an edge as the union of several explicit vectors.
        The range type of the forest must be booleans.
        The forest must not be for relations.
        @param  vlist Array of vectors. Each vector has dimension equal
                      to one plus the largest variable handle in the domain.
                      A vector \a x indicates a set of variable assignments,
                      where x[vh] less than 0 means
                      "don't care for variable vh", otherwise x[vh] gives
                      the variable assignment for vh.
        @param  N     Number of vectors (dimension of \a vlist).
        @param  e     returns a handle to a node in the forest, such that
                      f(v_1, ..., v_n) = 1, iff there is a vector
                      x in \a vlist corresponding to the variable assignments
                      v_1, ..., v_n; f(v_1, ..., v_n) = 0, otherwise.
        @return       An appropriate error code.
    */
    virtual error createEdge(const int* const* vlist, int N, dd_edge &e) = 0;

    /** Create an edge as the union of several vectors and return values.
        The range type of the forest must be integers.
        The forest must not be for relations.
        @param  vlist Array of vectors. Each vector has dimension equal
                      to one plus the largest variable handle in the domain.
                      A vector \a x indicates a set of variable assignments,
                      where x[vh] less than 0 means
                      "don't care for variable vh", otherwise x[vh] gives
                      the variable assignment for vh.
        @param  terms Array of return values, same dimension as \a vlist.
        @param  N     Number of vectors (dimension of \a vlist).
        @param  e     returns a handle to a node in the forest that encodes
                      the function f: f(v_1, ..., v_n) = terms[j]
                      iff j is the smallest integer such that vector vlist[j]
                      corresponds to the variable assignments v_1, ..., v_n;
                      f(v_1, ..., v_n) = 0, otherwise.
        @return       An appropriate error code.
    */
    virtual error createEdge(const int* const* vlist, const int* terms, int N,
        dd_edge &e) = 0;

    /** Create an edge as the union of several vectors and return values.
        The range type of the forest must be reals.
        The forest must not be for relations.
        @param  vlist Array of vectors. Each vector has dimension equal
                      to one plus the largest variable handle in the domain.
                      A vector \a x indicates a set of variable assignments,
                      where x[vh] less than 0 means
                      "don't care for variable vh", otherwise x[vh] gives
                      the variable assignment for vh.
        @param  terms Array of return values, same dimension as \a vlist.
        @param  N     Number of vectors (dimension of \a vlist).
        @param  e     returns a handle to a node in the forest that encodes
                      the function f: f(v_1, ..., v_n) = terms[j]
                      iff j is the smallest integer such that vector vlist[j]
                      corresponds to the variable assignments v_1, ..., v_n;
                      f(v_1, ..., v_n) = 0, otherwise.
        @return       An appropriate error code.
    */
    virtual error createEdge(const int* const* vlist, const float* terms,
        int N, dd_edge &e) = 0;


    /** Create an edge as the union of several explicit matrices.
        The range type of the forest must be booleans.
        The forest must be for relations.
        @param  vlist   Array of vectors. Each vector has dimension equal to
                        one plus the largest variable handle in the domain.
                        A vector \a x indicates a set of unprimed variable
                        assignments, where x[vh] equal to -1 means
                        "don't care for unprimed variable vh", x[vh] equal
                        to -2 means "don't change for variable vh",
                        otherwise x[vh] gives the variable assignment for
                        unprimed variable vh.
        @param  vplist  Array of vectors, same dimension as \a vlist.
                        A vector \a x = vplist[i] indicates a set of primed
                        variable assignments, where x[vh] less than 0 means
                        "don't change for variable vh" iff we are not
                        changing the unprimed variable, and
                        "don't care for primed variable vh" otherwise.
                        Otherwise x[vh] gives the variable assignment for
                        primed variable vh.
        @param  N       Number of vectors (dimension of \a vlist).
        @param  e       returns a handle to a node in the forest, such that
                        f(v_1, v'_1, ..., v_n, v'_n) = 1, iff there is a
                        vector vlist[i] corresponding to the variable
                        assignments v_1, ..., v_n, and vplist[i] corresponds
                        to variable assignments v'_1, ..., v'_n.
                        f(v_1, v'_1, ..., v_n, v'_n) = 0, otherwise.
        @return         An appropriate error code.
    */
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge &e) = 0;


    /** Create an edge as the union of several explicit matrices.
        The range type of the forest must be integers.
        The forest must be for relations.
        @param  vlist   Array of vectors. Each vector has dimension equal to
                        one plus the largest variable handle in the domain.
                        A vector \a x indicates a set of unprimed variable
                        assignments, where x[vh] equal to -1 means
                        "don't care for unprimed variable vh", x[vh] equal
                        to -2 means "don't change for variable vh",
                        otherwise x[vh] gives the variable assignment for
                        unprimed variable vh.
        @param  vplist  Array of vectors, same dimension as \a vlist.
                        A vector \a x = vplist[i] indicates a set of primed
                        variable assignments, where x[vh] less than 0 means
                        "don't change for variable vh" iff we are not
                        changing the unprimed variable, and
                        "don't care for primed variable vh" otherwise.
                        Otherwise x[vh] gives the variable assignment for
                        primed variable vh.
        @param  terms   Array of return values, same dimension as \a vlist.
        @param  N       Number of vectors (dimension of \a vlist).
        @param  e       returns a handle to a node in the forest, such that
                        f(v_1, v'_1, ..., v_n, v'_n) = term[i], iff there is
                        a vector vlist[i] corresponding to the variable
                        assignments v_1, ..., v_n, and vplist[i] corresponds
                        to variable assignments v'_1, ..., v'_n.
                        f(v_1, v'_1, ..., v_n, v'_n) = 0, otherwise.
        @return         An appropriate error code.
    */
    virtual error createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e) = 0;


    /** Create an edge as the union of several explicit matrices.
        The range type of the forest must be reals.
        The forest must be for relations.
        @param  vlist   Array of vectors. Each vector has dimension equal to
                        one plus the largest variable handle in the domain.
                        A vector \a x indicates a set of unprimed variable
                        assignments, where x[vh] equal to -1 means
                        "don't care for unprimed variable vh", x[vh] equal
                        to -2 means "don't change for variable vh",
                        otherwise x[vh] gives the variable assignment for
                        unprimed variable vh.
        @param  vplist  Array of vectors, same dimension as \a vlist.
                        A vector \a x = vplist[i] indicates a set of primed
                        variable assignments, where x[vh] less than 0 means
                        "don't change for variable vh" iff we are not
                        changing the unprimed variable, and
                        "don't care for primed variable vh" otherwise.
                        Otherwise x[vh] gives the variable assignment for
                        primed variable vh.
        @param  terms   Array of return values, same dimension as \a vlist.
        @param  N       Number of vectors (dimension of \a vlist).
        @param  e       returns a handle to a node in the forest, such that
                        f(v_1, v'_1, ..., v_n, v'_n) = term[i], iff there is
                        a vector vlist[i] corresponding to the variable
                        assignments v_1, ..., v_n, and vplist[i] corresponds
                        to variable assignments v'_1, ..., v'_n.
                        f(v_1, v'_1, ..., v_n, v'_n) = 0, otherwise.
        @return         An appropriate error code.
    */
    virtual error createEdge(const int* const* vlist, const int* const* vplist, 
        const float* terms, int N, dd_edge &e) = 0;


    /** Create an edge for a boolean constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.
        @return       An appropriate error code.
    */
    virtual error createEdge(bool val, dd_edge &e) = 0;

    /** Create an edge for an integer constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.
        @return       An appropriate error code.
    */
    virtual error createEdge(int val, dd_edge &e) = 0;

    /** Create an edge for an integer constant.
        @param  val   Requested constant.
        @param  e     returns a handle to a node in the forest for
                      function f = \a val.
        @return       An appropriate error code.
    */
    virtual error createEdge(float val, dd_edge &e) = 0;

    /** Evaluate the function encoded by an edge.
        The forest must not be a relation, and must have range type of boolean.
        @param  f     Edge (function) to evaluate.
        @param  vlist List of variable assignments, of dimension one higher
                      than the largest variable handle.
        @param  term  Output parameter, will be set to
                      f(v1 = vlist[v1], ..., vn = vlist[vn]).
        @return       An appropriate error code.
    */
    virtual error evaluate(const dd_edge &f, const int* vlist, bool &term)
      const = 0;

    /** Evaluate the function encoded by an edge.
        The forest must not be a relation, and must have range type of integer.
        @param  f     Edge (function) to evaluate.
        @param  vlist List of variable assignments, of dimension one higher
                      than the largest variable handle.
        @param  term  Output parameter, will be set to
                      f(v1 = vlist[v1], ..., vn = vlist[vn]).
        @return       An appropriate error code.
    */
    virtual error evaluate(const dd_edge &f, const int* vlist, int &term)
      const = 0;

    /** Evaluate the function encoded by an edge.
        The forest must not be a relation, and must have range type of real.
        @param  f     Edge (function) to evaluate.
        @param  vlist List of variable assignments, of dimension one higher
                      than the largest variable handle.
        @param  term  Output parameter, will be set to
                      f(v1 = vlist[v1], ..., vn = vlist[vn]).
        @return       An appropriate error code.
    */
    virtual error evaluate(const dd_edge &f, const int* vlist, float &term)
      const = 0;

    /** Evaluate the function encoded by an edge.
        The forest must be a relation, and must have range type of boolean.
        @param  f       Edge (function) to evaluate.
        @param  vlist   List of variable assignments for unprimed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  vplist  List of variable assignments for primed variables, 
                        of dimension one higher than the largest variable
                        handle.
        @param  term    Output parameter, will be set to
                        f(v1 = vlist[v1], v'1 = vplist[v1], ...,
                        vn = vlist[vn], v'n = vplist[vn]).
        @return         An appropriate error code.
    */
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const = 0;

    /** Evaluate the function encoded by an edge.
        The forest must be a relation, and must have range type of integer.
        @param  f       Edge (function) to evaluate.
        @param  vlist   List of variable assignments for unprimed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  vplist  List of variable assignments for primed variables, 
                        of dimension one higher than the largest variable
                        handle.
        @param  term    Output parameter, will be set to
                        f(v1 = vlist[v1], v'1 = vplist[v1], ...,
                        vn = vlist[vn], v'n = vplist[vn]).
        @return         An appropriate error code.
    */
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, int &term) const = 0;

    /** Evaluate the function encoded by an edge.
        The forest must be a relation, and must have range type of real.
        @param  f       Edge (function) to evaluate.
        @param  vlist   List of variable assignments for unprimed variables,
                        of dimension one higher than the largest variable
                        handle.
        @param  vplist  List of variable assignments for primed variables, 
                        of dimension one higher than the largest variable
                        handle.
        @param  term    Output parameter, will be set to
                        f(v1 = vlist[v1], v'1 = vplist[v1], ...,
                        vn = vlist[vn], v'n = vplist[vn]).
        @return         An appropriate error code.
    */
    virtual error evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, float &term) const = 0;


    /** Returns a state from the MDD represented by \a f.
        @param  f       Edge.
        @param  vlist   Output parameter used to return a state from \a f.
        @return         INVALID_OPERATION if forest is not a MDD.
                        SUCCESS if a state is found,
                        INVALID_ASSIGNMENT otherwise.
    */
    virtual error findFirstElement(const dd_edge& f, int* vlist) const = 0;


    /** Returns a transition from the MXD represented by \a f.
        @param  f       Edge.
        @param  vlist   Output parameter used to return a transition from \a f.
        @param  vplist  Output parameter used to return a transition from \a f.
        @return         INVALID_OPERATION if forest is not a MXD.
                        SUCCESS if a transition is found,
                        INVALID_ASSIGNMENT otherwise.
    */
    virtual error findFirstElement(const dd_edge& f, int* vlist, int* vplist)
        const = 0;


    /** Display all active (i.e., connected) nodes in the forest.
        This is primarily for aid in debugging.
        @param  strm      File stream to write to.
        @param  verbosity How much information to display.
                          0 : just statistics.
                          1 : all forest nodes + statistics.
                          2 : internal forest + statistics.
    */
    virtual void showInfo(FILE* strm, int verbosity=0) = 0;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                          domain class                          *
// *                                                                *
// *                                                                *
// ******************************************************************


/** Domain class.
    Abstract base class.
    A domain is an ordered collection of variables,
    along with a rich set of operations for adding and removing variables.
    Variables are described by "handles" which have the following properties:
    - variable handles are non-negative integers.
    - when a new variable is created, it is given the lowest unused handle.
    - there is a special variable, TERMINALS, that refers to terminal nodes.
  
    When a domain is destroyed, all of its forests are destroyed.
*/
class domain {
  public:
    /// Special variable handle for the terminal nodes.
    static const int TERMINALS = 0;

    /** Error codes for domain methods.
        "SUCCESS" is fixed to the code 0, so we can
        test for errors using:  if (err) { HandleError(); }
    */
    enum error {
      /// Returned whenever an operation is successful.
      SUCCESS=0,
      /// Returned whenever the requested operation is not yet implemented!
      NOT_IMPLEMENTED,
      /// Returned whenever an operation fails due to lack of memory.
      INSUFFICIENT_MEMORY,
      /// We expected an empty domain, but it wasn't
      DOMAIN_NOT_EMPTY,
      /// Returned if a provided variable handle is erroneous.
      INVALID_VARIABLE,
      /// Returned if a provided variable bound is out of range.
      INVALID_BOUND
    };

    /** Create all variables at once, from the bottom up.
        Requires the domain to be "empty" (containing no variables or forests).
        When successful, variable handle \a N will refer to the top-most 
        variable, and variable handle 1 will refer to the bottom-most variable.

        @param  bounds  variable bounds.
                        bounds[i-1] gives the bound for variable i.
        @param  N       Number of variables.
        @return         An appropriate error code.
    */
    virtual error createVariablesBottomUp(const int* bounds, int N) = 0;

    /** Create a forest in this domain.
        Conceptually, a forest is a structure used to represent a
        collection of functions over a common domain. For simplicity
        (although it is a slight abuse of notation) a forest may represent
        "vectors" or "sets" over a domain, or "matrices" or "relations".
        In case of matrices / relations, the forest uses primed and unprimed
        versions of every variable in the domain. 

        @param  rel     Is this a relation / matrix, versus a set / vector.
        @param  t       Range type of the functions, namely,
                        booleans, integers, or reals.
        @param  ev      Edge labeling mechanism, i.e., should this be a
                        Multi-terminal decision diagram forest,
                        edge-valued with plus/times decision diagram forest.
        @return 0       if an error occurs, a new forest otherwise.
    */
    virtual forest* createForest(bool rel, forest::range_type t,
        forest::edge_labeling ev) = 0;

    /// Get the number of variables in this domain.
    /// Note that TERMINALS is not counted as a level.
    /// So, if the domain has N variables in it (excluding the
    /// TERMINAL level which is always present), N is returned.
    virtual int getNumVariables() const = 0;

    /** Get the specified bound of a variable.
      @param  vh      Variable handle.
      @param  prime   If prime is true, get the bound for the primed variable.
      @return         The bound set for variable \a vh.
                      If \a vh is TERMINALS, returns 0.
                      If \a vh is invalid, returns -1.
     */
    virtual int getVariableBound(int vh, bool prime = false) const = 0;

    /** Get the topmost variable.
        @return         The variable handle for the top-most variable.
                        If there are no variables, returns TERMINALS.
    */
    virtual int getTopVariable() const = 0;

    /** Get the variable immediately above this one.
        @param  vh      Any variable handle.
        @return         The variable appearing on top of this one. If \a vh
                        is already the top-most variable, returns -1.
    */
    virtual int getVariableAbove(int vh) const = 0;

    /** Get the variable immediately below this one.
        @param  vh      Any variable handle.
        @return         The variable appearing below this one. If \a vh is the
                        bottom-most variable, returns \a TERMINALS. If \a vh is
                        \a TERMINALS, returns -1.
    */
    virtual int getVariableBelow(int vh) const = 0;

    /// Get the number of forests associated with this domain.
    virtual int getNumForests() const = 0;

    /** Display lots of information about the domain.
        This is primarily for aid in debugging.
        @param  strm    File stream to write to.
    */
    virtual void showInfo(FILE* strm) = 0;

    /** Produce a human-readable name for an error code.
        @param  e   Error code.
        @return     A string that describes the error, suitable for printing.
    */
    static const char* getErrorCodeName(error e);

    /// Destructor.
    virtual ~domain();

  protected:
    /// Constructor, creates an empty domain.
    domain();
};

/** Front-end function to create a domain.
    This is required because domain is an abstract class.
*/
domain* MEDDLY_createDomain();


// ******************************************************************
// *                                                                *
// *                                                                *
// *                         dd_edge  class                         *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Structure for handles to edges in forests.

    A dd_edge is a handle for user manipulation of functions stored in forests.

    There are a few useful operations that can be applied directly
    to a dd_edge; all the rest are done either through the "parent" forest,
    or through operations in the compute_manager class.  These include:
   
    - Deletion of a dd_edge.  This will cause the parent forest
      to recycle nodes as appropriate.
    
    - Checking for equality of two dd_edges, using the method equals().
*/
class expert_forest;

class dd_edge {
  public:
    /** Constructor.
        Creates an empty edge in forest \a p.
        @param  p     forest to which this dd_edge will belong to.
    */
    dd_edge(forest* p);

    /// Destructor.  Will notify parent as appropriate.
    ~dd_edge();

    /** Copy Constructor.
        @param  e       dd_edge to copy.
    */
    dd_edge(const dd_edge &e);

    /** Clears the contents of this edge. It will belong to the same
        forest as before.
    */
    void clear();

    /** Obtain a modifiable copy of the forest owning this edge.
        @return         the forest to which this edge belongs to.
    */
    forest* getForest() const;

    /** Get this dd_edge's node handle.
        @return         the node handle.
    */
    int getNode() const;

    /** Get this dd_edge's edge value (only valid for edge-valued MDDs).
        Note: EV+MDDs use Integer edge-values while EV*MDDs use
        Real edge-value, so use the appropriate method.
        @return         edge value (if applicable).
    */
    void getEdgeValue(int& edgeValue) const;
    void getEdgeValue(float& edgeValue) const;

    /** Get this dd_edge's level handle.
        @return         the level handle.
    */
    int getLevel() const;

    /** Get node cardinality.
        Provided for backward compatibility.
        Use apply(CARDINALITY, ...) instead.
        @return         the cardinality of the node.
    */
    double getCardinality() const;

    /** Counts the number of unique nodes in this decision diagram.
        @return         the number of unique nodes starting at the root node
                        of this dd_edge.
    */
    unsigned getNodeCount() const;

    /** Counts the number of unique edges in this decision diagram.
        @param  countZeroes
                        if true, the stored zero edges are also counted
                        (sparse nodes do not store zero edges, so this
                        does not effect them; truncated nodes do store
                        some zero edges, so those edges will be counted).
        @return         the number of unique edges starting at the root node
                        of this dd_edge.
    */
    unsigned getEdgeCount(bool countZeroes = false) const;

    /** Modifies the dd_edge fields.
        The dd_edge is cleared (it will still belong to the same forest),
        and the dd_edge data is filled with the data provided.
        @param  node    node handle.
        @param  value   value of edge coming into the node (only useful
                        for edge-valued MDDs)
        @param  level   level handle.
    */
    void set(int node, int value, int level);
    void set(int node, float value, int level);

    class iterator {
      public:
        enum iter_type {
          DEFAULT=0,
          ROW,
          COLUMN
        };
        iterator();
        ~iterator();
        iterator(const iterator& iter);
        iterator& operator=(const iterator& iter);
        operator bool() const;
        void operator--();
        void operator++();
        bool operator!=(const iterator& iter) const;
        bool operator==(const iterator& iter) const;
        const int* getAssignments() const;
        const int* getPrimedAssignments() const;
        // Highest level at which the current minterm differs
        // from the previous minterm.
        int getLevel() const;
        void getValue(int& edgeValue) const;
        void getValue(float& edgeValue) const;
        
        iterator(dd_edge* e, iter_type t, const int* minterm);
        bool findFirstColumn(int height, int node);
        bool findNextColumn(int height);
        bool findFirstRow(int height, int node);
        bool findNextRow(int height);

      private:
        friend class dd_edge;
        void incrNonRelation();
        void incrRelation();
        void incrNonIdentRelation();
        void incrRow();
        void incrColumn();
        void incrNonIdentRow();
        void incrNonIdentColumn();

        dd_edge*  e;
        unsigned  size;
        int*      element;
        int*      nodes;
        int*      pelement;
        int*      pnodes;
        iter_type type;
        int       foundPathAtLevel;
    };

    typedef iterator const_iterator;

    /** Returns an iterator to the first element of the dd_edge.
        The iterator can be used to visit the elements in the DD in
        lexicographic order.
        @return         an iterator pointing to the first element.
    */
    const_iterator begin();

    /** Returns an iterator to the first element of the dd_edge
        with minterm as the "from" component of the element
        (the entire element being from->to).

        This iterator can then be used to visit the elements
        in the DD with the same "from" component.
        The elements are visited in lexicographic order of the
        "to" component.
        
        This is only valid for DD that store relations.

        If the relation is thought of as a matrix with the Y-axis
        representing the "from" components and the X-axis representing
        the "to" components, this iterator is useful for visiting all
        the "to"s that correspond to a single "from".

        @return         an iterator pointing to the first element.
    */
    const_iterator beginRow(const int* minterm);

    /** Returns an iterator to the first element of the dd_edge
        with minterm as the "to" component of the element
        (the entire element being to->from).

        This iterator can then be used to visit the elements
        in the DD with the same "to" component.
        The elements are visited in lexicographic order of the
        "from" component.
        
        This is only valid for DD that store relations.

        If the relation is thought of as a matrix with the Y-axis
        representing the "from" components and the X-axis representing
        the "to" components, this iterator is useful for visiting all
        the "from"s that correspond to a single "to".

        @return         an iterator pointing to the first element.
    */
    const_iterator beginColumn(const int* minterm);

    /** Assignment operator.
        @param  e       dd_edge to copy.
        @return         the new dd_edge.
    */
    dd_edge& operator=(const dd_edge &e);

    /** Check for equality.
        @return true    iff this edge has the same parent and refers to
                        the same edge as \a e.
    */
    bool operator==(const dd_edge& e) const;

    /** Check for inequality.
        @return true    iff this edge does not refer to the same edge as \a e.
    */
    bool operator!=(const dd_edge& e) const;

    /** Plus operator.
        BOOLEAN forests: Union; INTEGER/REAL forests: Addition.
        @param  e       dd_edge to Union/Add with this dd_edge.
        @return         \a this + \a e.
    */
    const dd_edge operator+(const dd_edge& e) const;

    /** Compound Plus operator.
        BOOLEAN forests: Union; INTEGER/REAL forests: Addition.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge to Union/Add with this dd_edge.
        @return         \a this + \a e.
    */
    dd_edge& operator+=(const dd_edge &e);

    /** Star operator.
        BOOLEAN forests: Intersection; INTEGER/REAL forests: Multiplication.
        @param  e       dd_edge to Intersection/Multiply with this dd_edge.
        @return         \a this * \a e.
    */
    const dd_edge operator*(const dd_edge& e) const;

    /** Compound Star operator.
        BOOLEAN forests: Intersection; INTEGER/REAL forests: Multiplication.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge to Intersection/Multiply with this dd_edge.
        @return         \a this * \a e.
    */
    dd_edge& operator*=(const dd_edge &e);

    /** Minus operator.
        BOOLEAN forests: Difference; INTEGER/REAL forests: Subtraction.
        @param  e       dd_edge for difference/subtract.
        @return         \a this - \a e.
    */
    const dd_edge operator-(const dd_edge& e) const;

    /** Compound Minus operator.
        BOOLEAN forests: Difference; INTEGER/REAL forests: Subtraction.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge for difference/subtract.
        @return         \a this - \a e.
    */
    dd_edge& operator-=(const dd_edge &e);

    /** Divide operator.
        BOOLEAN forests: INVALID; INTEGER/REAL forests: Division.
        @param  e       dd_edge for division.
        @return         \a this / \a e.
    */
    const dd_edge operator/(const dd_edge& e) const;

    /** Compound Divide operator.
        BOOLEAN forests: INVALID; INTEGER/REAL forests: Division.
        This edge is overwritten with the result of the operation.
        @param  e       dd_edge for division.
        @return         \a this / \a e.
    */
    dd_edge& operator/=(const dd_edge &e);

#if 0
    //  Not implemented
    /** Less-than operator.
        BOOLEAN forests: INVALID; INTEGER/REAL forests: Less-than.
        This returns an MTMDD containing all elements whose terminal value
        in \a this dd_edge is less than the terminal value in dd_edge \a e.

        Example:
        Let's assume A and B are dd_edges belonging to the same forest with
        INTEGERs for terminal values. Let there be two variables in this
        forest and each variable can take up values 0 and 1. An element is
        represented as (2nd variable, 1st variable, terminal value)

        Assuming,
        A = {(0, 0, 1), (0, 1, 4), (1, 0, 0), (1, 1, -1)}
        B = {(0, 0, 2), (0, 1, 0), (1, 0, 3), (1, 1, 0)}
        A < B = {(0, 0, 1), (0, 1, 0), (1, 0, 1), (1, 1, 1)}
        where 0 and 1 represent False and True respectively.

        @param  e       dd_edge for less-than.
        @return         dd_edge representing \a this < \a e.
    */
    const dd_edge operator<(const dd_edge& e) const;
#endif

    /** Display the edge information.
        This is primarily for aid in debugging.
        Note that cardinality only works for MDDs, MTMDDs, MXDs and MTMXDs.
        @param  strm      File stream to write to.
        @param  verbosity 0: default
                          1: default + cardinality
                          2: default + displays graph rooted at this node.
                          3: default + cardinality + graph.
    */
    void show(FILE* strm, int verbosity = 0) const;

  private:
    friend class expert_forest;
    friend class iterator;

    void setIndex(int index);
    int getIndex() const;

    forest *parent;
    int node;
    int value;
    int level;
    int index;

    op_info* opPlus;
    op_info* opStar;
    op_info* opMinus;
    op_info* opDivide;

    void updateIterators();

    bool            updateNeeded;
    const_iterator* beginIterator;
};


// ******************************************************************
// *                                                                *
// *                                                                *
// *                     compute_manager  class                     *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Class that handles operations on functions in forests.
    Abstract base class.

    Substantial operations on decision diagrams are here.
    I.e., this class manages operations on dd_edges
    in (common or different) forests that produce new dd_edges.
    
  
    TODO: do we need a compute table cleanup interface here,
    or compute table policy settings here?

*/
class compute_manager {
  public:
    /// built-in operators.
    enum op_code {
      /** Create a copy of a dd_edge.
          The copy may be stored in any forest as long as it belongs to the
          same domain as the original and the transformation is valid.
          Copying is valid with the following:
          MDD to MTMDD, MTMDD to MDD, MXD to MTMXD, MTMXD to MXD.
      */
      COPY=0,

      /// Unary operation.  Return the number of variable assignments 
      /// so that the function evaluates to non-zero.
      CARDINALITY,
      /// Set operation for forests with range_type of BOOLEAN. All operands
      /// must belong to the same forest.
      UNION,
      /// Set operation for forests with range_type of BOOLEAN. All operands
      /// must belong to the same forest.
      INTERSECTION,
      /// Set operation for forests with range_type of BOOLEAN. All operands
      /// must belong to the same forest.
      DIFFERENCE,
      /// Set operation for forests with range_type of BOOLEAN. All operands
      /// must belong to the same forest.
      COMPLEMENT,

      /// Binary operation.  Combines two functions into a single one,
      /// where the operands are MDDs and the result is an MXD.
      /// Specifically, for MDD operands f and g, produces MXD h where
      /// h(xn, x'n, ..., x1, x'1) = f(xn, ..., x1) * g(x'n, ..., x'1)
      /// Works for BOOLEAN forests.
      CROSS,

      /// Unary operation for forests with range_type of INTEGER or REAL.
      /// Return the largest value returned by the function.
      MAX_RANGE,
      /// Unary operation for forests with range_type of INTEGER or REAL.
      /// Return the smallest value returned by the function.
      MIN_RANGE,

      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      MIN,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      MAX,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      PLUS,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      MINUS,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      MULTIPLY,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      DIVIDE,

      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      EQUAL,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      NOT_EQUAL,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      LESS_THAN,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      LESS_THAN_EQUAL,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      GREATER_THAN,
      /// For forests with range_type of INTEGER and REAL. All operands must
      /// belong to the same forest.
      GREATER_THAN_EQUAL,

      /** Image operations on a set-of-states with a transition function.
          The first operand must be the set-of-states and the second operand
          must be the transition function. The result is a set-of-states that
          must be stored in the same forest as the first operand.
          
          Applies to:
          PRE_IMAGE, POST_IMAGE, REACHABLE_STATES_DFS, REACHABLE_STATES_BFS,
          REVERSE_REACHABLE_DFS, REVERSE_REACHABLE_BFS.
      */
      PRE_IMAGE,
      POST_IMAGE,
      REACHABLE_STATES_DFS,
      REACHABLE_STATES_BFS,
      REVERSE_REACHABLE_DFS,
      REVERSE_REACHABLE_BFS,

      /// Convert MDD to EV+MDD index set
      CONVERT_TO_INDEX_SET

      // there will be many more codes...
    };

    /** Error codes for compute_manager methods.
        "SUCCESS" is fixed to the code 0, so we can
        test for errors using:  if (err) { HandleError(); }
    */
    enum error {
      /// Returned whenever an operation is successful.
      SUCCESS=0,
      /// Returned whenever the requested operation is not yet implemented!
      NOT_IMPLEMENTED,
      /// Returned whenever an operation fails due to lack of memory.
      INSUFFICIENT_MEMORY,
      /// Unknown operation (bad operation handle).
      UNKNOWN_OPERATION,
      /// Requested operation requires same forest, it wasn't.
      FOREST_MISMATCH,
      /// Requested operation not supported for forest, operand, or result type.
      TYPE_MISMATCH,
      /// Requested operation requires different number of operands.
      WRONG_NUMBER,
      /// If we are counting something and it won't fit in an integer / float.
      OVERFLOW
    };

  public:
    compute_manager();
    virtual ~compute_manager();

    /** Removes all cached computation results from the compute table.
    */
    virtual void clearComputeTable() = 0;

    /** Change the type of hash table being used for the compute cache.
        The compute cache by default uses a hash table with chaining (i.e. a
        list of entries at each hash location). This can be changed to a
        hash-table without chaining. The maximum size of the hash table can
        also be fixed (to limit the amount of memory it uses).
        This operation is valid as long as nothing has been added to the cache
        previously.
        @param  chaining
                      enables/disables chaining in hash table.
        @param  size  size of cache (default 16777216).
        @return       An appropriate error code
    */
    virtual error setHashTablePolicy(bool chaining, unsigned size = 16777216u)
        = 0;

    /** Produce a human-readable name for an error code.
        @param  e     Error code.
        @return       A string that describes the error, suitable for printing.
    */
    static const char* getErrorCodeName(error e);

    /** Obtain a human-readable name for an operation.
        @param  op    Operator handle.
        @return       A human-readable string constant (e.g., "Union").
    */
    virtual const char* getOperationName(op_code op) const = 0;

    /** Apply a unary operator.
        The operand and the result are not necessarily in the same forest,
        but they must belong to forests that share the same domain.
        This is useful, for instance, for copying a function to a new forest.
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Output parameter: the result, where \a c = \a op \a a.
        @return       An appropriate error code.
    */
    virtual error apply(op_code op, const dd_edge &a, dd_edge &c) = 0;

    /** Apply a unary operator.
        For operators whose result is an integer.
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Output parameter: the result, where \a c = \a op \a a.
        @return       An appropriate error code.
    */
    virtual error apply(op_code op, const dd_edge &a, long &c) = 0;

    /** Apply a unary operator.
        For operators whose result is a real.
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Output parameter: the result, where \a c = \a op \a a.
        @return       An appropriate error code.
    */
    virtual error apply(op_code op, const dd_edge &a, double &c) = 0;

    virtual error apply(op_code op, const dd_edge &a, ct_object &c) = 0;

#ifdef __GMP_H__
    /** Apply a unary operator.
        For operators whose result is an arbitrary-precision integer
        (as supplied by the GNU MP library).
        @param  op    Operator handle.
        @param  a     Operand.
        @param  c     Input: an initialized MP integer.
                      Output parameter: the result, where \a c = \a op \a a.
        @return       An appropriate error code.
    */
    inline error apply(op_code op, const dd_edge &a, mpz_t &c) {
      ct_object& x = MEDDLY_get_mpz_wrapper();
      error e = apply(op, a, x);
      MEDDLY_unwrap(x, c);
      return e; 
    }
#endif


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
        @return       An appropriate error code.
    */ 
    virtual error apply(op_code op, const dd_edge &a, const dd_edge &b,
        dd_edge &c) = 0;


    /**
        Computes y = y + xA.
        x and y are vectors, stored explicitly, and A is a matrix.

        @param  y       Vector; dimension must be enough for largest y index.
        @param  y_ind   Function to determine how minterms are mapped
                        to indexes for vector y.  A value of infinity
                        can be used to ignore minterms.  Should be an
                        EV+MDD.

        @param  x       Vector; dimension must be enough for largest x index.
        @param  x_ind   Function to determine how minterms are mapped
                        to indexes for vector x.  A value of infinity
                        can be used to ignore minterms.  Should be an
                        EV+MDD.

        @param  A       Real-valued matrix, as an MTMxD or EV*MxD with
                        same domain as y_ind and x_ind.

        @return         An appropriate error code.  
                        TBD: list the possible errors.
    */
    virtual error vectorMatrixMultiply(double* y, const dd_edge &y_ind,
                  const double* x, const dd_edge &x_ind, const dd_edge &A) = 0;

    /**
        Computes y = y + Ax.
        x and y are vectors, stored explicitly, and A is a matrix.

        @param  y       Vector; dimension must be enough for largest y index.
        @param  y_ind   Function to determine how minterms are mapped
                        to indexes for vector y.  A value of infinity
                        can be used to ignore minterms.  Should be an
                        EV+MDD.

        @param  A       Real-valued matrix, as an MTMxD or EV*MxD with
                        same domain as y_ind and x_ind.

        @param  x       Vector; dimension must be enough for largest x index.
        @param  x_ind   Function to determine how minterms are mapped
                        to indexes for vector x.  A value of infinity
                        can be used to ignore minterms.  Should be an
                        EV+MDD.

        @return         An appropriate error code.  
                        TBD: list the possible errors.
    */
    virtual error matrixVectorMultiply(double* y, const dd_edge &y_ind,
                  const dd_edge &A, const double* x, const dd_edge &x_ind) = 0;

    /** Display compute table information.
        This is primarily for aid in debugging.
        @param  strm  File stream to write to.
    */
    virtual void showComputeTable(FILE* strm) const = 0;

    /** Get the number of entries in the compute table.
        @return       The number of entries in the compute table.
    */
    virtual long getNumCacheEntries() const = 0;
};

/** Function to build and initialize the compute manager.
    Built-in operations are initialized here.
    Multiple calls will return the same compute_manager.
*/
compute_manager* MEDDLY_getComputeManager();






























// **********************************************************************
//
//                    Inlined methods for dd_edge
//
// **********************************************************************

inline
forest* dd_edge::getForest() const
{
  return parent;
}

inline
int dd_edge::getNode() const
{
  return node;
}

inline
void dd_edge::getEdgeValue(int& ev) const
{
  ev = value;
}

inline
int dd_edge::getLevel() const
{
  return level;
}

inline
void dd_edge::setIndex(int index)
{
  this->index = index;
}

inline
int dd_edge::getIndex() const
{
  return index;
}

inline
void dd_edge::clear()
{
  assert(index != -1);
  set(0, 0, 0);
  updateNeeded = true;
}

inline
bool dd_edge::operator==(const dd_edge& e) const
{
  return (this == &e) ||
         (parent == e.parent && node == e.node &&
          value == e.value && level == e.level);
}

inline
bool dd_edge::operator!=(const dd_edge& e) const
{
  return !(*this == e);
}

inline
const dd_edge dd_edge::operator+(const dd_edge& e) const
{
  return dd_edge(*this) += e;
}

inline
const dd_edge dd_edge::operator*(const dd_edge& e) const
{
  return dd_edge(*this) *= e;
}

inline
const dd_edge dd_edge::operator-(const dd_edge& e) const
{
  return dd_edge(*this) -= e;
}

inline
dd_edge::iterator::operator bool() const
{
  return nodes != 0 && nodes[0] != 0;
}

inline
int dd_edge::iterator::getLevel() const
{
  return foundPathAtLevel;
}

#endif
