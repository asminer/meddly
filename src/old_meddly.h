
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
#include "defines.h"

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

  // ******************************************************************
  // *                     Named unary operations                     *
  // ******************************************************************

  /// Unary operation.  Return the number of variable assignments
  /// so that the function evaluates to non-zero.
  // extern const unary_opname* CARDINALITY;

  /// For BOOLEAN forests, flip the return values.
  // extern const unary_opname* COMPLEMENT;

  /// Find the largest value returned by the function.
  // extern const unary_opname* MAX_RANGE;

  /// Find the smallest value returned by the function.
  // extern const unary_opname* MIN_RANGE;

  /// Convert MDD to EV+MDD index set.  A special case of COPY, really.
  // extern const unary_opname* CONVERT_TO_INDEX_SET;

  /// Extract cycles (EV+MDD) from transitive closure (EV+MxD)
  // extern const unary_opname* CYCLE;

  /// Randomly select one state from a set of states
  // extern const unary_opname* SELECT;

  // ******************************************************************
  // *                    Named  binary operations                    *
  // ******************************************************************

#if 0
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
#endif

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
  // *                  object destruction functions                  *
  // ******************************************************************


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

#endif
