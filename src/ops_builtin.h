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

#ifndef MEDDLY_OPS_BUILTIN_H
#define MEDDLY_OPS_BUILTIN_H

#include "defines.h"

namespace MEDDLY {
    class opname;
    class unary_opname;
    class binary_opname;
    class numerical_opname;

    class satpregen_opname;
    class satotf_opname;
    class satimpl_opname;
    class sathyb_opname;
    class constrained_opname;

    class initializer_list;

    // ******************************************************************
    // *                                                                *
    // *                    unary  operation handles                    *
    // *                                                                *
    // ******************************************************************

    /// Unary operation.  Return the number of variable assignments
    /// so that the function evaluates to non-zero.
    unary_opname* CARDINALITY();

    /// For BOOLEAN forests, flip the return values.
    unary_opname* COMPLEMENT();

    /// Convert MDD to EV+MDD index set.  A special case of COPY, really.
    unary_opname* CONVERT_TO_INDEX_SET();

    /// Copy a function across forests, with the same domain.
    unary_opname* COPY();

    /// Extract cycles (EV+MDD) from transitive closure (EV+MxD)
    unary_opname* CYCLE();

    unary_opname* EQUANT();
    /// Find the largest value returned by the function.
    unary_opname* MAX_RANGE();

    /// Find the smallest value returned by the function.
    unary_opname* MIN_RANGE();

    /// Randomly select one state from a set of states
    unary_opname* SELECT();

    unary_opname* EXTRACTFROM();

    unary_opname* EXTRACTCOVERED();

    // ******************************************************************
    // *                                                                *
    // *                    binary operation handles                    *
    // *                                                                *
    // ******************************************************************

    /// Set union operation for forests with range_type of BOOLEAN
    binary_opname* UNION();

    /// Set intersection operation for forests with range_type of BOOLEAN
    binary_opname* INTERSECTION();

    /// Set difference operation for forests with range_type of BOOLEAN
    binary_opname* DIFFERENCE();

    /// Combine two functions into a single one, where the operands are MDDs
    /// and the result is an MXD.  Specifically, for MDD operands f and g,
    /// produces MXD h where
    ///     h(xn, x'n, ..., x1, x'1) = f(xn, ..., x1) * g(x'n, ..., x'1).
    /// Works for BOOLEAN forests.
    binary_opname* CROSS();

    /// The minimum of two functions, with range_type INTEGER or REAL
    binary_opname* MINIMUM();

    /// The maximum of two functions, with range_type INTEGER or REAL
    binary_opname* MAXIMUM();

    /// Add two functions, with range type INTEGER and REAL
    binary_opname* PLUS();

    /// Subtract two functions, with range type INTEGER and REAL
    binary_opname* MINUS();

    /// Multiply two functions, with range type INTEGER and REAL
    binary_opname* MULTIPLY();

    /// Divide two functions, with range type INTEGER and REAL
    binary_opname* DIVIDE();

    /// Take the remainder of two functions, with range type INTEGER
    binary_opname* MODULO();

    /// Compare for equality, two functions with range type INTEGER or REAL
    binary_opname* EQUAL();

    /// Compare for inequality, two functions with range type INTEGER or REAL
    binary_opname* NOT_EQUAL();

    /// Compare for <, two functions with range type INTEGER or REAL
    binary_opname* LESS_THAN();

    /// Compare for <=, two functions with range type INTEGER or REAL
    binary_opname* LESS_THAN_EQUAL();

    /// Compare for >, two functions with range type INTEGER or REAL
    binary_opname* GREATER_THAN();

    /// Compare for >=, two functions with range type INTEGER or REAL
    binary_opname* GREATER_THAN_EQUAL();

    /** Plus operation used to compute transitive closure and further
        minimum witness. The first operand must be an EV+MxD and the second
        operand must be an EV+MDD. The result is an EV+MxD.
    */
    binary_opname* PRE_PLUS();
    binary_opname* POST_PLUS();

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
    binary_opname* PRE_IMAGE();
    binary_opname* POST_IMAGE();
    binary_opname* TC_POST_IMAGE();
    binary_opname* MRC_POST_IMAGE();
    binary_opname* REACHABLE_STATES_DFS();
    binary_opname* REACHABLE_STATES_BFS();
    binary_opname* REVERSE_REACHABLE_DFS();
    binary_opname* REVERSE_REACHABLE_BFS();
    binary_opname* COV_TC();

    /** Vector matrix multiply, where the first argument is vector (MDD),
        the second argument is a matrix (MXD),
        and the result is a vector (MDD).
    */
    binary_opname* VM_MULTIPLY();

    /** Matrix vector multiply, where the first argument is a matrix (MXD),
        the second argument is a vector (MDD),
        and the result is a vector (MDD).
    */
    binary_opname* MV_MULTIPLY();

    /** Matrix multiplication, where the first argument is a matrix (MXD),
        the second argument is a matrix (MXD),
        and the result is a matrix (MXD),
        such that, C[m][n] += A[m][i] * B[i][n], for all m, n and i.
    */
    binary_opname* MM_MULTIPLY();

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
    numerical_opname* EXPLVECT_MATR_MULT();
    // extern const numerical_opname* VECT_MATR_MULT; // renamed!

    /** Computes y = y + Ax.
        x and y are vectors, stored explicitly, and A is a matrix.
        x_ind and y_ind specify how minterms are mapped to indexes
        for vectors x and y, respectively.
    */
    numerical_opname* MATR_EXPLVECT_MULT();
    // extern const numerical_opname* MATR_VECT_MULT; // renamed!

    // ******************************************************************
    // *                                                                *
    // *                  Named saturation operations                   *
    // *                                                                *
    // ******************************************************************

    /** Forward reachability using saturation.
        Transition relation is already known.
    */
    satpregen_opname* SATURATION_FORWARD();

    /** Backward reachability using saturation.
        Transition relation is already known.
    */
    satpregen_opname* SATURATION_BACKWARD();

    /** Forward reachability using saturation.
        Transition relation is not completely known,
        will be built along with reachability set.
    */
    satotf_opname* SATURATION_OTF_FORWARD();

    satotf_opname* BFS_OTF_FORWARD();

    satotf_opname* COV();


    /** Forward reachability using saturation.
        Transition relation is specified implicitly.
    */
    satimpl_opname* SATURATION_IMPL_FORWARD();

    /** Forward reachability using saturation.
        Allows hybrid representation of transition relation.
    */
    sathyb_opname* SATURATION_HYB_FORWARD();

    /** Minimum-witness operations.
    */
    constrained_opname* CONSTRAINED_BACKWARD_BFS();
    constrained_opname* CONSTRAINED_FORWARD_DFS();
    constrained_opname* CONSTRAINED_BACKWARD_DFS();
    constrained_opname* TRANSITIVE_CLOSURE_DFS();

    // ******************************************************************
    // *                                                                *
    // *                   Initialization  front-end                    *
    // *                                                                *
    // ******************************************************************

    /// Build the initializer for builtins.
    initializer_list* makeBuiltinInitializer(initializer_list* prev);
};

#endif // #include guard
