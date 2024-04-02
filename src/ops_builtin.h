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

#include "oper.h"


namespace MEDDLY {

    // ******************************************************************
    // *                                                                *
    // *                  "built-in"  unary operations                  *
    // *                                                                *
    // ******************************************************************

    class unary_operation;

    /// Return the number of variable assignments
    /// so that the function evaluates to non-zero.
    unary_operation* CARDINALITY(forest* arg, opnd_type res);

    /// For BOOLEAN forests, flip the return values.
    unary_operation* COMPLEMENT(forest* arg, forest* res);

    /// Convert MDD to EV+MDD index set.  A special case of COPY, really.
    unary_operation* CONVERT_TO_INDEX_SET(forest* arg, forest* res);

    /// Copy a function across forests, with the same domain.
    unary_operation* COPY(forest* arg, forest* res);

    /// Extract cycles (EV+MDD) from transitive closure (EV+MxD)
    unary_operation* CYCLE(forest* arg, forest* res);

    /// Find the largest value returned by the function.
    unary_operation* MAX_RANGE(forest* arg, opnd_type res);

    /// Find the smallest value returned by the function.
    unary_operation* MIN_RANGE(forest* arg, opnd_type res);

    /// Randomly select one state from a set of states
    unary_operation* SELECT(forest* arg, forest* res);

    // ******************************************************************
    // *                                                                *
    // *                  "built-in" binary operations                  *
    // *                                                                *
    // ******************************************************************

    class binary_operation;

    /// Set union operation for forests with range_type of BOOLEAN
    binary_operation* UNION(forest* a, forest* b, forest* c);

    /// Set intersection operation for forests with range_type of BOOLEAN
    binary_operation* INTERSECTION(forest* a, forest* b, forest* c);

    /// Set difference operation for forests with range_type of BOOLEAN
    binary_operation* DIFFERENCE(forest* a, forest* b, forest* c);

    /// Combine two functions into a single one, where the operands are MDDs
    /// and the result is an MXD.  Specifically, for MDD operands f and g,
    /// produces MXD h where
    ///     h(xn, x'n, ..., x1, x'1) = f(xn, ..., x1) * g(x'n, ..., x'1).
    /// Works for BOOLEAN forests.
    binary_operation* CROSS(forest* a, forest* b, forest* c);

    /// The minimum of two functions, with range_type INTEGER or REAL
    binary_operation* MINIMUM(forest* a, forest* b, forest* c);

    /// The maximum of two functions, with range_type INTEGER or REAL
    binary_operation* MAXIMUM(forest* a, forest* b, forest* c);

    /// Add two functions, with range type INTEGER and REAL
    binary_operation* PLUS(forest* a, forest* b, forest* c);

    /// Subtract two functions, with range type INTEGER and REAL
    binary_operation* MINUS(forest* a, forest* b, forest* c);

    /// Multiply two functions, with range type INTEGER and REAL
    binary_operation* MULTIPLY(forest* a, forest* b, forest* c);

    /// Divide two functions, with range type INTEGER and REAL
    binary_operation* DIVIDE(forest* a, forest* b, forest* c);

    /// Take the remainder of two functions, with range type INTEGER
    binary_operation* MODULO(forest* a, forest* b, forest* c);

    /// Compare for equality, two functions with range type INTEGER or REAL
    binary_operation* EQUAL(forest* a, forest* b, forest* c);

    /// Compare for inequality, two functions with range type INTEGER or REAL
    binary_operation* NOT_EQUAL(forest* a, forest* b, forest* c);

    /// Compare for <, two functions with range type INTEGER or REAL
    binary_operation* LESS_THAN(forest* a, forest* b, forest* c);

    /// Compare for <=, two functions with range type INTEGER or REAL
    binary_operation* LESS_THAN_EQUAL(forest* a, forest* b, forest* c);

    /// Compare for >, two functions with range type INTEGER or REAL
    binary_operation* GREATER_THAN(forest* a, forest* b, forest* c);

    /// Compare for >=, two functions with range type INTEGER or REAL
    binary_operation* GREATER_THAN_EQUAL(forest* a, forest* b, forest* c);

    /** Plus operation used to compute transitive closure and further
        minimum witness. The first operand must be an EV+MxD and the second
        operand must be an EV+MDD. The result is an EV+MxD.
    */
    binary_operation* PRE_PLUS(forest* a, forest* b, forest* c);
    binary_operation* POST_PLUS(forest* a, forest* b, forest* c);

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
    binary_operation* PRE_IMAGE(forest* a, forest* b, forest* c);
    binary_operation* POST_IMAGE(forest* a, forest* b, forest* c);
    binary_operation* TC_POST_IMAGE(forest* a, forest* b, forest* c);
    binary_operation* REACHABLE_STATES_DFS(forest* a, forest* b, forest* c);
    binary_operation* REACHABLE_STATES_BFS(forest* a, forest* b, forest* c);
    binary_operation* REVERSE_REACHABLE_DFS(forest* a, forest* b, forest* c);
    binary_operation* REVERSE_REACHABLE_BFS(forest* a, forest* b, forest* c);

    /** Vector matrix multiply, where the first argument is vector (MDD),
        the second argument is a matrix (MXD),
        and the result is a vector (MDD).
    */
    binary_operation* VM_MULTIPLY(forest* a, forest* b, forest* c);

    /** Matrix vector multiply, where the first argument is a matrix (MXD),
        the second argument is a vector (MDD),
        and the result is a vector (MDD).
    */
    binary_operation* MV_MULTIPLY(forest* a, forest* b, forest* c);

    /** Matrix multiplication, where the first argument is a matrix (MXD),
        the second argument is a matrix (MXD),
        and the result is a matrix (MXD),
        such that, C[m][n] += A[m][i] * B[i][n], for all m, n and i.
    */
    binary_operation* MM_MULTIPLY(forest* a, forest* b, forest* c);


    // ******************************************************************
    // *                                                                *
    // *                "built-in"  numerical operations                *
    // *                                                                *
    // ******************************************************************

    class numerical_operation;
    class dd_edge;

    numerical_operation* EXPLVECT_MATR_MULT(const dd_edge &xind,
            const dd_edge &A, const dd_edge &yind);

    numerical_operation* MATR_EXPLVECT_MULT(const dd_edge &xind,
            const dd_edge &A, const dd_edge &yind);

    // ******************************************************************
    // *                                                                *
    // *                  Named saturation operations                   *
    // *                                                                *
    // ******************************************************************

    class saturation_operation;
    class pregen_relation;
    class otf_relation;
    class implicit_relation;
    class hybrid_relation;

    saturation_operation* SATURATION_FORWARD(forest* inF,
            pregen_relation* nsf, forest* outf);
    saturation_operation* SATURATION_BACKWARD(forest* inF,
            pregen_relation* nsf, forest* outf);

    saturation_operation* SATURATION_OTF_FORWARD(forest* inF,
            otf_relation* nsf, forest* outf);
    saturation_operation* SATURATION_OTF_BACKWARD(forest* inF,
            otf_relation* nsf, forest* outf);

    saturation_operation* SATURATION_IMPL_FORWARD(forest* inF,
            implicit_relation* nsf, forest* outf);

    saturation_operation* SATURATION_HYB_FORWARD(forest* inF,
            hybrid_relation* nsf, forest* outf);

    // OLD BELOW

    // class satpregen_opname;
    class satotf_opname;
    class satimpl_opname;
    class sathyb_opname;
    class constrained_opname;

    /** Forward reachability using saturation.
        Transition relation is already known.
    */
    // satpregen_opname* SATURATION_FORWARD();

    /** Backward reachability using saturation.
        Transition relation is already known.
    */
    // satpregen_opname* SATURATION_BACKWARD();

    /** Forward reachability using saturation.
        Transition relation is not completely known,
        will be built along with reachability set.
    */
//     satotf_opname* SATURATION_OTF_FORWARD();

    /** Forward reachability using saturation.
        Transition relation is specified implicitly.
    */
    // satimpl_opname* SATURATION_IMPL_FORWARD();

    /** Forward reachability using saturation.
        Allows hybrid representation of transition relation.
    */
    // sathyb_opname* SATURATION_HYB_FORWARD();

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

    class initializer_list;

    /// Build the initializer for builtins.
    initializer_list* makeBuiltinInitializer(initializer_list* prev);
};

#endif // #include guard
