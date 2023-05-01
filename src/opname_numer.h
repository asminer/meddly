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

#ifndef MEDDLY_OPNAME_NUMER_H
#define MEDDLY_OPNAME_NUMER_H

#include "opname.h"

namespace MEDDLY {
    class numerical_opname;
    class dd_edge;

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
// *                inlined numerical_opname methods                *
// *                                                                *
// ******************************************************************

inline MEDDLY::specialized_operation*
MEDDLY::numerical_opname::buildOperation(const dd_edge &x_ind,
    const dd_edge &A, const dd_edge &y_ind) const
{
  numerical_args na(x_ind, A, y_ind);
  na.setAutoDestroy(false); // na will be destroyed when we return
  return buildOperation(&na);
}



#endif // #include guard
