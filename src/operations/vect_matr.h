
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

#ifndef VECT_MATR_H
#define VECT_MATR_H

// meddly.h must be included before this one.

class op_param;

/** Vector matrix multiplication for ev+mdd indexes, mtmxd matrices.
    Computes y += xA.

      @param  pt      Parameter types; holds forests for y_ind, x_ind, A.
      @param  ht      Height; outermost call should be number of levels.
      @param  y       y vector.
      @param  y_ind   EV+MDD node for y indexes.
      @param  x       x vector.
      @param  x_ind   EV+MDD node for x indexes.
      @param  A       MTMDD node for A matrix.

      @return         An appropriate error code.
*/
compute_manager::error vectorMatrixMult_evplus_mt(
  const op_param* pt, int ht, double* y, int y_ind, double* x, int x_ind, int A
);

/** Vector matrix multiplication for ev+mdd indexes, ev*mxd matrices.
    Computes y += xA.

      @param  pt      Parameter types; holds forests for y_ind, x_ind, A.
      @param  ht      Height; outermost call should be number of levels.
      @param  y       y vector.
      @param  y_ind   EV+MDD node for y indexes.
      @param  x       x vector.
      @param  x_ind   EV+MDD node for x indexes.
      @param  A       MTMDD node for A matrix.

      @return         An appropriate error code.
*/
compute_manager::error vectorMatrixMult_evplus_evtimes(
  const op_param* pt, int ht, double* y, int y_ind, double* x, int x_ind, int A
);




/** Matrix vector multiplication for ev+mdd indexes, mtmxd matrices.
    Computes y += Ax.

      @param  pt      Parameter types; holds forests for y_ind, x_ind, A.
      @param  ht      Height; outermost call should be number of levels.
      @param  y       y vector.
      @param  y_ind   EV+MDD node for y indexes.
      @param  A       MTMDD node for A matrix.
      @param  x       x vector.
      @param  x_ind   EV+MDD node for x indexes.

      @return         An appropriate error code.
*/
compute_manager::error matrixVectorMult_evplus_mt(
  const op_param* pt, int ht, double* y, int y_ind, int A, double* x, int x_ind
);

/** Matrix vector multiplication for ev+mdd indexes, ev*mxd matrices.
    Computes y += Ax.

      @param  pt      Parameter types; holds forests for y_ind, x_ind, A.
      @param  ht      Height; outermost call should be number of levels.
      @param  y       y vector.
      @param  y_ind   EV+MDD node for y indexes.
      @param  A       MTMDD node for A matrix.
      @param  x       x vector.
      @param  x_ind   EV+MDD node for x indexes.

      @return         An appropriate error code.
*/
compute_manager::error matrixVectorMult_evplus_evtimes(
  const op_param* pt, int ht, double* y, int y_ind, int A, double* x, int x_ind
);

#endif

