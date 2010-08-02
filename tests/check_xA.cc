
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

/*
    Tests vector matrix multiplication.
*/


#include <stdio.h>

#include "meddly.h"

const int vars[] = {0, 10, 10, 10};

compute_manager* CM;

bool build_oz(forest* indf, forest* mxd, dd_edge &ss, dd_edge &P)
{
  const int R[] = {0, 3, 1, 4};
  const int N[] = {0, 1, 4, 9};
  const int S[] = {0, 2, 6, 5};

  // Build state indexes

  const int* sslist[] = { R, N, S };
  const int indexes[] = { 0, 1, 2 };

  forest::error fe;
  fe = indf->createEdge(sslist, indexes, 3, ss);
  if (fe) {
    printf("Couldn't build indexes: %s\n", indf->getErrorCodeName(fe));
    return false;
  }

  // Build matrix elements

  const int* fromlist[] = {
    R,    R,    R,    N,    N,    S,    S,    S
  };
  const int* tolist[] = {
    R,    N,    S,    R,    S,    R,    N,    S
  };
  const float problist[] = {
    0.5,  0.25, 0.25, 0.5,  0.5,  0.25, 0.25, 0.5
  };

  
  fe = mxd->createEdge(fromlist, tolist, problist, 8, P);
  if (fe) {
    printf("Couldn't build matrix: %s\n", mxd->getErrorCodeName(fe));
    return false;
  }

  return true;
}

void by_hand_oz_xA(double* y, const double* const x)
{
  y[0] += x[0]/2 + x[1]/2 + x[2]/4;
  y[1] += x[0]/4 +          x[2]/4;
  y[2] += x[0]/4 + x[1]/2 + x[2]/2;
}

void by_hand_oz_Ax(double* y, const double* const x)
{
  y[0] += x[0]/2 + x[1]/4 + x[2]/4;
  y[1] += x[0]/2 +          x[2]/2;
  y[2] += x[0]/4 + x[1]/4 + x[2]/2;
}

bool equals(double *xgot, double *xans, int N)
{
  for (int i=0; i<N; i++) {
    double delta = xgot[i] - xans[i];
    if (xgot[i]) delta /= xgot[i];
    if (delta > -1e-5 && delta < 1e-5) continue;
    // not equal, print an error
    printf("Got     : [%lf", xgot[0]);
    for (int j=1; j<N; j++) printf(", %lf", xgot[j]);
    printf("]\nExpected: [%lf", xans[0]);
    for (int j=1; j<N; j++) printf(", %lf", xans[j]);
    printf("]\n");
    return false;
  }
  return true;
}

bool xA_check(dd_edge &ss, dd_edge &P)
{
  int i;
  double p[3];
  double q[3];
  double q_alt[3];
  p[0] = 0; p[1] = 1; p[2] = 0;
  printf("xA multiplications:\n");
  for (i=0; i<9; i++) {
    printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
    q[0] = q[1] = q[2] = 0;
    compute_manager::error ce = CM->vectorMatrixMultiply(q, ss, p, ss, P);
    if (ce) {
      printf("Couldn't multiply: %s\n", CM->getErrorCodeName(ce));
      return false;
    }
    q_alt[0] = q_alt[1] = q_alt[2] = 0;
    by_hand_oz_xA(q_alt, p);
    if (!equals(q, q_alt, 3)) return false;
    p[0] = q[0];
    p[1] = q[1];
    p[2] = q[2];
  }
  printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
  return true;
}

bool Ax_check(dd_edge &ss, dd_edge &P)
{
  int i;
  double p[3];
  double q[3];
  double q_alt[3];
  p[0] = 0; p[1] = 1; p[2] = 0;
  printf("Ax multiplications:\n");
  for (i=0; i<9; i++) {
    printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
    q[0] = q[1] = q[2] = 0;
    compute_manager::error ce = CM->matrixVectorMultiply(q, ss, P, p, ss);
    if (ce) {
      printf("Couldn't multiply: %s\n", CM->getErrorCodeName(ce));
      return false;
    }
    q_alt[0] = q_alt[1] = q_alt[2] = 0;
    by_hand_oz_Ax(q_alt, p);
    if (!equals(q, q_alt, 3)) return false;
    p[0] = q[0];
    p[1] = q[1];
    p[2] = q[2];
  }
  printf("p%d: [%lf, %lf, %lf]\n", i, p[0], p[1], p[2]);
  return true;
}

int main(int argc, const char** argv)
{
  CM = MEDDLY_getComputeManager();
  assert(CM);

  domain* ozd = MEDDLY_createDomain();
  assert(ozd);
  assert(domain::SUCCESS == ozd->createVariablesBottomUp(vars, 3));
  forest* evpmdds = ozd->createForest(0, forest::INTEGER, forest::EVPLUS);
  assert(evpmdds);
  forest* mtmxds = ozd->createForest(1, forest::REAL, forest::MULTI_TERMINAL);
  assert(mtmxds);

  dd_edge ss(evpmdds);
  dd_edge P(mtmxds);

  if (!build_oz(evpmdds, mtmxds, ss, P)) return 1;
  if (!xA_check(ss, P)) return 1;
  if (!Ax_check(ss, P)) return 1;

  return 0;
}
