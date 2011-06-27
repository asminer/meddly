
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

#include "../defines.h"
#include "vect_matr.h"

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::numerical_opname* MEDDLY::initVectorMatrixMult(const settings &s)
{
  // TBD
  return 0;
}

MEDDLY::numerical_opname* MEDDLY::initMatrixVectorMult(const settings &s)
{
  // TBD
  return 0;
}

#if 0

using namespace MEDDLY;

void
vmprimed_evplus_mt(
  const op_param* pt, int ht, double* y, int y_ind, const double* x, int x_ind, int A
)
{
  // Handles the primed levels of A
  const expert_forest* fA = pt[4].readForest();
  assert(fA);
  if (0==ht) {
    y[0] += x[0] * fA->getReal(A);
    return;
  }
  const expert_forest* fy = pt[1].readForest();
  assert(fy);

  // It should be impossible for an indexing function to skip levels, right?   
  assert(fy->getNodeHeight(y_ind) == ht);
  int yv, ys;

  //
  // Sparse y_ind
  //
  if (fy->isSparseNode(y_ind)) {
    ys = fy->getSparseNodeSize(y_ind);

    if (fA->isSparseNode(A)) {
      // A is sparse, y_ind is sparse; intersect them
      int As = fA->getSparseNodeSize(A);
      int yp = 0;
      int yj = fy->getSparseNodeIndex(y_ind, yp);
      int Ap = 0;
      int Aj = fA->getSparseNodeIndex(A, Ap);
      for (;;) {
        if (Aj < yj) {
          Ap++;
          if (Ap >= As) break;
          Aj = fA->getSparseNodeIndex(A, Ap);
          continue;
        }
        if (yj < Aj) {
          yp++;
          if (yp >= ys) break;
          yj = fy->getSparseNodeIndex(y_ind, yp);
          continue;
        }
        fy->getSparseNodeEdgeValue(y_ind, yp, yv);
        vectorMatrixMult_evplus_mt(
          pt, ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
          x, x_ind, fA->getSparseNodeDownPtr(A, Ap)
        );
        yp++;
        if (yp >= ys) break;
        yj = fy->getSparseNodeIndex(y_ind, yp);
        Ap++;
        if (Ap >= As) break;
        Aj = fA->getSparseNodeIndex(A, Ap);
      } // loop
      return;
    } 
    // A is full, y_ind is sparse; intersect them
    int As = fA->getFullNodeSize(A);
    for (int yp=0; yp<ys; yp++) {
      int i = fy->getSparseNodeIndex(y_ind, yp);
      if (i >= As) break;
      int Ad = fA->getFullNodeDownPtr(A, i);
      if (0==Ad) continue;
      fy->getSparseNodeEdgeValue(y_ind, yp, yv);
      vectorMatrixMult_evplus_mt(
        pt, ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), x, x_ind, Ad
      );
    } // for yp
    return;
  } // y_ind is sparse

  //
  // Full y_ind
  //
  assert(fy->isFullNode(y_ind));
  ys = fy->getFullNodeSize(y_ind);

  if (fA->isSparseNode(A)) {
    // A is sparse, y_ind is full; intersect them
    int As = fA->getSparseNodeSize(A);
    for (int Ap=0; Ap<As; Ap++) {
      int i = fA->getSparseNodeIndex(A, Ap);
      if (i >= ys) break;
      int yd = fy->getFullNodeDownPtr(y_ind, i);
      if (0==yd) continue;
      fy->getFullNodeEdgeValue(y_ind, i, yv);
      vectorMatrixMult_evplus_mt(
        pt, ht-1, y+yv, yd, x, x_ind, fA->getSparseNodeDownPtr(A, Ap)
      );
    } // for Ap
    return;
  }

  // A is full, y_ind is full; intersect them
  int As = fA->getFullNodeSize(A);
  int sz = MIN(ys, As);
  for (int i=0; i<sz; i++) {
    int yd = fy->getFullNodeDownPtr(y_ind, i);
    if (0==yd) continue;
    int Ad = fA->getFullNodeDownPtr(A, i);
    if (0==Ad) continue;
    fy->getFullNodeEdgeValue(y_ind, i, yv);
    vectorMatrixMult_evplus_mt(pt, ht-1, y+yv, yd, x, x_ind, Ad);
  } // for i
}

void 
MEDDLY::vectorMatrixMult_evplus_mt(
  const op_param* pt, int ht, double* y, int y_ind, const double* x, int x_ind, int A
)
{
  // Handles the unprimed levels of A
  const expert_forest* fA = pt[4].readForest();
  assert(fA);

  if (0==ht) {
    y[0] += x[0] * fA->getReal(A);
    return;
  }

  const expert_forest* fx = pt[3].readForest();
  assert(fx);

  // It should be impossible for an indexing function to skip levels, right?   
  assert(fx->getNodeHeight(x_ind) == ht);
  int xv, xs;

  //
  // A is identity matrix times a constant; exploit that if we can
  //
  if (0==fA->getNodeHeight(A) && (x_ind == y_ind)) {
    const expert_forest* fy = pt[1].readForest();
    assert(fy);
    if (fx == fy) {
      assert(fx->isIndexSet(x_ind));
      // yes we can
      float v = fA->getReal(A);
      for (long i = fx->getIndexSetCardinality(x_ind)-1; i>=0; i--) {
        y[i] += x[i] * v;
      }
      return;
    }
  }

  //
  // Identity node in A
  //
  if (fA->getNodeHeight(A) < ht) {
    const expert_forest* fy = pt[1].readForest();
    assert(fy);
    int yv, ys;

    if (fx->isSparseNode(x_ind)) {
      xs = fx->getSparseNodeSize(x_ind);

      if (fy->isSparseNode(y_ind)) {
        // x_ind and y_ind are sparse; intersect them
        int ys = fy->getSparseNodeSize(y_ind);
        int xp = 0;
        int xi = fx->getSparseNodeIndex(x_ind, xp);
        int yp = 0;
        int yi = fy->getSparseNodeIndex(y_ind, yp);
        for (;;) {
          if (yi < xi) {
            yp++;
            if (yp >= ys) break;
            yi = fy->getSparseNodeIndex(y_ind, yp);
            continue;
          }
          if (xi < yi) {
            xp++;
            if (xp >= xs) break;
            xi = fx->getSparseNodeIndex(x_ind, xp);
            continue;
          }
          fx->getSparseNodeEdgeValue(x_ind, xp, xv);
          fy->getSparseNodeEdgeValue(y_ind, yp, yv);
          vectorMatrixMult_evplus_mt(
            pt, ht-1,
            y+yv, fy->getSparseNodeDownPtr(y_ind, yp),
            x+xv, fx->getSparseNodeDownPtr(x_ind, xp), 
            A
          );
          xp++;
          if (xp >= xs) break;
          xi = fx->getSparseNodeIndex(x_ind, xp);
          yp++;
          if (yp >= ys) break;
          yi = fy->getSparseNodeIndex(y_ind, yp);
        } // loop
        return;
      } // if y_ind is sparse

      assert(fy->isFullNode(y_ind));
      // x_ind is sparse, y_ind is full; intersect them
      ys = fy->getFullNodeSize(y_ind);
      for (int xp=0; xp<xs; xp++) {
        int i = fx->getSparseNodeIndex(x_ind, xp);
        if (i >= ys) break;
        int yd = fy->getFullNodeDownPtr(y_ind, i);
        if (0==yd) continue;
        fx->getSparseNodeEdgeValue(x_ind, xp, xv);
        fy->getFullNodeEdgeValue(y_ind, i, yv);
        vectorMatrixMult_evplus_mt(
          pt, ht-1, y+yv, yd, x+xv, fx->getSparseNodeDownPtr(x_ind, xp), A
        );
      } // for xp
      return;
    } // if x_ind is sparse

    assert(fx->isFullNode(x_ind));
    xs = fx->getFullNodeSize(x_ind);

    if (fy->isSparseNode(y_ind)) {
      // x_ind is full, y_ind is sparse; intersect them
      ys = fy->getSparseNodeSize(y_ind);
      for (int yp=0; yp<ys; yp++) {
        int i = fy->getSparseNodeIndex(y_ind, yp);
        if (i >= xs) break;
        int xd = fx->getFullNodeDownPtr(x_ind, i);
        if (0==xd) continue;
        fx->getFullNodeEdgeValue(x_ind, i, xv);
        fy->getSparseNodeEdgeValue(y_ind, yp, yv);
        vectorMatrixMult_evplus_mt(
          pt, ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), x+xv, xd, A
        );
      } // for xp
      return;
    } // if y_ind is sparse

    assert(fy->isFullNode(y_ind));
    // x_ind is full, y_ind is full; intersect them
    ys = fy->getFullNodeSize(y_ind);
    int sz = MIN(xs, ys);
    for (int i=0; i<sz; i++) {
      int xd = fx->getFullNodeDownPtr(x_ind, i);
      if (0==xd) continue;
      int yd = fy->getFullNodeDownPtr(y_ind, i);
      if (0==xd) continue;
      fx->getFullNodeEdgeValue(x_ind, i, xv);
      fy->getFullNodeEdgeValue(y_ind, i, yv);
      vectorMatrixMult_evplus_mt(pt, ht-1, y+yv, yd, x+xv, xd, A);
    } // for i
    return;

  } // if A is identity

  //
  // Sparse x_ind
  //
  if (fx->isSparseNode(x_ind)) {
    xs = fx->getSparseNodeSize(x_ind);

    if (fA->isSparseNode(A)) {
      // A is sparse, x_ind is sparse; intersect them
      int As = fA->getSparseNodeSize(A);
      int xp = 0;
      int xi = fx->getSparseNodeIndex(x_ind, xp);
      int Ap = 0;
      int Ai = fA->getSparseNodeIndex(A, Ap);
      for (;;) {
        if (Ai < xi) {
          Ap++;
          if (Ap >= As) break;
          Ai = fA->getSparseNodeIndex(A, Ap);
          continue;
        }
        if (xi < Ai) {
          xp++;
          if (xp >= xs) break;
          xi = fx->getSparseNodeIndex(x_ind, xp);
          continue;
        }
        fx->getSparseNodeEdgeValue(x_ind, xp, xv);
        vmprimed_evplus_mt(
          pt, ht, y, y_ind,
          x+xv, fx->getSparseNodeDownPtr(x_ind, xp), 
          fA->getSparseNodeDownPtr(A, Ap)
        );
        xp++;
        if (xp >= xs) break;
        xi = fx->getSparseNodeIndex(x_ind, xp);
        Ap++;
        if (Ap >= As) break;
        Ai = fA->getSparseNodeIndex(A, Ap);
      } // loop
      return;
    } 
    // A is full, x_ind is sparse; intersect them
    int As = fA->getFullNodeSize(A);
    for (int xp=0; xp<xs; xp++) {
      int i = fx->getSparseNodeIndex(x_ind, xp);
      if (i >= As) break;
      int Ad = fA->getFullNodeDownPtr(A, i);
      if (0==Ad) continue;
      fx->getSparseNodeEdgeValue(x_ind, xp, xv);
      vmprimed_evplus_mt(
        pt, ht, y, y_ind, x+xv, fx->getSparseNodeDownPtr(x_ind, xp), Ad
      );
    } // for xp
    return;
  } // x_ind is sparse

  //
  // Full x_ind
  //
  assert(fx->isFullNode(x_ind));
  xs = fx->getFullNodeSize(x_ind);

  if (fA->isSparseNode(A)) {
    // A is sparse, x_ind is full; intersect them
    int As = fA->getSparseNodeSize(A);
    for (int Ap=0; Ap<As; Ap++) {
      int i = fA->getSparseNodeIndex(A, Ap);
      if (i >= xs) break;
      int xd = fx->getFullNodeDownPtr(x_ind, i);
      if (0==xd) continue;
      fx->getFullNodeEdgeValue(x_ind, i, xv);
      vmprimed_evplus_mt(
        pt, ht, y, y_ind, x+xv, xd, fA->getSparseNodeDownPtr(A, Ap)
      );
    } // for Ap
    return;
  }

  // A is full, x_ind is full; intersect them
  int As = fA->getFullNodeSize(A);
  int sz = MIN(xs, As);
  for (int i=0; i<sz; i++) {
    int xd = fx->getFullNodeDownPtr(x_ind, i);
    if (0==xd) continue;
    int Ad = fA->getFullNodeDownPtr(A, i);
    if (0==Ad) continue;
    fx->getFullNodeEdgeValue(x_ind, i, xv);
    vmprimed_evplus_mt(pt, ht, y, y_ind, x+xv, xd, Ad);
  } // for i

}


// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void 
MEDDLY::vectorMatrixMult_evplus_evtimes(
  const op_param* pt, int ht, double* y, int y_ind, const double* x, int x_ind, int A
)
{
  // TBD
  throw error(error::NOT_IMPLEMENTED);
}


// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void
mvprimed_evplus_mt(
  const op_param* pt, int ht, double* y, int y_ind, int A, const double* x, int x_ind
)
{
  // Handles the primed levels of A
  const expert_forest* fA = pt[2].readForest();
  assert(fA);
  if (0==ht) {
    y[0] += x[0] * fA->getReal(A);
    return;
  }
  const expert_forest* fx = pt[4].readForest();
  assert(fx);

  // It should be impossible for an indexing function to skip levels, right?   
  assert(fx->getNodeHeight(x_ind) == ht);
  int xv, xs;

  //
  // Sparse x_ind
  //
  if (fx->isSparseNode(x_ind)) {
    xs = fx->getSparseNodeSize(x_ind);

    if (fA->isSparseNode(A)) {
      // A is sparse, x_ind is sparse; intersect them
      int As = fA->getSparseNodeSize(A);
      int xp = 0;
      int xj = fx->getSparseNodeIndex(x_ind, xp);
      int Ap = 0;
      int Aj = fA->getSparseNodeIndex(A, Ap);
      for (;;) {
        if (Aj < xj) {
          Ap++;
          if (Ap >= As) break;
          Aj = fA->getSparseNodeIndex(A, Ap);
          continue;
        }
        if (xj < Aj) {
          xp++;
          if (xp >= xs) break;
          xj = fx->getSparseNodeIndex(x_ind, xp);
          continue;
        }
        fx->getSparseNodeEdgeValue(x_ind, xp, xv);
        matrixVectorMult_evplus_mt(
          pt, ht-1, y, y_ind,
          fA->getSparseNodeDownPtr(A, Ap),
          x+xv, fx->getSparseNodeDownPtr(x_ind, xp)
        );
        xp++;
        if (xp >= xs) break;
        xj = fx->getSparseNodeIndex(x_ind, xp);
        Ap++;
        if (Ap >= As) break;
        Aj = fA->getSparseNodeIndex(A, Ap);
      } // loop
      return;
    } 
    // A is full, x_ind is sparse; intersect them
    int As = fA->getFullNodeSize(A);
    for (int xp=0; xp<xs; xp++) {
      int i = fx->getSparseNodeIndex(x_ind, xp);
      if (i >= As) break;
      int Ad = fA->getFullNodeDownPtr(A, i);
      if (0==Ad) continue;
      fx->getSparseNodeEdgeValue(x_ind, xp, xv);
      matrixVectorMult_evplus_mt(
        pt, ht-1, y, y_ind, Ad,
        x+xv, fx->getSparseNodeDownPtr(x_ind, xp)
      );
    } // for xp
    return;
  } // x_ind is sparse

  //
  // Full x_ind
  //
  assert(fx->isFullNode(x_ind));
  xs = fx->getFullNodeSize(x_ind);

  if (fA->isSparseNode(A)) {
    // A is sparse, x_ind is full; intersect them
    int As = fA->getSparseNodeSize(A);
    for (int Ap=0; Ap<As; Ap++) {
      int i = fA->getSparseNodeIndex(A, Ap);
      if (i >= xs) break;
      int xd = fx->getFullNodeDownPtr(x_ind, i);
      if (0==xd) continue;
      fx->getFullNodeEdgeValue(x_ind, i, xv);
      matrixVectorMult_evplus_mt(
        pt, ht-1, y, y_ind, fA->getSparseNodeDownPtr(A, Ap), x+xv, xd
      );
    } // for Ap
    return;
  }

  // A is full, x_ind is full; intersect them
  int As = fA->getFullNodeSize(A);
  int sz = MIN(xs, As);
  for (int i=0; i<sz; i++) {
    int xd = fx->getFullNodeDownPtr(x_ind, i);
    if (0==xd) continue;
    int Ad = fA->getFullNodeDownPtr(A, i);
    if (0==Ad) continue;
    fx->getFullNodeEdgeValue(x_ind, i, xv);
    matrixVectorMult_evplus_mt(pt, ht-1, y, y_ind, Ad, x+xv, xd);
  } // for i

}



void 
MEDDLY::matrixVectorMult_evplus_mt(
  const op_param* pt, int ht, double* y, int y_ind, int A, const double* x, int X_ind
)
{
  // Handles the unprimed levels of A
  const expert_forest* fA = pt[2].readForest();
  assert(fA);
  if (0==ht) {
    y[0] += fA->getReal(A) * x[0];
    return;
  }
  const expert_forest* fy = pt[1].readForest();
  assert(fy);

  // It should be impossible for an indexing function to skip levels, right?   
  assert(fy->getNodeHeight(y_ind) == ht);
  int yv, ys;

  //
  // A is identity matrix times a constant; exploit that if we can
  //
  if (0==fA->getNodeHeight(A) && (X_ind == y_ind)) {
    const expert_forest* fx = pt[4].readForest();
    assert(fx);
    if (fx == fy) {
      assert(fy->isIndexSet(y_ind));
      // yes we can
      float v = fA->getReal(A);
      for (long i = fy->getIndexSetCardinality(y_ind)-1; i>=0; i--) {
        y[i] += x[i] * v;
      }
      return;
    }
  }

  //
  // Identity node in A
  //
  if (fA->getNodeHeight(A) < ht) {
    const expert_forest* fX = pt[4].readForest();
    assert(fX);
    int Xv, Xs;

    if (fy->isSparseNode(y_ind)) {
      ys = fy->getSparseNodeSize(y_ind);

      if (fX->isSparseNode(X_ind)) {
        // X_ind and y_ind are sparse; intersect them
        int Xs = fX->getSparseNodeSize(X_ind);
        int yp = 0;
        int yi = fy->getSparseNodeIndex(y_ind, yp);
        int Xp = 0;
        int Xi = fX->getSparseNodeIndex(X_ind, Xp);
        for (;;) {
          if (Xi < yi) {
            Xp++;
            if (Xp >= Xs) break;
            Xi = fX->getSparseNodeIndex(X_ind, Xp);
            continue;
          }
          if (yi < Xi) {
            yp++;
            if (yp >= ys) break;
            yi = fy->getSparseNodeIndex(y_ind, yp);
            continue;
          }
          fy->getSparseNodeEdgeValue(y_ind, yp, yv);
          fX->getSparseNodeEdgeValue(X_ind, Xp, Xv);
          matrixVectorMult_evplus_mt(
            pt, ht-1,
            y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
            A,
            x+Xv, fX->getSparseNodeDownPtr(X_ind, Xp)
          );
          yp++;
          if (yp >= ys) break;
          yi = fy->getSparseNodeIndex(y_ind, yp);
          Xp++;
          if (Xp >= Xs) break;
          Xi = fX->getSparseNodeIndex(X_ind, Xp);
        } // loop
        return;
      } // if X_ind is sparse

      assert(fX->isFullNode(X_ind));
      // y_ind is sparse, X_ind is full; intersect them
      Xs = fX->getFullNodeSize(X_ind);
      for (int yp=0; yp<ys; yp++) {
        int i = fy->getSparseNodeIndex(y_ind, yp);
        if (i >= Xs) break;
        int Xd = fX->getFullNodeDownPtr(X_ind, i);
        if (0==Xd) continue;
        fy->getSparseNodeEdgeValue(y_ind, yp, yv);
        fX->getFullNodeEdgeValue(X_ind, i, Xv);
        matrixVectorMult_evplus_mt(
          pt, ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), A, x+Xv, Xd
        );
      } // for yp
      return;
    } // if y_ind is sparse

    assert(fy->isFullNode(y_ind));
    ys = fy->getFullNodeSize(y_ind);

    if (fX->isSparseNode(X_ind)) {
      // y_ind is full, X_ind is sparse; intersect them
      Xs = fX->getSparseNodeSize(X_ind);
      for (int Xp=0; Xp<Xs; Xp++) {
        int i = fX->getSparseNodeIndex(X_ind, Xp);
        if (i >= ys) break;
        int yd = fy->getFullNodeDownPtr(y_ind, i);
        if (0==yd) continue;
        fy->getFullNodeEdgeValue(y_ind, i, yv);
        fX->getSparseNodeEdgeValue(X_ind, Xp, Xv);
        matrixVectorMult_evplus_mt(
          pt, ht-1, y+yv, yd, A, x+Xv, fX->getSparseNodeDownPtr(X_ind, Xp)
        );
      } // for yp
      return;
    } // if X_ind is sparse

    assert(fX->isFullNode(X_ind));
    // y_ind is full, X_ind is full; intersect them
    Xs = fX->getFullNodeSize(X_ind);
    int sz = MIN(ys, Xs);
    for (int i=0; i<sz; i++) {
      int yd = fy->getFullNodeDownPtr(y_ind, i);
      if (0==yd) continue;
      int Xd = fX->getFullNodeDownPtr(X_ind, i);
      if (0==yd) continue;
      fy->getFullNodeEdgeValue(y_ind, i, yv);
      fX->getFullNodeEdgeValue(X_ind, i, Xv);
      matrixVectorMult_evplus_mt(pt, ht-1, y+yv, yd, A, x+Xv, Xd);
    } // for i
    return;

  } // if A is identity

  //
  // Sparse y_ind
  //
  if (fy->isSparseNode(y_ind)) {
    ys = fy->getSparseNodeSize(y_ind);

    if (fA->isSparseNode(A)) {
      // A is sparse, y_ind is sparse; intersect them
      int As = fA->getSparseNodeSize(A);
      int yp = 0;
      int yi = fy->getSparseNodeIndex(y_ind, yp);
      int Ap = 0;
      int Ai = fA->getSparseNodeIndex(A, Ap);
      for (;;) {
        if (Ai < yi) {
          Ap++;
          if (Ap >= As) break;
          Ai = fA->getSparseNodeIndex(A, Ap);
          continue;
        }
        if (yi < Ai) {
          yp++;
          if (yp >= ys) break;
          yi = fy->getSparseNodeIndex(y_ind, yp);
          continue;
        }
        fy->getSparseNodeEdgeValue(y_ind, yp, yv);
        mvprimed_evplus_mt(
          pt, ht, 
          y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
          fA->getSparseNodeDownPtr(A, Ap),
          x, X_ind
        );
        yp++;
        if (yp >= ys) break;
        yi = fy->getSparseNodeIndex(y_ind, yp);
        Ap++;
        if (Ap >= As) break;
        Ai = fA->getSparseNodeIndex(A, Ap);
      } // loop
      return;
    } 
    // A is full, y_ind is sparse; intersect them
    int As = fA->getFullNodeSize(A);
    for (int yp=0; yp<ys; yp++) {
      int i = fy->getSparseNodeIndex(y_ind, yp);
      if (i >= As) break;
      int Ad = fA->getFullNodeDownPtr(A, i);
      if (0==Ad) continue;
      fy->getSparseNodeEdgeValue(y_ind, yp, yv);
      mvprimed_evplus_mt(
        pt, ht, 
        y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
        Ad, x, X_ind
      );
    } // for yp
    return;
  } // y_ind is sparse

  //
  // Full y_ind
  //
  assert(fy->isFullNode(y_ind));
  ys = fy->getFullNodeSize(y_ind);

  if (fA->isSparseNode(A)) {
    // A is sparse, y_ind is full; intersect them
    int As = fA->getSparseNodeSize(A);
    for (int Ap=0; Ap<As; Ap++) {
      int i = fA->getSparseNodeIndex(A, Ap);
      if (i >= ys) break;
      int yd = fy->getFullNodeDownPtr(y_ind, i);
      if (0==yd) continue;
      fy->getFullNodeEdgeValue(y_ind, i, yv);
      mvprimed_evplus_mt(
        pt, ht, y+yv, yd, fA->getSparseNodeDownPtr(A, Ap), x, X_ind
      );
    } // for Ap
    return;
  }

  // A is full, y_ind is full; intersect them
  int As = fA->getFullNodeSize(A);
  int sz = MIN(ys, As);
  for (int i=0; i<sz; i++) {
    int yd = fy->getFullNodeDownPtr(y_ind, i);
    if (0==yd) continue;
    int Ad = fA->getFullNodeDownPtr(A, i);
    if (0==Ad) continue;
    fy->getFullNodeEdgeValue(y_ind, i, yv);
    mvprimed_evplus_mt(pt, ht, y+yv, yd, Ad, x, X_ind);
  } // for i

}


// ----------------------------------------------------------------------
// ----------------------------------------------------------------------


void 
MEDDLY::matrixVectorMult_evplus_evtimes(
  const op_param* pt, int ht, double* y, int y_ind, int A, const double* x, int x_ind
)
{
  // TBD
  throw error(error::NOT_IMPLEMENTED);
}

#endif

