
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

namespace MEDDLY {
  class base_evplus_mt;

  class VM_evplus_mt;

  class MV_evplus_mt;

  class VM_opname;
  class MV_opname;
};

// ******************************************************************
// *                                                                *
// *                      base_evplus_mt class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::base_evplus_mt : public numerical_operation {
  public:
    base_evplus_mt(const numerical_opname* code, const dd_edge &x_ind,
      const dd_edge& A, const dd_edge &y_ind);

    virtual ~base_evplus_mt();

    virtual void compute(double* y, const double* x);

    virtual void compute(int ht, double* y, int y_ind, const double* x, 
      int x_ind, int A) = 0;

    virtual bool isStaleEntry(const int*) {
      throw error(error::MISCELLANEOUS);
    }
    virtual void discardEntry(const int*) {
      throw error(error::MISCELLANEOUS);
    }
    virtual void showEntry(FILE*, const int*) const {
      throw error(error::MISCELLANEOUS);
    }

  protected:
    const expert_forest* fx;
    const expert_forest* fA;
    const expert_forest* fy;
    int x_root;
    int A_root;
    int y_root;
    int L;
};

MEDDLY::base_evplus_mt::base_evplus_mt(const numerical_opname* code, 
  const dd_edge &x_ind, const dd_edge& A, const dd_edge &y_ind)
 : numerical_operation(code)
{
  fx = (const expert_forest*) x_ind.getForest();
  fA = (const expert_forest*) A.getForest();
  fy = (const expert_forest*) y_ind.getForest();
  MEDDLY_DCASSERT(fx);
  MEDDLY_DCASSERT(fA);
  MEDDLY_DCASSERT(fy);
  x_root = x_ind.getNode();
  A_root = A.getNode();
  y_root = y_ind.getNode();
  L = fx->getDomain()->getNumVariables();
}

MEDDLY::base_evplus_mt::~base_evplus_mt()
{
}

void MEDDLY::base_evplus_mt::compute(double* y, const double* x)
{
  compute(L, y, y_root, x, x_root, A_root);
}

// ******************************************************************
// *                                                                *
// *                       VM_evplus_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::VM_evplus_mt : public base_evplus_mt {
  public:
    VM_evplus_mt(const numerical_opname* code, const dd_edge &x_ind,
      const dd_edge& A, const dd_edge &y_ind);

    virtual void compute(int ht, double* y, int y_ind, const double* x, 
      int x_ind, int A);

    void comp_pr(int ht, double* y, int y_ind, const double* x, 
      int x_ind, int A);

};

MEDDLY::VM_evplus_mt::VM_evplus_mt(const numerical_opname* code, 
  const dd_edge &x_ind, const dd_edge& A, const dd_edge &y_ind)
  : base_evplus_mt(code, x_ind, A, y_ind)
{
}

void MEDDLY::VM_evplus_mt::compute(int ht, double* y, int y_ind, 
  const double* x, int x_ind, int A)
{
  // Handles the unprimed levels of A
  if (0==ht) {
    y[0] += x[0] * fA->getReal(A);
    return;
  }

  // It should be impossible for an indexing function to skip levels, right?   
  MEDDLY_DCASSERT(fx->getNodeHeight(x_ind) == ht);
  int xv, xs;

  //
  // A is identity matrix times a constant; exploit that if we can
  //
  if (0==fA->getNodeHeight(A) && (x_ind == y_ind)) {
    if (fx == fy) {
      MEDDLY_DCASSERT(fx->isIndexSet(x_ind));
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
          compute(
            ht-1,
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

      MEDDLY_DCASSERT(fy->isFullNode(y_ind));
      // x_ind is sparse, y_ind is full; intersect them
      ys = fy->getFullNodeSize(y_ind);
      for (int xp=0; xp<xs; xp++) {
        int i = fx->getSparseNodeIndex(x_ind, xp);
        if (i >= ys) break;
        int yd = fy->getFullNodeDownPtr(y_ind, i);
        if (0==yd) continue;
        fx->getSparseNodeEdgeValue(x_ind, xp, xv);
        fy->getFullNodeEdgeValue(y_ind, i, yv);
        compute(ht-1, y+yv, yd, x+xv, fx->getSparseNodeDownPtr(x_ind, xp), A);
      } // for xp
      return;
    } // if x_ind is sparse

    MEDDLY_DCASSERT(fx->isFullNode(x_ind));
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
        compute(ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), x+xv, xd, A);
      } // for xp
      return;
    } // if y_ind is sparse

    MEDDLY_DCASSERT(fy->isFullNode(y_ind));
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
      compute(ht-1, y+yv, yd, x+xv, xd, A);
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
        comp_pr(
          ht, y, y_ind,
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
      comp_pr(
        ht, y, y_ind, x+xv, fx->getSparseNodeDownPtr(x_ind, xp), Ad
      );
    } // for xp
    return;
  } // x_ind is sparse

  //
  // Full x_ind
  //
  MEDDLY_DCASSERT(fx->isFullNode(x_ind));
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
      comp_pr(
        ht, y, y_ind, x+xv, xd, fA->getSparseNodeDownPtr(A, Ap)
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
    comp_pr(ht, y, y_ind, x+xv, xd, Ad);
  } // for i

}

void MEDDLY::VM_evplus_mt::comp_pr(int ht, double* y, int y_ind, 
  const double* x, int x_ind, int A)
{
  // Handles the primed levels of A
  if (0==ht) {
    y[0] += x[0] * fA->getReal(A);
    return;
  }

  // It should be impossible for an indexing function to skip levels, right?   
  MEDDLY_DCASSERT(fy->getNodeHeight(y_ind) == ht);
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
        compute(
          ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
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
      compute(
        ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), x, x_ind, Ad
      );
    } // for yp
    return;
  } // y_ind is sparse

  //
  // Full y_ind
  //
  MEDDLY_DCASSERT(fy->isFullNode(y_ind));
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
      compute(
        ht-1, y+yv, yd, x, x_ind, fA->getSparseNodeDownPtr(A, Ap)
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
    compute(ht-1, y+yv, yd, x, x_ind, Ad);
  } // for i
}




// ******************************************************************
// *                                                                *
// *                       MV_evplus_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::MV_evplus_mt : public base_evplus_mt {
  public:
    MV_evplus_mt(const numerical_opname* code, const dd_edge &x_ind,
      const dd_edge& A, const dd_edge &y_ind);

    virtual void compute(int ht, double* y, int y_ind, const double* x, 
      int x_ind, int A);

    void comp_pr(int ht, double* y, int y_ind, const double* x, 
      int x_ind, int A);

};

MEDDLY::MV_evplus_mt::MV_evplus_mt(const numerical_opname* code, 
  const dd_edge &x_ind, const dd_edge& A, const dd_edge &y_ind)
  : base_evplus_mt(code, x_ind, A, y_ind)
{
}

void MEDDLY::MV_evplus_mt::compute(int ht, double* y, int y_ind, 
  const double* x, int X_ind, int A)
{
  // Handles the unprimed levels of A
  if (0==ht) {
    y[0] += fA->getReal(A) * x[0];
    return;
  }

  // It should be impossible for an indexing function to skip levels, right?   
  MEDDLY_DCASSERT(fy->getNodeHeight(y_ind) == ht);
  int yv, ys;

  //
  // A is identity matrix times a constant; exploit that if we can
  //
  if (0==fA->getNodeHeight(A) && (X_ind == y_ind)) {
    if (fx == fy) {
      MEDDLY_DCASSERT(fy->isIndexSet(y_ind));
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
    int Xv, Xs;

    if (fy->isSparseNode(y_ind)) {
      ys = fy->getSparseNodeSize(y_ind);

      if (fx->isSparseNode(X_ind)) {
        // X_ind and y_ind are sparse; intersect them
        int Xs = fx->getSparseNodeSize(X_ind);
        int yp = 0;
        int yi = fy->getSparseNodeIndex(y_ind, yp);
        int Xp = 0;
        int Xi = fx->getSparseNodeIndex(X_ind, Xp);
        for (;;) {
          if (Xi < yi) {
            Xp++;
            if (Xp >= Xs) break;
            Xi = fx->getSparseNodeIndex(X_ind, Xp);
            continue;
          }
          if (yi < Xi) {
            yp++;
            if (yp >= ys) break;
            yi = fy->getSparseNodeIndex(y_ind, yp);
            continue;
          }
          fy->getSparseNodeEdgeValue(y_ind, yp, yv);
          fx->getSparseNodeEdgeValue(X_ind, Xp, Xv);
          compute(
            ht-1,
            y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
            x+Xv, fx->getSparseNodeDownPtr(X_ind, Xp),
            A
          );
          yp++;
          if (yp >= ys) break;
          yi = fy->getSparseNodeIndex(y_ind, yp);
          Xp++;
          if (Xp >= Xs) break;
          Xi = fx->getSparseNodeIndex(X_ind, Xp);
        } // loop
        return;
      } // if X_ind is sparse

      MEDDLY_DCASSERT(fx->isFullNode(X_ind));
      // y_ind is sparse, X_ind is full; intersect them
      Xs = fx->getFullNodeSize(X_ind);
      for (int yp=0; yp<ys; yp++) {
        int i = fy->getSparseNodeIndex(y_ind, yp);
        if (i >= Xs) break;
        int Xd = fx->getFullNodeDownPtr(X_ind, i);
        if (0==Xd) continue;
        fy->getSparseNodeEdgeValue(y_ind, yp, yv);
        fx->getFullNodeEdgeValue(X_ind, i, Xv);
        compute(ht-1, y+yv, fy->getSparseNodeDownPtr(y_ind, yp), x+Xv, Xd, A);
      } // for yp
      return;
    } // if y_ind is sparse

    MEDDLY_DCASSERT(fy->isFullNode(y_ind));
    ys = fy->getFullNodeSize(y_ind);

    if (fx->isSparseNode(X_ind)) {
      // y_ind is full, X_ind is sparse; intersect them
      Xs = fx->getSparseNodeSize(X_ind);
      for (int Xp=0; Xp<Xs; Xp++) {
        int i = fx->getSparseNodeIndex(X_ind, Xp);
        if (i >= ys) break;
        int yd = fy->getFullNodeDownPtr(y_ind, i);
        if (0==yd) continue;
        fy->getFullNodeEdgeValue(y_ind, i, yv);
        fx->getSparseNodeEdgeValue(X_ind, Xp, Xv);
        compute(ht-1, y+yv, yd, x+Xv, fx->getSparseNodeDownPtr(X_ind, Xp), A);
      } // for yp
      return;
    } // if X_ind is sparse

    MEDDLY_DCASSERT(fx->isFullNode(X_ind));
    // y_ind is full, X_ind is full; intersect them
    Xs = fx->getFullNodeSize(X_ind);
    int sz = MIN(ys, Xs);
    for (int i=0; i<sz; i++) {
      int yd = fy->getFullNodeDownPtr(y_ind, i);
      if (0==yd) continue;
      int Xd = fx->getFullNodeDownPtr(X_ind, i);
      if (0==yd) continue;
      fy->getFullNodeEdgeValue(y_ind, i, yv);
      fx->getFullNodeEdgeValue(X_ind, i, Xv);
      compute(ht-1, y+yv, yd, x+Xv, Xd, A);
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
        comp_pr(
          ht, 
          y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
          x, X_ind,
          fA->getSparseNodeDownPtr(A, Ap)
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
      comp_pr(
        ht, 
        y+yv, fy->getSparseNodeDownPtr(y_ind, yp), 
        x, X_ind, Ad
      );
    } // for yp
    return;
  } // y_ind is sparse

  //
  // Full y_ind
  //
  MEDDLY_DCASSERT(fy->isFullNode(y_ind));
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
      comp_pr(ht, y+yv, yd, x, X_ind, fA->getSparseNodeDownPtr(A, Ap));
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
    comp_pr(ht, y+yv, yd, x, X_ind, Ad);
  } // for i

}

void MEDDLY::MV_evplus_mt::comp_pr(int ht, double* y, int y_ind, 
  const double* x, int x_ind, int A)
{
  // Handles the primed levels of A
  if (0==ht) {
    y[0] += x[0] * fA->getReal(A);
    return;
  }

  // It should be impossible for an indexing function to skip levels, right?   
  MEDDLY_DCASSERT(fx->getNodeHeight(x_ind) == ht);
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
        compute(
          ht-1, y, y_ind,
          x+xv, fx->getSparseNodeDownPtr(x_ind, xp),
          fA->getSparseNodeDownPtr(A, Ap)
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
      compute(
        ht-1, y, y_ind, x+xv, fx->getSparseNodeDownPtr(x_ind, xp), Ad
      );
    } // for xp
    return;
  } // x_ind is sparse

  //
  // Full x_ind
  //
  MEDDLY_DCASSERT(fx->isFullNode(x_ind));
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
      compute(
        ht-1, y, y_ind, x+xv, xd, fA->getSparseNodeDownPtr(A, Ap)
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
    compute(ht-1, y, y_ind, x+xv, xd, Ad);
  } // for i

}



// ******************************************************************
// *                                                                *
// *                        VM_opname  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::VM_opname : public numerical_opname {
  public:
    VM_opname();
    virtual numerical_operation* buildOperation(const dd_edge &x_ind,
      const dd_edge& A, const dd_edge& y_ind) const;
};

MEDDLY::VM_opname::VM_opname() : numerical_opname("VectMatrMult")
{
}

MEDDLY::numerical_operation* 
MEDDLY::VM_opname::buildOperation(const dd_edge &x_ind, const dd_edge& A,
  const dd_edge& y_ind) const
{
  const expert_forest* fx = (const expert_forest*) x_ind.getForest();
  const expert_forest* fA = (const expert_forest*) A.getForest();
  const expert_forest* fy = (const expert_forest*) y_ind.getForest();

  // everyone must use the same domain
  if (      (fx->getDomain() != fy->getDomain()) 
        ||  (fx->getDomain() != fA->getDomain())  )
  {
    throw error(error::DOMAIN_MISMATCH);
  }

  // Check edge types
  if (
           (fy->getRangeType() != forest::INTEGER) 
        || (fy->isForRelations())
        || (fx->getRangeType() != forest::INTEGER)
        || (fx->isForRelations())
        || (fA->getRangeType() != forest::REAL)
        || (!fA->isForRelations())
      ) 
  {
    throw error(error::TYPE_MISMATCH);
  }

  // A can't be fully reduced.
  if (fA->isFullyReduced()) {
    throw error(error::TYPE_MISMATCH);
  }

  // For now, fy and fx must be EV+MDDs.
  if (     (fy->getEdgeLabeling() != forest::EVPLUS) 
        || (fx->getEdgeLabeling() != forest::EVPLUS) )
  {
    throw error(error::NOT_IMPLEMENTED);
  }

  switch (fA->getEdgeLabeling()) {
    case forest::MULTI_TERMINAL:
      return new VM_evplus_mt(this, x_ind, A, y_ind);

    case forest::EVTIMES:
      throw error(error::NOT_IMPLEMENTED);

    default:
      throw error(error::TYPE_MISMATCH);
  };
}

// ******************************************************************
// *                                                                *
// *                        MV_opname  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::MV_opname : public numerical_opname {
  public:
    MV_opname();
    virtual numerical_operation* buildOperation(const dd_edge &x_ind,
      const dd_edge& A, const dd_edge& y_ind) const;
};

MEDDLY::MV_opname::MV_opname() : numerical_opname("MatrVectMult")
{
}

MEDDLY::numerical_operation* 
MEDDLY::MV_opname::buildOperation(const dd_edge &x_ind, const dd_edge& A,
  const dd_edge& y_ind) const
{
  const expert_forest* fx = (const expert_forest*) x_ind.getForest();
  const expert_forest* fA = (const expert_forest*) A.getForest();
  const expert_forest* fy = (const expert_forest*) y_ind.getForest();

  // everyone must use the same domain
  if (      (fx->getDomain() != fy->getDomain()) 
        ||  (fx->getDomain() != fA->getDomain())  )
  {
    throw error(error::DOMAIN_MISMATCH);
  }

  // Check edge types
  if (
           (fy->getRangeType() != forest::INTEGER) 
        || (fy->isForRelations())
        || (fx->getRangeType() != forest::INTEGER)
        || (fx->isForRelations())
        || (fA->getRangeType() != forest::REAL)
        || (!fA->isForRelations())
      ) 
  {
    throw error(error::TYPE_MISMATCH);
  }

  // A can't be fully reduced.
  if (fA->isFullyReduced()) {
    throw error(error::TYPE_MISMATCH);
  }

  // For now, fy and fx must be EV+MDDs.
  if (     (fy->getEdgeLabeling() != forest::EVPLUS) 
        || (fx->getEdgeLabeling() != forest::EVPLUS) )
  {
    throw error(error::NOT_IMPLEMENTED);
  }

  switch (fA->getEdgeLabeling()) {
    case forest::MULTI_TERMINAL:
      return new MV_evplus_mt(this, x_ind, A, y_ind);

    case forest::EVTIMES:
      throw error(error::NOT_IMPLEMENTED);

    default:
      throw error(error::TYPE_MISMATCH);
  };
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::numerical_opname* MEDDLY::initVectorMatrixMult(const settings &s)
{
  return new VM_opname;
}

MEDDLY::numerical_opname* MEDDLY::initMatrixVectorMult(const settings &s)
{
  return new MV_opname;
}

