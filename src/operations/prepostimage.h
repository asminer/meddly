
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



#ifndef PREPOSTIMAGE_H
#define PREPOSTIMAGE_H

#include "operation_ext.h"

/* Pre-image and Post-image operations */

// -------------------------- MDD MXD Image operations ----------------------


class mdd_mxd_image_operation : public mdd_apply_operation {
  public:
    mdd_mxd_image_operation();
    virtual ~mdd_mxd_image_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual const char* getName() const { return "Mdd-Mxd Image Operation"; }
    virtual bool isCommutative() const { return false; }

    virtual int compute(op_info* owner, int a, int b) = 0;

#if 0
  protected:
    virtual bool findResult(op_info* owner, int a, int b, int& c);
    virtual void saveResult(op_info* owner, int a, int b, int c);
#endif

  private:
    // Disabling this function
    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mdd_post_image : public mdd_mxd_image_operation {
  public:
    static mdd_post_image* getInstance();
    virtual int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Post Image"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_post_image();
    mdd_post_image(const mdd_post_image& copy);
    mdd_post_image& operator=(const mdd_post_image& copy);
    virtual ~mdd_post_image();

    virtual int compute(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int fullFull(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int fullSparse(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int sparseFull(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int sparseSparse(op_info* owner, op_info* unionOp,
        int mdd, int mxd);
    virtual int expandMdd(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int expandMxd(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual void expandMxdHelper (op_info* owner, op_info* unionOp,
        int mdd, int iMxd, int result);
};


class mdd_pre_image : public mdd_post_image {
  public:
    static mdd_pre_image* getInstance();
    virtual const char* getName() const { return "Mdd-Mxd Pre Image"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_pre_image();
    mdd_pre_image(const mdd_pre_image& copy);
    mdd_pre_image& operator=(const mdd_pre_image& copy);
    virtual ~mdd_pre_image();

    virtual int compute(op_info* owner, op_info* unionOp, int mdd, int mxd);

    virtual int expandMxd (op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int expandMxdByOneLevel(op_info* owner, op_info* unionOp,
        int mdd, int iMxd); 

    virtual int expand(op_info* owner, op_info* unionOp, int mdd, int mxd);
    virtual int expandByOneLevel(op_info* owner, op_info* unionOp,
        int mdd, int iMxd);

    // No need for mdd_pre_image::expandMdd()
    // since we are reusing mdd_post_image::expandMdd()
    // This works correctly because when a mxd level is skipped
    // due to an identity reduced node:
    // result[j] += pre_image(mdd[i], mxd[j][i])
    // is only valid when i == j (those are the only non-zero entries for
    // an identity reduced node).
    // Therefore, the above becomes
    // result[i] += pre_image(mdd[i], mxd[i][i])
    // which is the same expansion used in mdd_post_image::expandMdd().
};


// And these?

class mtmdd_post_image : public mdd_post_image {
  public:
    static mtmdd_post_image* getInstance();
    virtual const char* getName() const { return "MtMdd Post-Image"; }
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mtmdd_post_image();
    mtmdd_post_image(const mtmdd_post_image& copy);
    mtmdd_post_image& operator=(const mtmdd_post_image& copy);
    virtual ~mtmdd_post_image();

    virtual int compute(op_info* owner, op_info* plusOp, int mdd, int mxd);
};


class mtmdd_pre_image : public mdd_pre_image {
  public:
    static mtmdd_pre_image* getInstance();
    virtual const char* getName() const { return "MtMdd Pre-Image"; }
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mtmdd_pre_image();
    mtmdd_pre_image(const mtmdd_pre_image& copy);
    mtmdd_pre_image& operator=(const mtmdd_pre_image& copy);
    virtual ~mtmdd_pre_image();

    virtual int compute(op_info* owner, op_info* plusOp, int mdd, int mxd);
};


#endif

