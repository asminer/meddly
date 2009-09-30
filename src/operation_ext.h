
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



#ifndef OPERATION_EXTENSIONS_H
#define OPERATION_EXTENSIONS_H

#include "../include/meddly_expert.h"

#include <vector>


/** MDD element-wise operations

    Abstract class for element-wise MDD operations.
*/
class mdd_apply_operation : public operation {
  public:
    mdd_apply_operation(const char *name, bool commutative);
    virtual ~mdd_apply_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Calls compute(op_info*, dd_edge, dd_edge, dd_edge)
    virtual compute_manager::error compute(op_info* owner, dd_edge** operands);

    // Returns an error
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        dd_edge& b);

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual int compute(op_info* owner, int a, int b);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;

    virtual bool findResult(op_info* owner, int a, int b, int& c);
    virtual void saveResult(op_info* owner, int a, int b, int c);

    virtual void expandA(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void expandB(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void fullFull(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void fullSparse(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void sparseFull(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void sparseSparse(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);

  private:
    bool commutative;
};


class mdd_union : public mdd_apply_operation {
  public:
    static mdd_union* getInstance();
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mdd_union(const char* name);
    mdd_union(const mdd_union& copy);
    mdd_union& operator=(const mdd_union& copy);
    ~mdd_union();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);

#if 1
    virtual void expandA(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    // use expandA(b, a) instead of expandB(a, b)
    virtual void fullFull(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void sparseSparse(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void fullSparse(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    // use fullSparse(b, a) instead of sparseFull(a, b)
#endif
};


class mdd_intersection : public mdd_apply_operation {
  public:
    static mdd_intersection* getInstance();
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mdd_intersection(const char* name);
    mdd_intersection(const mdd_intersection& copy);
    mdd_intersection& operator=(const mdd_intersection& copy);
    ~mdd_intersection();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
#if 1
    virtual void expandA(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void expandB(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void fullFull(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void fullSparse(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void sparseFull(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
    virtual void sparseSparse(op_info* owner, expert_forest* expertForest,
        int a, int b, int result, int resultSize);
#endif
};


class mdd_difference : public mdd_apply_operation {
  public:
    static mdd_difference* getInstance();

  protected:
    mdd_difference(const char* name);
    mdd_difference(const mdd_difference& copy);
    mdd_difference& operator=(const mdd_difference& copy);
    ~mdd_difference();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// ------------------ MXD operations -------------------------------

// HERE: Test!!

class mxd_apply_operation : public mdd_apply_operation {
  public:
    mxd_apply_operation(const char *name, bool commutative);
    virtual ~mxd_apply_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;

    virtual void expandA(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
    virtual void singleExpandA(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
    virtual void expandB(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
    virtual void singleExpandB(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
};


class mxd_alt_apply_operation : public operation {
  public:
    mxd_alt_apply_operation(const char *name, bool commutative);
    virtual ~mxd_alt_apply_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Calls compute(op_info*, dd_edge, dd_edge, dd_edge)
    virtual compute_manager::error compute(op_info* owner, dd_edge** operands);

    // Returns an error
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        dd_edge& b);

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual int compute(op_info* owner, int level, int a, int b);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;

    virtual bool findResult(op_info* owner, int k, int a, int b, int& c);
    virtual void saveResult(op_info* owner, int k, int a, int b, int c);

    virtual void expandSkippedLevel(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
    virtual void expandA(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
    virtual void singleExpandA(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
    virtual void expandB(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);
    virtual void singleExpandB(op_info* owner, int a, int b,
        expert_forest* expertForest,
        int result, int resultLevel, int resultSize);

  private:
    bool commutative;
};


class mxd_union : public mxd_apply_operation {
  public:
    static mxd_union* getInstance();

  protected:
    mxd_union(const char* name);
    mxd_union(const mxd_union& copy);
    mxd_union& operator=(const mxd_union& copy);
    ~mxd_union();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mxd_intersection : public mxd_apply_operation {
  public:
    static mxd_intersection* getInstance();

  protected:
    mxd_intersection(const char* name);
    mxd_intersection(const mxd_intersection& copy);
    mxd_intersection& operator=(const mxd_intersection& copy);
    ~mxd_intersection();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mxd_difference : public mxd_apply_operation {
  public:
    static mxd_difference* getInstance();

  protected:
    mxd_difference(const char* name);
    mxd_difference(const mxd_difference& copy);
    mxd_difference& operator=(const mxd_difference& copy);
    ~mxd_difference();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// ------------------------ MTMXD Apply operations ---------------------------

/** MTXDD element-wise operations

    Abstract class for element-wise MTMXD operations.
*/
// TODO: TEST
class mtmxd_apply_operation : public mxd_apply_operation {
  public:
    mtmxd_apply_operation(const char *name, bool commutative);
    virtual ~mtmxd_apply_operation();
    virtual compute_manager::error typeCheck(const op_info* owner);
  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;
};


class mtmxd_plus : public mtmxd_apply_operation {
  public:
    static mtmxd_plus* getInstance();

  protected:
    mtmxd_plus(const char* name);
    mtmxd_plus(const mtmxd_plus& copy);
    mtmxd_plus& operator=(const mtmxd_plus& copy);
    ~mtmxd_plus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_minus : public mtmxd_apply_operation {
  public:
    static mtmxd_minus* getInstance();

  protected:
    mtmxd_minus(const char* name);
    mtmxd_minus(const mtmxd_minus& copy);
    mtmxd_minus& operator=(const mtmxd_minus& copy);
    ~mtmxd_minus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_multiply : public mtmxd_apply_operation {
  public:
    static mtmxd_multiply* getInstance();

  protected:
    mtmxd_multiply(const char* name);
    mtmxd_multiply(const mtmxd_multiply& copy);
    mtmxd_multiply& operator=(const mtmxd_multiply& copy);
    ~mtmxd_multiply();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_min : public mtmxd_apply_operation {
  public:
    static mtmxd_min* getInstance();

  protected:
    mtmxd_min(const char* name);
    mtmxd_min(const mtmxd_min& copy);
    mtmxd_min& operator=(const mtmxd_min& copy);
    ~mtmxd_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_max : public mtmxd_apply_operation {
  public:
    static mtmxd_max* getInstance();

  protected:
    mtmxd_max(const char* name);
    mtmxd_max(const mtmxd_max& copy);
    mtmxd_max& operator=(const mtmxd_max& copy);
    ~mtmxd_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_not_equal : public mtmxd_apply_operation {
  public:
    static mtmxd_not_equal* getInstance();

  protected:
    mtmxd_not_equal(const char* name);
    mtmxd_not_equal(const mtmxd_not_equal& copy);
    mtmxd_not_equal& operator=(const mtmxd_not_equal& copy);
    ~mtmxd_not_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_less_than : public mtmxd_apply_operation {
  public:
    static mtmxd_less_than* getInstance();

  protected:
    mtmxd_less_than(const char* name);
    mtmxd_less_than(const mtmxd_less_than& copy);
    mtmxd_less_than& operator=(const mtmxd_less_than& copy);
    ~mtmxd_less_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_greater_than : public mtmxd_apply_operation {
  public:
    static mtmxd_greater_than* getInstance();

  protected:
    mtmxd_greater_than(const char* name);
    mtmxd_greater_than(const mtmxd_greater_than& copy);
    mtmxd_greater_than& operator=(const mtmxd_greater_than& copy);
    ~mtmxd_greater_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// same as mtmxd_apply_operation, except this is for mtmxd operations
// with op(0,0)!=0 (such as /, <=, >=, ==)
class mtmxd_alt_apply_operation : public mxd_alt_apply_operation {
  public:
    mtmxd_alt_apply_operation(const char *name, bool commutative);
    virtual ~mtmxd_alt_apply_operation();
    virtual compute_manager::error typeCheck(const op_info* owner);
  protected:
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;
};


class mtmxd_divide : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_divide* getInstance();

  protected:
    mtmxd_divide(const char* name);
    mtmxd_divide(const mtmxd_divide& copy);
    mtmxd_divide& operator=(const mtmxd_divide& copy);
    ~mtmxd_divide();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_less_than_equal : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_less_than_equal* getInstance();

  protected:
    mtmxd_less_than_equal(const char* name);
    mtmxd_less_than_equal(const mtmxd_less_than_equal& copy);
    mtmxd_less_than_equal& operator=(const mtmxd_less_than_equal& copy);
    ~mtmxd_less_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_greater_than_equal : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_greater_than_equal* getInstance();

  protected:
    mtmxd_greater_than_equal(const char* name);
    mtmxd_greater_than_equal(const mtmxd_greater_than_equal& copy);
    mtmxd_greater_than_equal& operator=(const mtmxd_greater_than_equal& copy);
    ~mtmxd_greater_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_equal : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_equal* getInstance();

  protected:
    mtmxd_equal(const char* name);
    mtmxd_equal(const mtmxd_equal& copy);
    mtmxd_equal& operator=(const mtmxd_equal& copy);
    ~mtmxd_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// -------------------------- MDD MXD Image operations ----------------------


class mdd_mxd_image_operation : public mdd_apply_operation {
  public:
    mdd_mxd_image_operation(const char *name);
    virtual ~mdd_mxd_image_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);

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

  protected:
    mdd_post_image(const char* name);
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

  protected:
    mdd_pre_image(const char* name);
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


// Traditional reachability analysis
class mdd_reachability_bfs : public mdd_mxd_image_operation {
  public:
    static mdd_reachability_bfs* getInstance();
    int compute(op_info* owner, int a, int b);

  protected:
    mdd_reachability_bfs(const char* name);
    mdd_reachability_bfs(const mdd_reachability_bfs& copy);
    mdd_reachability_bfs& operator=(const mdd_reachability_bfs& copy);
    virtual ~mdd_reachability_bfs();
};


// Reachability via "saturation" algorithm
class mdd_reachability_dfs : public mdd_mxd_image_operation {
  public:
    static mdd_reachability_dfs* getInstance();
    int compute(op_info* owner, int a, int b);

  protected:
    mdd_reachability_dfs(const char* name);
    mdd_reachability_dfs(const mdd_reachability_dfs& copy);
    mdd_reachability_dfs& operator=(const mdd_reachability_dfs& copy);
    virtual ~mdd_reachability_dfs();

    void initialize(op_info* owner);
    void clear();

    void splitMxd(int mxd);
    int saturate(int mdd);
    void saturateHelper(int& mdd);
    int recFire(int mdd, int mxd);

  private:
    op_info*       owner;         // pointer to dfs reachability operation
    expert_forest* ddf;           // MDD forest
    expert_forest* xdf;           // MXD forest
    expert_domain* ed;            // domain
    expert_compute_manager* ecm;  // compute manager

    // Next-state function is split and stored here (see Saturation algorithm).
    std::vector<int> splits;

    // scratch.size () == number of variable handles in the domain
    // scratch[level_handle].size() == level_bound(level_handle)
    std::vector< std::vector<int> > scratch;
    std::vector< std::vector<bool> > curr, next;

    op_info*          mddUnionOp;
    op_info*          mxdIntersectionOp;
    op_info*          mxdDifferenceOp;

    mdd_union*        mddUnion;
    mxd_intersection* mxdIntersection;
    mxd_difference*   mxdDifference;
};


// ------------------------ MTMDD Apply operations ---------------------------

/** MTMDD element-wise operations

    Abstract class for element-wise MTMDD operations.
*/
class mtmdd_apply_operation : public mdd_apply_operation {
  public:
    mtmdd_apply_operation(const char *name, bool commutative);
    virtual ~mtmdd_apply_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;
};


#if 0
class mtmdd_plus : public mtmdd_apply_operation {
#else
class mtmdd_plus : public mdd_union {
#endif
  public:
    static mtmdd_plus* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);

  protected:
    mtmdd_plus(const char* name);
    mtmdd_plus(const mtmdd_plus& copy);
    mtmdd_plus& operator=(const mtmdd_plus& copy);
    ~mtmdd_plus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_minus : public mtmdd_apply_operation {
  public:
    static mtmdd_minus* getInstance();

  protected:
    mtmdd_minus(const char* name);
    mtmdd_minus(const mtmdd_minus& copy);
    mtmdd_minus& operator=(const mtmdd_minus& copy);
    ~mtmdd_minus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_multiply : public mtmdd_apply_operation {
  public:
    static mtmdd_multiply* getInstance();

  protected:
    mtmdd_multiply(const char* name);
    mtmdd_multiply(const mtmdd_multiply& copy);
    mtmdd_multiply& operator=(const mtmdd_multiply& copy);
    ~mtmdd_multiply();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_divide : public mtmdd_apply_operation {
  public:
    static mtmdd_divide* getInstance();

  protected:
    mtmdd_divide(const char* name);
    mtmdd_divide(const mtmdd_divide& copy);
    mtmdd_divide& operator=(const mtmdd_divide& copy);
    ~mtmdd_divide();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_min : public mtmdd_apply_operation {
  public:
    static mtmdd_min* getInstance();

  protected:
    mtmdd_min(const char* name);
    mtmdd_min(const mtmdd_min& copy);
    mtmdd_min& operator=(const mtmdd_min& copy);
    ~mtmdd_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_max : public mtmdd_apply_operation {
  public:
    static mtmdd_max* getInstance();

  protected:
    mtmdd_max(const char* name);
    mtmdd_max(const mtmdd_max& copy);
    mtmdd_max& operator=(const mtmdd_max& copy);
    ~mtmdd_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_or_min : public mtmdd_apply_operation {
  public:
    static mtmdd_or_min* getInstance();

  protected:
    mtmdd_or_min(const char* name);
    mtmdd_or_min(const mtmdd_or_min& copy);
    mtmdd_or_min& operator=(const mtmdd_or_min& copy);
    ~mtmdd_or_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_or_max : public mtmdd_apply_operation {
  public:
    static mtmdd_or_max* getInstance();

  protected:
    mtmdd_or_max(const char* name);
    mtmdd_or_max(const mtmdd_or_max& copy);
    mtmdd_or_max& operator=(const mtmdd_or_max& copy);
    ~mtmdd_or_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_and_min : public mtmdd_apply_operation {
  public:
    static mtmdd_and_min* getInstance();

  protected:
    mtmdd_and_min(const char* name);
    mtmdd_and_min(const mtmdd_and_min& copy);
    mtmdd_and_min& operator=(const mtmdd_and_min& copy);
    ~mtmdd_and_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_and_max : public mtmdd_apply_operation {
  public:
    static mtmdd_and_max* getInstance();

  protected:
    mtmdd_and_max(const char* name);
    mtmdd_and_max(const mtmdd_and_max& copy);
    mtmdd_and_max& operator=(const mtmdd_and_max& copy);
    ~mtmdd_and_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_less_than : public mtmdd_apply_operation {
  public:
    static mtmdd_less_than* getInstance();

  protected:
    mtmdd_less_than(const char* name);
    mtmdd_less_than(const mtmdd_less_than& copy);
    mtmdd_less_than& operator=(const mtmdd_less_than& copy);
    ~mtmdd_less_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_less_than_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_less_than_equal* getInstance();

  protected:
    mtmdd_less_than_equal(const char* name);
    mtmdd_less_than_equal(const mtmdd_less_than_equal& copy);
    mtmdd_less_than_equal& operator=(const mtmdd_less_than_equal& copy);
    ~mtmdd_less_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_greater_than : public mtmdd_apply_operation {
  public:
    static mtmdd_greater_than* getInstance();

  protected:
    mtmdd_greater_than(const char* name);
    mtmdd_greater_than(const mtmdd_greater_than& copy);
    mtmdd_greater_than& operator=(const mtmdd_greater_than& copy);
    ~mtmdd_greater_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_greater_than_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_greater_than_equal* getInstance();

  protected:
    mtmdd_greater_than_equal(const char* name);
    mtmdd_greater_than_equal(const mtmdd_greater_than_equal& copy);
    mtmdd_greater_than_equal& operator=(const mtmdd_greater_than_equal& copy);
    ~mtmdd_greater_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_equal* getInstance();

  protected:
    mtmdd_equal(const char* name);
    mtmdd_equal(const mtmdd_equal& copy);
    mtmdd_equal& operator=(const mtmdd_equal& copy);
    ~mtmdd_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_not_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_not_equal* getInstance();

  protected:
    mtmdd_not_equal(const char* name);
    mtmdd_not_equal(const mtmdd_not_equal& copy);
    mtmdd_not_equal& operator=(const mtmdd_not_equal& copy);
    ~mtmdd_not_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// ----------------------- MTMDD-To-MDD apply operations ---------------------

/** MTMDD element-wise operations that result in MDDs

    Abstract class.
*/
class mtmdd_to_mdd_apply_operation : public mtmdd_apply_operation {
  public:
    mtmdd_to_mdd_apply_operation(const char *name, bool commutative);
    virtual ~mtmdd_to_mdd_apply_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);
    virtual int computeHelper(op_info* owner, int a, int b);

  protected:
    /// Scratch is a 2-D array.
    /// 
    /// The purpose of this array is to prevent allocation/deallocation of
    /// temporary int[] that may be needed in APPLY-like operations.
    ///
    /// It has the same number of levels as the domain. So variable j in
    /// the domain is mapped to scratch[j]. Each scratch[j] is an int[] of
    /// size == domain->getVariableBound(j)
    void buildScratch(const expert_domain* d);
    /// Checks if this scratch is compatible with domain.
    /// If it is not compatible, it rebuilds it.
    void initScratch(const expert_domain* d);
    /// Delete the scratch
    void deleteScratch();
    /// Returns a int[] of size == domain->getVariableBound(j).
    int* getScratch(int level);
    /// Assigns scratch[level][j] == val for all valid j.
    void setScratch(int level, int val);

  private:
    int** scratch;
    unsigned int nLevels;
    int* levelSizes;
};


class mtmdd_to_mdd_less_than : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_less_than* getInstance();

  protected:
    mtmdd_to_mdd_less_than(const char* name);
    mtmdd_to_mdd_less_than(const mtmdd_to_mdd_less_than& copy);
    mtmdd_to_mdd_less_than& operator=(const mtmdd_to_mdd_less_than& copy);
    ~mtmdd_to_mdd_less_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_less_than_equal : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_less_than_equal* getInstance();

  protected:
    mtmdd_to_mdd_less_than_equal(const char* name);
    mtmdd_to_mdd_less_than_equal(const mtmdd_to_mdd_less_than_equal& copy);
    mtmdd_to_mdd_less_than_equal&
      operator=(const mtmdd_to_mdd_less_than_equal& copy);
    ~mtmdd_to_mdd_less_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_greater_than : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_greater_than* getInstance();

  protected:
    mtmdd_to_mdd_greater_than(const char* name);
    mtmdd_to_mdd_greater_than(const mtmdd_to_mdd_greater_than& copy);
    mtmdd_to_mdd_greater_than&
      operator=(const mtmdd_to_mdd_greater_than& copy);
    ~mtmdd_to_mdd_greater_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_greater_than_equal : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_greater_than_equal* getInstance();

  protected:
    mtmdd_to_mdd_greater_than_equal(const char* name);
    mtmdd_to_mdd_greater_than_equal(
        const mtmdd_to_mdd_greater_than_equal& copy);
    mtmdd_to_mdd_greater_than_equal&
      operator=(const mtmdd_to_mdd_greater_than_equal& copy);
    ~mtmdd_to_mdd_greater_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_equal : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_equal* getInstance();

  protected:
    mtmdd_to_mdd_equal(const char* name);
    mtmdd_to_mdd_equal(const mtmdd_to_mdd_equal& copy);
    mtmdd_to_mdd_equal& operator=(const mtmdd_to_mdd_equal& copy);
    ~mtmdd_to_mdd_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_not_equal : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_not_equal* getInstance();

  protected:
    mtmdd_to_mdd_not_equal(const char* name);
    mtmdd_to_mdd_not_equal(const mtmdd_to_mdd_not_equal& copy);
    mtmdd_to_mdd_not_equal& operator=(const mtmdd_to_mdd_not_equal& copy);
    ~mtmdd_to_mdd_not_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


/** Conversion operations (Unary)

    Abstract class for element-wise operations that convert elements from
    one type of forest to another. For example from MTMDD to MDD and
    vice-versa.
*/
class conversion_operation : public operation {
  public:
    conversion_operation(const char *name);
    virtual ~conversion_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Calls compute(op_info*, dd_edge, dd_edge)
    virtual compute_manager::error compute(op_info* owner, dd_edge** operands);

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        dd_edge& b);

    // Returns an error
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual int compute(op_info* owner, int a);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int& b) = 0;

    virtual bool findResult(op_info* owner, int a, int& b);
    virtual void saveResult(op_info* owner, int a, int c);
};


class mtmdd_to_mdd : public conversion_operation {
  public:
    static mtmdd_to_mdd* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);

  protected:
    mtmdd_to_mdd(const char* name);
    mtmdd_to_mdd(const mtmdd_to_mdd& copy);
    mtmdd_to_mdd& operator=(const mtmdd_to_mdd& copy);
    ~mtmdd_to_mdd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

class mdd_to_mtmdd : public conversion_operation {
  public:
    static mdd_to_mtmdd* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);

  protected:
    mdd_to_mtmdd(const char* name);
    mdd_to_mtmdd(const mdd_to_mtmdd& copy);
    mdd_to_mtmdd& operator=(const mdd_to_mtmdd& copy);
    ~mdd_to_mtmdd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

class mtmxd_to_mxd : public conversion_operation {
  public:
    static mtmxd_to_mxd* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);

  protected:
    mtmxd_to_mxd(const char* name);
    mtmxd_to_mxd(const mtmxd_to_mxd& copy);
    mtmxd_to_mxd& operator=(const mtmxd_to_mxd& copy);
    ~mtmxd_to_mxd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

class mxd_to_mtmxd : public conversion_operation {
  public:
    static mxd_to_mtmxd* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);

  protected:
    mxd_to_mtmxd(const char* name);
    mxd_to_mtmxd(const mxd_to_mtmxd& copy);
    mxd_to_mtmxd& operator=(const mxd_to_mtmxd& copy);
    ~mxd_to_mtmxd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

/// Convert MTMDDs to EVMDDs
class mtmdd_to_evmdd : public operation {
  public:
    static mtmdd_to_evmdd* getInstance();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;
    // Calls compute(op_info*, dd_edge, dd_edge)
    virtual compute_manager::error compute(op_info* owner, dd_edge** operands);
    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        dd_edge& b);
    // Returns an error
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

  protected:
    mtmdd_to_evmdd(const char* name);
    mtmdd_to_evmdd(const mtmdd_to_evmdd& copy);
    mtmdd_to_evmdd& operator=(const mtmdd_to_evmdd& copy);
    ~mtmdd_to_evmdd();

    virtual void compute(op_info* owner, int a, int& b, int &ev);
    virtual bool checkTerminals(op_info* op, int a, int& b, int &ev);
    virtual bool findResult(op_info* owner, int a, int& b, int &ev);
    virtual void saveResult(op_info* owner, int a, int b, int ev);

    // for later: once evmdds can use real edge-values
    virtual bool checkTerminals(op_info* op, int a, int& b, float &ev);
    virtual bool findResult(op_info* owner, int a, int& b, float &ev);
    virtual void saveResult(op_info* owner, int a, int b, float ev);
};


class mtmdd_post_image : public mdd_post_image {
  public:
    static mtmdd_post_image* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mtmdd_post_image(const char* name);
    mtmdd_post_image(const mtmdd_post_image& copy);
    mtmdd_post_image& operator=(const mtmdd_post_image& copy);
    virtual ~mtmdd_post_image();

    virtual int compute(op_info* owner, op_info* plusOp, int mdd, int mxd);
};


class mtmdd_pre_image : public mdd_pre_image {
  public:
    static mtmdd_pre_image* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);

  protected:
    mtmdd_pre_image(const char* name);
    mtmdd_pre_image(const mtmdd_pre_image& copy);
    mtmdd_pre_image& operator=(const mtmdd_pre_image& copy);
    virtual ~mtmdd_pre_image();

    virtual int compute(op_info* owner, op_info* plusOp, int mdd, int mxd);
};


#endif


