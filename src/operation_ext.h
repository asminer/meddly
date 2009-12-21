
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


//#define ALT_SATURATE_HELPER


/** MDD element-wise operations

    Abstract class for element-wise MDD operations.
*/
class mdd_apply_operation : public operation {
  public:
    mdd_apply_operation();
    virtual ~mdd_apply_operation();

    virtual const char* getName() const { return "Mdd Apply"; }
    virtual int getKeyLength() const { return 2; }
    virtual int getAnsLength() const { return 1; }
    virtual int getCacheEntryLength() const { return 3; }

    virtual int getKeyLengthInBytes() const { return 8; }
    virtual int getAnsLengthInBytes() const { return 4; }
    virtual int getCacheEntryLengthInBytes() const { return 12; }

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

    // Calls compute(op_info*, int, int) and stores result in c.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    // Implements APPLY. Returns the result of (a op b).
    // Calls checkTerminals() for termination condition of recursion.
    virtual int compute(op_info* owner, int a, int b);

  protected:
    // If terminal condition is reached, returns true and the result in c.
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;

    // Returns true if the operation is commutative (i.e. A op B == B op A)
    virtual bool isCommutative() const = 0;

    // Searches the compute table for (a op b). If the result is found,
    // returns true and the result via c (c's incount is incremented by 1).
    // Otherwise, returns false.
    virtual bool findResult(op_info* owner, int a, int b, int& c);

    // Stores the tuple (a op b == c) in the compute table.
    // Cachecounts for a, b, and c are incremented by 1.
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
};


class mdd_union : public mdd_apply_operation {
  public:
    static mdd_union* getInstance();
#if 0
    virtual int compute(op_info* owner, int a, int b);
#endif

    virtual const char* getName() const { return "Mdd Union"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mdd_union();
    mdd_union(const mdd_union& copy);
    mdd_union& operator=(const mdd_union& copy);
    ~mdd_union();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);

#if 0
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
#if 0
    virtual int compute(op_info* owner, int a, int b);
#endif

    virtual const char* getName() const { return "Mdd Intersection"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mdd_intersection();
    mdd_intersection(const mdd_intersection& copy);
    mdd_intersection& operator=(const mdd_intersection& copy);
    ~mdd_intersection();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
#if 0
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

    virtual const char* getName() const { return "Mdd Difference"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_difference();
    mdd_difference(const mdd_difference& copy);
    mdd_difference& operator=(const mdd_difference& copy);
    ~mdd_difference();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// ------------------ MXD operations -------------------------------

// HERE: Test!!

class mxd_apply_operation : public mdd_apply_operation {
  public:
    mxd_apply_operation();
    virtual ~mxd_apply_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mxd Apply"; }

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
    mxd_alt_apply_operation();
    virtual ~mxd_alt_apply_operation();
    virtual const char* getName() const { return "Mxd Alternative Apply"; }

    virtual int getKeyLength() const { return 3; }
    virtual int getAnsLength() const { return 1; }
    virtual int getCacheEntryLength() const { return 4; }

    virtual int getKeyLengthInBytes() const { return 12; }
    virtual int getAnsLengthInBytes() const { return 4; }
    virtual int getCacheEntryLengthInBytes() const { return 16; }

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

    // Returns true if the operation is commutative (i.e. A op B == B op A)
    virtual bool isCommutative() const = 0;

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
};


class mxd_union : public mxd_apply_operation {
  public:
    static mxd_union* getInstance();
    virtual const char* getName() const { return "Mxd Union"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mxd_union();
    mxd_union(const mxd_union& copy);
    mxd_union& operator=(const mxd_union& copy);
    ~mxd_union();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mxd_intersection : public mxd_apply_operation {
  public:
    static mxd_intersection* getInstance();
    virtual const char* getName() const { return "Mxd Intersection"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mxd_intersection();
    mxd_intersection(const mxd_intersection& copy);
    mxd_intersection& operator=(const mxd_intersection& copy);
    ~mxd_intersection();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mxd_difference : public mxd_apply_operation {
  public:
    static mxd_difference* getInstance();
    virtual const char* getName() const { return "Mxd Difference"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mxd_difference();
    mxd_difference(const mxd_difference& copy);
    mxd_difference& operator=(const mxd_difference& copy);
    ~mxd_difference();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// ------------------------ MTMXD Apply operations ---------------------------

// HERE --------

/** MTXDD element-wise operations

    Abstract class for element-wise MTMXD operations.
*/
// TODO: TEST
class mtmxd_apply_operation : public mxd_apply_operation {
  public:
    mtmxd_apply_operation();
    virtual ~mtmxd_apply_operation();
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual const char* getName() const { return "MtMxd Apply"; }

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;
};


class mtmxd_plus : public mtmxd_apply_operation {
  public:
    static mtmxd_plus* getInstance();
    virtual const char* getName() const { return "MtMxd Plus"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmxd_plus();
    mtmxd_plus(const mtmxd_plus& copy);
    mtmxd_plus& operator=(const mtmxd_plus& copy);
    ~mtmxd_plus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_minus : public mtmxd_apply_operation {
  public:
    static mtmxd_minus* getInstance();
    virtual const char* getName() const { return "MtMxd Minus"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmxd_minus();
    mtmxd_minus(const mtmxd_minus& copy);
    mtmxd_minus& operator=(const mtmxd_minus& copy);
    ~mtmxd_minus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_multiply : public mtmxd_apply_operation {
  public:
    static mtmxd_multiply* getInstance();
    virtual const char* getName() const { return "MtMxd Multiply"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmxd_multiply();
    mtmxd_multiply(const mtmxd_multiply& copy);
    mtmxd_multiply& operator=(const mtmxd_multiply& copy);
    ~mtmxd_multiply();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_min : public mtmxd_apply_operation {
  public:
    static mtmxd_min* getInstance();
    virtual const char* getName() const { return "MtMxd Min"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmxd_min();
    mtmxd_min(const mtmxd_min& copy);
    mtmxd_min& operator=(const mtmxd_min& copy);
    ~mtmxd_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_max : public mtmxd_apply_operation {
  public:
    static mtmxd_max* getInstance();
    virtual const char* getName() const { return "MtMxd Max"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmxd_max();
    mtmxd_max(const mtmxd_max& copy);
    mtmxd_max& operator=(const mtmxd_max& copy);
    ~mtmxd_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_not_equal : public mtmxd_apply_operation {
  public:
    static mtmxd_not_equal* getInstance();
    virtual const char* getName() const { return "MtMxd Not-Equal"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmxd_not_equal();
    mtmxd_not_equal(const mtmxd_not_equal& copy);
    mtmxd_not_equal& operator=(const mtmxd_not_equal& copy);
    ~mtmxd_not_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_less_than : public mtmxd_apply_operation {
  public:
    static mtmxd_less_than* getInstance();
    virtual const char* getName() const { return "MtMxd Less-Than"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmxd_less_than();
    mtmxd_less_than(const mtmxd_less_than& copy);
    mtmxd_less_than& operator=(const mtmxd_less_than& copy);
    ~mtmxd_less_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_greater_than : public mtmxd_apply_operation {
  public:
    static mtmxd_greater_than* getInstance();
    virtual const char* getName() const { return "MtMxd Greater-Than"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmxd_greater_than();
    mtmxd_greater_than(const mtmxd_greater_than& copy);
    mtmxd_greater_than& operator=(const mtmxd_greater_than& copy);
    ~mtmxd_greater_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


// same as mtmxd_apply_operation, except this is for mtmxd operations
// with op(0,0)!=0 (such as /, <=, >=, ==)
class mtmxd_alt_apply_operation : public mxd_alt_apply_operation {
  public:
    mtmxd_alt_apply_operation();
    virtual ~mtmxd_alt_apply_operation();
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual const char* getName() const { return "MtMxd Alternative Apply"; }
  protected:
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;
};


class mtmxd_divide : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_divide* getInstance();
    virtual const char* getName() const { return "MtMxd Divide"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmxd_divide();
    mtmxd_divide(const mtmxd_divide& copy);
    mtmxd_divide& operator=(const mtmxd_divide& copy);
    ~mtmxd_divide();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_less_than_equal : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_less_than_equal* getInstance();
    virtual const char* getName() const { return "MtMxd Less-Than-Or-Equal"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmxd_less_than_equal();
    mtmxd_less_than_equal(const mtmxd_less_than_equal& copy);
    mtmxd_less_than_equal& operator=(const mtmxd_less_than_equal& copy);
    ~mtmxd_less_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_greater_than_equal : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_greater_than_equal* getInstance();
    virtual const char* getName() const {
      return "MtMxd Greater-Than-Or-Equal";
    }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmxd_greater_than_equal();
    mtmxd_greater_than_equal(const mtmxd_greater_than_equal& copy);
    mtmxd_greater_than_equal& operator=(const mtmxd_greater_than_equal& copy);
    ~mtmxd_greater_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmxd_equal : public mtmxd_alt_apply_operation {
  public:
    static mtmxd_equal* getInstance();
    virtual const char* getName() const { return "MtMxd Equal"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmxd_equal();
    mtmxd_equal(const mtmxd_equal& copy);
    mtmxd_equal& operator=(const mtmxd_equal& copy);
    ~mtmxd_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


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


// Traditional reachability analysis
class mdd_reachability_bfs : public mdd_mxd_image_operation {
  public:
    static mdd_reachability_bfs* getInstance();
    int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Reachability BFS"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_reachability_bfs();
    mdd_reachability_bfs(const mdd_reachability_bfs& copy);
    mdd_reachability_bfs& operator=(const mdd_reachability_bfs& copy);
    virtual ~mdd_reachability_bfs();
};


// Reachability via "saturation" algorithm
class mdd_reachability_dfs : public mdd_mxd_image_operation {
  public:
    static mdd_reachability_dfs* getInstance();
    int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mdd-Mxd Reachability DFS"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mdd_reachability_dfs();
    mdd_reachability_dfs(const mdd_reachability_dfs& copy);
    mdd_reachability_dfs& operator=(const mdd_reachability_dfs& copy);
    virtual ~mdd_reachability_dfs();

    void initialize(op_info* owner);
    void clear();

    void splitMxd(int mxd);
    int saturate(int mdd);
    int recFire(int mdd, int mxd);
#ifdef ALT_SATURATE_HELPER
    void saturateHelper(int mdd);
    void saturateHelperUnPrimeFull(int mdd, int mxd);
    void saturateHelperUnPrimeSparse(int mdd, int mxd);
    void saturateHelperPrimeFull(int mdd, int i, int mxdI,
      std::vector<bool>& next);
    void saturateHelperPrimeSparse(int mdd, int i, int mxdI,
      std::vector<bool>& next);

    void recFireExpandMdd(int mdd, int mxd, int result);
    void recFireExpandMxd(int mdd, int mxd, int result);
    void recFireFF(int mdd, int mxd, int result);
    void recFireFS(int mdd, int mxd, int result);
    void recFireSF(int mdd, int mxd, int result);
    void recFireSS(int mdd, int mxd, int result);
    void recFirePrime(int mdd, int mxd, int result);
#else
    void saturateHelper(int mddLevel, std::vector<int>& mdd);
#endif

    int getMddUnion(int a, int b);
    int getMxdIntersection(int a, int b);
    int getMxdDifference(int a, int b);

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
    mtmdd_apply_operation();
    virtual ~mtmdd_apply_operation();

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual const char* getName() const { return "MtMdd Apply"; }

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
    virtual const char* getName() const { return "MtMdd Plus"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_plus();
    mtmdd_plus(const mtmdd_plus& copy);
    mtmdd_plus& operator=(const mtmdd_plus& copy);
    ~mtmdd_plus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_minus : public mtmdd_apply_operation {
  public:
    static mtmdd_minus* getInstance();
    virtual const char* getName() const { return "MtMdd Minus"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_minus();
    mtmdd_minus(const mtmdd_minus& copy);
    mtmdd_minus& operator=(const mtmdd_minus& copy);
    ~mtmdd_minus();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_multiply : public mtmdd_apply_operation {
  public:
    static mtmdd_multiply* getInstance();
    virtual const char* getName() const { return "MtMdd Multiply"; }
    virtual bool isCommutative() const { return true; }


  protected:
    mtmdd_multiply();
    mtmdd_multiply(const mtmdd_multiply& copy);
    mtmdd_multiply& operator=(const mtmdd_multiply& copy);
    ~mtmdd_multiply();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_divide : public mtmdd_apply_operation {
  public:
    static mtmdd_divide* getInstance();
    virtual const char* getName() const { return "MtMdd Divide"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_divide();
    mtmdd_divide(const mtmdd_divide& copy);
    mtmdd_divide& operator=(const mtmdd_divide& copy);
    ~mtmdd_divide();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_min : public mtmdd_apply_operation {
  public:
    static mtmdd_min* getInstance();
    virtual const char* getName() const { return "MtMdd Min"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_min();
    mtmdd_min(const mtmdd_min& copy);
    mtmdd_min& operator=(const mtmdd_min& copy);
    ~mtmdd_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_max : public mtmdd_apply_operation {
  public:
    static mtmdd_max* getInstance();
    virtual const char* getName() const { return "MtMdd Max"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_max();
    mtmdd_max(const mtmdd_max& copy);
    mtmdd_max& operator=(const mtmdd_max& copy);
    ~mtmdd_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_or_min : public mtmdd_apply_operation {
  public:
    static mtmdd_or_min* getInstance();
    virtual const char* getName() const { return "MtMdd Or-Min"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_or_min();
    mtmdd_or_min(const mtmdd_or_min& copy);
    mtmdd_or_min& operator=(const mtmdd_or_min& copy);
    ~mtmdd_or_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_or_max : public mtmdd_apply_operation {
  public:
    static mtmdd_or_max* getInstance();
    virtual const char* getName() const { return "MtMdd Or-Max"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_or_max();
    mtmdd_or_max(const mtmdd_or_max& copy);
    mtmdd_or_max& operator=(const mtmdd_or_max& copy);
    ~mtmdd_or_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_and_min : public mtmdd_apply_operation {
  public:
    static mtmdd_and_min* getInstance();
    virtual const char* getName() const { return "MtMdd And-Min"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_and_min();
    mtmdd_and_min(const mtmdd_and_min& copy);
    mtmdd_and_min& operator=(const mtmdd_and_min& copy);
    ~mtmdd_and_min();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_and_max : public mtmdd_apply_operation {
  public:
    static mtmdd_and_max* getInstance();
    virtual const char* getName() const { return "MtMdd And-Max"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_and_max();
    mtmdd_and_max(const mtmdd_and_max& copy);
    mtmdd_and_max& operator=(const mtmdd_and_max& copy);
    ~mtmdd_and_max();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_less_than : public mtmdd_apply_operation {
  public:
    static mtmdd_less_than* getInstance();
    virtual const char* getName() const { return "MtMdd Less-Than"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_less_than();
    mtmdd_less_than(const mtmdd_less_than& copy);
    mtmdd_less_than& operator=(const mtmdd_less_than& copy);
    ~mtmdd_less_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_less_than_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_less_than_equal* getInstance();
    virtual const char* getName() const { return "MtMdd Less-Than-Or-Equal"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_less_than_equal();
    mtmdd_less_than_equal(const mtmdd_less_than_equal& copy);
    mtmdd_less_than_equal& operator=(const mtmdd_less_than_equal& copy);
    ~mtmdd_less_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_greater_than : public mtmdd_apply_operation {
  public:
    static mtmdd_greater_than* getInstance();
    virtual const char* getName() const { return "MtMdd Greater-Than"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_greater_than();
    mtmdd_greater_than(const mtmdd_greater_than& copy);
    mtmdd_greater_than& operator=(const mtmdd_greater_than& copy);
    ~mtmdd_greater_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_greater_than_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_greater_than_equal* getInstance();
    virtual const char* getName() const {
      return "MtMdd Greater-Than-Or-Equal";
    }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_greater_than_equal();
    mtmdd_greater_than_equal(const mtmdd_greater_than_equal& copy);
    mtmdd_greater_than_equal& operator=(const mtmdd_greater_than_equal& copy);
    ~mtmdd_greater_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_equal* getInstance();
    virtual const char* getName() const { return "MtMdd Equal"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_equal();
    mtmdd_equal(const mtmdd_equal& copy);
    mtmdd_equal& operator=(const mtmdd_equal& copy);
    ~mtmdd_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_not_equal : public mtmdd_apply_operation {
  public:
    static mtmdd_not_equal* getInstance();
    virtual const char* getName() const { return "MtMdd Not-Equal"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_not_equal();
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
    mtmdd_to_mdd_apply_operation();
    virtual ~mtmdd_to_mdd_apply_operation();
    virtual const char* getName() const { return "MtMdd To Mdd Apply"; }

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
    virtual const char* getName() const { return "MtMdd To Mdd Less-Than"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_to_mdd_less_than();
    mtmdd_to_mdd_less_than(const mtmdd_to_mdd_less_than& copy);
    mtmdd_to_mdd_less_than& operator=(const mtmdd_to_mdd_less_than& copy);
    ~mtmdd_to_mdd_less_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_less_than_equal : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_less_than_equal* getInstance();
    virtual const char* getName() const {
      return "MtMdd To Mdd Less-Than-Or-Equal";
    }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_to_mdd_less_than_equal();
    mtmdd_to_mdd_less_than_equal(const mtmdd_to_mdd_less_than_equal& copy);
    mtmdd_to_mdd_less_than_equal&
      operator=(const mtmdd_to_mdd_less_than_equal& copy);
    ~mtmdd_to_mdd_less_than_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_greater_than : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_greater_than* getInstance();
    virtual const char* getName() const { return "MtMdd To Mdd Greater-Than"; }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_to_mdd_greater_than();
    mtmdd_to_mdd_greater_than(const mtmdd_to_mdd_greater_than& copy);
    mtmdd_to_mdd_greater_than&
      operator=(const mtmdd_to_mdd_greater_than& copy);
    ~mtmdd_to_mdd_greater_than();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_greater_than_equal : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_greater_than_equal* getInstance();
    virtual const char* getName() const {
      return "MtMdd To Mdd Greater-Than-Or-Equal";
    }
    virtual bool isCommutative() const { return false; }

  protected:
    mtmdd_to_mdd_greater_than_equal();
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
    virtual const char* getName() const { return "MtMdd To Mdd Equal"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_to_mdd_equal();
    mtmdd_to_mdd_equal(const mtmdd_to_mdd_equal& copy);
    mtmdd_to_mdd_equal& operator=(const mtmdd_to_mdd_equal& copy);
    ~mtmdd_to_mdd_equal();

    virtual bool checkTerminals(op_info* op, int a, int b, int& c);
};


class mtmdd_to_mdd_not_equal : public mtmdd_to_mdd_apply_operation {
  public:
    static mtmdd_to_mdd_not_equal* getInstance();
    virtual const char* getName() const { return "MtMdd To Mdd Not-Equal"; }
    virtual bool isCommutative() const { return true; }

  protected:
    mtmdd_to_mdd_not_equal();
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
    conversion_operation();
    virtual ~conversion_operation();

    virtual const char* getName() const { return "Conversion operation"; }

    virtual int getKeyLength() const { return 1; }
    virtual int getAnsLength() const { return 1; }
    virtual int getCacheEntryLength() const { return 2; }

    virtual int getKeyLengthInBytes() const { return 4; }
    virtual int getAnsLengthInBytes() const { return 4; }
    virtual int getCacheEntryLengthInBytes() const { return 8; }

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
    virtual const char* getName() const { return "Convert MtMdd to Mdd"; }

  protected:
    mtmdd_to_mdd();
    mtmdd_to_mdd(const mtmdd_to_mdd& copy);
    mtmdd_to_mdd& operator=(const mtmdd_to_mdd& copy);
    ~mtmdd_to_mdd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

class mdd_to_mtmdd : public conversion_operation {
  public:
    static mdd_to_mtmdd* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual const char* getName() const { return "Convert Mdd to MtMdd"; }

  protected:
    mdd_to_mtmdd();
    mdd_to_mtmdd(const mdd_to_mtmdd& copy);
    mdd_to_mtmdd& operator=(const mdd_to_mtmdd& copy);
    ~mdd_to_mtmdd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

class mtmxd_to_mxd : public conversion_operation {
  public:
    static mtmxd_to_mxd* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual const char* getName() const { return "Convert MtMxd to Mxd"; }

  protected:
    mtmxd_to_mxd();
    mtmxd_to_mxd(const mtmxd_to_mxd& copy);
    mtmxd_to_mxd& operator=(const mtmxd_to_mxd& copy);
    ~mtmxd_to_mxd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

class mxd_to_mtmxd : public conversion_operation {
  public:
    static mxd_to_mtmxd* getInstance();
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual const char* getName() const { return "Convert Mxd to MtMxd"; }

  protected:
    mxd_to_mtmxd();
    mxd_to_mtmxd(const mxd_to_mtmxd& copy);
    mxd_to_mtmxd& operator=(const mxd_to_mtmxd& copy);
    ~mxd_to_mtmxd();

    virtual bool checkTerminals(op_info* op, int a, int& b);
};

/// Convert MTMDDs to EVMDDs
class mtmdd_to_evmdd : public operation {
  public:
    static mtmdd_to_evmdd* getInstance();

    virtual const char* getName() const { return "Convert MtMdd to EvMdd"; }

    virtual int getKeyLength() const { return 1; }
    virtual int getAnsLength() const { return 2; }
    virtual int getCacheEntryLength() const { return 3; }

    virtual int getKeyLengthInBytes() const { return 4; }
    virtual int getAnsLengthInBytes() const { return 8; }
    virtual int getCacheEntryLengthInBytes() const { return 12; }

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
    mtmdd_to_evmdd();
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


class mdd_to_evplusmdd_index_set : public mtmdd_to_evmdd {
  public:
    static mdd_to_evplusmdd_index_set* getInstance();
    virtual const char* getName() const { return "Convert Mdd to Index Set"; }
    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual compute_manager::error compute(op_info* owner,
        const dd_edge& a, dd_edge& b);

  protected:
    mdd_to_evplusmdd_index_set();
    mdd_to_evplusmdd_index_set(const mdd_to_evplusmdd_index_set& copy);
    mdd_to_evplusmdd_index_set& operator=(const
        mdd_to_evplusmdd_index_set& copy);
    ~mdd_to_evplusmdd_index_set();

    virtual void compute(op_info* owner, int a, int height, int& b, int& bev);
};


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


class evmdd_apply_operation : public operation {
  public:
    evmdd_apply_operation();
    virtual ~evmdd_apply_operation();

    virtual const char* getName() const { return "EVMDD Apply"; }

    virtual int getKeyLength() const { return 4; }
    virtual int getAnsLength() const { return 2; }
    virtual int getCacheEntryLength() const { return 6; }

    virtual int getKeyLengthInBytes() const { return 16; }
    virtual int getAnsLengthInBytes() const { return 8; }
    virtual int getCacheEntryLengthInBytes() const { return 24; }

    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);

    // Calls compute(op_info*, dd_edge, dd_edge, dd_edge)
    virtual compute_manager::error compute(op_info* owner, dd_edge** operands);

    // Returns an error
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        dd_edge& b);

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c) = 0;
};


class evplusmdd_apply_operation : public evmdd_apply_operation {
  public:
    evplusmdd_apply_operation();
    virtual ~evplusmdd_apply_operation();

    virtual const char* getName() const { return "EV+MDD Apply"; }

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual compute_manager::error compute(op_info* owner, int a, int aev,
        int b, int bev, int& c, int& cev);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int aev, int b, int bev,
        int& c, int& cev) = 0;

    // Returns true if the operation is commutative (i.e. A op B == B op A)
    virtual bool isCommutative() const = 0;

    virtual bool findResult(op_info* owner, int a, int aev,
        int b, int bev, int& c, int& cev);
    virtual void saveResult(op_info* owner, int a, int aev,
        int b, int bev, int c, int cev);
};


class evtimesmdd_apply_operation : public evmdd_apply_operation {
  public:
    evtimesmdd_apply_operation();
    virtual ~evtimesmdd_apply_operation();

    virtual const char* getName() const { return "EV*MDD Apply"; }

    virtual compute_manager::error typeCheck(const op_info* owner);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual compute_manager::error compute(op_info* owner, int a, float aev,
        int b, float bev, int& c, float& cev);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, float aev, int b, float bev,
        int& c, float& cev) = 0;

    // Returns true if the operation is commutative (i.e. A op B == B op A)
    virtual bool isCommutative() const = 0;

    virtual bool findResult(op_info* owner, int a, float aev,
        int b, float bev, int& c, float& cev);
    virtual void saveResult(op_info* owner, int a, float aev,
        int b, float bev, int c, float cev);
};


class evplusmdd_plus : public evplusmdd_apply_operation {
  public:
    virtual const char* getName() const { return "EV+MDD Plus"; }
    virtual bool isCommutative() const { return true; }
    virtual bool checkTerminals(op_info* op, int a, int aev, int b, int bev,
        int& c, int& cev);
    static evplusmdd_plus* getInstance();

  protected:
    evplusmdd_plus() {}
    evplusmdd_plus(const evplusmdd_plus& copy);
    evplusmdd_plus& operator=(const evplusmdd_plus& copy);
    virtual ~evplusmdd_plus() {}
};


class evplusmdd_multiply : public evplusmdd_apply_operation {
  public:
    virtual const char* getName() const { return "EV+MDD Multiply"; }
    virtual bool isCommutative() const { return true; }
    virtual bool checkTerminals(op_info* op, int a, int aev, int b, int bev,
        int& c, int& cev);
    static evplusmdd_multiply* getInstance();

  protected:
    evplusmdd_multiply() {}
    evplusmdd_multiply(const evplusmdd_multiply& copy);
    evplusmdd_multiply& operator=(const evplusmdd_multiply& copy);
    virtual ~evplusmdd_multiply() {}
};


class evplusmdd_minus : public evplusmdd_apply_operation {
  public:
    virtual const char* getName() const { return "EV+MDD Minus"; }
    virtual bool isCommutative() const { return false; }
    virtual bool checkTerminals(op_info* op, int a, int aev, int b, int bev,
        int& c, int& cev);
    static evplusmdd_minus* getInstance();

  protected:
    evplusmdd_minus() {}
    evplusmdd_minus(const evplusmdd_minus& copy);
    evplusmdd_minus& operator=(const evplusmdd_minus& copy);
    virtual ~evplusmdd_minus() {}
};


class evtimesmdd_plus : public evtimesmdd_apply_operation {
  public:
    virtual const char* getName() const { return "EV*MDD Plus"; }
    virtual bool isCommutative() const { return true; }
    virtual bool checkTerminals(op_info* op, int a, float aev, int b, float bev,
        int& c, float& cev);
    static evtimesmdd_plus* getInstance();

  protected:
    evtimesmdd_plus() {}
    evtimesmdd_plus(const evtimesmdd_plus& copy);
    evtimesmdd_plus& operator=(const evtimesmdd_plus& copy);
    virtual ~evtimesmdd_plus() {}
};


class evtimesmdd_multiply : public evtimesmdd_apply_operation {
  public:
    virtual const char* getName() const { return "EV*MDD Multiply"; }
    virtual bool isCommutative() const { return true; }
    virtual bool checkTerminals(op_info* op, int a, float aev, int b, float bev,
        int& c, float& cev);
    static evtimesmdd_multiply* getInstance();

  protected:
    evtimesmdd_multiply() {}
    evtimesmdd_multiply(const evtimesmdd_multiply& copy);
    evtimesmdd_multiply& operator=(const evtimesmdd_multiply& copy);
    virtual ~evtimesmdd_multiply() {}
};


class evtimesmdd_minus : public evtimesmdd_apply_operation {
  public:
    virtual const char* getName() const { return "EV*MDD Minus"; }
    virtual bool isCommutative() const { return false; }
    virtual bool checkTerminals(op_info* op, int a, float aev, int b, float bev,
        int& c, float& cev);
    static evtimesmdd_minus* getInstance();

  protected:
    evtimesmdd_minus() {}
    evtimesmdd_minus(const evtimesmdd_minus& copy);
    evtimesmdd_minus& operator=(const evtimesmdd_minus& copy);
    virtual ~evtimesmdd_minus() {}
};


#endif


