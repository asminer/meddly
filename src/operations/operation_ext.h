
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

#include "../defines.h"
#include <vector>

using namespace MEDDLY;

/** MDD element-wise operations

    Abstract class for element-wise MDD operations.
*/
class mdd_apply_operation : public old_operation {
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

    virtual void typeCheck(const op_info* owner);
    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Calls compute(op_info*, dd_edge, dd_edge, dd_edge)
    virtual void compute(op_info* owner, dd_edge** operands);

    // Returns an error
    virtual void compute(op_info* owner, const dd_edge& a,
        dd_edge& b);

    // Calls compute(op_info*, int, int) and stores result in c.
    virtual void compute(op_info* owner, const dd_edge& a,
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

    virtual void typeCheck(const op_info* owner);
    virtual int compute(op_info* owner, int a, int b);
    virtual const char* getName() const { return "Mxd Apply"; }

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;

    virtual int computeIdent(op_info* owner, int a, int b);
    virtual int computeNonIdent(op_info* owner, int a, int b);

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


class mxd_alt_apply_operation : public old_operation {
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

    virtual void typeCheck(const op_info* owner);
    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Calls compute(op_info*, dd_edge, dd_edge, dd_edge)
    virtual void compute(op_info* owner, dd_edge** operands);

    // Returns an error
    virtual void compute(op_info* owner, const dd_edge& a,
        dd_edge& b);

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual void compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual int compute(op_info* owner, int level, int a, int b);

  protected:
    /// To be implemented by derived classes
    virtual bool checkTerminals(op_info* op, int a, int b, int& c) = 0;

    virtual int computeIdent(op_info* owner, int level, int a, int b);
    virtual int computeNonIdent(op_info* owner, int level, int a, int b);

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

    virtual int computeIdent(op_info* owner, int a, int b);
    virtual int computeIdentExpandA(op_info* owner, int a, int b);
    virtual int computeIdentExpandOneLevel(op_info* owner, int a, int b);

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

/** MTXDD element-wise operations

    Abstract class for element-wise MTMXD operations.
*/
// TODO: TEST
class mtmxd_apply_operation : public mxd_apply_operation {
  public:
    mtmxd_apply_operation();
    virtual ~mtmxd_apply_operation();
    virtual void typeCheck(const op_info* owner);
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


// same as mtmxd_apply_operation, except this is for mtmxd operations
// with op(0,0)!=0 (such as /, <=, >=, ==)
class mtmxd_alt_apply_operation : public mxd_alt_apply_operation {
  public:
    mtmxd_alt_apply_operation();
    virtual ~mtmxd_alt_apply_operation();
    virtual void typeCheck(const op_info* owner);
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



// ------------------------ MTMDD Apply operations ---------------------------

/** MTMDD element-wise operations

    Abstract class for element-wise MTMDD operations.
*/
class mtmdd_apply_operation : public mdd_apply_operation {
  public:
    mtmdd_apply_operation();
    virtual ~mtmdd_apply_operation();

    virtual void typeCheck(const op_info* owner);
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
    virtual void typeCheck(const op_info* owner);
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







// --------------------------- EVMDD operations ------------------------------

class evmdd_apply_operation : public old_operation {
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
    virtual void compute(op_info* owner, dd_edge** operands);

    // Returns an error
    virtual void compute(op_info* owner, const dd_edge& a,
        dd_edge& b);

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual void compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c) = 0;
};


class evplusmdd_apply_operation : public evmdd_apply_operation {
  public:
    evplusmdd_apply_operation();
    virtual ~evplusmdd_apply_operation();

    virtual const char* getName() const { return "EV+MDD Apply"; }

    virtual void typeCheck(const op_info* owner);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual void compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual void compute(op_info* owner, int a, int aev,
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

    virtual void typeCheck(const op_info* owner);
    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    // Implements APPLY operation -- calls checkTerminals to compute
    // result for terminal nodes.
    virtual void compute(op_info* owner, const dd_edge& a,
        const dd_edge& b, dd_edge& c);

    virtual void compute(op_info* owner, int a, float aev,
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
    static evplusmdd_plus* getInstance();
    virtual const char* getName() const { return "EV+MDD Plus"; }
    virtual bool isCommutative() const { return true; }
    virtual bool checkTerminals(op_info* op, int a, int aev, int b, int bev,
        int& c, int& cev);
    virtual void compute(op_info* owner, int a, int aev,
        int b, int bev, int& c, int& cev);

  protected:
    evplusmdd_plus() {}
    evplusmdd_plus(const evplusmdd_plus& copy);
    evplusmdd_plus& operator=(const evplusmdd_plus& copy);
    virtual ~evplusmdd_plus() {}
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


inline expert_forest* getExpertForest(op_info* op, int index) {
  return op->p[index].getForest();
  // return smart_cast<expert_forest*>(op->f[index]);
}

inline const expert_forest* getExpertForest(const op_info* op, int index) {
  return op->p[index].readForest();
  // return smart_cast<const expert_forest*>(op->f[index]);
}


#endif


