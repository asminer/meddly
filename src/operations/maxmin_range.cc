
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
#include "maxmin_range.h"
#include "../compute_cache.h"

// ******************************************************************
// *                                                                *
// *                      maxminrange_op class                      *
// *                                                                *
// ******************************************************************

/** Abstract base class for max, min range operations.
*/
class maxminrange_op : public operation {
  public:
    maxminrange_op()                        { }
    virtual ~maxminrange_op()               { }

    virtual int getKeyLength() const        { return 1; }
    virtual int getKeyLengthInBytes() const { return sizeof(int); }

    virtual bool isEntryStale(const op_info* owner, const int* entryData);
    virtual void discardEntry(op_info* owner, const int* entryData);

  protected:
    inline compute_manager::error 
    type_check(const op_info* owner, op_param::type p1type) 
    {
        if (owner == 0)
          return compute_manager::UNKNOWN_OPERATION;
        if (owner->op == 0 || owner->p == 0 || owner->cc == 0)
          return compute_manager::TYPE_MISMATCH;
        if (owner->nParams != 2)
          return compute_manager::WRONG_NUMBER;
        if (owner->p[1].getType() != p1type)
          return compute_manager::TYPE_MISMATCH;
        return compute_manager::SUCCESS;
    }
};

bool maxminrange_op::
isEntryStale(const op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] <--> DD node to find max/min range of
  // data[1] <--> result value (not a DD node)
  DCASSERT(owner->nParams == 2);
  return owner->p[0].getForest()->isStale(data[0]);
}

void maxminrange_op::
discardEntry(op_info* owner, const int* data)
{
  // data[] is of size owner.nParams
  // data[0] <--> DD node to find max/min range of
  // data[1] <--> result value (not a DD node)
  DCASSERT(owner->nParams == 2);
  owner->p[0].getForest()->uncacheNode(data[0]);
}


// ******************************************************************
// *                                                                *
// *                    int_maxminrange_op class                    *
// *                                                                *
// ******************************************************************

/** Abstract base class for max, min range operations with integer range.
*/
class int_maxminrange_op : public maxminrange_op {
  public:
    int_maxminrange_op()                    { }
    virtual ~int_maxminrange_op()           { }

    virtual int getAnsLength() const  { 
      return sizeof(long) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(long) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(long); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(long);
    }

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, op_param::INTEGER);
    }

    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
      long& b);

    // for derived classes
    virtual long compute(op_info* owner, int a) = 0;
};

void
int_maxminrange_op::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): %ld(L)]",
      owner->op->getName(), data[0], ((const long*)(data+1))[0]
  );
}

compute_manager::error 
int_maxminrange_op::
compute(op_info* owner, const dd_edge& a, long& b)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  b = compute(owner, a.getNode());
  return compute_manager::SUCCESS;
}


// ******************************************************************
// *                                                                *
// *                     int_max_range_op class                     *
// *                                                                *
// ******************************************************************

class int_max_range_op : public int_maxminrange_op {
  public:
    int_max_range_op()                    { }
    virtual ~int_max_range_op()           { }
    virtual const char* getName() const   { return "Max_range"; }

    static int_max_range_op* getInstance();

    virtual long compute(op_info* owner, int a);
};

int_max_range_op* 
int_max_range_op::
getInstance()
{
  static int_max_range_op instance;
  return &instance;
}

long int_max_range_op::compute(op_info* owner, int a)
{
  expert_forest* f = owner->p[0].getForest();

  // Terminal case
  if (f->isTerminalNode(a)) return f->getInteger(a);
  
  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
  }

  // recurse
  long max;
  if (f->isFullNode(a)) {
    // Full node
    int asize = f->getFullNodeSize(a);
    max = compute(owner, f->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(owner, f->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = f->getSparseNodeSize(a);
    max = compute(owner, f->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(owner, f->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((long*)(ansEntry+1))[0] = max;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
  return max;
}


// ******************************************************************
// *                                                                *
// *                     int_min_range_op class                     *
// *                                                                *
// ******************************************************************

class int_min_range_op : public int_maxminrange_op {
  public:
    int_min_range_op()                    { }
    virtual ~int_min_range_op()           { }
    virtual const char* getName() const   { return "Min_range"; }

    static int_min_range_op* getInstance();

    virtual long compute(op_info* owner, int a);
};

int_min_range_op* 
int_min_range_op::
getInstance()
{
  static int_min_range_op instance;
  return &instance;
}

long int_min_range_op::compute(op_info* owner, int a)
{
  expert_forest* f = owner->p[0].getForest();

  // Terminal case
  if (f->isTerminalNode(a)) return f->getInteger(a);
  
  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
  }

  // recurse
  long min;
  if (f->isFullNode(a)) {
    // Full node
    int asize = f->getFullNodeSize(a);
    min = compute(owner, f->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(owner, f->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = f->getSparseNodeSize(a);
    min = compute(owner, f->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(owner, f->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((long*)(ansEntry+1))[0] = min;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
  return min;
}


// ******************************************************************
// *                                                                *
// *                   real_maxminrange_op  class                   *
// *                                                                *
// ******************************************************************

/** Abstract base class for max, min range operations with real range.
*/
class real_maxminrange_op : public maxminrange_op {
  public:
    real_maxminrange_op()                    { }
    virtual ~real_maxminrange_op()           { }

    virtual int getAnsLength() const  { 
      return sizeof(double) / sizeof(int); 
    }
    virtual int getCacheEntryLength() const { 
      return 1 + sizeof(double) / sizeof(int);
    }

    virtual int getAnsLengthInBytes() const { 
      return sizeof(double); 
    }
    virtual int getCacheEntryLengthInBytes() const { 
      return sizeof(int) + sizeof(double);
    }

    virtual compute_manager::error typeCheck(const op_info* owner) {
      return type_check(owner, op_param::REAL);
    }

    virtual void showEntry(const op_info* owner, FILE* strm,
        const int *entryData) const;

    virtual compute_manager::error compute(op_info* owner, const dd_edge& a,
      double& b);

    // for derived classes
    virtual double compute(op_info* owner, int a) = 0;
};

void
real_maxminrange_op::
showEntry(const op_info* owner, FILE* strm, const int *data) const
{
  DCASSERT(owner->nParams == 2);
  fprintf(strm, "[%s(%d): %le]",
      owner->op->getName(), data[0], ((const double*)(data+1))[0]
  );
}

compute_manager::error 
real_maxminrange_op::
compute(op_info* owner, const dd_edge& a, double& b)
{
  if (0==owner) return compute_manager::TYPE_MISMATCH;
  b = compute(owner, a.getNode());
  return compute_manager::SUCCESS;
}


// ******************************************************************
// *                                                                *
// *                    real_max_range_op  class                    *
// *                                                                *
// ******************************************************************

class real_max_range_op : public real_maxminrange_op {
  public:
    real_max_range_op()                   { }
    virtual ~real_max_range_op()          { }
    virtual const char* getName() const   { return "Max_range"; }

    static real_max_range_op* getInstance();

    virtual double compute(op_info* owner, int a);
};

real_max_range_op* 
real_max_range_op::
getInstance()
{
  static real_max_range_op instance;
  return &instance;
}

double real_max_range_op::compute(op_info* owner, int a)
{
  expert_forest* f = owner->p[0].getForest();

  // Terminal case
  if (f->isTerminalNode(a)) return f->getReal(a);
  
  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  double max;
  if (f->isFullNode(a)) {
    // Full node
    int asize = f->getFullNodeSize(a);
    max = compute(owner, f->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(owner, f->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = f->getSparseNodeSize(a);
    max = compute(owner, f->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(owner, f->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(double)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((double*)(ansEntry+1))[0] = max;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
  return max;
}


// ******************************************************************
// *                                                                *
// *                    real_min_range_op  class                    *
// *                                                                *
// ******************************************************************

class real_min_range_op : public real_maxminrange_op {
  public:
    real_min_range_op()                   { }
    virtual ~real_min_range_op()          { }
    virtual const char* getName() const   { return "Min_range"; }

    static real_min_range_op* getInstance();

    virtual double compute(op_info* owner, int a);
};

real_min_range_op* 
real_min_range_op::
getInstance()
{
  static real_min_range_op instance;
  return &instance;
}

double real_min_range_op::compute(op_info* owner, int a)
{
  expert_forest* f = owner->p[0].getForest();

  // Terminal case
  if (f->isTerminalNode(a)) return f->getReal(a);
  
  // Check compute table
  const int* cacheEntry = owner->cc->find(owner, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  double min;
  if (f->isFullNode(a)) {
    // Full node
    int asize = f->getFullNodeSize(a);
    min = compute(owner, f->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(owner, f->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = f->getSparseNodeSize(a);
    min = compute(owner, f->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(owner, f->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(double)/sizeof(int)];
  owner->p[0].getForest()->cacheNode(a);
  ansEntry[0] = a;
  ((double*)(ansEntry+1))[0] = min;

  owner->cc->add(owner, const_cast<const int*>(ansEntry));
  return min;
}



// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

operation* getMaxRangeOperation(const op_param &ft, const op_param &rt)
{
  switch (rt.getType()) {
    case op_param::INTEGER:
      if (!ft.isIntForest())    return 0;
      if (ft.isMT())            return int_max_range_op::getInstance();
      return 0;

    case op_param::REAL:
      if (!ft.isRealForest())   return 0;
      if (ft.isMT())            return real_max_range_op::getInstance();
      return 0;

    default:
      return 0;
  } // switch

  return 0;
}

operation* getMinRangeOperation(const op_param &ft, const op_param &rt)
{
  switch (rt.getType()) {
    case op_param::INTEGER:
      if (!ft.isIntForest())    return 0;
      if (ft.isMT())            return int_min_range_op::getInstance();
      return 0;

    case op_param::REAL:
      if (!ft.isRealForest())   return 0;
      if (ft.isMT())            return real_min_range_op::getInstance();
      return 0;

    default:
      return 0;
  } // switch

  return 0;
}

