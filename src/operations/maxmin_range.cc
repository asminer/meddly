
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "maxmin_range.h"
#include "../compute_table.h"

namespace MEDDLY {

  class range_int;
  class range_real;

  class maxrange_int;
  class minrange_int;

  class maxrange_real;
  class minrange_real;

  class maxrange_opname;
  class minrange_opname;
};

// ******************************************************************
// *                                                                *
// *                        range_int  class                        *
// *                                                                *
// ******************************************************************

/// Abstract base class: max or min range that returns an integer.
class MEDDLY::range_int : public unary_operation {
  public:
    range_int(const unary_opname* oc, expert_forest* arg);

    // common
    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
};

MEDDLY::range_int::range_int(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, arg, INTEGER)
{
  key_length = 1;
  ans_length = sizeof(long) / sizeof(int); 
}

bool MEDDLY::range_int::isEntryStale(const int* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::range_int::discardEntry(const int* data)
{
  argF->uncacheNode(data[0]);
}

void MEDDLY::range_int::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d): %ld(L)]",
      getName(), data[0], ((const long*)(data+1))[0]
  );
}


// ******************************************************************
// *                                                                *
// *                        range_real class                        *
// *                                                                *
// ******************************************************************

/// Abstract base class: max or min range that returns a real.
class MEDDLY::range_real : public unary_operation {
  public:
    range_real(const unary_opname* oc, expert_forest* arg);

    // common
    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
};

MEDDLY::range_real::range_real(const unary_opname* oc, expert_forest* arg)
 : unary_operation(oc, arg, INTEGER)
{
  key_length = 1;
  ans_length = sizeof(double) / sizeof(int); 
}

bool MEDDLY::range_real::isEntryStale(const int* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::range_real::discardEntry(const int* data)
{
  argF->uncacheNode(data[0]);
}

void MEDDLY::range_real::showEntry(FILE* strm, const int *data) const
{
  fprintf(strm, "[%s(%d): %le]",
      getName(), data[0], ((const double*)(data+1))[0]
  );
}


// ******************************************************************
// *                                                                *
// *                       maxrange_int class                       *
// *                                                                *
// ******************************************************************

/// Max range, returns an integer
class MEDDLY::maxrange_int : public range_int {
public:
  maxrange_int(const unary_opname* oc, expert_forest* arg)
    : range_int(oc, arg) { }
  virtual void compute(const dd_edge &arg, long &res) {
    res = compute(arg.getNode());
  }
  long compute(int a);
};

long MEDDLY::maxrange_int::compute(int a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) return argF->getInteger(a);
  
  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
  }

  // recurse
  long max;
  if (argF->isFullNode(a)) {
    // Full node
    int asize = argF->getFullNodeSize(a);
    max = compute(argF->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(argF->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = argF->getSparseNodeSize(a);
    max = compute(argF->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(argF->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((long*)(ansEntry+1))[0] = max;

  CT->add(this, const_cast<const int*>(ansEntry));
  return max;
}



// ******************************************************************
// *                                                                *
// *                       minrange_int class                       *
// *                                                                *
// ******************************************************************

/// Min range, returns an integer
class MEDDLY::minrange_int : public range_int {
public:
  minrange_int(const unary_opname* oc, expert_forest* arg)
    : range_int(oc, arg) { }
  virtual void compute(const dd_edge &arg, long &res) {
    res = compute(arg.getNode());
  }
  long compute(int a);
};

long MEDDLY::minrange_int::compute(int a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) return argF->getInteger(a);
  
  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const long*)(cacheEntry+1))[0];
  }

  // recurse
  long min;
  if (argF->isFullNode(a)) {
    // Full node
    int asize = argF->getFullNodeSize(a);
    min = compute(argF->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(argF->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = argF->getSparseNodeSize(a);
    min = compute(argF->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(argF->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(long)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((long*)(ansEntry+1))[0] = min;

  CT->add(this, const_cast<const int*>(ansEntry));
  return min;
}



// ******************************************************************
// *                                                                *
// *                      maxrange_real  class                      *
// *                                                                *
// ******************************************************************

/// Max range, returns a real
class MEDDLY::maxrange_real : public range_real {
public:
  maxrange_real(const unary_opname* oc, expert_forest* arg)
    : range_real(oc, arg) { }
  virtual void compute(const dd_edge &arg, double &res) {
    res = compute(arg.getNode());
  }
  double compute(int a);
};

double MEDDLY::maxrange_real::compute(int a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) return argF->getReal(a);
  
  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  double max;
  if (argF->isFullNode(a)) {
    // Full node
    int asize = argF->getFullNodeSize(a);
    max = compute(argF->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(argF->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = argF->getSparseNodeSize(a);
    max = compute(argF->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      max = MAX(max, compute(argF->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(double)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((double*)(ansEntry+1))[0] = max;

  CT->add(this, const_cast<const int*>(ansEntry));
  return max;
}



// ******************************************************************
// *                                                                *
// *                      minrange_real  class                      *
// *                                                                *
// ******************************************************************

/// Min range, returns a real
class MEDDLY::minrange_real : public range_real {
public:
  minrange_real(const unary_opname* oc, expert_forest* arg)
    : range_real(oc, arg) { }
  virtual void compute(const dd_edge &arg, double &res) {
    res = compute(arg.getNode());
  }
  double compute(int a);
};

double MEDDLY::minrange_real::compute(int a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) return argF->getReal(a);
  
  // Check compute table
  const int* cacheEntry = CT->find(this, &a);
  if (cacheEntry) {
    return ((const double*)(cacheEntry+1))[0];
  }

  // recurse
  double min;
  if (argF->isFullNode(a)) {
    // Full node
    int asize = argF->getFullNodeSize(a);
    min = compute(argF->getFullNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(argF->getFullNodeDownPtr(a, i)));
    } // for i
  } else {
    // Sparse node
    int asize = argF->getSparseNodeSize(a);
    min = compute(argF->getSparseNodeDownPtr(a, 0));
    for (int i = 1; i < asize; ++i) {
      min = MIN(min, compute(argF->getSparseNodeDownPtr(a, i)));
    } // for i
  } 

  // Add entry to compute table
  static int ansEntry[1+sizeof(double)/sizeof(int)];
  ansEntry[0] = argF->cacheNode(a);
  ((double*)(ansEntry+1))[0] = min;

  CT->add(this, const_cast<const int*>(ansEntry));
  return min;
}



// ******************************************************************
// *                                                                *
// *                     maxrange_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::maxrange_opname : public unary_opname {
  public:
    maxrange_opname();
    virtual unary_operation*
      buildOperation(const forest* ar, opnd_type res) const;
};

MEDDLY::maxrange_opname::maxrange_opname() : unary_opname("Max_range")
{
}

MEDDLY::unary_operation*
MEDDLY::maxrange_opname::buildOperation(const forest* ar, opnd_type res) const
{
  if (0==ar) return 0;

  if (ar->getEdgeLabeling() != forest::MULTI_TERMINAL)
    throw error(error::NOT_IMPLEMENTED);

  switch (res) {
    case INTEGER:
      if (forest::INTEGER != ar->getRangeType())
        throw error(error::TYPE_MISMATCH);
      return new maxrange_int(this, (expert_forest*) ar);

    case REAL:
      if (forest::REAL != ar->getRangeType())
        throw error(error::TYPE_MISMATCH);
      return new maxrange_real(this, (expert_forest*) ar);

    default:
      throw error(error::TYPE_MISMATCH);
  } // switch

  throw error(error::MISCELLANEOUS);
}

// ******************************************************************
// *                                                                *
// *                     minrange_opname  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::minrange_opname : public unary_opname {
  public:
    minrange_opname();
    virtual unary_operation*
      buildOperation(const forest* ar, opnd_type res) const;
};

MEDDLY::minrange_opname::minrange_opname() : unary_opname("Min_range")
{
}

MEDDLY::unary_operation*
MEDDLY::minrange_opname::buildOperation(const forest* ar, opnd_type res) const
{
  if (0==ar) return 0;

  if (ar->getEdgeLabeling() != forest::MULTI_TERMINAL)
    throw error(error::NOT_IMPLEMENTED);

  switch (res) {
    case INTEGER:
      if (forest::INTEGER != ar->getRangeType())
        throw error(error::TYPE_MISMATCH);
      return new minrange_int(this, (expert_forest*) ar);

    case REAL:
      if (forest::REAL != ar->getRangeType())
        throw error(error::TYPE_MISMATCH);
      return new minrange_real(this, (expert_forest*) ar);

    default:
      throw error(error::TYPE_MISMATCH);
  } // switch

  throw error(error::MISCELLANEOUS);
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeMaxRange(const settings &s)
{
  return new maxrange_opname;
}

MEDDLY::unary_opname* MEDDLY::initializeMinRange(const settings &s)
{
  return new minrange_opname;
}

