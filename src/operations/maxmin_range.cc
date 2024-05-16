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

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"

namespace MEDDLY {
    class range_int;
    class range_real;

    class maxrange_int;
    class minrange_int;

    class maxrange_real;
    class minrange_real;

    unary_list MAXRANGE_cache;
    unary_list MINRANGE_cache;
};

// ******************************************************************
// *                                                                *
// *                        range_int  class                        *
// *                                                                *
// ******************************************************************

/// Abstract base class: max or min range that returns an integer.
class MEDDLY::range_int : public unary_operation {
  public:
    range_int(unary_list& oc, forest* arg);

  protected:
    inline ct_entry_key*
    findResult(node_handle a, int &b)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = CTresult[0].readI();
      CT0->recycle(CTsrch);
      return 0;
    }
    inline long saveResult(ct_entry_key* Key,
      node_handle a, int &b)
    {
      CTresult[0].reset();
      CTresult[0].writeI(b);
      CT0->addEntry(Key, CTresult[0]);
      return b;
    }
};

MEDDLY::range_int::range_int(unary_list& oc, forest* arg)
 : unary_operation(oc, 1, arg, opnd_type::INTEGER)
{
  ct_entry_type* et = new ct_entry_type(oc.getName(), "N:I");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
}

// ******************************************************************
// *                                                                *
// *                        range_real class                        *
// *                                                                *
// ******************************************************************

/// Abstract base class: max or min range that returns a real.
class MEDDLY::range_real : public unary_operation {
  public:
    range_real(unary_list& oc, forest* arg);

  protected:
    inline ct_entry_key* findResult(node_handle a, float &b) {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = CTresult[0].readF();
      CT0->recycle(CTsrch);
      return 0;
    }
    inline float saveResult(ct_entry_key* Key,
      node_handle a, float &b)
    {
      CTresult[0].reset();
      CTresult[0].writeF(b);
      CT0->addEntry(Key, CTresult[0]);
      return b;
    }
};

MEDDLY::range_real::range_real(unary_list& oc, forest* arg)
 : unary_operation(oc, 1, arg, opnd_type::REAL)
{
  ct_entry_type* et = new ct_entry_type(oc.getName(), "N:F");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
}

// ******************************************************************
// *                                                                *
// *                       maxrange_int class                       *
// *                                                                *
// ******************************************************************

/// Max range, returns an integer
class MEDDLY::maxrange_int : public range_int {
    public:
        maxrange_int(forest* arg) : range_int(MAXRANGE_cache, arg) { }
        virtual void compute(const dd_edge &arg, long &res) {
            res = compute_r(arg.getNode());
        }
        int compute_r(node_handle a);
};

int MEDDLY::maxrange_int::compute_r(node_handle a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) {
      terminal t;
      t.setFromHandle(terminal_type::INTEGER, a);
      return t.getInteger();
      // return int_Tencoder::handle2value(a);
  }

  // Check compute table
  int max;
  ct_entry_key* Key = findResult(a, max);
  if (0==Key) return max;

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // recurse
  max = compute_r(A->down(0));
  for (unsigned i=1; i<A->getSize(); i++) {
    max = MAX(max, compute_r(A->down(i)));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  return saveResult(Key, a, max);
}

// ******************************************************************
// *                                                                *
// *                       minrange_int class                       *
// *                                                                *
// ******************************************************************

/// Min range, returns an integer
class MEDDLY::minrange_int : public range_int {
    public:
        minrange_int(forest* arg) : range_int(MINRANGE_cache, arg) { }
        virtual void compute(const dd_edge &arg, long &res) {
            res = compute_r(arg.getNode());
        }
        int compute_r(node_handle a);
};

int MEDDLY::minrange_int::compute_r(node_handle a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) {
      terminal t;
      t.setFromHandle(terminal_type::INTEGER, a);
      return t.getInteger();
      // return int_Tencoder::handle2value(a);
  }

  // Check compute table
  int min;
  ct_entry_key* Key = findResult(a, min);
  if (0==Key) return min;

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // recurse
  min = compute_r(A->down(0));
  for (unsigned i=1; i<A->getSize(); i++) {
    min = MIN(min, compute_r(A->down(i)));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  return saveResult(Key, a, min);
}



// ******************************************************************
// *                                                                *
// *                      maxrange_real  class                      *
// *                                                                *
// ******************************************************************

/// Max range, returns a real
class MEDDLY::maxrange_real : public range_real {
    public:
        maxrange_real(forest* arg) : range_real(MAXRANGE_cache, arg) { }
        virtual void compute(const dd_edge &arg, double &res) {
            res = compute_r(arg.getNode());
        }
        float compute_r(node_handle a);
};

float MEDDLY::maxrange_real::compute_r(node_handle a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) {
      terminal t;
      t.setFromHandle(terminal_type::REAL, a);
      return t.getReal();
      // return float_Tencoder::handle2value(a);
  }

  // Check compute table
  float max;
  ct_entry_key* Key = findResult(a, max);
  if (0==Key) return max;

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // recurse
  max = compute_r(A->down(0));
  for (unsigned i=1; i<A->getSize(); i++) {
    max = MAX(max, compute_r(A->down(i)));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  return saveResult(Key, a, max);
}



// ******************************************************************
// *                                                                *
// *                      minrange_real  class                      *
// *                                                                *
// ******************************************************************

/// Min range, returns a real
class MEDDLY::minrange_real : public range_real {
    public:
        minrange_real(forest* arg) : range_real(MINRANGE_cache, arg) { }
        virtual void compute(const dd_edge &arg, double &res) {
            res = compute_r(arg.getNode());
        }
        float compute_r(node_handle a);
};

float MEDDLY::minrange_real::compute_r(node_handle a)
{
  // Terminal case
  if (argF->isTerminalNode(a)) {
      terminal t;
      t.setFromHandle(terminal_type::REAL, a);
      return t.getReal();
      // return float_Tencoder::handle2value(a);
  }

  // Check compute table
  float min;
  ct_entry_key* Key = findResult(a, min);
  if (0==Key) return min;

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // recurse
  min = compute_r(A->down(0));
  for (unsigned i=1; i<A->getSize(); i++) {
    min = MIN(min, compute_r(A->down(i)));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Add entry to compute table
  return saveResult(Key, a, min);
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_operation* MEDDLY::MAX_RANGE(forest* arg, opnd_type res)
{
    if (!arg) return nullptr;
    unary_operation* uop = MAXRANGE_cache.find(arg, res);
    if (uop) {
        return uop;
    }

    if (arg->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

    switch (res) {
        case opnd_type::INTEGER:
            if (range_type::INTEGER != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return MAXRANGE_cache.add(new maxrange_int(arg));

        case opnd_type::REAL:
            if (range_type::REAL != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return MAXRANGE_cache.add(new maxrange_real(arg));

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    } // switch

    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
}

void MEDDLY::MAX_RANGE_init()
{
    MAXRANGE_cache.reset("Max_range");
}

void MEDDLY::MAX_RANGE_done()
{
    MEDDLY_DCASSERT(MAXRANGE_cache.isEmpty());
}

// ******************************************************************

MEDDLY::unary_operation* MEDDLY::MIN_RANGE(forest* arg, opnd_type res)
{
    if (!arg) return nullptr;
    unary_operation* uop = MINRANGE_cache.find(arg, res);
    if (uop) {
        return uop;
    }

    if (arg->getEdgeLabeling() != edge_labeling::MULTI_TERMINAL)
        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

    switch (res) {
        case opnd_type::INTEGER:
            if (range_type::INTEGER != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return MINRANGE_cache.add(new minrange_int(arg));

        case opnd_type::REAL:
            if (range_type::REAL != arg->getRangeType())
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            return MINRANGE_cache.add(new minrange_real(arg));

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    } // switch

    throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
}

void MEDDLY::MIN_RANGE_init()
{
    MINRANGE_cache.reset("Min_range");
}

void MEDDLY::MIN_RANGE_done()
{
    MEDDLY_DCASSERT(MINRANGE_cache.isEmpty());
}

