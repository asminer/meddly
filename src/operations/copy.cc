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
#include "copy.h"

#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_unary.h"

// #define DEBUG_COPY_COMPUTE_ALL

namespace MEDDLY {
    class copy_MT;
};

// ******************************************************************
// *                                                                *
// *                         copy_MT  class                         *
// *                                                                *
// ******************************************************************

/// Abstract base class for copying between multi-terminal DDs.
class MEDDLY::copy_MT : public unary_operation {
    protected:
        copy_MT(unary_list &cache, forest* arg, forest* res);
        virtual node_handle compute_r(node_handle a) = 0;

        inline ct_entry_key*
        findResult(node_handle a, node_handle &b)
        {
            ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
            MEDDLY_DCASSERT(CTsrch);
            CTsrch->writeN(a);
            CT0->find(CTsrch, CTresult[0]);
            if (!CTresult[0]) return CTsrch;
            b = resF->linkNode(CTresult[0].readN());
            CT0->recycle(CTsrch);
            return 0;
        }
        inline node_handle saveResult(ct_entry_key* Key,
            node_handle a, node_handle b)
        {
            CTresult[0].reset();
            CTresult[0].writeN(b);
            CT0->addEntry(Key, CTresult[0]);
            return b;
        }
    public:
        virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);


};


MEDDLY::copy_MT
:: copy_MT(unary_list &cache, forest* arg, forest* res)
 : unary_operation(cache, 1, arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct_entry_type* et = new ct_entry_type(cache.getName(), "N:N");
    et->setForestForSlot(0, arg);
    et->setForestForSlot(2, res);
    registerEntryType(0, et);
    buildCTs();
}

void MEDDLY::copy_MT::computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
{
  node_handle result = compute_r(arg.getNode());
  res.set(result);
}

// ******************************************************************
// *                                                                *
// *                       copy_MT_tmpl class                       *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

    template <typename RESULT>
    class copy_MT_tmpl : public copy_MT {
        protected:
            copy_MT_tmpl(forest* A, forest* R)
                : copy_MT(cache, A, R) { }

            virtual node_handle compute_r(node_handle a)
            {
                if (argF->getReductionRule() == resF->getReductionRule()) {
                    return computeSkip(-1, a);  // same skipping rule, ok
                } else {
                    // need to visit every level...
                    return computeAll(-1, resF->getMaxLevelIndex(), a);
                }
            }

            node_handle computeSkip(int in, node_handle a);
            node_handle computeAll(int in, int k, node_handle a);

        public:
            static unary_list cache;

            inline static unary_operation* build(forest* a, forest* r)
            {
                unary_operation* uop =  cache.findOperation(a, r);
                if (uop) {
                    return uop;
                }
                return cache.addOperation(new copy_MT_tmpl<RESULT>(a, r));
            }
    };

};  // namespace MEDDLY

template <typename RESULT>
MEDDLY::unary_list MEDDLY::copy_MT_tmpl<RESULT>::cache;

template <typename RESULT>
MEDDLY::node_handle MEDDLY::copy_MT_tmpl<RESULT>::computeSkip(int in, node_handle a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    RESULT aTerm;
    argF->getValueFromHandle(a, aTerm);
    return resF->handleForValue(aTerm);
  }

  // Check compute table
  node_handle b;
  ct_entry_key* Key = findResult(a, b);
  if (0==Key) return b;

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
      nb->setSparse(z, A->index(z), computeSkip(int(A->index(z)), A->down(z)));
  }

  // Handle extensible edge, if any
  if (A->isExtensible()) nb->markAsExtensible();

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  b = resF->createReducedNode(in, nb);

  // Add to compute table
  return saveResult(Key, a, b);
}


template <typename RESULT>
MEDDLY::node_handle MEDDLY::copy_MT_tmpl<RESULT>::computeAll(int in, int k, node_handle a)
{
  // Check terminals
  if (0==k) {
    RESULT aTerm;
    argF->getValueFromHandle(a, aTerm);
    return resF->handleForValue(aTerm);
  }

  // 0 is 0 is 0, I think...
  if (0==a) {
    return 0;
  }

#ifdef DEBUG_COPY_COMPUTE_ALL
  fprintf(stderr, "copy(%d, %d, %d)\n", in, k, a);
#endif

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  node_handle b;
  ct_entry_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = findResult(a, b);
    if (0==Key) return b;
  }
  int nextk;
  if (resF->isForRelations()) {
    nextk = resF->downLevel(k);
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::New(argF);
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A->initIdentity(argF, k, unsigned(in), a, SPARSE_ONLY);
    } else {
      A->initRedundant(argF, k, a, SPARSE_ONLY);
    }
  } else {
    argF->unpackNode(A, a, SPARSE_ONLY);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
      nb->setSparse(z, A->index(z),
              computeAll(int(A->index(z)), nextk, A->down(z)));
  }

  // Handle extensible edge, if any
  if (A->isExtensible()) nb->markAsExtensible();

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  b = resF->createReducedNode(in, nb);

  // Add to compute table
  if (Key) saveResult(Key, a, b);
  return b;
}



// ******************************************************************
// *                                                                *
// *                        copy_MT2EV class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  template <typename TYPE>
  class copy_MT2EV : public unary_operation {
      protected:
        copy_MT2EV(forest* arg, forest* res, const char* pattern)
            : unary_operation(cache, 1, arg, res)
        {
            //
            // Pattern should be of the form "N:xN" where x is the EV Type.
            //
            // entry[0]: mt node
            // entry[1]: EV value (output)
            // entry[2]: EV node (output)
            //
            ct_entry_type* et = new ct_entry_type(cache.getName(), pattern);
            et->setForestForSlot(0, arg);
            et->setForestForSlot(3, res);
            registerEntryType(0, et);
            buildCTs();
        }

        virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag)
        {
            node_handle b;
            TYPE bev = Inf<TYPE>();
            if (argF->getReductionRule() == resF->getReductionRule()) {
                computeSkip(-1, arg.getNode(), b, bev);  // same skipping rule, ok
            } else {
                // need to visit every level...
                computeAll(-1, resF->getMaxLevelIndex(), arg.getNode(), b, bev);
            }
            res.set(b, bev);
        }
        void computeSkip(int in, node_handle a, node_handle &b, TYPE &bev);
        void computeAll(int in, int k, node_handle a, node_handle &b, TYPE &bev);

    protected:
        inline ct_entry_key*
        inCache(node_handle a, node_handle &b, TYPE &bev)
        {
            ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
            MEDDLY_DCASSERT(CTsrch);
            CTsrch->writeN(a);
            CT0->find(CTsrch, CTresult[0]);
            if (CTresult[0]) {
                CTresult[0].read_ev(bev);
                b = resF->linkNode(CTresult[0].readN());
                CT0->recycle(CTsrch);
                return 0;
            }
            return CTsrch;
        }

        inline void addToCache(ct_entry_key* Key,
            node_handle a, node_handle b, long bev)
        {
            MEDDLY_DCASSERT(bev != Inf<long>());
            CTresult[0].reset();
            CTresult[0].writeL(bev);
            CTresult[0].writeN(b);
            CT0->addEntry(Key, CTresult[0]);
        }

        inline void addToCache(ct_entry_key* Key,
            node_handle a, node_handle b, float bev)
        {
            MEDDLY_DCASSERT(bev != Inf<float>());
            CTresult[0].reset();
            CTresult[0].writeF(bev);
            CTresult[0].writeN(b);
            CT0->addEntry(Key, CTresult[0]);
        }

        public:
            static unary_list cache;

            inline static unary_operation* build(forest* a, forest* r,
                    const char* p)
            {
                unary_operation* uop =  cache.findOperation(a, r);
                if (uop) {
                    return uop;
                }
                return cache.addOperation(new copy_MT2EV<TYPE>(a, r, p));
            }
  };

};  // namespace MEDDLY

template <typename TYPE>
MEDDLY::unary_list MEDDLY::copy_MT2EV<TYPE>::cache;


template <typename TYPE>
void MEDDLY::copy_MT2EV<TYPE>
::computeSkip(int in, node_handle a, node_handle &b, TYPE &bev)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    MEDDLY_DCASSERT(a != argF->getTransparentNode());
    if (argF->getRangeType() == range_type::BOOLEAN) {
      bev = 0;
    }
    else {
      argF->getValueFromHandle(a, bev);
    }
    terminal t(true);
    b = t.getHandle();
    return;
  }

  // Check compute table
  ct_entry_key* Key = inCache(a, b, bev);
  if (0==Key) return;

  // Initialize sparse node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
    node_handle d;
    TYPE dev = Inf<TYPE>();
    computeSkip(int(A->index(z)), A->down(z), d, dev);
    nb->setSparse(z, A->index(z), edge_value(dev), d);
    // nb->i_ref(z) = A->index(z);
    // nb->d_ref(z) = d;
    // nb->setEdge(z, dev);
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  resF->createReducedNode(in, nb, bev, b);

  // Add to compute table
  addToCache(Key, a, b, bev);
}

template <typename TYPE>
void MEDDLY::copy_MT2EV<TYPE>
::computeAll(int in, int k, node_handle a, node_handle &b, TYPE &bev)
{
  // Check terminals
  if (0==k) {
    MEDDLY_DCASSERT(a != argF->getTransparentNode());
    if (argF->getRangeType() == range_type::BOOLEAN) {
      bev = 0;
    }
    else {
      argF->getValueFromHandle(a, bev);
    }
    terminal t(true);
    b = t.getHandle();
    return;
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  ct_entry_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = inCache(a, b, bev);
    if (0==Key) return;
  }

  // What's below?
  int nextk;
  if (resF->isForRelations()) {
    nextk = resF->downLevel(k);
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::New(argF);
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A->initIdentity(argF, k, unsigned(in), a, SPARSE_ONLY);
    } else {
      A->initRedundant(argF, k, a, SPARSE_ONLY);
    }
  } else {
    argF->unpackNode(A, a, SPARSE_ONLY);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
    TYPE dev;
    node_handle dn;
    computeAll(int(A->index(z)), nextk, A->down(z), dn, dev);
    nb->setSparse(z, A->index(z), edge_value(dev), dn);

    // nb->i_ref(z) = A->index(z);
    // computeAll(int(A->index(z)), nextk, A->down(z), nb->d_ref(z), dev);
    // nb->setEdge(z, dev);
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  resF->createReducedNode(in, nb, bev, b);

  // Add to compute table
  if (Key) addToCache(Key, a, b, bev);
}

// ******************************************************************
// *                                                                *
// *                        copy_EV2MT class                        *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  template <typename TYPE, class OP>
  class copy_EV2MT : public unary_operation {
    protected:
      copy_EV2MT(forest* arg, forest* res, const char* pattern)
          : unary_operation(cache, 1, arg, res)
      {
        //
        // Pattern should be of the form "xN:N"
        //
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: mt node (output)
        //
        ct_entry_type* et = new ct_entry_type(cache.getName(), pattern);
        et->setForestForSlot(1, arg);
        et->setForestForSlot(3, res);
        registerEntryType(0, et);
        buildCTs();
      }

      virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag) {
        TYPE ev;
        node_handle b;
        arg.getEdgeValue(ev);
        if (argF->getReductionRule() == resF->getReductionRule()) {
          b = computeSkip(-1, ev, arg.getNode());  // same skipping rule, ok
        } else {
          // need to visit every level...
          b = computeAll(-1, resF->getMaxLevelIndex(), ev, arg.getNode());
        }
        res.set(b);
      }
      node_handle computeSkip(int in, TYPE ev, node_handle a);
      node_handle computeAll(int in, int k, TYPE ev, node_handle a);

      inline ct_entry_key*
      inCache(TYPE ev, node_handle a, node_handle &b)
      {
        ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->write_ev(ev);
        CTsrch->writeN(a);
        CT0->find(CTsrch, CTresult[0]);
        if (CTresult[0]) {
          b = resF->linkNode(CTresult[0].readN());
          CT0->recycle(CTsrch);
          return 0;
        }
        return CTsrch;
      }

      inline void addToCache(ct_entry_key* Key,
        TYPE ev, node_handle a, node_handle b)
      {
        CTresult[0].reset();
        CTresult[0].writeN(b);
        CT0->addEntry(Key, CTresult[0]);
      }

        public:
            static unary_list cache;

            inline static unary_operation* build(forest* a, forest* r,
                    const char* p)
            {
                unary_operation* uop =  cache.findOperation(a, r);
                if (uop) {
                    return uop;
                }
                return cache.addOperation(new copy_EV2MT<TYPE,OP>(a, r, p));
            }

  };

};  // namespace MEDDLY

template <typename TYPE, class OP>
MEDDLY::unary_list MEDDLY::copy_EV2MT<TYPE,OP>::cache;

template <typename TYPE, class OP>
MEDDLY::node_handle  MEDDLY::copy_EV2MT<TYPE,OP>
::computeSkip(int in, TYPE ev, node_handle a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    MEDDLY_DCASSERT(a != argF->getTransparentNode());
    if (resF->getRangeType() == range_type::BOOLEAN) {
      terminal t(true);
      return t.getHandle();
    }
    else {
      return resF->handleForValue(ev);
    }
  }

  // Check compute table
  node_handle b;
  ct_entry_key* Key = inCache(ev, a, b);
  if (0==Key) return b;

  // Initialize sparse node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
    TYPE aev;
    A->edgeval(z).get(aev);
    nb->setSparse(z, A->index(z),
            computeSkip(int(A->index(z)), OP::apply(ev, aev), A->down(z))
    );
    // A->getEdge(z, aev);
    // nb->i_ref(z) = A->index(z);
    // nb->d_ref(z) = computeSkip(int(A->index(z)), OP::apply(ev, aev), A->down(z));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  b = resF->createReducedNode(in, nb);

  // Add to compute table
  addToCache(Key, ev, a, b);
  return b;
}

template <typename TYPE, class OP>
MEDDLY::node_handle  MEDDLY::copy_EV2MT<TYPE,OP>
::computeAll(int in, int k, TYPE ev, node_handle a)
{
  // Check terminals
  if (0==k) {
    MEDDLY_DCASSERT(a != argF->getTransparentNode());
    if (resF->getRangeType() == range_type::BOOLEAN) {
      terminal t(true);
      return t.getHandle();
    }
    else {
      return resF->handleForValue(ev);
    }
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  node_handle b;
  ct_entry_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = inCache(ev, a, b);
    if (0==Key) return b;
  }

  // What's below?
  int nextk;
  if (resF->isForRelations()) {
    nextk = resF->downLevel(k);
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::New(argF);
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      TYPE rev;
      OP::redundant(rev);
      A->initIdentity(argF, k, unsigned(in), rev, a, SPARSE_ONLY);
    } else {
      TYPE rev;
      OP::redundant(rev);
      A->initRedundant(argF, k, rev, a, SPARSE_ONLY);
    }
  } else {
    argF->unpackNode(A, a, SPARSE_ONLY);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
    TYPE aev;
    A->edgeval(z).get(aev);
    nb->setSparse(z, A->index(z),
        computeAll(int(A->index(z)), nextk, OP::apply(ev, aev), A->down(z))
    );
    // A->getEdge(z, aev);
    // nb->i_ref(z) = A->index(z);
    // nb->d_ref(z) = computeAll(int(A->index(z)), nextk, OP::apply(ev, aev), A->down(z));
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  b = resF->createReducedNode(in, nb);

  // Add to compute table
  if (Key) addToCache(Key, ev, a, b);
  return b;
}

// ******************************************************************
// *                                                                *
// *                     copy_EV2EV_fast  class                     *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  // 1-1 mapping between input edges and output edges
  template <typename INTYPE, typename OUTTYPE>
  class copy_EV2EV_fast : public unary_operation {
    protected:
      copy_EV2EV_fast(forest* arg, forest* res)
        : unary_operation(cache, 1, arg, res)
      {
        // entry[0]: EV node
        // entry[1]: EV node
        ct_entry_type* et = new ct_entry_type(cache.getName(), "N:N");
        et->setForestForSlot(0, arg);
        et->setForestForSlot(2, res);
        registerEntryType(0, et);
        buildCTs();
      }

      virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag) {
        INTYPE av;
        node_handle bn;
        arg.getEdgeValue(av);
        bn = computeSkip(-1, arg.getNode());
        OUTTYPE bv = av;
        res.set(bn, bv);
      }

      node_handle computeSkip(int in, node_handle a);

      inline ct_entry_key*
      findResult(node_handle a, node_handle &b)
      {
        ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->writeN(a);
        CT0->find(CTsrch, CTresult[0]);
        if (!CTresult[0]) return CTsrch;
        b = resF->linkNode(CTresult[0].readN());
        CT0->recycle(CTsrch);
        return 0;
      }
      inline node_handle saveResult(ct_entry_key* Key,
        node_handle a, node_handle b)
      {
        CTresult[0].reset();
        CTresult[0].writeN(b);
        CT0->addEntry(Key, CTresult[0]);
        return b;
      }

        public:
            static unary_list cache;

            inline static unary_operation* build(forest* a, forest* r)
            {
                unary_operation* uop =  cache.findOperation(a, r);
                if (uop) {
                    return uop;
                }
                return cache.addOperation(
                        new copy_EV2EV_fast<INTYPE,OUTTYPE>(a, r)
                );
            }

  };

};  // namespace MEDDLY

template <typename INTYPE, typename OUTTYPE>
MEDDLY::unary_list MEDDLY::copy_EV2EV_fast<INTYPE,OUTTYPE>::cache;


template <typename INTYPE, typename OUTTYPE>
MEDDLY::node_handle
MEDDLY::copy_EV2EV_fast<INTYPE,OUTTYPE>::computeSkip(int in, node_handle a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
      terminal t(a != 0);
      return t.getHandle();
  }

  // Check compute table
  node_handle b;
  ct_entry_key* Key = findResult(a, b);
  if (0==Key) return b;

  // Initialize node reader
  unpacked_node* A = argF->newUnpacked(a, SPARSE_ONLY);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
    INTYPE av;
    A->edgeval(z).get(av);
    OUTTYPE bv = av;
    nb->setSparse(z, A->index(z), edge_value(bv),
        computeSkip(int(A->index(z)), A->down(z))
    );

    //nb->i_ref(z) = A->index(z);
    //nb->d_ref(z) = computeSkip(int(A->index(z)), A->down(z));
    //nb->setEdge(z, bv);
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  OUTTYPE bv;
  resF->createReducedNode(in, nb, bv, b);
  // bv should be the redundant/identity value

  // Add to compute table
  return saveResult(Key, a, b);
}


// ******************************************************************
// *                                                                *
// *                     copy_EV2EV_slow  class                     *
// *                                                                *
// ******************************************************************

namespace MEDDLY {

  template <typename INTYPE, class INOP, typename OUTTYPE>
  class copy_EV2EV_slow : public unary_operation {
    protected:
      copy_EV2EV_slow(forest* arg, forest* res, const char* pattern)
        : unary_operation(cache, 1, arg, res)
      {
        //
        // Pattern is of the form "xN:yN"
        //
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: EV value
        // entry[3]: EV node
        //
        ct_entry_type* et = new ct_entry_type(cache.getName(), pattern);
        et->setForestForSlot(1, arg);
        et->setForestForSlot(4, res);
        registerEntryType(0, et);
        buildCTs();
      }

      virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag) {
        INTYPE av;
        node_handle an, bn;
        OUTTYPE bv;
        an = arg.getNode();

        if (argF->isTransparentEdge(an, arg.getEdgeValue())
                && argF->getTransparentNode() == resF->getTransparentNode())
        {
            edge_value bev;
            resF->getTransparentEdge(bn, bev);
            bev.get(bv);
        }
        else {
          // need to visit every level...
          arg.getEdgeValue(av);
          computeAll(-1, resF->getMaxLevelIndex(), av, an, bv, bn);
        }
        res.set(bn, bv);
      }
      void computeAll(int in, int k, INTYPE av, node_handle an,
        OUTTYPE &bv, node_handle &bn);

      inline ct_entry_key*
      inCache(INTYPE av, node_handle an, OUTTYPE &bv, node_handle &bn)
      {
        ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->write_ev(av);
        CTsrch->writeN(an);
        CT0->find(CTsrch, CTresult[0]);
        if (CTresult[0]) {
          CTresult[0].read_ev(bv);
          bn = resF->linkNode(CTresult[0].readN());
          CT0->recycle(CTsrch);
          return 0;
        }
        return CTsrch;
      }

      inline void addToCache(ct_entry_key* Key,
        INTYPE av, node_handle an, long bv, node_handle &bn)
      {
        CTresult[0].reset();
        CTresult[0].writeL(bv);
        CTresult[0].writeN(bn);
        CT0->addEntry(Key, CTresult[0]);
      }

      inline void addToCache(ct_entry_key* Key,
        INTYPE av, node_handle a, float bv, node_handle &bn)
      {
        CTresult[0].reset();
        CTresult[0].writeF(bv);
        CTresult[0].writeN(bn);
        CT0->addEntry(Key, CTresult[0]);
      }

        public:
            static unary_list cache;

            inline static unary_operation* build(forest* a, forest* r,
                    const char* p)
            {
                unary_operation* uop =  cache.findOperation(a, r);
                if (uop) {
                    return uop;
                }
                return cache.addOperation(
                        new copy_EV2EV_slow<INTYPE,INOP,OUTTYPE>(a, r, p)
                );
            }

  };

};  // namespace MEDDLY

template <typename INTYPE, class INOP, typename OUTTYPE>
MEDDLY::unary_list MEDDLY::copy_EV2EV_slow<INTYPE,INOP,OUTTYPE>::cache;


template <typename INTYPE, class INOP, typename OUTTYPE>
void MEDDLY::copy_EV2EV_slow<INTYPE,INOP,OUTTYPE>
::computeAll(int in, int k, INTYPE av, node_handle an,
  OUTTYPE &bv, node_handle &bn)
{
  // Check terminals
  if (0==k) {
    bv = av;
    terminal t(an != 0);
    bn = t.getHandle();
    return;
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(an);

  // Check compute table
  ct_entry_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = inCache(av, an, bv, bn);
    if (0==Key) return;
  }

  // What's below?
  int nextk;
  if (resF->isForRelations()) {
    nextk = resF->downLevel(k);
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  unpacked_node* A = unpacked_node::New(argF);
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      INTYPE rev;
      INOP::redundant(rev);
      A->initIdentity(argF, k, unsigned(in), rev, an, SPARSE_ONLY);
    } else {
      INTYPE rev;
      INOP::redundant(rev);
      A->initRedundant(argF, k, rev, an, SPARSE_ONLY);
    }
  } else {
    argF->unpackNode(A, an, SPARSE_ONLY);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getSize());

  // recurse
  for (unsigned z=0; z<A->getSize(); z++) {
    INTYPE adv;
    A->edgeval(z).get(adv);
    // A->getEdge(z, adv);
    // nb->i_ref(z) = A->index(z);
    OUTTYPE bdv;
    node_handle bdn;
    computeAll(int(A->index(z)), nextk, INOP::apply(av, adv), A->down(z), bdv, bdn);
    nb->setSparse(z, A->index(z), edge_value(bdv), bdn);
  }

  // Cleanup
  unpacked_node::Recycle(A);

  // Reduce
  resF->createReducedNode(in, nb, bv, bn);

  // Add to compute table
  if (Key) addToCache(Key, av, an, bv, bn);
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    class PLUS {
        public:
            template <typename T>
            static inline T apply(T a, T b) { return a+b; }
            template <typename T>
            static inline void redundant(T &ev) { ev = 0; }
    };
    class TIMES {
        public:
            template <typename T>
            static inline T apply(T a, T b) { return a*b; }
            template <typename T>
            static inline void redundant(T &ev) { ev = 1; }
    };
};

MEDDLY::unary_operation* MEDDLY::COPY(forest* arg, forest* res)
{
    if (!arg || !res) {
        return nullptr;
    }

    if (arg->isForRelations() != res->isForRelations()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    if (arg->isMultiTerminal() && res->isMultiTerminal())
    {
        //
        // MT copies, handled by the new template class!
        //
        switch (res->getRangeType()) {
            case range_type::BOOLEAN:
                return copy_MT_tmpl<bool>::build(arg, res);

            case range_type::INTEGER:
                return copy_MT_tmpl<int>::build(arg, res);

            case range_type::REAL:
                return copy_MT_tmpl<float>::build(arg, res);

            default:  // any other types?
                throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        };
    }

    //
    // Must be at least one EV forest
    //

    if (arg->isMultiTerminal() &&
        (res->isEVPlus() || res->isEVTimes()))
    {
        //
        // MT to EV conversion
        //
        switch (res->getRangeType()) {
            case range_type::INTEGER:
                return copy_MT2EV<long>::build(arg, res, "N:LN");

            case range_type::REAL:
                return copy_MT2EV<float>::build(arg, res, "N:FN");

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        };
    }

    if (res->isMultiTerminal() && (arg->isEVPlus() || arg->isIndexSet()))
    {
        //
        // EV+ to MT conversion
        //
        switch (arg->getRangeType()) {
            case range_type::INTEGER:
                return copy_EV2MT<long,PLUS>::build(arg, res, "LN:N");

            case range_type::REAL:
                return copy_EV2MT<float,PLUS>::build(arg, res, "FN:N");

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        };
    }

    if (res->isMultiTerminal() && arg->isEVTimes())
    {
        //
        // EV* to MT conversion  (untested!)
        //
        switch (arg->getRangeType()) {
            case range_type::INTEGER:
                return copy_EV2MT<long,TIMES>::build(arg, res, "LN:N");

            case range_type::REAL:
                return copy_EV2MT<float,TIMES>::build(arg, res, "FN:N");

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        };
    }

    //
    // That's it for any MT arguments
    //
    if (res->isMultiTerminal() || arg->isMultiTerminal()) {
        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    }


    //
    // "fast" EV to EV copies (same RR, same operation, no info loss)
    //
    if (arg->getReductionRule() == res->getReductionRule()) {

        if ( ((arg->isEVPlus() || arg->isIndexSet()) && res->isEVPlus())
            || (arg->isEVTimes() && res->isEVTimes()) )
        {

            switch (arg->getRangeType()) {
                case range_type::INTEGER:
                    switch (res->getRangeType()) {
                        case range_type::INTEGER:
                            return copy_EV2EV_fast<long,long>::build(arg, res);
                        case range_type::REAL:
                            return copy_EV2EV_fast<long,float>::build(arg, res);
                        default:
                            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
                    };
                    break;    // in case anything falls through

                case range_type::REAL:
                    switch (res->getRangeType()) {
                        case range_type::INTEGER:
                            break;
                            // not safe to go from real -> integer c way
                        case range_type::REAL:
                            return copy_EV2EV_fast<float,float>::build(arg, res);
                        default:
                            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
                    };
                    break;    // things may fall through

                default:
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            };
        }

    }

    //
    // Generic EV+ to EV, may be slow
    //
    if (arg->isEVPlus() || arg->isIndexSet()) {

        switch (arg->getRangeType()) {

            case range_type::INTEGER:
                switch (res->getRangeType()) {
                    case range_type::INTEGER:
                        return copy_EV2EV_slow<long,PLUS,long>::build(arg, res, "LN:LN");
                    case range_type::REAL:
                        return copy_EV2EV_slow<long,PLUS,float>::build(arg, res, "LN:FN");
                    default:
                        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
                };
                break;

            case range_type::REAL:
                switch (res->getRangeType()) {
                    case range_type::INTEGER:
                        return copy_EV2EV_slow<float,PLUS,long>::build(arg, res, "FN:LN");
                    case range_type::REAL:
                        return copy_EV2EV_slow<float,PLUS,float>::build(arg, res, "FN:FN");
                    default:
                        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
                };
                break;

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        }
    }

    //
    // Generic EV* to EV, may be slow
    //
    if (arg->isEVTimes()) {

        switch (arg->getRangeType()) {

            case range_type::INTEGER:
                switch (res->getRangeType()) {
                    case range_type::INTEGER:
                        return copy_EV2EV_slow<long,TIMES,long>::build(arg, res, "LN:LN");
                    case range_type::REAL:
                        return copy_EV2EV_slow<long,TIMES,float>::build(arg, res, "LN:FN");
                    default:
                        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
                };
                break;

            case range_type::REAL:
                switch (res->getRangeType()) {
                    case range_type::INTEGER:
                        return copy_EV2EV_slow<float,TIMES,long>::build(arg, res, "FN:LN");
                    case range_type::REAL:
                        return copy_EV2EV_slow<float,TIMES,float>::build(arg, res, "FN:FN");
                    default:
                        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
                };
                break;

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        }
    }

    //
    // Catch all for any other cases
    //
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::COPY_init()
{
    copy_MT_tmpl<bool>::cache.reset("Copy");
    copy_MT_tmpl<int>::cache.reset("Copy");
    copy_MT_tmpl<float>::cache.reset("Copy");

    copy_MT2EV<long>::cache.reset("Copy");
    copy_MT2EV<float>::cache.reset("Copy");

    copy_EV2MT<long,PLUS>::cache.reset("Copy");
    copy_EV2MT<float,PLUS>::cache.reset("Copy");
    copy_EV2MT<long,TIMES>::cache.reset("Copy");
    copy_EV2MT<float,TIMES>::cache.reset("Copy");

    copy_EV2EV_fast<long,long>::cache.reset("Copy");
    copy_EV2EV_fast<long,float>::cache.reset("Copy");
    copy_EV2EV_fast<float,float>::cache.reset("Copy");

    copy_EV2EV_slow<long,PLUS,long>::cache.reset("Copy");
    copy_EV2EV_slow<long,PLUS,float>::cache.reset("Copy");
    copy_EV2EV_slow<float,PLUS,long>::cache.reset("Copy");
    copy_EV2EV_slow<float,PLUS,float>::cache.reset("Copy");

    copy_EV2EV_slow<long,TIMES,long>::cache.reset("Copy");
    copy_EV2EV_slow<long,TIMES,float>::cache.reset("Copy");
    copy_EV2EV_slow<float,TIMES,long>::cache.reset("Copy");
    copy_EV2EV_slow<float,TIMES,float>::cache.reset("Copy");
}

void MEDDLY::COPY_done()
{
    MEDDLY_DCASSERT(copy_MT_tmpl<bool>::cache.isEmpty());
    MEDDLY_DCASSERT(copy_MT_tmpl<int>::cache.isEmpty());
    MEDDLY_DCASSERT(copy_MT_tmpl<float>::cache.isEmpty());

    MEDDLY_DCASSERT(copy_MT2EV<long>::cache.isEmpty());
    MEDDLY_DCASSERT(copy_MT2EV<float>::cache.isEmpty());

    MEDDLY_DCASSERT( (copy_EV2MT<long,PLUS>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2MT<float,PLUS>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2MT<long,TIMES>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2MT<float,TIMES>::cache.isEmpty()) );

    MEDDLY_DCASSERT( (copy_EV2EV_fast<long,long>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_fast<long,float>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_fast<float,float>::cache.isEmpty()) );

    MEDDLY_DCASSERT( (copy_EV2EV_slow<long,PLUS,long>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_slow<long,PLUS,float>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_slow<float,PLUS,long>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_slow<float,PLUS,float>::cache.isEmpty()) );

    MEDDLY_DCASSERT( (copy_EV2EV_slow<long,TIMES,long>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_slow<long,TIMES,float>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_slow<float,TIMES,long>::cache.isEmpty()) );
    MEDDLY_DCASSERT( (copy_EV2EV_slow<float,TIMES,float>::cache.isEmpty()) );
}

