
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
#include "copy.h"

// #define DEBUG_COPY_COMPUTE_ALL

namespace MEDDLY {
  class copy_MT;

  class copy_opname;
};

// ******************************************************************
// *                                                                *
// *                         copy_MT  class                         *
// *                                                                *
// ******************************************************************

/// Abstract base class for copying between multi-terminal DDs.
class MEDDLY::copy_MT : public unary_operation {
  public:
    copy_MT(const unary_opname* oc, expert_forest* arg, expert_forest* res);

#ifdef OLD_OP_CT
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle*);
#endif

    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle* entryData, bool key_only) const;
#endif

    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res);
  protected:
    virtual node_handle compute_r(node_handle a) = 0;

    inline compute_table::entry_key* 
    findResult(node_handle a, node_handle &b) 
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeN(a);
#ifdef OLD_OP_CT
      compute_table::entry_result& cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      b = resF->linkNode(cacheFind.readN());
#else
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      b = resF->linkNode(CTresult[0].readN());
#endif
      CT0->recycle(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::entry_key* Key, 
      node_handle a, node_handle b) 
    {
#ifdef OLD_OP_CT
      argF->cacheNode(a);
      resF->cacheNode(b);
      static compute_table::entry_result result(1);
      result.reset();
      result.writeN(b);
      CT0->addEntry(Key, result);
#else
      CTresult[0].reset();
      CTresult[0].writeN(b);
      CT0->addEntry(Key, CTresult[0]);
#endif
      return b;
    }
};

#ifdef OLD_OP_CT
MEDDLY::copy_MT
:: copy_MT(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 1, 1, arg, res)
{
  // mtres = res;
}
#else
MEDDLY::copy_MT
:: copy_MT(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 1, arg, res)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:N");
  et->setForestForSlot(0, arg);
  et->setForestForSlot(2, res);
  registerEntryType(0, et);
  buildCTs();
}
#endif

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
bool MEDDLY::copy_MT::isStaleEntry(const node_handle* entryData)
{
  return 
    argF->isStale(entryData[0]) ||
    resF->isStale(entryData[1]);
}
#else
MEDDLY::forest::node_status
MEDDLY::copy_MT::getStatusOfEntry(const node_handle* data)
{
  MEDDLY::forest::node_status a = argF->getNodeStatus(data[0]);
  MEDDLY::forest::node_status c = resF->getNodeStatus(data[1]);

  if (a == MEDDLY::forest::DEAD ||
      c == MEDDLY::forest::DEAD)
    return MEDDLY::forest::DEAD;
  else if (a == MEDDLY::forest::RECOVERABLE ||
      c == MEDDLY::forest::RECOVERABLE)
    return MEDDLY::forest::RECOVERABLE;
  else
    return MEDDLY::forest::ACTIVE;
}
#endif

void MEDDLY::copy_MT::discardEntry(const node_handle* entryData)
{
  argF->uncacheNode(entryData[0]);
  resF->uncacheNode(entryData[1]);
}

void MEDDLY::copy_MT::showEntry(output &strm, const node_handle* eD, bool key_only) const
{
  strm << "[" << getName() << "(" << long(eD[0]) << "): ";
  if (key_only) {
    strm << "?";
  } else {
    strm << long(eD[1]);
  }
  strm << "]";
}

#endif // OLD_OP_CT

void MEDDLY::copy_MT::computeDDEdge(const dd_edge &arg, dd_edge &res)
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
    public:
      copy_MT_tmpl(const unary_opname* N, expert_forest* A, expert_forest* R)
        : copy_MT(N, A, R) { }
    protected:
      virtual node_handle compute_r(node_handle a) {
        if (argF->getReductionRule() == resF->getReductionRule()) {
          return computeSkip(-1, a);  // same skipping rule, ok
        } else {
          // need to visit every level...
          return computeAll(-1, resF->getNumVariables(), a);
        }
      }
  
      node_handle computeSkip(int in, node_handle a);
      node_handle computeAll(int in, int k, node_handle a);
  };

};  // namespace MEDDLY

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
  compute_table::entry_key* Key = findResult(a, b);
  if (0==Key) return b;

  // Initialize node reader
  unpacked_node* A = unpacked_node::useUnpackedNode();
  A->initFromNode(argF, a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    nb->i_ref(z) = A->i(z);
    nb->d_ref(z) = computeSkip(A->i(z), A->d(z));
  }

  // Handle extensible edge, if any
  if (A->isExtensible()) nb->markAsExtensible();

  // Cleanup
  unpacked_node::recycle(A);

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
  compute_table::entry_key* Key = 0;
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
  unpacked_node* A = unpacked_node::useUnpackedNode();
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A->initIdentity(argF, k, in, a, false);
    } else {
      A->initRedundant(argF, k, a, false);
    }
  } else {
    A->initFromNode(argF, a, false);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    nb->i_ref(z) = A->i(z);
    nb->d_ref(z) = computeAll(A->i(z), nextk, A->d(z));
  }

  // Handle extensible edge, if any
  if (A->isExtensible()) nb->markAsExtensible();

  // Cleanup
  unpacked_node::recycle(A);

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
    public:
#ifdef OLD_OP_CT
      copy_MT2EV(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc,
          sizeof(node_handle) / sizeof(node_handle),
          (sizeof(TYPE) + sizeof(node_handle)) / sizeof(node_handle),
          arg, res)
      {
        // entry[0]: mt node 
        // entry[1]: EV value (output)
        // entry[2]: EV node (output)
      }
#else
      copy_MT2EV(const unary_opname* oc, expert_forest* arg, expert_forest* res, 
        const char* pattern) : unary_operation(oc, 1, arg, res)
      {
        //
        // Pattern should be of the form "N:xN" where x is the EV Type.
        //
        // entry[0]: mt node 
        // entry[1]: EV value (output)
        // entry[2]: EV node (output)
        //
        compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), pattern);
        et->setForestForSlot(0, arg);
        et->setForestForSlot(3, res);
        registerEntryType(0, et);
        buildCTs();
      }
#endif

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[0]) ||
          resF->isStale(entryData[(sizeof(node_handle) + sizeof(TYPE)) / sizeof(node_handle)]);
      }
#else
      virtual MEDDLY::forest::node_status
      getStatusOfEntry(const node_handle* data)
      {
        MEDDLY::forest::node_status a = argF->getNodeStatus(data[0]);
        MEDDLY::forest::node_status c = resF->getNodeStatus(data[2]);

        if (a == MEDDLY::forest::DEAD ||
            c == MEDDLY::forest::DEAD)
          return MEDDLY::forest::DEAD;
        else if (a == MEDDLY::forest::RECOVERABLE ||
            c == MEDDLY::forest::RECOVERABLE)
          return MEDDLY::forest::RECOVERABLE;
        else
          return MEDDLY::forest::ACTIVE;
      }
#endif

      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[0]);
        resF->uncacheNode(entryData[(sizeof(node_handle) + sizeof(TYPE)) / sizeof(node_handle)]);
      }
      virtual void showEntry(output &strm, const node_handle* entryData, bool key_only) const {
        TYPE ev;
        compute_table::readEV(entryData + sizeof(node_handle) / sizeof(node_handle), ev);
        strm << "[" << getName()
          << "(" << long(entryData[0])
          << "): ";
        if (key_only) {
          strm << "?";
        } else {
          strm << "<" << ev << ", " << long(entryData[(sizeof(node_handle) + sizeof(TYPE)) / sizeof(node_handle)]) << ">";
        } 
        strm << "]";
      }

#endif // OLD_OP_CT

      virtual void computeDDEdge(const dd_edge &arg, dd_edge &res) {
        node_handle b;
        TYPE bev = Inf<TYPE>();
        if (argF->getReductionRule() == resF->getReductionRule()) {
          computeSkip(-1, arg.getNode(), b, bev);  // same skipping rule, ok
        } else {
          // need to visit every level...
          computeAll(-1, resF->getNumVariables(), arg.getNode(), b, bev);
        }
        res.set(b, bev);
      }
      void computeSkip(int in, node_handle a, node_handle &b, TYPE &bev);
      void computeAll(int in, int k, node_handle a, node_handle &b, TYPE &bev);

    protected:
      inline compute_table::entry_key* 
      inCache(node_handle a, node_handle &b, TYPE &bev) 
      {
#ifdef OLD_OP_CT
        compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
        compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->writeN(a);
#ifdef OLD_OP_CT
        compute_table::entry_result& cacheFind = CT0->find(CTsrch);
        if (cacheFind) {
          cacheFind.read_ev(bev);
          b = resF->linkNode(cacheFind.readN());
          CT0->recycle(CTsrch);
          return 0;
        }
#else
        CT0->find(CTsrch, CTresult[0]);
        if (CTresult[0]) {
          CTresult[0].read_ev(bev);
          b = resF->linkNode(CTresult[0].readN());
          CT0->recycle(CTsrch);
          return 0;
        }
#endif
        return CTsrch;
      }

      inline void addToCache(compute_table::entry_key* Key,
        node_handle a, node_handle b, long bev) 
      {
        MEDDLY_DCASSERT(bev != Inf<long>());

#ifdef OLD_OP_CT
        argF->cacheNode(a);
        resF->cacheNode(b);
        static compute_table::entry_result result(1 + sizeof(long)/sizeof(node_handle));
        result.reset();
        result.writeL(bev);
        result.writeN(b);
        CT0->addEntry(Key, result);
#else
        CTresult[0].reset();
        CTresult[0].writeL(bev);
        CTresult[0].writeN(b);
        CT0->addEntry(Key, CTresult[0]);
#endif
      }

      inline void addToCache(compute_table::entry_key* Key,
        node_handle a, node_handle b, float bev) 
      {
        MEDDLY_DCASSERT(bev != Inf<float>());

#ifdef OLD_OP_CT
        argF->cacheNode(a);
        resF->cacheNode(b);
        static compute_table::entry_result result(2);
        result.reset();
        result.writeF(bev);
        result.writeN(b);
        CT0->addEntry(Key, result);
#else
        CTresult[0].reset();
        CTresult[0].writeF(bev);
        CTresult[0].writeN(b);
        CT0->addEntry(Key, CTresult[0]);
#endif
      }

  };

};  // namespace MEDDLY

template <typename TYPE>
void MEDDLY::copy_MT2EV<TYPE>
::computeSkip(int in, node_handle a, node_handle &b, TYPE &bev)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    MEDDLY_DCASSERT(a != argF->getTransparentNode());
    if (argF->getRangeType() == forest::BOOLEAN) {
      bev = 0;
    }
    else {
      argF->getValueFromHandle(a, bev);
    }
    b = expert_forest::bool_Tencoder::value2handle(true);
    return;
  }

  // Check compute table
  compute_table::entry_key* Key = inCache(a, b, bev);
  if (0==Key) return;

  // Initialize sparse node reader
  unpacked_node* A = unpacked_node::useUnpackedNode();
  A->initFromNode(argF, a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    node_handle d;
    TYPE dev = Inf<TYPE>();
    computeSkip(A->i(z), A->d(z), d, dev);
    nb->i_ref(z) = A->i(z);
    nb->d_ref(z) = d;
    nb->setEdge(z, dev);
  }

  // Cleanup
  unpacked_node::recycle(A);

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
    if (argF->getRangeType() == forest::BOOLEAN) {
      bev = 0;
    }
    else {
      argF->getValueFromHandle(a, bev);
    }
    b = expert_forest::bool_Tencoder::value2handle(true);
    return;
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  compute_table::entry_key* Key = 0;
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
  unpacked_node* A = unpacked_node::useUnpackedNode();
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A->initIdentity(argF, k, in, a, false);
    } else {
      A->initRedundant(argF, k, a, false);
    }
  } else {
    A->initFromNode(argF, a, false);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    TYPE dev;
    nb->i_ref(z) = A->i(z);
    computeAll(A->i(z), nextk, A->d(z), nb->d_ref(z), dev);
    nb->setEdge(z, dev);
  }

  // Cleanup
  unpacked_node::recycle(A);

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
    public:
#ifdef OLD_OP_CT
      copy_EV2MT(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc,
          (sizeof(TYPE) + sizeof(node_handle)) / sizeof(node_handle),
          sizeof(node_handle) / sizeof(node_handle),
          arg, res)
      {
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: mt node (output)
      }
#else
      copy_EV2MT(const unary_opname* oc, expert_forest* arg, expert_forest* res, 
        const char* pattern) : unary_operation(oc, 1, arg, res)
      {
        //
        // Pattern should be of the form "xN:N"
        //
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: mt node (output)
        //
        compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), pattern);
        et->setForestForSlot(1, arg);
        et->setForestForSlot(3, res);
        registerEntryType(0, et);
        buildCTs();
      }
#endif

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[sizeof(TYPE) / sizeof(node_handle)]) ||
          resF->isStale(entryData[(sizeof(TYPE) + sizeof(node_handle)) / sizeof(node_handle)]);
      }
#else
      virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* data) {
        MEDDLY::forest::node_status a = argF->getNodeStatus(data[1]);
        MEDDLY::forest::node_status c = resF->getNodeStatus(data[2]);

        if (a == MEDDLY::forest::DEAD ||
            c == MEDDLY::forest::DEAD)
          return MEDDLY::forest::DEAD;
        else if (a == MEDDLY::forest::RECOVERABLE ||
            c == MEDDLY::forest::RECOVERABLE)
          return MEDDLY::forest::RECOVERABLE;
        else
          return MEDDLY::forest::ACTIVE;
      }
#endif
      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[sizeof(TYPE) / sizeof(node_handle)]);
        resF->uncacheNode(entryData[(sizeof(TYPE) + sizeof(node_handle)) / sizeof(node_handle)]);
      }
      virtual void showEntry(output &strm, const node_handle* entryData, bool key_only) const {
        TYPE ev;
        compute_table::readEV(entryData, ev);
        strm << "[" << getName()
          << "(<" << ev << "," << long(entryData[sizeof(TYPE) / sizeof(node_handle)])
          << ">): ";
        if (key_only) {
          strm << "?";
        } else {
          strm << long(entryData[(sizeof(TYPE) + sizeof(node_handle)) / sizeof(node_handle)]);
        }
        strm << "]";
      }

#endif // OLD_OP_CT

      virtual void computeDDEdge(const dd_edge &arg, dd_edge &res) {
        TYPE ev;
        node_handle b;
        arg.getEdgeValue(ev);
        if (argF->getReductionRule() == resF->getReductionRule()) {
          b = computeSkip(-1, ev, arg.getNode());  // same skipping rule, ok
        } else {
          // need to visit every level...
          b = computeAll(-1, resF->getNumVariables(), ev, arg.getNode());
        }
        res.set(b);
      }
      node_handle computeSkip(int in, TYPE ev, node_handle a);
      node_handle computeAll(int in, int k, TYPE ev, node_handle a);

    protected:
      inline compute_table::entry_key* 
      inCache(TYPE ev, node_handle a, node_handle &b) 
      {
#ifdef OLD_OP_CT
        compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
        compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->write_ev(ev);
        CTsrch->writeN(a);
#ifdef OLD_OP_CT
        compute_table::entry_result& cacheFind = CT0->find(CTsrch);
        if (cacheFind) {
          b = resF->linkNode(cacheFind.readN());
          CT0->recycle(CTsrch);
          return 0;
        }
#else
        CT0->find(CTsrch, CTresult[0]);
        if (CTresult[0]) {
          b = resF->linkNode(CTresult[0].readN());
          CT0->recycle(CTsrch);
          return 0;
        }
#endif
        return CTsrch;
      }

      inline void addToCache(compute_table::entry_key* Key, 
        TYPE ev, node_handle a, node_handle b) 
      {
#ifdef OLD_OP_CT
        argF->cacheNode(a);
        resF->cacheNode(b);
        static compute_table::entry_result result(1);
        result.reset();
        result.writeN(b);
        CT0->addEntry(Key, result);
#else
        CTresult[0].reset();
        CTresult[0].writeN(b);
        CT0->addEntry(Key, CTresult[0]);
#endif
      }

  };

};  // namespace MEDDLY

template <typename TYPE, class OP>
MEDDLY::node_handle  MEDDLY::copy_EV2MT<TYPE,OP>
::computeSkip(int in, TYPE ev, node_handle a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    MEDDLY_DCASSERT(a != argF->getTransparentNode());
    if (resF->getRangeType() == forest::BOOLEAN) {
      return expert_forest::bool_Tencoder::value2handle(true);
    }
    else {
      return resF->handleForValue(ev);
    }
  }

  // Check compute table
  node_handle b;
  compute_table::entry_key* Key = inCache(ev, a, b);
  if (0==Key) return b;

  // Initialize sparse node reader
  unpacked_node* A = unpacked_node::useUnpackedNode();
  A->initFromNode(argF, a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    TYPE aev;
    A->getEdge(z, aev);
    nb->i_ref(z) = A->i(z);
    nb->d_ref(z) = computeSkip(A->i(z), OP::apply(ev, aev), A->d(z));
  }

  // Cleanup
  unpacked_node::recycle(A);

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
    if (resF->getRangeType() == forest::BOOLEAN) {
      return expert_forest::bool_Tencoder::value2handle(true);
    }
    else {
      return resF->handleForValue(ev);
    }
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  node_handle b;
  compute_table::entry_key* Key = 0;
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
  unpacked_node* A = unpacked_node::useUnpackedNode();
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      TYPE rev;
      OP::redundant(rev);
      A->initIdentity(argF, k, in, rev, a, false);
    } else {
      TYPE rev;
      OP::redundant(rev);
      A->initRedundant(argF, k, rev, a, false);
    }
  } else {
    A->initFromNode(argF, a, false);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    TYPE aev;
    A->getEdge(z, aev);
    nb->i_ref(z) = A->i(z);
    nb->d_ref(z) = computeAll(A->i(z), nextk, OP::apply(ev, aev), A->d(z));
  }

  // Cleanup
  unpacked_node::recycle(A);

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
    public:
#ifdef OLD_OP_CT
      copy_EV2EV_fast(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc,
          sizeof(node_handle) / sizeof(node_handle),
          sizeof(node_handle) / sizeof(node_handle),
          arg, res)
      {
        // entry[0]: EV node
        // entry[1]: EV node 
      }
#else
      copy_EV2EV_fast(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc, 1, arg, res)
      {
        // entry[0]: EV node
        // entry[1]: EV node 
        compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:N");
        et->setForestForSlot(0, arg);
        et->setForestForSlot(2, res);
        registerEntryType(0, et);
        buildCTs();
      }
#endif

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[0]) ||
          resF->isStale(entryData[1]);
      }
#else
      virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* data) {
        MEDDLY::forest::node_status a = argF->getNodeStatus(data[0]);
        MEDDLY::forest::node_status c = resF->getNodeStatus(data[1]);

        if (a == MEDDLY::forest::DEAD ||
            c == MEDDLY::forest::DEAD)
          return MEDDLY::forest::DEAD;
        else if (a == MEDDLY::forest::RECOVERABLE ||
            c == MEDDLY::forest::RECOVERABLE)
          return MEDDLY::forest::RECOVERABLE;
        else
          return MEDDLY::forest::ACTIVE;
      }
#endif
      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[0]);
        resF->uncacheNode(entryData[1]);
      }
      virtual void showEntry(output &strm, const node_handle* entryData, bool key_only) const {
        strm << "[" << getName()
          << "(<?," << long(entryData[0])
          << ">): ";
        if (key_only) {
          strm << "?";
        } else {
          strm << "<?," << long(entryData[1]) << ">";
        }
        strm << "]";
      }

#endif // OLD_OP_CT

      virtual void computeDDEdge(const dd_edge &arg, dd_edge &res) {
        INTYPE av;
        node_handle bn;
        arg.getEdgeValue(av);
        bn = computeSkip(-1, arg.getNode());
        OUTTYPE bv = av;
        res.set(bn, bv);
      }

      node_handle computeSkip(int in, node_handle a);

    protected:
      inline compute_table::entry_key* 
      findResult(node_handle a, node_handle &b) 
      {
#ifdef OLD_OP_CT
        compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
        compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->writeN(a);
#ifdef OLD_OP_CT
        compute_table::entry_result& cacheFind = CT0->find(CTsrch);
        if (!cacheFind) return CTsrch;
        b = resF->linkNode(cacheFind.readN());
#else
        CT0->find(CTsrch, CTresult[0]);
        if (!CTresult[0]) return CTsrch;
        b = resF->linkNode(CTresult[0].readN());
#endif
        CT0->recycle(CTsrch);
        return 0;
      }
      inline node_handle saveResult(compute_table::entry_key* Key,
        node_handle a, node_handle b) 
      {
#ifdef OLD_OP_CT
        argF->cacheNode(a);
        resF->cacheNode(b);
        static compute_table::entry_result result(1);
        result.reset();
        result.writeN(b);
        CT0->addEntry(Key, result);
#else
        CTresult[0].reset();
        CTresult[0].writeN(b);
        CT0->addEntry(Key, CTresult[0]);
#endif
        return b;
      }
  };

};  // namespace MEDDLY

template <typename INTYPE, typename OUTTYPE>
MEDDLY::node_handle 
MEDDLY::copy_EV2EV_fast<INTYPE,OUTTYPE>::computeSkip(int in, node_handle a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    return expert_forest::bool_Tencoder::value2handle(a != 0);
  }

  // Check compute table
  node_handle b;
  compute_table::entry_key* Key = findResult(a, b);
  if (0==Key) return b;

  // Initialize node reader
  unpacked_node* A = unpacked_node::useUnpackedNode();
  A->initFromNode(argF, a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  unpacked_node* nb = unpacked_node::newSparse(resF, level, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    nb->i_ref(z) = A->i(z);
    nb->d_ref(z) = computeSkip(A->i(z), A->d(z));
    INTYPE av;
    OUTTYPE bv;
    A->getEdge(z, av);
    bv = av;
    nb->setEdge(z, bv);
  }

  // Cleanup
  unpacked_node::recycle(A);

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
    public:
#ifdef OLD_OP_CT
      copy_EV2EV_slow(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc,
          (sizeof(INTYPE) + sizeof(node_handle)) / sizeof(node_handle),
          (sizeof(OUTTYPE) + sizeof(node_handle)) / sizeof(node_handle),
          arg, res)
      {
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: EV value
        // entry[3]: EV node 
      }
#else
      copy_EV2EV_slow(const unary_opname* oc, expert_forest* arg, expert_forest* res, 
        const char* pattern) : unary_operation(oc, 1, arg, res)
      {
        //
        // Pattern is of the form "xN:yN"
        //
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: EV value
        // entry[3]: EV node 
        //
        compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), pattern);
        et->setForestForSlot(1, arg);
        et->setForestForSlot(4, res);
        registerEntryType(0, et);
        buildCTs();
      }
#endif

#ifdef OLD_OP_CT

#ifndef USE_NODE_STATUS
      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[sizeof(INTYPE) / sizeof(node_handle)]) ||
          resF->isStale(entryData[(sizeof(INTYPE) + sizeof(node_handle) + sizeof(OUTTYPE)) / sizeof(node_handle)]);
      }
#else
      virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* data) {
        MEDDLY::forest::node_status a = argF->getNodeStatus(data[1]);
        MEDDLY::forest::node_status c = resF->getNodeStatus(data[3]);

        if (a == MEDDLY::forest::DEAD ||
            c == MEDDLY::forest::DEAD)
          return MEDDLY::forest::DEAD;
        else if (a == MEDDLY::forest::RECOVERABLE ||
            c == MEDDLY::forest::RECOVERABLE)
          return MEDDLY::forest::RECOVERABLE;
        else
          return MEDDLY::forest::ACTIVE;
      }
#endif
      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[sizeof(INTYPE) / sizeof(node_handle)]);
        resF->uncacheNode(entryData[(sizeof(INTYPE) + sizeof(node_handle) + sizeof(OUTTYPE)) / sizeof(node_handle)]);
      }
      virtual void showEntry(output &strm, const node_handle* entryData, bool key_only) const {
        INTYPE ev1;
        compute_table::readEV(entryData, ev1);
        node_handle n1 = entryData[sizeof(INTYPE) / sizeof(node_handle)];
        OUTTYPE ev2;
        compute_table::readEV(entryData + (sizeof(INTYPE) + sizeof(node_handle)) / sizeof(node_handle), ev2);
        node_handle n2 = entryData[(sizeof(INTYPE) + sizeof(node_handle) + sizeof(OUTTYPE)) / sizeof(node_handle)];

        strm << "[" << getName()
          << "(<" << ev1 << "," << n1
          << ">): ";
        if (key_only) {
          strm << "?";
        } else {
          strm << "<" << ev2 << "," << n2 << ">";
        }
        strm << "]";
      }

#endif // OLD_OP_CT

      virtual void computeDDEdge(const dd_edge &arg, dd_edge &res) {
        INTYPE av;
        node_handle an, bn;
        OUTTYPE bv;
        arg.getEdgeValue(av);
        an = arg.getNode();

        if (argF->isTransparentEdge(an, &av) && argF->getTransparentNode() == resF->getTransparentNode()) {
          resF->getTransparentEdge(bn, &bv);
        }
        else {
          // need to visit every level...
          computeAll(-1, resF->getNumVariables(), av, an, bv, bn);
        }
        res.set(bn, bv);
      }
      void computeAll(int in, int k, INTYPE av, node_handle an,
        OUTTYPE &bv, node_handle &bn);

    protected:
      inline compute_table::entry_key* 
      inCache(INTYPE av, node_handle an, OUTTYPE &bv, node_handle &bn) 
      {
#ifdef OLD_OP_CT
        compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
        compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->write_ev(av);
        CTsrch->writeN(an);
#ifdef OLD_OP_CT
        compute_table::entry_result& cacheFind = CT0->find(CTsrch);
        if (cacheFind) {
          cacheFind.read_ev(bv);
          bn = resF->linkNode(cacheFind.readN());
          CT0->recycle(CTsrch);
          return 0;
        }
#else
        CT0->find(CTsrch, CTresult[0]);
        if (CTresult[0]) {
          CTresult[0].read_ev(bv);
          bn = resF->linkNode(CTresult[0].readN());
          CT0->recycle(CTsrch);
          return 0;
        }
#endif
        return CTsrch;
      }

      inline void addToCache(compute_table::entry_key* Key,
        INTYPE av, node_handle an, long bv, node_handle &bn)
      {
#ifdef OLD_OP_CT
        argF->cacheNode(an);
        resF->cacheNode(bn);
        static compute_table::entry_result result(1 + sizeof(long)/sizeof(node_handle));
        result.reset();
        result.writeL(bv);
        result.writeN(bn);
        CT0->addEntry(Key, result);
#else
        CTresult[0].reset();
        CTresult[0].writeL(bv);
        CTresult[0].writeN(bn);
        CT0->addEntry(Key, CTresult[0]);
#endif
      }

      inline void addToCache(compute_table::entry_key* Key,
        INTYPE av, node_handle a, float bv, node_handle &bn)
      {
#ifdef OLD_OP_CT
        argF->cacheNode(a);
        resF->cacheNode(bn);
        static compute_table::entry_result result(2);
        result.reset();
        result.writeF(bv);
        result.writeN(bn);
        CT0->addEntry(Key, result);
#else
        CTresult[0].reset();
        CTresult[0].writeF(bv);
        CTresult[0].writeN(bn);
        CT0->addEntry(Key, CTresult[0]);
#endif
      }


  };

};  // namespace MEDDLY

template <typename INTYPE, class INOP, typename OUTTYPE>
void MEDDLY::copy_EV2EV_slow<INTYPE,INOP,OUTTYPE>
::computeAll(int in, int k, INTYPE av, node_handle an, 
  OUTTYPE &bv, node_handle &bn)
{
  // Check terminals
  if (0==k) {
    bv = av;
    bn = expert_forest::bool_Tencoder::value2handle(an != 0);
    return;
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(an);

  // Check compute table
  compute_table::entry_key* Key = 0;
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
  unpacked_node* A = unpacked_node::useUnpackedNode();
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      INTYPE rev;
      INOP::redundant(rev);
      A->initIdentity(argF, k, in, rev, an, false);
    } else {
      INTYPE rev;
      INOP::redundant(rev);
      A->initRedundant(argF, k, rev, an, false);
    }
  } else {
    A->initFromNode(argF, an, false);
  }

  // Initialize node builder
  unpacked_node* nb = unpacked_node::newSparse(resF, k, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    INTYPE adv;
    OUTTYPE bdv;
    A->getEdge(z, adv);
    nb->i_ref(z) = A->i(z);
    computeAll(A->i(z), nextk, INOP::apply(av, adv), A->d(z), bdv, nb->d_ref(z));
    nb->setEdge(z, bdv);
  }

  // Cleanup
  unpacked_node::recycle(A);

  // Reduce
  resF->createReducedNode(in, nb, bv, bn);

  // Add to compute table
  if (Key) addToCache(Key, av, an, bv, bn);
}

// ******************************************************************
// *                                                                *
// *                       copy_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::copy_opname : public unary_opname {
  public:
    copy_opname();
    virtual unary_operation* 
      buildOperation(expert_forest* ar, expert_forest* res) const;


  private:
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

MEDDLY::copy_opname::copy_opname()
 : unary_opname("Copy")
{
}

MEDDLY::unary_operation*
MEDDLY::copy_opname
::buildOperation(expert_forest* arg, expert_forest* res) const
{
  if (0==arg || 0==res) return 0;

  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);

  if (arg->isForRelations() != res->isForRelations())
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

  if (arg->isMultiTerminal() && res->isMultiTerminal())
  {
    // 
    // MT copies, handled by the new template class!
    //
    switch (res->getRangeType()) {
      case forest::BOOLEAN:
        return new copy_MT_tmpl<bool>(this, arg, res);

      case forest::INTEGER:
        return new copy_MT_tmpl<int>(this, arg, res);

      case forest::REAL:
        return new copy_MT_tmpl<float>(this, arg, res);


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
#ifdef OLD_OP_CT
      case forest::INTEGER:
        return new copy_MT2EV<long>(this, arg, res);

      case forest::REAL:
        return new copy_MT2EV<float>(this, arg, res);
#else
      case forest::INTEGER:
        return new copy_MT2EV<long>(this, arg, res, "N:LN");

      case forest::REAL:
        return new copy_MT2EV<float>(this, arg, res, "N:FN");
#endif

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
#ifdef OLD_OP_CT
      case forest::INTEGER:
        return new copy_EV2MT<long,PLUS>(this, arg, res);

      case forest::REAL:
        return new copy_EV2MT<float,PLUS>(this, arg, res);
#else
      case forest::INTEGER:
        return new copy_EV2MT<long,PLUS>(this, arg, res, "LN:N");

      case forest::REAL:
        return new copy_EV2MT<float,PLUS>(this, arg, res, "FN:N");
#endif

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
#ifdef OLD_OP_CT
      case forest::INTEGER:
        return new copy_EV2MT<long,TIMES>(this, arg, res);

      case forest::REAL:
        return new copy_EV2MT<float,TIMES>(this, arg, res);
#else
      case forest::INTEGER:
        return new copy_EV2MT<long,TIMES>(this, arg, res, "LN:N");

      case forest::REAL:
        return new copy_EV2MT<float,TIMES>(this, arg, res, "FN:N");
#endif

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
        case forest::INTEGER:
            switch (res->getRangeType()) {
                case forest::INTEGER:
                    return new copy_EV2EV_fast<long,long>(this, arg, res);
                case forest::REAL:
                    return new copy_EV2EV_fast<long,float>(this, arg, res);
                default:
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            };
            break;    // in case anything falls through

        case forest::REAL:
            switch (res->getRangeType()) {
                case forest::INTEGER:
                    break;    // not safe to go from real -> integer this way
                case forest::REAL:
                    return new copy_EV2EV_fast<float,float>(this, arg, res);
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

      case forest::INTEGER:
          switch (res->getRangeType()) {
#ifdef OLD_OP_CT
            case forest::INTEGER:
                return new copy_EV2EV_slow<long,PLUS,long>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<long,PLUS,float>(this, arg, res);
#else
            case forest::INTEGER:
                return new copy_EV2EV_slow<long,PLUS,long>(this, arg, res, "LN:LN");
            case forest::REAL:
                return new copy_EV2EV_slow<long,PLUS,float>(this, arg, res, "LN:FN");
#endif
            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
          };
        
      case forest::REAL:
          switch (res->getRangeType()) {
#ifdef OLD_OP_CT
            case forest::INTEGER:
                return new copy_EV2EV_slow<float,PLUS,long>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<float,PLUS,float>(this, arg, res);
#else
            case forest::INTEGER:
                return new copy_EV2EV_slow<float,PLUS,long>(this, arg, res, "FN:LN");
            case forest::REAL:
                return new copy_EV2EV_slow<float,PLUS,float>(this, arg, res, "FN:FN");
#endif
            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
          };

      default:
          throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
  }

  //
  // Generic EV* to EV, may be slow
  //
  if (arg->isEVTimes()) {

    switch (arg->getRangeType()) {

      case forest::INTEGER:
          switch (res->getRangeType()) {
#ifdef OLD_OP_CT
            case forest::INTEGER:
                return new copy_EV2EV_slow<long,TIMES,long>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<long,TIMES,float>(this, arg, res);
#else
            case forest::INTEGER:
                return new copy_EV2EV_slow<long,TIMES,long>(this, arg, res, "LN:LN");
            case forest::REAL:
                return new copy_EV2EV_slow<long,TIMES,float>(this, arg, res, "LN:FN");
#endif
            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
          };
        
      case forest::REAL:
          switch (res->getRangeType()) {
#ifdef OLD_OP_CT
            case forest::INTEGER:
                return new copy_EV2EV_slow<float,TIMES,long>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<float,TIMES,float>(this, arg, res);
#else
            case forest::INTEGER:
                return new copy_EV2EV_slow<float,TIMES,long>(this, arg, res, "FN:LN");
            case forest::REAL:
                return new copy_EV2EV_slow<float,TIMES,float>(this, arg, res, "FN:FN");
#endif
            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
          };

      default:
          throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }
  }

  //
  // Catch all for any other cases
  //
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeCopy()
{
  return new copy_opname;
}

