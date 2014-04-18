
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

    virtual bool isStaleEntry(const node_handle* entryData);
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(FILE* strm, const node_handle* entryData) const;
    virtual void compute(const dd_edge &arg, dd_edge &res);
  protected:
    virtual node_handle compute(node_handle a) = 0;

    inline compute_table::search_key* 
    findResult(node_handle a, node_handle &b) 
    {
      compute_table::search_key* CTsrch = useCTkey();
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->reset();
      CTsrch->writeNH(a);
      compute_table::search_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      b = resF->linkNode(cacheFind.readNH());
      doneCTkey(CTsrch);
      return 0;
    }
    inline node_handle saveResult(compute_table::search_key* Key, 
      node_handle a, node_handle b) 
    {
      argF->cacheNode(a);
      compute_table::entry_builder &entry = CT->startNewEntry(Key);
      // entry.writeKeyNH(argF->cacheNode(a));
      entry.writeResultNH(resF->cacheNode(b));
      CT->addEntry();
      return b;
    }
};

MEDDLY::copy_MT
:: copy_MT(const unary_opname* oc, expert_forest* arg, expert_forest* res)
 : unary_operation(oc, 1, 1, arg, res)
{
  // mtres = res;
}

bool MEDDLY::copy_MT::isStaleEntry(const node_handle* entryData)
{
  return 
    argF->isStale(entryData[0]) ||
    resF->isStale(entryData[1]);
}

void MEDDLY::copy_MT::discardEntry(const node_handle* entryData)
{
  argF->uncacheNode(entryData[0]);
  resF->uncacheNode(entryData[1]);
}

void MEDDLY::copy_MT::showEntry(FILE* strm, const node_handle* entryData) const
{
  fprintf(strm, "[%s(%d) %d]", getName(), entryData[0], entryData[1]);
}

void MEDDLY::copy_MT::compute(const dd_edge &arg, dd_edge &res)
{
  node_handle result = compute(arg.getNode());
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
      virtual node_handle compute(node_handle a) {
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
  compute_table::search_key* Key = findResult(a, b);
  if (0==Key) return b;

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  node_builder& nb = resF->useSparseBuilder(level, A->getNNZs());


  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    nb.i(z) = A->i(z);
    nb.d(z) = computeSkip(A->i(z), A->d(z));
  }

  // Cleanup
  node_reader::recycle(A);

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

#ifdef DEBUG_COPY_COMPUTE_ALL
  fprintf(stderr, "copy(%d, %d, %d)\n", in, k, a);
#endif

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  node_handle b;
  compute_table::search_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = findResult(a, b);
    if (0==Key) return b;
  }
  int nextk;
  if (resF->isForRelations()) {
    nextk = (k>0) ? -k : -k-1;
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  node_reader* A;
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A = argF->initIdentityReader(k, in, a, false);
    } else {
      A = argF->initRedundantReader(k, a, false);
    }
  } else {
    A = argF->initNodeReader(a, false);
  }

  // Initialize node builder
  node_builder& nb = resF->useSparseBuilder(k, A->getNNZs());


  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    nb.i(z) = A->i(z);
    nb.d(z) = computeAll(A->i(z), nextk, A->d(z));
  }

  // Cleanup
  node_reader::recycle(A);

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
      copy_MT2EV(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc, 1, 2, arg, res)
      {
        // entry[0]: mt node 
        // entry[1]: EV value (output)
        // entry[2]: EV node (output)
      }

      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[0]) ||
          resF->isStale(entryData[2]);
      }
      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[0]);
        resF->uncacheNode(entryData[2]);
      }
      virtual void showEntry(FILE* strm, const node_handle* entryData) const {
        fprintf(strm, "[%s(%d) <", getName(), entryData[0]);
        TYPE ev;
        compute_table::readEV(entryData+1, ev);
        show(strm, ev);
        fprintf(strm, ", %d>]", entryData[2]);
      }
      virtual void compute(const dd_edge &arg, dd_edge &res) {
        node_handle b;
        TYPE bev;
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
      inline compute_table::search_key* 
      inCache(node_handle a, node_handle &b, TYPE &bev) 
      {
        compute_table::search_key* CTsrch = useCTkey();
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->reset();
        CTsrch->writeNH(a);
        compute_table::search_result &cacheFind = CT->find(CTsrch);
        if (cacheFind) {
          cacheFind.read(bev);
          b = resF->linkNode(cacheFind.readNH());
          doneCTkey(CTsrch);
          return 0;
        }
        return CTsrch;
      }

      inline void addToCache(compute_table::search_key* Key,
        node_handle a, node_handle b, TYPE bev) 
      {
        argF->cacheNode(a);
        compute_table::entry_builder &entry = CT->startNewEntry(Key);
        // entry.writeKeyNH(argF->cacheNode(a));
        entry.writeResult(bev);
        entry.writeResultNH(resF->cacheNode(b));
        CT->addEntry();
      }

    private:
      static inline void show(FILE* strm, int ev)    { fprintf(strm, "%d", ev); }
      static inline void show(FILE* strm, float ev)  { fprintf(strm, "%f", ev); }
  };

};  // namespace MEDDLY

template <typename TYPE>
void MEDDLY::copy_MT2EV<TYPE>
::computeSkip(int in, node_handle a, node_handle &b, TYPE &bev)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    argF->getValueFromHandle(a, bev);
    b = expert_forest::bool_Tencoder::value2handle(true);
    return;
  }

  // Check compute table
  compute_table::search_key* Key = inCache(a, b, bev);
  if (0==Key) return;

  // Initialize sparse node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  node_builder& nb = resF->useSparseBuilder(level, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    node_handle d;
    TYPE dev;
    computeSkip(A->i(z), A->d(z), d, dev);
    nb.i(z) = A->i(z);
    nb.d(z) = d;
    nb.setEdge(z, dev);
  }

  // Cleanup
  node_reader::recycle(A);

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
    argF->getValueFromHandle(a, bev);
    b = expert_forest::bool_Tencoder::value2handle(true);
    return;
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  compute_table::search_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = inCache(a, b, bev);
    if (0==Key) return;
  }

  // What's below?
  int nextk;
  if (resF->isForRelations()) {
    nextk = (k>0) ? -k : -k-1;
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  node_reader* A;
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A = argF->initIdentityReader(k, in, a, false);
    } else {
      A = argF->initRedundantReader(k, a, false);
    }
  } else {
    A = argF->initNodeReader(a, false);
  }

  // Initialize node builder
  node_builder& nb = resF->useSparseBuilder(k, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    TYPE dev;
    nb.i(z) = A->i(z);
    computeAll(A->i(z), nextk, A->d(z), nb.d(z), dev);
    nb.setEdge(z, dev);
  }

  // Cleanup
  node_reader::recycle(A);

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
      copy_EV2MT(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc, 2, 1, arg, res)
      {
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: mt node (output)
      }

      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[1]) ||
          resF->isStale(entryData[2]);
      }
      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[1]);
        resF->uncacheNode(entryData[2]);
      }
      virtual void showEntry(FILE* strm, const node_handle* entryData) const {
        fprintf(strm, "[%s(<", getName());
        TYPE ev;
        compute_table::readEV(entryData, ev);
        show(strm, ev);
        fprintf(strm, ",%d> %d]", entryData[1], entryData[2]);
      }
      virtual void compute(const dd_edge &arg, dd_edge &res) {
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
      inline compute_table::search_key* 
      inCache(TYPE ev, node_handle a, node_handle &b) 
      {
        compute_table::search_key* CTsrch = useCTkey();
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->reset();
        CTsrch->write(ev);
        CTsrch->writeNH(a);
        compute_table::search_result &cacheFind = CT->find(CTsrch);
        if (cacheFind) {
          b = resF->linkNode(cacheFind.readNH());
          doneCTkey(CTsrch);
          return 0;
        }
        return CTsrch;
      }

      inline void addToCache(compute_table::search_key* Key, 
        TYPE ev, node_handle a, node_handle b) 
      {
        argF->cacheNode(a);
        compute_table::entry_builder &entry = CT->startNewEntry(Key);
        /*
        entry.writeKey(ev);
        entry.writeKeyNH(argF->cacheNode(a));
        */
        entry.writeResultNH(resF->cacheNode(b));
        CT->addEntry();
      }

    private:
      static inline void show(FILE* strm, int ev)    { fprintf(strm, "%d", ev); }
      static inline void show(FILE* strm, float ev)  { fprintf(strm, "%f", ev); }
  };

};  // namespace MEDDLY

template <typename TYPE, class OP>
MEDDLY::node_handle  MEDDLY::copy_EV2MT<TYPE,OP>
::computeSkip(int in, TYPE ev, node_handle a)
{
  // Check terminals
  if (argF->isTerminalNode(a)) {
    return resF->handleForValue(ev);
  }

  // Check compute table
  node_handle b;
  compute_table::search_key* Key = inCache(ev, a, b);
  if (0==Key) return b;

  // Initialize sparse node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  node_builder& nb = resF->useSparseBuilder(level, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    TYPE aev;
    A->getEdge(z, aev);
    nb.i(z) = A->i(z);
    nb.d(z) = computeSkip(A->i(z), OP::apply(ev, aev), A->d(z));
  }

  // Cleanup
  node_reader::recycle(A);

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
    return resF->handleForValue(ev);
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(a);

  // Check compute table
  node_handle b;
  compute_table::search_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = inCache(ev, a, b);
    if (0==Key) return b;
  }

  // What's below?
  int nextk;
  if (resF->isForRelations()) {
    nextk = (k>0) ? -k : -k-1;
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  node_reader* A;
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A = argF->initIdentityReader(k, in, a, false);
    } else {
      TYPE rev;
      OP::redundant(rev);
      A = argF->initRedundantReader(k, rev, a, false);
    }
  } else {
    A = argF->initNodeReader(a, false);
  }

  // Initialize node builder
  node_builder& nb = resF->useSparseBuilder(k, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    TYPE aev;
    A->getEdge(z, aev);
    nb.i(z) = A->i(z);
    nb.d(z) = computeAll(A->i(z), nextk, OP::apply(ev, aev), A->d(z));
  }

  // Cleanup
  node_reader::recycle(A);

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
      copy_EV2EV_fast(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc, 1, 1, arg, res)
      {
        // entry[0]: EV node
        // entry[1]: EV node 
      }
      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[0]) ||
          resF->isStale(entryData[1]);
      }
      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[0]);
        resF->uncacheNode(entryData[1]);
      }
      virtual void showEntry(FILE* strm, const node_handle* entryData) const {
        fprintf(strm, "[%s(<?,%d>) <?,%d>]", getName(), entryData[0], entryData[1]);
      }
      virtual void compute(const dd_edge &arg, dd_edge &res) {
        INTYPE av;
        node_handle bn;
        arg.getEdgeValue(av);
        bn = computeSkip(-1, arg.getNode());
        OUTTYPE bv = av;
        res.set(bn, bv);
      }

      node_handle computeSkip(int in, node_handle a);

    protected:
      inline compute_table::search_key* 
      findResult(node_handle a, node_handle &b) 
      {
        compute_table::search_key* CTsrch = useCTkey();
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->reset();
        CTsrch->writeNH(a);
        compute_table::search_result &cacheFind = CT->find(CTsrch);
        if (!cacheFind) return CTsrch;
        b = resF->linkNode(cacheFind.readNH());
        doneCTkey(CTsrch);
        return 0;
      }
      inline node_handle saveResult(compute_table::search_key* Key,
        node_handle a, node_handle b) 
      {
        argF->cacheNode(a);
        compute_table::entry_builder &entry = CT->startNewEntry(Key);
        /*
        entry.writeKeyNH(argF->cacheNode(a));
        */
        entry.writeResultNH(resF->cacheNode(b));
        CT->addEntry();
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
    return expert_forest::bool_Tencoder::value2handle(true);
  }

  // Check compute table
  node_handle b;
  compute_table::search_key* Key = findResult(a, b);
  if (0==Key) return b;

  // Initialize node reader
  node_reader* A = argF->initNodeReader(a, false);

  // Initialize node builder
  const int level = argF->getNodeLevel(a);
  node_builder& nb = resF->useSparseBuilder(level, A->getNNZs());


  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    nb.i(z) = A->i(z);
    nb.d(z) = computeSkip(A->i(z), A->d(z));
    INTYPE av;
    OUTTYPE bv;
    A->getEdge(z, av);
    bv = av;
    nb.setEdge(z, bv);
  }

  // Cleanup
  node_reader::recycle(A);

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
      copy_EV2EV_slow(const unary_opname* oc, expert_forest* arg, 
        expert_forest* res) : unary_operation(oc, 2, 2, arg, res)
      {
        // entry[0]: EV value
        // entry[1]: EV node
        // entry[2]: EV value
        // entry[3]: EV node 
      }

      virtual bool isStaleEntry(const node_handle* entryData) {
        return 
          argF->isStale(entryData[1]) ||
          resF->isStale(entryData[3]);
      }
      virtual void discardEntry(const node_handle* entryData) {
        argF->uncacheNode(entryData[1]);
        resF->uncacheNode(entryData[3]);
      }
      virtual void showEntry(FILE* strm, const node_handle* entryData) const {
        fprintf(strm, "[%s(<", getName());
        INTYPE ev1;
        compute_table::readEV(entryData, ev1);
        show(strm, ev1);
        fprintf(strm, ",%d> <", entryData[1]);
        OUTTYPE ev2; 
        compute_table::readEV(entryData+2, ev2);
        show(strm, ev2);
        fprintf(strm, ",%d>]", entryData[3]);
      }
      virtual void compute(const dd_edge &arg, dd_edge &res) {
        INTYPE av;
        node_handle an, bn;
        OUTTYPE bv;
        arg.getEdgeValue(av);
        an = arg.getNode();
        // need to visit every level...
        computeAll(-1, resF->getNumVariables(), av, an, bv, bn);
        res.set(bn, bv);
      }
      void computeAll(int in, int k, INTYPE av, node_handle an,
        OUTTYPE &bv, node_handle &bn);

    protected:
      inline compute_table::search_key* 
      inCache(INTYPE av, node_handle an, OUTTYPE &bv, node_handle &bn) 
      {
        compute_table::search_key* CTsrch = useCTkey();
        MEDDLY_DCASSERT(CTsrch);
        CTsrch->reset();
        CTsrch->write(av);
        CTsrch->writeNH(an);
        compute_table::search_result &cacheFind = CT->find(CTsrch);
        if (cacheFind) {
          cacheFind.read(bv);
          bn = resF->linkNode(cacheFind.readNH());
          doneCTkey(CTsrch);
          return 0;
        }
        return CTsrch;
      }

      inline void addToCache(compute_table::search_key* Key, 
        INTYPE av, node_handle an, OUTTYPE bv, node_handle &bn) 
      {
        argF->cacheNode(an);
        compute_table::entry_builder &entry = CT->startNewEntry(Key);
        /*
        entry.writeKey(av);
        entry.writeKey(argF->cacheNode(an));
        */
        entry.writeResult(bv);
        entry.writeResultNH(resF->cacheNode(bn));
        CT->addEntry();
      }

    private:
      static inline void show(FILE* strm, int ev)    { fprintf(strm, "%d", ev); }
      static inline void show(FILE* strm, float ev)  { fprintf(strm, "%f", ev); }
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
    bn = expert_forest::bool_Tencoder::value2handle(true);
    return;
  }

  // Get level number
  const int aLevel = argF->getNodeLevel(an);

  // Check compute table
  compute_table::search_key* Key = 0;
  if (k == aLevel && k>0) {
    Key = inCache(av, an, bv, bn);
    if (0==Key) return;
  }

  // What's below?
  int nextk;
  if (resF->isForRelations()) {
    nextk = (k>0) ? -k : -k-1;
  } else {
    nextk = k-1;
  }

  // Initialize node reader
  node_reader* A;
  if (isLevelAbove(k, aLevel)) {
    if (k<0 && argF->isIdentityReduced()) {
      A = argF->initIdentityReader(k, in, an, false);
    } else {
      INTYPE rev;
      INOP::redundant(rev);
      A = argF->initRedundantReader(k, rev, an, false);
    }
  } else {
    A = argF->initNodeReader(an, false);
  }

  // Initialize node builder
  node_builder& nb = resF->useSparseBuilder(k, A->getNNZs());

  // recurse
  for (int z=0; z<A->getNNZs(); z++) {
    INTYPE adv;
    OUTTYPE bdv;
    A->getEdge(z, adv);
    nb.i(z) = A->i(z);
    computeAll(A->i(z), nextk, INOP::apply(av, adv), A->d(z), bdv, nb.d(z));
    nb.setEdge(z, bdv);
  }

  // Cleanup
  node_reader::recycle(A);

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
    throw error(error::DOMAIN_MISMATCH);

  if (arg->isForRelations() != res->isForRelations())
    throw error(error::TYPE_MISMATCH);

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
        throw error(error::NOT_IMPLEMENTED);
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
      case forest::INTEGER:
        return new copy_MT2EV<int>(this, arg, res);

      case forest::REAL:
        return new copy_MT2EV<float>(this, arg, res);

      default:
        throw error(error::TYPE_MISMATCH);
    };
  }

  if (res->isMultiTerminal() && (arg->isEVPlus() || arg->isIndexSet()))
  {
    //
    // EV+ to MT conversion
    //
    switch (arg->getRangeType()) {
      case forest::INTEGER:
        return new copy_EV2MT<int,PLUS>(this, arg, res);

      case forest::REAL:
        return new copy_EV2MT<float,PLUS>(this, arg, res);

      default:
        throw error(error::TYPE_MISMATCH);
    };
  }

  if (res->isMultiTerminal() && arg->isEVTimes())
  {
    //
    // EV* to MT conversion  (untested!)
    //
    switch (arg->getRangeType()) {
      case forest::INTEGER:
        return new copy_EV2MT<int,TIMES>(this, arg, res);

      case forest::REAL:
        return new copy_EV2MT<float,TIMES>(this, arg, res);

      default:
        throw error(error::TYPE_MISMATCH);
    };
  }

  //
  // That's it for any MT arguments
  //
  if (res->isMultiTerminal() || arg->isMultiTerminal()) {
    throw error(error::NOT_IMPLEMENTED);
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
                    return new copy_EV2EV_fast<int,int>(this, arg, res);
                case forest::REAL:
                    return new copy_EV2EV_fast<int,float>(this, arg, res);
                default:
                    throw error(error::TYPE_MISMATCH);
            };
            break;    // in case anything falls through

        case forest::REAL:
            switch (res->getRangeType()) {
                case forest::INTEGER:
                    break;    // not safe to go from real -> integer this way
                case forest::REAL:
                    return new copy_EV2EV_fast<float,float>(this, arg, res);
                default:
                    throw error(error::TYPE_MISMATCH);
            };
            break;    // things may fall through

        default:
            throw error(error::TYPE_MISMATCH);
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
            case forest::INTEGER:
                return new copy_EV2EV_slow<int,PLUS,int>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<int,PLUS,float>(this, arg, res);
            default:
                throw error(error::TYPE_MISMATCH);
          };
        
      case forest::REAL:
          switch (res->getRangeType()) {
            case forest::INTEGER:
                return new copy_EV2EV_slow<float,PLUS,int>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<float,PLUS,float>(this, arg, res);
            default:
                throw error(error::TYPE_MISMATCH);
          };

      default:
          throw error(error::TYPE_MISMATCH);
    }
  }

  //
  // Generic EV* to EV, may be slow
  //
  if (arg->isEVTimes()) {

    switch (arg->getRangeType()) {

      case forest::INTEGER:
          switch (res->getRangeType()) {
            case forest::INTEGER:
                return new copy_EV2EV_slow<int,TIMES,int>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<int,TIMES,float>(this, arg, res);
            default:
                throw error(error::TYPE_MISMATCH);
          };
        
      case forest::REAL:
          switch (res->getRangeType()) {
            case forest::INTEGER:
                return new copy_EV2EV_slow<float,TIMES,int>(this, arg, res);
            case forest::REAL:
                return new copy_EV2EV_slow<float,TIMES,float>(this, arg, res);
            default:
                throw error(error::TYPE_MISMATCH);
          };

      default:
          throw error(error::TYPE_MISMATCH);
    }
  }

  //
  // Catch all for any other cases
  //
  throw error(error::NOT_IMPLEMENTED);

}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeCopy(const settings &s)
{
  return new copy_opname;
}

