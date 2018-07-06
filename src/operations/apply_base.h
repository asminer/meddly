
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

#ifndef APPLY_BASE_H
#define APPLY_BASE_H

/*
    Useful base classes for binary apply operations.
*/

namespace MEDDLY {
  class generic_binary_mdd;
  class generic_binary_mxd;
  class generic_binbylevel_mxd;
  class generic_binary_ev;
  class generic_binary_evplus;
  class generic_binary_evplus_mxd;
  class generic_binary_evtimes;
}

// ******************************************************************

class MEDDLY::generic_binary_mdd : public binary_operation {
  public:
    generic_binary_mdd(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_mdd();

  public:
#ifdef OLD_OP_CT
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle *entryData) const;
#endif
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(node_handle a, node_handle b);
    virtual node_handle compute_normal(node_handle a, node_handle b);
    virtual node_handle compute_ext(node_handle a, node_handle b);

  protected:
#ifdef OLD_OP_CT
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif
#endif

    inline compute_table::entry_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeN(b);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
      }
      compute_table::entry_result& cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* K, 
      node_handle a, node_handle b, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
#ifdef OLD_OP_CT
      static compute_table::entry_result result(1);
#else
      static compute_table::entry_result result(etype[0]);
#endif
      result.reset();
      result.writeN(resF->cacheNode(c));
      CT0->addEntry(K, result);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c) = 0;
    
};

// ******************************************************************

class MEDDLY::generic_binary_mxd : public binary_operation {
  public:
    generic_binary_mxd(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_mxd();

  public:
#ifdef OLD_OP_CT
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle *entryData) const;
#endif
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(node_handle a, node_handle b);
    virtual node_handle compute_normal(node_handle a, node_handle b);
    virtual node_handle compute_ext(node_handle a, node_handle b);

  protected:
    node_handle compute_r(int i, int k, node_handle a, node_handle b);
    node_handle compute_r_normal(int i, int k, node_handle a, node_handle b);
    node_handle compute_r_ext(int i, int k, node_handle a, node_handle b);

  protected:
#ifdef OLD_OP_CT
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif
#endif

    inline compute_table::entry_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeN(b);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
      }
      compute_table::entry_result &cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* Key, 
      node_handle a, node_handle b, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
#ifdef OLD_OP_CT
      static compute_table::entry_result result(1);
#else
      static compute_table::entry_result result(etype[0]);
#endif
      result.reset();
      result.writeN(resF->cacheNode(c));
      CT0->addEntry(Key, result);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c) = 0;
};

// ******************************************************************

class MEDDLY::generic_binbylevel_mxd : public binary_operation {
  public:
    generic_binbylevel_mxd(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binbylevel_mxd();

  public:
#ifdef OLD_OP_CT
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle *entryData) const;
#endif
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(int level, node_handle a, node_handle b);

  protected:
#ifdef OLD_OP_CT
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif
#endif

    inline compute_table::entry_key* 
    findResult(int k, node_handle a, node_handle b, node_handle &c) 
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeI(k);
      if (can_commute && a > b) {
        CTsrch->writeN(b);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
      }
      compute_table::entry_result &cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* Key,
      int k, node_handle a, node_handle b, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
#ifdef OLD_OP_CT
      static compute_table::entry_result result(1);
#else
      static compute_table::entry_result result(etype[0]);
#endif
      result.reset();
      result.writeN(resF->cacheNode(c));
      CT0->addEntry(Key, result);
    }

    node_handle compute_r(int i, int level, node_handle a, node_handle b);

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c) = 0;
};
    

// ******************************************************************

class MEDDLY::generic_binary_ev : public binary_operation {
  public:
    generic_binary_ev(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_ev();

#ifdef OLD_OP_CT
  public:
    virtual void discardEntry(const node_handle* entryData);

  protected:
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif
#endif
};

// ******************************************************************

class MEDDLY::generic_binary_evplus : public generic_binary_ev {
  public:
    generic_binary_evplus(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evplus();

  public:
#ifdef OLD_OP_CT
    virtual void showEntry(output &strm, const node_handle *entryData) const;
#endif
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    virtual compute_table::entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c)
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeL(bev);
        CTsrch->writeN(b);
        CTsrch->writeL(aev);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeL(aev);
        CTsrch->writeN(a);
        CTsrch->writeL(bev);
        CTsrch->writeN(b);
      }
      compute_table::entry_result &cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cev = cacheFind.readL();
      c = resF->linkNode(cacheFind.readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    virtual void saveResult(compute_table::entry_key* Key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
#ifdef OLD_OP_CT
      static compute_table::entry_result result(1 + sizeof(long)/sizeof(node_handle));
#else
      static compute_table::entry_result result(etype[0]);
#endif
      result.reset();
      result.writeL(cev);
      result.writeN(resF->cacheNode(c));
      CT0->addEntry(Key, result);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long &cev, node_handle &c) = 0;
};

// ******************************************************************

class MEDDLY::generic_binary_evplus_mxd : public generic_binary_ev {
  public:
    generic_binary_evplus_mxd(const binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evplus_mxd();

  public:
#ifdef OLD_OP_CT
    virtual void showEntry(output &strm, const node_handle *entryData) const;
#endif
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    void compute_r(int in, int level, long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    virtual compute_table::entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c)
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeL(bev);
        CTsrch->writeN(b);
        CTsrch->writeL(aev);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeL(aev);
        CTsrch->writeN(a);
        CTsrch->writeL(bev);
        CTsrch->writeN(b);
      }
      compute_table::entry_result &cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cev = cacheFind.readL();
      c = resF->linkNode(cacheFind.readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    virtual void saveResult(compute_table::entry_key* Key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
#ifdef OLD_OP_CT
      static compute_table::entry_result result(1 + sizeof(long)/sizeof(node_handle));
#else
      static compute_table::entry_result result(etype[0]);
#endif
      result.reset();
      result.writeL(cev);
      result.writeN(resF->cacheNode(c));
      CT0->addEntry(Key, result);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(long aev, node_handle a, long bev, node_handle b,
      long &cev, node_handle &c) = 0;
};

// ******************************************************************

class MEDDLY::generic_binary_evtimes : public generic_binary_ev {
  public:
    generic_binary_evtimes(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evtimes();

  public:
#ifdef OLD_OP_CT
    virtual void showEntry(output &strm, const node_handle *entryData) const;
#endif
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(float aev, node_handle a, float bev, node_handle b, 
      float& cev, node_handle &c);

    virtual void compute_k(int in, int k, float aev, node_handle a,
      float bev, node_handle b, float& cev, node_handle& c);

  protected:
    inline compute_table::entry_key* findResult(float aev, node_handle a, 
      float bev, node_handle b, float& cev, node_handle &c) 
    {
#ifdef OLD_OP_CT
      compute_table::entry_key* CTsrch = CT0->useEntryKey(this);
#else
      compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
#endif
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeF(bev);
        CTsrch->writeN(b);
        CTsrch->writeF(aev);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeF(aev);
        CTsrch->writeN(a);
        CTsrch->writeF(bev);
        CTsrch->writeN(b);
      }
      compute_table::entry_result &cacheFind = CT0->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cev = cacheFind.readF();
      c = resF->linkNode(cacheFind.readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* Key, float aev, 
      node_handle a, float bev, node_handle b, float cev, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
#ifdef OLD_OP_CT
      static compute_table::entry_result result(2);
#else
      static compute_table::entry_result result(etype[0]);
#endif
      result.reset();
      result.writeF(cev);
      result.writeN(resF->cacheNode(c));
      CT0->addEntry(Key, result);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(float aev, node_handle a,
      float bev, node_handle b, float &cev, node_handle &c) = 0;
};


#endif

