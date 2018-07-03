
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
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle *entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(node_handle a, node_handle b);
    virtual node_handle compute_normal(node_handle a, node_handle b);
    virtual node_handle compute_ext(node_handle a, node_handle b);

  protected:
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif

    inline compute_table::entry_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
      compute_table::entry_key* CTsrch = CT->useEntryKey(this);
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeNH(b);
        CTsrch->writeNH(a);
      } else {
        CTsrch->writeNH(a);
        CTsrch->writeNH(b);
      }
      compute_table::entry_result& cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readNH());
      CT->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* K, 
      node_handle a, node_handle b, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      static compute_table::entry_result result(1);
      result.reset();
      result.writeN(resF->cacheNode(c));
      CT->addEntry(K, result);
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
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle *entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(node_handle a, node_handle b);
    virtual node_handle compute_normal(node_handle a, node_handle b);
    virtual node_handle compute_ext(node_handle a, node_handle b);

  protected:
    node_handle compute_r(int i, int k, node_handle a, node_handle b);
    node_handle compute_r_normal(int i, int k, node_handle a, node_handle b);
    node_handle compute_r_ext(int i, int k, node_handle a, node_handle b);

  protected:
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif

    inline compute_table::entry_key* 
    findResult(node_handle a, node_handle b, node_handle &c) 
    {
      compute_table::entry_key* CTsrch = CT->useEntryKey(this);
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeNH(b);
        CTsrch->writeNH(a);
      } else {
        CTsrch->writeNH(a);
        CTsrch->writeNH(b);
      }
      compute_table::entry_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readNH());
      CT->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* Key, 
      node_handle a, node_handle b, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      static compute_table::entry_result result(1);
      result.reset();
      result.writeN(resF->cacheNode(c));
      CT->addEntry(Key, result);
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
    virtual void discardEntry(const node_handle* entryData);
    virtual void showEntry(output &strm, const node_handle *entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(int level, node_handle a, node_handle b);

  protected:
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
#endif

    inline compute_table::entry_key* 
    findResult(int k, node_handle a, node_handle b, node_handle &c) 
    {
      compute_table::entry_key* CTsrch = CT->useEntryKey(this);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->write(k);
      if (can_commute && a > b) {
        CTsrch->writeNH(b);
        CTsrch->writeNH(a);
      } else {
        CTsrch->writeNH(a);
        CTsrch->writeNH(b);
      }
      compute_table::entry_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      c = resF->linkNode(cacheFind.readNH());
      CT->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* Key,
      int k, node_handle a, node_handle b, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      static compute_table::entry_result result(1);
      result.reset();
      result.writeN(resF->cacheNode(c));
      CT->addEntry(Key, result);
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

  public:
    virtual void discardEntry(const node_handle* entryData);

  protected:
#ifndef USE_NODE_STATUS
    virtual bool isStaleEntry(const node_handle* entryData);
#else
    virtual MEDDLY::forest::node_status getStatusOfEntry(const node_handle* entryData);
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
    virtual void showEntry(output &strm, const node_handle *entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    virtual compute_table::entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c)
    {
      compute_table::entry_key* CTsrch = CT->useEntryKey(this);
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->write(bev);
        CTsrch->writeNH(b);
        CTsrch->write(aev);
        CTsrch->writeNH(a);
      } else {
        CTsrch->write(aev);
        CTsrch->writeNH(a);
        CTsrch->write(bev);
        CTsrch->writeNH(b);
      }
      compute_table::entry_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cacheFind.read(cev);
      c = resF->linkNode(cacheFind.readNH());
      CT->recycle(CTsrch);
      return 0;
    }

    virtual void saveResult(compute_table::entry_key* Key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      static compute_table::entry_result result(1 + sizeof(long)/sizeof(node_handle));
      result.reset();
      result.writeL(cev);
      result.writeN(resF->cacheNode(c));
      CT->addEntry(Key, result);
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
    virtual void showEntry(output &strm, const node_handle *entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    void compute_r(int in, int level, long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    virtual compute_table::entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c)
    {
      compute_table::entry_key* CTsrch = CT->useEntryKey(this);
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->write(bev);
        CTsrch->writeNH(b);
        CTsrch->write(aev);
        CTsrch->writeNH(a);
      } else {
        CTsrch->write(aev);
        CTsrch->writeNH(a);
        CTsrch->write(bev);
        CTsrch->writeNH(b);
      }
      compute_table::entry_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cacheFind.read(cev);
      c = resF->linkNode(cacheFind.readNH());
      CT->recycle(CTsrch);
      return 0;
    }

    virtual void saveResult(compute_table::entry_key* Key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      static compute_table::entry_result result(1 + sizeof(long)/sizeof(node_handle));
      result.reset();
      result.writeL(cev);
      result.writeN(resF->cacheNode(c));
      CT->addEntry(Key, result);
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
    virtual void showEntry(output &strm, const node_handle *entryData) const;
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(float aev, node_handle a, float bev, node_handle b, 
      float& cev, node_handle &c);

    virtual void compute_k(int in, int k, float aev, node_handle a,
      float bev, node_handle b, float& cev, node_handle& c);

  protected:
    inline compute_table::entry_key* findResult(float aev, node_handle a, 
      float bev, node_handle b, float& cev, node_handle &c) 
    {
      compute_table::entry_key* CTsrch = CT->useEntryKey(this);
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->write(bev);
        CTsrch->writeNH(b);
        CTsrch->write(aev);
        CTsrch->writeNH(a);
      } else {
        CTsrch->write(aev);
        CTsrch->writeNH(a);
        CTsrch->write(bev);
        CTsrch->writeNH(b);
      }
      compute_table::entry_result &cacheFind = CT->find(CTsrch);
      if (!cacheFind) return CTsrch;
      cacheFind.read(cev);
      c = resF->linkNode(cacheFind.readNH());
      CT->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(compute_table::entry_key* Key, float aev, 
      node_handle a, float bev, node_handle b, float cev, node_handle c) 
    {
      arg1F->cacheNode(a);
      arg2F->cacheNode(b);
      static compute_table::entry_result result(2);
      result.reset();
      result.writeF(cev);
      result.writeN(resF->cacheNode(c));
      CT->addEntry(Key, result);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(float aev, node_handle a,
      float bev, node_handle b, float &cev, node_handle &c) = 0;
};


#endif

