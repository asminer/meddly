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

#ifndef MEDDLY_APPLY_BASE_H
#define MEDDLY_APPLY_BASE_H

#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../compute_table.h"
#include "../oper_binary.h"

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
    generic_binary_mdd(binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_mdd();

  public:
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

    virtual node_handle compute(node_handle a, node_handle b);
    virtual node_handle compute_normal(node_handle a, node_handle b);
    virtual node_handle compute_ext(node_handle a, node_handle b);

  protected:
    inline ct_entry_key*
    findResult(node_handle a, node_handle b, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeN(b);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
      }
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(ct_entry_key* K,
      node_handle a, node_handle b, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(K, CTresult[0]);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c) = 0;

};

// ******************************************************************

class MEDDLY::generic_binary_mxd : public binary_operation {
  public:
    generic_binary_mxd(binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_mxd();

  public:
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

    virtual node_handle compute(node_handle a, node_handle b);
    virtual node_handle compute_normal(node_handle a, node_handle b);
    virtual node_handle compute_ext(node_handle a, node_handle b);

  protected:
    node_handle compute_r(int i, int k, node_handle a, node_handle b);
    node_handle compute_r_normal(int i, int k, node_handle a, node_handle b);
    node_handle compute_r_ext(int i, int k, node_handle a, node_handle b);

  protected:
    inline ct_entry_key*
    findResult(node_handle a, node_handle b, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->writeN(b);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
      }
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(ct_entry_key* Key,
      node_handle a, node_handle b, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(node_handle a, node_handle b, node_handle& c) = 0;
};

// ******************************************************************

class MEDDLY::generic_binbylevel_mxd : public binary_operation {
  public:
    generic_binbylevel_mxd(binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binbylevel_mxd();

  public:
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

    virtual node_handle compute(int level, node_handle a, node_handle b);

  protected:
    inline ct_entry_key*
    findResult(int k, node_handle a, node_handle b, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->writeI(k);
      if (can_commute && a > b) {
        CTsrch->writeN(b);
        CTsrch->writeN(a);
      } else {
        CTsrch->writeN(a);
        CTsrch->writeN(b);
      }
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(ct_entry_key* Key,
      int k, node_handle a, node_handle b, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
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
    generic_binary_ev(binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_ev();

};

// ******************************************************************

class MEDDLY::generic_binary_evplus : public generic_binary_ev {
  public:
    generic_binary_evplus(binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evplus();

  public:
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    virtual ct_entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
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
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      cev = CTresult[0].readL();
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    virtual void saveResult(ct_entry_key* Key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeL(cev);
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
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
    generic_binary_evplus_mxd(binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evplus_mxd();

  public:
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

    virtual void compute(long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    void compute_r(int in, int level, long aev, node_handle a, long bev, node_handle b, long& cev, node_handle &c);

  protected:
    virtual ct_entry_key* findResult(long aev, node_handle a,
      long bev, node_handle b, long& cev, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
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
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      cev = CTresult[0].readL();
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    virtual void saveResult(ct_entry_key* Key,
      long aev, node_handle a, long bev, node_handle b, long cev, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeL(cev);
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
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
    generic_binary_evtimes(binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evtimes();

  public:
    virtual void computeDDEdge(const dd_edge& a, const dd_edge& b, dd_edge &c, bool userFlag);

    virtual void compute(float aev, node_handle a, float bev, node_handle b,
      float& cev, node_handle &c);

    virtual void compute_k(int in, int k, float aev, node_handle a,
      float bev, node_handle b, float& cev, node_handle& c);

  protected:
    inline ct_entry_key* findResult(float aev, node_handle a,
      float bev, node_handle b, float& cev, node_handle &c)
    {
      ct_entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
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
      CT0->find(CTsrch, CTresult[0]);
      if (!CTresult[0]) return CTsrch;
      cev = CTresult[0].readF();
      c = resF->linkNode(CTresult[0].readN());
      CT0->recycle(CTsrch);
      return 0;
    }

    inline void saveResult(ct_entry_key* Key, float aev,
      node_handle a, float bev, node_handle b, float cev, node_handle c)
    {
      CTresult[0].reset();
      CTresult[0].writeF(cev);
      CTresult[0].writeN(c);
      CT0->addEntry(Key, CTresult[0]);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(float aev, node_handle a,
      float bev, node_handle b, float &cev, node_handle &c) = 0;
};


#endif

