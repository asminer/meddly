
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
    virtual void showEntry(FILE* strm, const node_handle *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(node_handle a, node_handle b);

  protected:
    virtual bool isStaleEntry(const node_handle* entryData);

    inline bool findResult(node_handle a, node_handle b, node_handle &c) {
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->key(0) = b;
        CTsrch->key(1) = a;
      } else {
        CTsrch->key(0) = a;
        CTsrch->key(1) = b;
      }
      const node_handle* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[2]);
      return true;
    }

    inline void saveResult(node_handle a, node_handle b, node_handle c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      if (can_commute && a > b) {
        entry.key(0) = arg2F->cacheNode(b);
        entry.key(1) = arg1F->cacheNode(a);
      } else {
        entry.key(0) = arg1F->cacheNode(a);
        entry.key(1) = arg2F->cacheNode(b);
      }
      entry.result(0) = resF->cacheNode(c);
      CT->addEntry();
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
    virtual void showEntry(FILE* strm, const node_handle *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(node_handle a, node_handle b);

  protected:
    node_handle compute(int i, int k, node_handle a, node_handle b);

  protected:
    virtual bool isStaleEntry(const node_handle* entryData);

    inline bool findResult(int in, node_handle a, node_handle b, node_handle &c) {
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->key(0) = in;
      if (can_commute && a > b) {
        CTsrch->key(1) = b;
        CTsrch->key(2) = a;
      } else {
        CTsrch->key(1) = a;
        CTsrch->key(2) = b;
      }
      const node_handle* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[3]);
      return true;
    }

    inline void saveResult(int in, node_handle a, node_handle b, node_handle c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      entry.key(0) = in;
      if (can_commute && a > b) {
        entry.key(1) = arg2F->cacheNode(b);
        entry.key(2) = arg1F->cacheNode(a);
      } else {
        entry.key(1) = arg1F->cacheNode(a);
        entry.key(2) = arg2F->cacheNode(b);
      }
      entry.result(0) = resF->cacheNode(c);
      CT->addEntry();
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
    virtual void showEntry(FILE* strm, const node_handle *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual node_handle compute(int level, node_handle a, node_handle b);

  protected:
    virtual bool isStaleEntry(const node_handle* entryData);

    inline bool findResult(int k, node_handle a, node_handle b, node_handle &c) {
      MEDDLY_DCASSERT(CTsrch);
      CTsrch->key(0) = k;
      if (can_commute && a > b) {
        CTsrch->key(1) = b;
        CTsrch->key(2) = a;
      } else {
        CTsrch->key(1) = a;
        CTsrch->key(2) = b;
      }
      const node_handle* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[3]);
      return true;
    }

    inline void saveResult(int k, node_handle a, node_handle b, node_handle c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      entry.key(0) = k;
      if (can_commute && a > b) {
        entry.key(1) = arg2F->cacheNode(b);
        entry.key(2) = arg1F->cacheNode(a);
      } else {
        entry.key(1) = arg1F->cacheNode(a);
        entry.key(2) = arg2F->cacheNode(b);
      }
      entry.result(0) = resF->cacheNode(c);
      CT->addEntry();
    }

    node_handle compute(int i, int level, node_handle a, node_handle b);

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
    virtual bool isStaleEntry(const node_handle* entryData);
};

// ******************************************************************

class MEDDLY::generic_binary_evplus : public generic_binary_ev {
  public:
    generic_binary_evplus(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evplus();

  public:
    virtual void showEntry(FILE* strm, const node_handle *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(int aev, node_handle a, int bev, node_handle b, int& cev, node_handle &c);

  protected:
    inline bool findResult(int aev, node_handle a, int bev, node_handle b, int& cev, node_handle &c) {
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->setKeyEV(0, bev);
        CTsrch->key(1) = b;
        CTsrch->setKeyEV(2, aev);
        CTsrch->key(3) = a;
      } else {
        CTsrch->setKeyEV(0, aev);
        CTsrch->key(1) = a;
        CTsrch->setKeyEV(2, bev);
        CTsrch->key(3) = b;
      }
      const node_handle* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      compute_table::readEV(cacheFind+4, cev);
      c = resF->linkNode(cacheFind[5]);
      return true;
    }

    inline void saveResult(int aev, node_handle a, int bev, node_handle b, int cev, node_handle c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      if (can_commute && a > b) {
        entry.setKeyEV(0, bev);
        entry.key(1) = arg2F->cacheNode(b);
        entry.setKeyEV(2, aev);
        entry.key(3) = arg1F->cacheNode(a);
      } else {
        entry.setKeyEV(0, aev);
        entry.key(1) = arg1F->cacheNode(a);
        entry.setKeyEV(2, bev);
        entry.key(3) = arg2F->cacheNode(b);
      }
      entry.setResultEV(0, cev);
      entry.result(1) = resF->cacheNode(c);
      CT->addEntry();
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(int aev, node_handle a, int bev, node_handle b, 
      int &cev, node_handle &c) = 0;
};

// ******************************************************************

class MEDDLY::generic_binary_evtimes : public generic_binary_ev {
  public:
    generic_binary_evtimes(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evtimes();

  public:
    virtual void showEntry(FILE* strm, const node_handle *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(float aev, node_handle a, float bev, node_handle b, 
      float& cev, node_handle &c);

  protected:
    inline bool findResult(float aev, node_handle a, float bev, node_handle b, 
      float& cev, node_handle &c) 
    {
      MEDDLY_DCASSERT(CTsrch);
      if (can_commute && a > b) {
        CTsrch->setKeyEV(0, bev);
        CTsrch->key(1) = b;
        CTsrch->setKeyEV(2, bev);
        CTsrch->key(3) = a;
      } else {
        CTsrch->setKeyEV(0, aev);
        CTsrch->key(1) = a;
        CTsrch->setKeyEV(2, bev);
        CTsrch->key(3) = b;
      }
      const node_handle* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      compute_table::readEV(cacheFind+4, cev);
      c = resF->linkNode(cacheFind[5]);
      return true;
    }

    inline void saveResult(float aev, node_handle a, float bev, node_handle b, float cev, node_handle c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      if (can_commute && a > b) {
        entry.setKeyEV(0, bev);
        entry.key(1) = arg2F->cacheNode(b);
        entry.setKeyEV(2, aev);
        entry.key(3) = arg1F->cacheNode(a);
      } else {
        entry.setKeyEV(0, aev);
        entry.key(1) = arg1F->cacheNode(a);
        entry.setKeyEV(2, bev);
        entry.key(3) = arg2F->cacheNode(b);
      }
      entry.setResultEV(0, cev);
      entry.result(1) = resF->cacheNode(c);
      CT->addEntry();
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(float aev, node_handle a, float bev, node_handle b, 
      float &cev, node_handle &c) = 0;
};


#endif

