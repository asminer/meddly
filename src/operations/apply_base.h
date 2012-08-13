
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
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual long compute(long a, long b);

  protected:
    virtual bool isStaleEntry(const int* entryData);

    inline bool findResult(long a, long b, long &c) {
      if (can_commute && a > b) {
        CTsrch.key(0) = b;
        CTsrch.key(1) = a;
      } else {
        CTsrch.key(0) = a;
        CTsrch.key(1) = b;
      }
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[2]);
      return true;
    }

    inline void saveResult(long a, long b, long c) {
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
    virtual bool checkTerminals(long a, long b, long& c) = 0;
    
};

// ******************************************************************

class MEDDLY::generic_binary_mxd : public binary_operation {
  public:
    generic_binary_mxd(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_mxd();

  public:
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual long compute(long a, long b);

  protected:
    long compute(int i, int k, long a, long b);

  protected:
    virtual bool isStaleEntry(const int* entryData);

    inline bool findResult(int in, long a, long b, long &c) {
      CTsrch.key(0) = in;
      if (can_commute && a > b) {
        CTsrch.key(1) = b;
        CTsrch.key(2) = a;
      } else {
        CTsrch.key(1) = a;
        CTsrch.key(2) = b;
      }
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[3]);
      return true;
    }

    inline void saveResult(int in, long a, long b, long c) {
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
    virtual bool checkTerminals(long a, long b, long& c) = 0;
};

// ******************************************************************

class MEDDLY::generic_binbylevel_mxd : public binary_operation {
  public:
    generic_binbylevel_mxd(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binbylevel_mxd();

  public:
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual long compute(int level, long a, long b);

  protected:
    virtual bool isStaleEntry(const int* entryData);

    inline bool findResult(int k, long a, long b, long &c) {
      CTsrch.key(0) = k;
      if (can_commute && a > b) {
        CTsrch.key(1) = b;
        CTsrch.key(2) = a;
      } else {
        CTsrch.key(1) = a;
        CTsrch.key(2) = b;
      }
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[3]);
      return true;
    }

    inline void saveResult(int k, long a, long b, long c) {
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

    long compute(int i, int level, long a, long b);

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(long a, long b, long& c) = 0;
};
    

// ******************************************************************

class MEDDLY::generic_binary_ev : public binary_operation {
  public:
    generic_binary_ev(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_ev();

  public:
    virtual void discardEntry(const int* entryData);

  protected:
    virtual bool isStaleEntry(const int* entryData);
};

// ******************************************************************

class MEDDLY::generic_binary_evplus : public generic_binary_ev {
  public:
    generic_binary_evplus(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evplus();

  public:
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(int aev, long a, int bev, long b, int& cev, long &c);

  protected:
    inline bool findResult(int aev, long a, int bev, long b, int& cev, long &c) {
      if (can_commute && a > b) {
        CTsrch.key(0) = bev;
        CTsrch.key(1) = b;
        CTsrch.key(2) = aev;
        CTsrch.key(3) = a;
      } else {
        CTsrch.key(0) = aev;
        CTsrch.key(1) = a;
        CTsrch.key(2) = bev;
        CTsrch.key(3) = b;
      }
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      cev = cacheFind[4];
      c = resF->linkNode(cacheFind[5]);
      return true;
    }

    inline bool saveResult(int aev, long a, int bev, long b, int cev, long c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      if (can_commute && a > b) {
        entry.key(0) = bev;
        entry.key(1) = arg2F->cacheNode(b);
        entry.key(2) = aev;
        entry.key(3) = arg1F->cacheNode(a);
      } else {
        entry.key(0) = aev;
        entry.key(1) = arg1F->cacheNode(a);
        entry.key(2) = bev;
        entry.key(3) = arg2F->cacheNode(b);
      }
      entry.result(0) = cev;
      entry.result(1) = resF->cacheNode(c);
      CT->addEntry();
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(int aev, long a, int bev, long b, 
      int &cev, long &c) = 0;
};

// ******************************************************************

class MEDDLY::generic_binary_evtimes : public generic_binary_ev {
  public:
    generic_binary_evtimes(const binary_opname* code, expert_forest* arg1, 
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_evtimes();

  public:
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual void compute(float aev, long a, float bev, long b, 
      float& cev, long &c);

  protected:
    inline bool findResult(float aev, long a, float bev, long b, 
      float& cev, long &c) 
    {
      if (can_commute && a > b) {
        CTsrch.key(0) = toInt(bev);
        CTsrch.key(1) = b;
        CTsrch.key(2) = toInt(aev);
        CTsrch.key(3) = a;
      } else {
        CTsrch.key(0) = toInt(aev);
        CTsrch.key(1) = a;
        CTsrch.key(2) = toInt(bev);
        CTsrch.key(3) = b;
      }
      const int* cacheFind = CT->find(CTsrch);
      if (0==cacheFind) return false;
      cev = toFloat(cacheFind[4]);
      c = resF->linkNode(cacheFind[5]);
      return true;
    }

    inline bool saveResult(float aev, long a, float bev, long b, float cev, long c) {
      compute_table::temp_entry &entry = CT->startNewEntry(this);
      if (can_commute && a > b) {
        entry.key(0) = toInt(bev);
        entry.key(1) = arg2F->cacheNode(b);
        entry.key(2) = toInt(aev);
        entry.key(3) = arg1F->cacheNode(a);
      } else {
        entry.key(0) = toInt(aev);
        entry.key(1) = arg1F->cacheNode(a);
        entry.key(2) = toInt(bev);
        entry.key(3) = arg2F->cacheNode(b);
      }
      entry.result(0) = toInt(cev);
      entry.result(1) = resF->cacheNode(c);
      CT->addEntry();
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(float aev, long a, float bev, long b, 
      float &cev, long &c) = 0;
};


#endif

