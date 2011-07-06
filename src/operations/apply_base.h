
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
    bool can_commute;
  protected:
    compute_table::search_key CTsrch;
  public:
    generic_binary_mdd(const binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_mdd();

  public:
    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual int compute(int a, int b);

  protected:
    inline void operationCommutes() {
      can_commute = (arg1F == arg2F);
    }

    inline bool findResult(int a, int b, int &c) {
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
      /*
      static int key[2];
      if (can_commute && a > b) {
        key[0] = b; key[1] = a;
      } else {
        key[0] = a; key[1] = b;
      }
      const int* cacheFind = CT->find(this, key);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[2]);
      return true;
      */
    }

    inline void saveResult(int a, int b, int c) {
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
      /*
      static int cacheEntry[3];
      if (can_commute && a > b) {
        cacheEntry[0] = arg1F->cacheNode(b); 
        cacheEntry[1] = arg2F->cacheNode(a);
      } else {
        cacheEntry[0] = arg1F->cacheNode(a); 
        cacheEntry[1] = arg2F->cacheNode(b);
      }
      cacheEntry[2] = resF->cacheNode(c);
      CT->add(this, cacheEntry);
      */
    }

    // not sure about these...

    virtual void expandA(int a, int b, int result, int resultSize);
    virtual void expandB(int a, int b, int result, int resultSize);
    virtual void fullFull(int a, int b, int result, int resultSize);
    virtual void fullSparse(int a, int b, int result, int resultSize);
    virtual void sparseFull(int a, int b, int result, int resultSize);
    virtual void sparseSparse(int a, int b, int result, int resultSize);

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(int a, int b, int& c) = 0;
    
};

// ******************************************************************

class MEDDLY::generic_binary_mxd : public generic_binary_mdd {
  public:
    generic_binary_mxd(const binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_mxd();

  public:
    virtual int compute(int a, int b);

  protected:
    virtual int computeIdent(int a, int b);
    virtual int computeNonIdent(int a, int b);

    virtual void expandA(int a, int b, int result, 
      int resultLevel, int resultSize);
    virtual void singleExpandA( int a, int b, int result, 
      int resultLevel, int resultSize);
    virtual void expandB(int a, int b, int result, 
      int resultLevel, int resultSize);
    virtual void singleExpandB(int a, int b, int result, 
      int resultLevel, int resultSize);
};

// ******************************************************************

class MEDDLY::generic_binbylevel_mxd : public binary_operation {
    bool can_commute;
  public:
    generic_binbylevel_mxd(const binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binbylevel_mxd();

  public:
    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);
    virtual void showEntry(FILE* strm, const int *entryData) const;
    virtual void compute(const dd_edge& a, const dd_edge& b, dd_edge &c);

    virtual int compute(int level, int a, int b);

  protected:
    inline void operationCommutes() {
      can_commute = (arg1F == arg2F);
    }

    inline bool findResult(int k, int a, int b, int &c) {
      static int key[3];
      key[0] = k;
      if (can_commute && a > b) {
        key[1] = b; key[2] = a;
      } else {
        key[1] = a; key[2] = b;
      }
      const int* cacheFind = CT->find(this, key);
      if (0==cacheFind) return false;
      c = resF->linkNode(cacheFind[3]);
      return true;
    }

    inline bool saveResult(int k, int a, int b, int c) {
      static int cacheEntry[4];
      cacheEntry[0] = k;
      if (can_commute && a > b) {
        cacheEntry[1] = arg1F->cacheNode(b); 
        cacheEntry[2] = arg2F->cacheNode(a);
      } else {
        cacheEntry[1] = arg1F->cacheNode(a); 
        cacheEntry[2] = arg2F->cacheNode(b);
      }
      cacheEntry[3] = resF->cacheNode(c);
      CT->add(this, cacheEntry);
    }

    virtual int computeIdent(int level, int a, int b);
    virtual int computeNonIdent(int level, int a, int b);

    virtual void expandSkippedLevel(int a, int b,
        int result, int resultLevel, int resultSize);
    virtual void expandA(int a, int b,
        int result, int resultLevel, int resultSize);
    virtual void singleExpandA(int a, int b,
        int result, int resultLevel, int resultSize);
    virtual void expandB(int a, int b,
        int result, int resultLevel, int resultSize);
    virtual void singleExpandB(int a, int b,
        int result, int resultLevel, int resultSize);

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(int a, int b, int& c) = 0;
};
    

// ******************************************************************

class MEDDLY::generic_binary_ev : public binary_operation {
  protected:
    bool can_commute;
  public:
    generic_binary_ev(const binary_opname* code, expert_forest* arg1,
      expert_forest* arg2, expert_forest* res);

  protected:
    virtual ~generic_binary_ev();

  public:
    virtual bool isEntryStale(const int* entryData);
    virtual void discardEntry(const int* entryData);

  protected:
    inline void operationCommutes() {
      can_commute = (arg1F == arg2F);
    }
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

    virtual void compute(int aev, int a, int bev, int b, int& cev, int &c);

  protected:
    inline bool findResult(int aev, int a, int bev, int b, int& cev, int &c) {
      static int key[4];
      if (can_commute && a > b) {
        key[0] = bev;   key[1] = b;
        key[2] = aev;   key[3] = a;
      } else {
        key[0] = aev;   key[1] = a;
        key[2] = bev;   key[3] = b;
      }
      const int* cacheFind = CT->find(this, key);
      if (0==cacheFind) return false;
      cev = cacheFind[4];
      c = resF->linkNode(cacheFind[5]);
      return true;
    }

    inline bool saveResult(int aev, int a, int bev, int b, int cev, int c) {
      static int cacheEntry[6];
      if (can_commute && a > b) {
        cacheEntry[0] = bev;  cacheEntry[1] = arg1F->cacheNode(b); 
        cacheEntry[2] = aev;  cacheEntry[3] = arg2F->cacheNode(a);
      } else {
        cacheEntry[0] = aev;  cacheEntry[1] = arg2F->cacheNode(a);
        cacheEntry[2] = bev;  cacheEntry[3] = arg1F->cacheNode(b); 
      }
      cacheEntry[4] = cev;  cacheEntry[5] = resF->cacheNode(c);
      CT->add(this, cacheEntry);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(int aev, int a, int bev, int b, 
      int &cev, int &c) = 0;
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

    virtual void compute(float aev, int a, float bev, int b, 
      float& cev, int &c);

  protected:
    inline bool findResult(float aev, int a, float bev, int b, 
      float& cev, int &c) 
    {
      static int key[4];
      if (can_commute && a > b) {
        key[0] = toInt(bev);  key[1] = b;
        key[2] = toInt(aev);  key[3] = a;
      } else {
        key[0] = toInt(aev);  key[1] = a;
        key[2] = toInt(bev);  key[3] = b;
      }
      const int* cacheFind = CT->find(this, key);
      if (0==cacheFind) return false;
      cev = toFloat(cacheFind[4]);
      c = resF->linkNode(cacheFind[5]);
      return true;
    }

    inline bool saveResult(int aev, int a, int bev, int b, int cev, int c) {
      static int cacheEntry[6];
      if (can_commute && a > b) {
        cacheEntry[0] = toInt(bev);   cacheEntry[1] = arg1F->cacheNode(b); 
        cacheEntry[2] = toInt(aev);   cacheEntry[3] = arg2F->cacheNode(a);
      } else {
        cacheEntry[0] = toInt(aev);   cacheEntry[1] = arg2F->cacheNode(a);
        cacheEntry[2] = toInt(bev);   cacheEntry[3] = arg1F->cacheNode(b); 
      }
      cacheEntry[4] = toInt(cev);   cacheEntry[5] = resF->cacheNode(c);
      CT->add(this, cacheEntry);
    }

  protected:
    // If terminal condition is reached, returns true and the result in c.
    // Must be provided in derived classes.
    virtual bool checkTerminals(float aev, int a, float bev, int b, 
      float &cev, int &c) = 0;
};


#endif

