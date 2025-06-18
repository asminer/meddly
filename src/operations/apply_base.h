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
    class generic_binary_ev;
    class generic_binary_evplus;
}

// ******************************************************************

class MEDDLY::generic_binary_ev : public binary_operation {
  public:
    generic_binary_ev(binary_list& code, forest* arg1,
      forest* arg2, forest* res);

  protected:
    virtual ~generic_binary_ev();

};

// ******************************************************************

class MEDDLY::generic_binary_evplus : public generic_binary_ev {
  public:
    generic_binary_evplus(binary_list& code, forest* arg1,
      forest* arg2, forest* res);

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
      if (canCommute() && a > b) {
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

#endif

