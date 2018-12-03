
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
#include "node_edge_count.h"

// #define DEBUG_CARD

namespace MEDDLY {
  class node_count_int;

  class node_count_opname;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        Node Count  operations                  *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                         node_count_int class                   *
// *                                                                *
// ******************************************************************

//  Abstract base class: node_count that returns an integer
class MEDDLY::node_count_int : public unary_operation {
public:
  node_count_int(const unary_opname* oc, expert_forest* arg);

  // common
  virtual bool isStaleEntry(const node_handle* entryData);
  virtual void discardEntry(const node_handle* entryData);
  virtual void showEntry(output &strm, const node_handle* entryData) const;

  virtual void compute(const dd_edge &arg, long &res) {
    res = compute_r(arg.getNode());
  }
  virtual void compute(int k, node_handle a, long &res) {
    res = compute_r(a);
  }
  long compute_r(node_handle a);
};

MEDDLY::node_count_int::node_count_int(const unary_opname* oc, expert_forest* arg)
#ifdef OLD_OP_CT
 : unary_operation(oc, 1, sizeof(long) / sizeof(node_handle), arg, INTEGER)
 {
 }

bool MEDDLY::node_count_int::isStaleEntry(const node_handle* data)
{
  return argF->isStale(data[0]);
}

void MEDDLY::node_count_int::discardEntry(const node_handle* data)
{
  argF->uncacheNode(data[0]);
}

void MEDDLY::node_count_int::showEntry(output &strm, const node_handle* data) const
{
  long answer;
  memcpy(&answer, data+1, sizeof(long));
  strm  << "[" << getName() << "(" << long(data[0]) 
        << "): " << answer << "(L)]";
}
#else
 : unary_operation(oc, 1, arg, INTEGER)
{
  compute_table::entry_type* et = new compute_table::entry_type(oc->getName(), "N:I");
  et->setForestForSlot(0, arg);
  registerEntryType(0, et);
  buildCTs();
}
#endif


long MEDDLY::node_count_int::compute_r(node_handle a)
{
  // Terminal cases
  if (argF->isTerminalNode(a)) {
    return a == 0? 0: 1;
  }

  // Check compute table
#ifdef OLD_OP_CT
  compute_table::search_key* CTsrch = useCTkey();
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->reset();
  CTsrch->writeNH(a); 
  compute_table::search_result &cacheEntry = CT->find(CTsrch);
  if (cacheEntry) {
    long answer;
    cacheEntry.read(answer);
    doneCTkey(CTsrch);
    return answer;
  }
#else
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (CTresult[0]) {
    long ans = CTresult[0].readI();
#ifdef DEBUG_MXD_COMPL
    fprintf(stderr, "\tin CT:   node_count(%ld) : %d\n", a, ans);
#endif
    CT0->recycle(CTsrch);
    return ans;
  }
#endif

  long node_count = argF->getNodeCount(a);

  // Add entry to compute table
#ifdef OLD_OP_CT
  argF->cacheNode(a);
  compute_table::entry_builder &entry = CT->startNewEntry(CTsrch);
  entry.writeResult(node_count);
  CT->addEntry();
#else
  CTresult[0].reset();
  CTresult[0].writeI(node_count);
  CT0->addEntry(CTsrch, CTresult[0]);
#endif

#ifdef DEBUG_CARD
  fprintf(stderr, "Node Count of node %d is %ld(L)\n", a, node_count);
#endif
  return node_count;
}


// ******************************************************************
// *                                                                *
// *                 node_count_opname  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::node_count_opname : public unary_opname {
  public:
    node_count_opname();
    virtual unary_operation* 
      buildOperation(expert_forest* ar, opnd_type res) const;
};

MEDDLY::node_count_opname::node_count_opname()
 : unary_opname("Node Count")
{
}

MEDDLY::unary_operation* 
MEDDLY::node_count_opname::buildOperation(expert_forest* arg, opnd_type res) const
{
  if (0==arg) return 0;
  switch (res) {
    case INTEGER:
      return new node_count_int(this, arg);

    default:
      throw error(error::TYPE_MISMATCH);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeNodeCount()
{
  return new node_count_opname;
}

#if 0
MEDDLY::unary_opname* MEDDLY::initializeEdgeCount()
{
  return new node_count_opname;
}

MEDDLY::unary_opname* MEDDLY::initializeNodeEdgeCount()
{
  return new node_count_opname;
}
#endif
