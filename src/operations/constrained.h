
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

#ifndef SAT_CONSTRAINED_H
#define SAT_CONSTRAINED_H

#include "meddly.h"
#include "meddly_expert.h"

namespace MEDDLY {
  class common_constrained;

  class constrained_bfs_opname;
  class constrained_bckwd_bfs;

  class constrained_dfs_opname;
  class constrained_bckwd_dfs;

  //class constraint_sat_opname;
  class constrained_saturation;

  constrained_opname* initConstrainedBFSBackward();
  constrained_opname* initConstrainedDFSBackward();
}

class MEDDLY::common_constrained: public specialized_operation
{
protected:
  expert_forest* consF;
  expert_forest* argF;
  expert_forest* transF;
  expert_forest* resF;

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

public:
  common_constrained(const constrained_opname* code, int kl, int al,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);
  ~common_constrained();
};

class MEDDLY::constrained_bfs_opname : public constrained_opname {
protected:
  bool forward;

public:
  constrained_bfs_opname(bool fwd);

  virtual specialized_operation* buildOperation(arguments* a) const;
};

class MEDDLY::constrained_bckwd_bfs: public common_constrained
{
protected:
  binary_operation* imageOp;
  binary_operation* plusOp;
  binary_operation* minOp;

  void iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);

  virtual bool isStaleEntry(const node_handle* entryData);
  virtual void discardEntry(const node_handle* entryData);
  virtual void showEntry(output &strm, const node_handle* entryData) const;

public:
  constrained_bckwd_bfs(const constrained_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual void compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res);
};

class MEDDLY::constrained_dfs_opname : public constrained_opname {
protected:
  bool forward;

public:
  constrained_dfs_opname(bool fwd);

  virtual specialized_operation* buildOperation(arguments* a) const;
};

class MEDDLY::constrained_bckwd_dfs: public common_constrained
{
protected:
  static const int NODE_INDICES_IN_KEY[4];

  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* minOp;

  node_handle* splits;

  compute_table::search_key* findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle c, long& dev, node_handle& d);
  void saveResult(compute_table::search_key* key,
    long aev, node_handle a, long bev, node_handle b, node_handle c, long dev, node_handle d);

  void splitMxd(node_handle mxd);
  void recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);

  virtual bool isStaleEntry(const node_handle* data);
  virtual void discardEntry(const node_handle* data);
  virtual void showEntry(output &strm, const node_handle *data) const;

public:
  constrained_bckwd_dfs(const constrained_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual void compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res);
  void compute(int aev, node_handle a, int bev, node_handle b, node_handle r, long& cev, node_handle& c);

  void saturateHelper(long aev, node_handle a, unpacked_node& nb);
};

class MEDDLY::constrained_saturation: public specialized_operation
{
protected:
  int NODE_INDICES_IN_KEY[3];

  constrained_bckwd_dfs* parent;

  expert_forest* consF;
  expert_forest* argF;
  expert_forest* resF;

  virtual ~constrained_saturation();

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  compute_table::search_key* findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c);
  void saveResult(compute_table::search_key* Key,
    long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c);

  virtual bool isStaleEntry(const node_handle* data);
  virtual void discardEntry(const node_handle* data);
  virtual void showEntry(output &strm, const node_handle *data) const;

public:
  constrained_saturation(constrained_bckwd_dfs* p,
    expert_forest* cons, expert_forest* arg, expert_forest* res);

  bool matches(const expert_forest* arg1, const expert_forest* arg2,
    const expert_forest* res) const;

  // high-level front-end
  void saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  void saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c);
};

#endif
