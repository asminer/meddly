
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

#ifndef MINIMUM_WITNESS_H
#define MINIMUM_WITNESS_H

#include "meddly.h"
#include "meddly_expert.h"

namespace MEDDLY {
  class common_constraint;

  class constraint_bfs_opname;
  class constraint_bckwd_bfs;

  class constraint_dfs_opname;
  class constraint_bckwd_dfs;

  class constraint_sat_opname;
  class constraint_saturation;

  minimum_witness_opname* initConstraintDFSBackward();
}

class MEDDLY::common_constraint: public specialized_operation
{
protected:
  expert_forest* consF;
  expert_forest* argF;
  expert_forest* transF;
  expert_forest* resF;

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

public:
  common_constraint(const minimum_witness_opname* code, int kl, int al,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);
  ~common_constraint();

  virtual dd_edge compute(const dd_edge& a, const dd_edge& b, const dd_edge& r) = 0;
};

class MEDDLY::constraint_bfs_opname : public minimum_witness_opname {
protected:
  bool forward;
public:
  constraint_bfs_opname(bool fwd);
  virtual specialized_operation* buildOperation(expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res) const;

  virtual specialized_operation* buildOperation(arguments* a) const
  {
    throw error::NOT_IMPLEMENTED;
  }
};

class MEDDLY::constraint_bckwd_bfs: public common_constraint
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
  constraint_bckwd_bfs(const minimum_witness_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual dd_edge compute(const dd_edge& a, const dd_edge& b, const dd_edge& r);
};

class MEDDLY::constraint_dfs_opname : public minimum_witness_opname {
protected:
  bool forward;
public:
  constraint_dfs_opname(bool fwd);
  virtual specialized_operation* buildOperation(expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res) const;

  virtual specialized_operation* buildOperation(arguments* a) const
  {
    throw error::NOT_IMPLEMENTED;
  }
};

class MEDDLY::constraint_bckwd_dfs: public common_constraint
{
protected:
  static const int NODE_INDICES_IN_KEY[4];

  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* plusOp;
  binary_operation* minOp;

  node_handle* splits;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, node_handle c, long& dev, node_handle& d);

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
  constraint_bckwd_dfs(const minimum_witness_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual dd_edge compute(const dd_edge& a, const dd_edge& b, const dd_edge& r);
  void compute(int aev, node_handle a, int bev, node_handle b, node_handle r, long& cev, node_handle& c);

  void saturateHelper(long aev, node_handle a, unpacked_node& nb);
};

class MEDDLY::constraint_saturation: public specialized_operation
{
protected:
  int NODE_INDICES_IN_KEY[3];

  constraint_bckwd_dfs* parent;

  expert_forest* consF;
  expert_forest* argF;
  expert_forest* resF;

  virtual ~constraint_saturation();

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
  constraint_saturation(constraint_bckwd_dfs* p,
    expert_forest* cons, expert_forest* arg, expert_forest* res);

  bool matches(const expert_forest* arg1, const expert_forest* arg2,
    const expert_forest* res) const;

  // high-level front-end
  void saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  void saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c);
};

#endif
