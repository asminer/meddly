
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

#ifndef TRANSITIVE_CLOSURE_H
#define TRANSITIVE_CLOSURE_H

#include "meddly.h"
#include "meddly_expert.h"

namespace MEDDLY {
  class common_transitive_closure;

  class transitive_closure_bfs_opname;
  class transitive_closure_forwd_bfs;

  class transitive_closure_dfs_opname;
  class transitive_closure_dfs;
  class transitive_closure_forwd_dfs;

  class transitive_closure_evplus;

  constrained_opname* initTransitiveClosureDFS();
}

class MEDDLY::common_transitive_closure: public specialized_operation
{
protected:
  expert_forest* consF;
  expert_forest* tcF;
  expert_forest* transF;
  expert_forest* resF;

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

public:
  common_transitive_closure(const constrained_opname* code, unsigned slots,
    expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res);
  ~common_transitive_closure();
};

class MEDDLY::transitive_closure_bfs_opname : public constrained_opname {
public:
  transitive_closure_bfs_opname();
  virtual specialized_operation* buildOperation(expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res) const;

  virtual specialized_operation* buildOperation(arguments* a) const
  {
    throw error::NOT_IMPLEMENTED;
  }
};


class MEDDLY::transitive_closure_forwd_bfs: public common_transitive_closure
{
protected:
  binary_operation* imageOp;
  binary_operation* plusOp;
  binary_operation* minOp;

  void iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);

public:
  transitive_closure_forwd_bfs(const constrained_opname* code,
    expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res);

  virtual void compute(const dd_edge &a, const dd_edge &b, const dd_edge &r, dd_edge &res);
};

class MEDDLY::transitive_closure_dfs_opname : public constrained_opname {
public:
  transitive_closure_dfs_opname();

  virtual specialized_operation* buildOperation(arguments* a) const;
};

class MEDDLY::transitive_closure_dfs: public common_transitive_closure
{
protected:
  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* minOp;

  node_handle* splits;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, node_handle c, long& dev, node_handle& d);

  compute_table::entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle c, long& dev, node_handle& d);
  void saveResult(compute_table::entry_key* key,
    long aev, node_handle a, long bev, node_handle b, node_handle c, long dev, node_handle d);

  void splitMxd(node_handle mxd);
  virtual void recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c) = 0;

public:
  transitive_closure_dfs(const constrained_opname* code,
    expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res);

  virtual void compute(const dd_edge &a, const dd_edge &b, const dd_edge &r, dd_edge &res);
  void _compute(int aev, node_handle a, int bev, node_handle b, node_handle r, long& cev, node_handle& c);

  virtual void saturateHelper(long aev, node_handle a, int in, unpacked_node& nb) = 0;
};

class MEDDLY::transitive_closure_forwd_dfs: public transitive_closure_dfs
{
protected:
  virtual void recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);

public:
  transitive_closure_forwd_dfs(const constrained_opname* code,
    expert_forest* cons, expert_forest* tc, expert_forest* trans, expert_forest* res);

  virtual void saturateHelper(long aev, node_handle a, int in, unpacked_node& nb);
};

class MEDDLY::transitive_closure_evplus: public specialized_operation
{
protected:
  transitive_closure_dfs* parent;

  expert_forest* tcF;
  expert_forest* consF;
  expert_forest* resF;

  virtual ~transitive_closure_evplus();

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  compute_table::entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c);
  void saveResult(compute_table::entry_key* Key,
    long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c);

public:
  transitive_closure_evplus(transitive_closure_dfs* p,
    expert_forest* cons, expert_forest* tc, expert_forest* res);

  // high-level front-end
  void saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  void saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c);
};

#endif
