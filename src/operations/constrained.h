
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

#ifndef MEDDLY_SAT_CONSTRAINED_H
#define MEDDLY_SAT_CONSTRAINED_H

#include "old_meddly.h"

#include "../opname_satur.h"
#include "../oper_special.h"
#include "../ct_entry_key.h"

namespace MEDDLY {
  class common_constrained;

  class constrained_bfs_opname;
  class constrained_bckwd_bfs_evplus;

  class constrained_dfs_opname;
  class constrained_dfs_mt;
  class constrained_forwd_dfs_mt;
  class constrained_bckwd_dfs_mt;
  class constrained_bckwd_dfs_evplus;

  class constrained_saturation_mt;
  class constrained_saturation_evplus;

  constrained_opname* initConstrainedBFSBackward();
  constrained_opname* initConstrainedDFSForward();
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
  common_constrained(constrained_opname* code, unsigned slots,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);
  ~common_constrained();
};

class MEDDLY::constrained_bfs_opname : public constrained_opname {
protected:
  bool forward;

public:
  constrained_bfs_opname(bool fwd);

  virtual specialized_operation* buildOperation(arguments* a);
};

class MEDDLY::constrained_bckwd_bfs_evplus: public common_constrained
{
protected:
  binary_operation* imageOp;
  binary_operation* plusOp;
  binary_operation* minOp;

  // void iterate(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);
  void iterate(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& c);

public:
  constrained_bckwd_bfs_evplus(constrained_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual void compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res);
};

class MEDDLY::constrained_dfs_opname : public constrained_opname {
protected:
  bool forward;

public:
  constrained_dfs_opname(bool fwd);

  virtual specialized_operation* buildOperation(arguments* a);
};

class MEDDLY::constrained_dfs_mt: public common_constrained
{
protected:
  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* unionOp;

  dd_edge* splits;

  ct_entry_key* findResult(node_handle a, node_handle b, node_handle r, node_handle& c);
  void saveResult(ct_entry_key* key,
    node_handle a, node_handle b, node_handle r, node_handle c);

  void splitMxd(const dd_edge& mxd);

public:
  constrained_dfs_mt(constrained_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual void compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res);
  void _compute(node_handle a, node_handle b, node_handle r, node_handle& c);

  virtual void saturateHelper(node_handle a, unpacked_node& nb) = 0;
};

class MEDDLY::constrained_forwd_dfs_mt: public constrained_dfs_mt
{
protected:
  void recFire(node_handle a, node_handle b, node_handle r, node_handle& c);

public:
  constrained_forwd_dfs_mt(constrained_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual void saturateHelper(node_handle a, unpacked_node& nb) override;
};

class MEDDLY::constrained_bckwd_dfs_mt: public constrained_dfs_mt
{
protected:
  void recFire(node_handle a, node_handle b, node_handle r, node_handle& c);

public:
  constrained_bckwd_dfs_mt(constrained_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual void saturateHelper(node_handle a, unpacked_node& nb) override;
};

class MEDDLY::constrained_saturation_mt: public specialized_operation
{
protected:
  constrained_dfs_mt* parent;

  expert_forest* consF;
  expert_forest* argF;
  expert_forest* resF;

  virtual ~constrained_saturation_mt();

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

  bool checkTerminals(node_handle a, node_handle b, node_handle& c);

  ct_entry_key* findResult(node_handle a, node_handle b, int level, node_handle &c);
  void saveResult(ct_entry_key* Key,
    node_handle a, node_handle b, int level, node_handle c);

public:
  constrained_saturation_mt(constrained_dfs_mt* p,
    expert_forest* cons, expert_forest* arg, expert_forest* res);

  bool matches(const expert_forest* arg1, const expert_forest* arg2,
    const expert_forest* res) const;

  // high-level front-end
  void saturate(node_handle a, node_handle b, node_handle& c);

  void saturate(node_handle a, node_handle b, int level, node_handle& c);
};

class MEDDLY::constrained_bckwd_dfs_evplus: public common_constrained
{
protected:
  binary_operation* mxdDifferenceOp;
  binary_operation* mxdIntersectionOp;
  binary_operation* minOp;

  dd_edge* splits;

  ct_entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, node_handle r, long& dev, node_handle& d);
  void saveResult(ct_entry_key* key,
    long aev, node_handle a, long bev, node_handle b, node_handle r, long dev, node_handle d);

  void splitMxd(const dd_edge& mxd);
  void recFire(long aev, node_handle a, long bev, node_handle b, node_handle r, long& cev, node_handle& c);

public:
  constrained_bckwd_dfs_evplus(constrained_opname* code,
    expert_forest* cons, expert_forest* arg, expert_forest* trans, expert_forest* res);

  virtual void compute(const dd_edge& a, const dd_edge& b, const dd_edge& r, dd_edge& res);
  void _compute(int aev, node_handle a, int bev, node_handle b, node_handle r, long& cev, node_handle& c);

  void saturateHelper(long aev, node_handle a, unpacked_node& nb);
};

class MEDDLY::constrained_saturation_evplus: public specialized_operation
{
protected:
  constrained_bckwd_dfs_evplus* parent;

  expert_forest* consF;
  expert_forest* argF;
  expert_forest* resF;

  virtual ~constrained_saturation_evplus();

  // Check if the variables orders of relevant forests are compatible
  virtual bool checkForestCompatibility() const;

  bool checkTerminals(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  ct_entry_key* findResult(long aev, node_handle a,
    long bev, node_handle b, int level, long& cev, node_handle &c);
  void saveResult(ct_entry_key* Key,
    long aev, node_handle a, long bev, node_handle b, int level, long cev, node_handle c);

public:
  constrained_saturation_evplus(constrained_bckwd_dfs_evplus* p,
    expert_forest* cons, expert_forest* arg, expert_forest* res);

  bool matches(const expert_forest* arg1, const expert_forest* arg2,
    const expert_forest* res) const;

  // high-level front-end
  void saturate(int aev, node_handle a, int bev, node_handle b, long& cev, node_handle& c);

  void saturate(int aev, node_handle a, int bev, node_handle b, int level, long& cev, node_handle& c);
};

#endif
