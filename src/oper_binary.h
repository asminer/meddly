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

#ifndef MEDDLY_OPER_BINARY_H
#define MEDDLY_OPER_BINARY_H

#include "oper.h"

namespace MEDDLY {
    class binary_operation;

    // ******************************************************************
    // *                      Operation management                      *
    // ******************************************************************

    /** Find, or build if necessary, a binary operation.
            @param  code    Operation we want
            @param  arg1    Argument 1 forest
            @param  arg2    Argument 2 forest
            @param  res     Result forest
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    binary_operation* getOperation(const binary_opname* code,
        expert_forest* arg1, expert_forest* arg2, expert_forest* res);

    /** Find, or build if necessary, a binary operation.
            @param  code    Operation we want
            @param  arg1    Argument 1 forest taken from this dd_edge
            @param  arg2    Argument 2 forest taken from this dd_edge
            @param  res     Result forest taken from this dd_edge
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    binary_operation* getOperation(const binary_opname* code,
        const dd_edge& arg1, const dd_edge& arg2, const dd_edge& res);

    /** Safely destroy the given binary operation.
        It should be unnecessary to call this directly.
    */
    void destroyOperation(binary_operation* &op);

};

// ******************************************************************
// *                                                                *
// *                     binary_operation class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a binary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::binary_operation : public operation {
  protected:
    bool can_commute;
    expert_forest* arg1F;
    expert_forest* arg2F;
    expert_forest* resF;
    opnd_type resultType;

    virtual ~binary_operation();
    void operationCommutes();

    // Check if the variables orders of relevant forests are compatible
    virtual bool checkForestCompatibility() const;

  public:
    binary_operation(const binary_opname* code, unsigned et_slots,
      expert_forest* arg1, expert_forest* arg2, expert_forest* res);

    bool matches(const expert_forest* arg1, const expert_forest* arg2,
      const expert_forest* res) const;

    // high-level front-end

    /**
      Checks forest comatability and then calls computeDDEdge().
    */
    void compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);
    void computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res);

    virtual void computeDDEdge(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res, bool userFlag)
      = 0;

    // low-level front ends.  TBD - REMOVE THESE BECAUSE THEY BREAK MARK & SWEEP

#ifdef KEEP_LL_COMPUTES

    /// Low-level compute on nodes a and b, return result.
    virtual node_handle compute(node_handle a, node_handle b);
    /// Low-level compute at level k on nodes a and b, return result.
    virtual node_handle compute(int k, node_handle a, node_handle b);

    /// Low-level compute on EV edges (av, ap) and (bv, bp), return result.
    virtual void compute(int av, node_handle ap, int bv, node_handle bp,
      int &cv, node_handle &cp);
    virtual void compute(long av, node_handle ap, long bv, node_handle bp,
      long &cv, node_handle &cp);
    virtual void compute(long av, node_handle ap, node_handle bp,
      long &cv, node_handle &cp);

    /// Low-level compute on EV edges (av, ap) and (bv, bp), return result.
    virtual void compute(float av, node_handle ap, float bv, node_handle bp,
      float &cv, node_handle &cp);

#endif

};

// ******************************************************************
// *                                                                *
// *                inlined binary_operation methods                *
// *                                                                *
// ******************************************************************


inline bool
MEDDLY::binary_operation::matches(const MEDDLY::expert_forest* arg1,
    const MEDDLY::expert_forest* arg2, const MEDDLY::expert_forest* res) const
{
  return (arg1 == arg1F && arg2 == arg2F && res == resF);
}

inline void
MEDDLY::binary_operation::operationCommutes()
{
  can_commute = (arg1F == arg2F);
}

inline bool
MEDDLY::binary_operation::checkForestCompatibility() const
{
  auto o1 = arg1F->variableOrder();
  auto o2 = arg2F->variableOrder();
  auto o3 = resF->variableOrder();
  return o1->is_compatible_with(*o2) && o1->is_compatible_with(*o3);
}

inline void
MEDDLY::binary_operation::compute(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(ar1, ar2, res, true);
}

inline void
MEDDLY::binary_operation::computeTemp(const dd_edge &ar1, const dd_edge &ar2, dd_edge &res)
{
  if (!checkForestCompatibility()) {
    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
  }
  computeDDEdge(ar1, ar2, res, false);
}

// ******************************************************************
// *                                                                *
// *                       inlined  functions                       *
// *                                                                *
// ******************************************************************


inline MEDDLY::binary_operation*
MEDDLY::getOperation(const binary_opname* code, const dd_edge& arg1,
    const dd_edge& arg2, const dd_edge& res)
{
  return getOperation(code, (MEDDLY::expert_forest*) arg1.getForest(),
      (MEDDLY::expert_forest*) arg2.getForest(), (MEDDLY::expert_forest*) res.getForest());
}


#endif // #include guard
