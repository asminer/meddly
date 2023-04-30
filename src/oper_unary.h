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

#ifndef MEDDLY_OPER_UNARY_H
#define MEDDLY_OPER_UNARY_H

#include "oper.h"

namespace MEDDLY {
    class dd_edge;
    class unary_operation;

    // ******************************************************************
    // *                      Operation management                      *
    // ******************************************************************

    /** Find, or build if necessary, a unary operation.
            @param  code    Operation we want
            @param  arg     Argument forest
            @param  res     Result forest
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    unary_operation* getOperation(const unary_opname* code,
        expert_forest* arg, expert_forest* res);

    /** Find, or build if necessary, a unary operation.
            @param  code    Operation we want
            @param  arg     Argument forest from this dd_edge
            @param  res     Result forest from this dd_edge
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    unary_operation* getOperation(const unary_opname* code,
        const dd_edge& arg, const dd_edge& res);

    /** Find, or build if necessary, a unary operation.
            @param  code    Operation we want
            @param  arg     Argument forest
            @param  res     Result type
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    unary_operation* getOperation(const unary_opname* code,
        expert_forest* arg, opnd_type result);

    /** Find, or build if necessary, a unary operation.
            @param  code    Operation we want
            @param  arg     Argument forest from this dd_edge
            @param  res     Result type
            @return         The matching operation, if it already exists;
                            a new operation, otherwise.
    */
    unary_operation* getOperation(const unary_opname* code,
        const dd_edge& arg, opnd_type result);


    /** Safely destroy the given unary operation.
        It should be unnecessary to call this directly.
    */
    void destroyOperation(unary_operation* &op);


    //
    //
    extern bool libraryRunning;
};


// ******************************************************************
// *                                                                *
// *                     unary_operation  class                     *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a unary operation in a specific forest.
    Specific operations will be derived from this class.
*/
class MEDDLY::unary_operation : public operation {
  protected:
    expert_forest* argF;
    expert_forest* resF;
    opnd_type resultType;

    virtual ~unary_operation();

    virtual bool checkForestCompatibility() const;

  public:
    unary_operation(const unary_opname* code, unsigned et_slots,
      expert_forest* arg, expert_forest* res);

    unary_operation(const unary_opname* code, unsigned et_slots,
      expert_forest* arg, opnd_type res);

    bool matches(const expert_forest* arg, const expert_forest* res)
      const;

    bool matches(const expert_forest* arg, opnd_type res) const;

    // high-level front-ends

    /**
      Checks forest comatability and then calls computeDDEdge().
    */
    void compute(const dd_edge &arg, dd_edge &res);
    void computeTemp(const dd_edge &arg, dd_edge &res);

    virtual void compute(const dd_edge &arg, long &res);
    virtual void compute(const dd_edge &arg, double &res);
    virtual void compute(const dd_edge &arg, ct_object &c);
    virtual void computeDDEdge(const dd_edge &arg, dd_edge &res, bool userFlag);
};

// ******************************************************************
// *                                                                *
// *                inlined  unary_operation methods                *
// *                                                                *
// ******************************************************************

inline bool
MEDDLY::unary_operation::matches(const MEDDLY::expert_forest* arg,
    const MEDDLY::expert_forest* res) const
{
  return (arg == argF && res == resF);
}

inline bool
MEDDLY::unary_operation::matches(const MEDDLY::expert_forest* arg, opnd_type res) const
{
  return (arg == argF && resultType == res);
}

#endif // #include guard
