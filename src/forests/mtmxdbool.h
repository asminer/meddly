
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

#ifndef MTMXDBOOL_H
#define MTMXDBOOL_H

#include "mtmxd.h"

namespace MEDDLY {
  class mt_mxd_bool;
};

// ******************************************************************

/** 
    Forest for multi-terminal, mxd, boolean range.
*/
class MEDDLY::mt_mxd_bool : public mtmxd_forest {
  // TODO: mxds can only be forest::IDENTITY_REDUCED
  public:

    mt_mxd_bool(int dsl, domain *d, const policies &p);
    ~mt_mxd_bool();

    using mtmxd_forest::createEdge;
    using mtmxd_forest::evaluate;

    virtual bool accumulate(int& tempNode, int* element, int* pelement);
    virtual void accumulate(int& a, int b);

    virtual int accumulate(int tempNode, bool cBM,
        int* element, int* pelement, int level);
    virtual int accumulateSkippedLevel(int tempNode,
        int* element, int* pelement, int level);

    virtual int accumulateMxd(int a, int b, bool cBM);
    virtual int accumulateMxdPrime(int a, int b, bool cBM);
    virtual int addPrimeReducedNodes(int a, int b);
    virtual int accumulateExpandA(int a, int b, bool cBM);
    virtual int buildQRIdentityNode(int node, int level);

    virtual void accumulateMxdHelper(int& a, int b, bool cBM,
        bool needsToMakeACopy,
        int (mt_mxd_bool::*function)(int, int, bool));

    // Refer to meddly.h
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        int N, dd_edge& e);
    virtual void createEdge(bool val, dd_edge &e);
    virtual void evaluate(const dd_edge& f, const int* vlist,
        const int* vplist, bool &term) const;

    // The following will either abort or return an error since they are not
    // applicable to this forest.
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const int* terms, int N, dd_edge &e);
    virtual void createEdge(const int* const* vlist, const int* const* vplist,
        const float* terms, int N, dd_edge &e);
    virtual void createEdge(int val, dd_edge &e);
    virtual void createEdge(float val, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist,
        const int* vplist, int &term) const;
    virtual void evaluate(const dd_edge &f, const int* vlist,
        const int* vplist, float &term) const;
};

#endif

