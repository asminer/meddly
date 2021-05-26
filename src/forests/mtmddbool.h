
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

#ifndef MTMDDBOOL_H
#define MTMDDBOOL_H
// #define DBG_MTMDD

#include "mtmdd.h"
#include <cfloat>
#include "../operations/mpz_object.h"
#include <gmp.h>
namespace MEDDLY {
  class mt_mdd_bool;
};

// ******************************************************************

/**
    Forest for multi-terminal, mdd, boolean range.
*/
class MEDDLY::mt_mdd_bool : public mtmdd_forest {
  public:

    mt_mdd_bool(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule=NULL, bool tv=false);
    ~mt_mdd_bool();

    virtual void createEdge(bool val, dd_edge &e);
    virtual void createEdge(const int* const* vlist, int N, dd_edge &e);
    virtual void createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a);
    virtual void evaluate(const dd_edge &f, const int* vlist, bool &term)
      const;

    virtual void showTerminal(output &s, node_handle tnode) const;
    virtual void writeTerminal(output &s, node_handle tnode) const;
    virtual node_handle readTerminal(input &s);
    virtual void underApproximate(dd_edge &e,long minThreshold, long maxThreshold,float desiredPercentage, int option);
    virtual void HeuristicUnderApproximate(dd_edge &e, long minThreshold, long maxThreshold,float desiredPercentage, int option,int deletedApproach, float rootStatePercentage);

  protected:
    virtual const char* codeChars() const;
private:
      void RemoveDuplicate(int lvl, std::map<int,int>map);
      int RemoveDuplicate2(int lvl, std::map<int,int>map,dd_edge &e,std::set<int>&RNA,std::set<int>&RNB);
      int maxid=0;
      void uniqueNodesforp(node_handle p,std::set<int> &result );
      void getNC(int lvl, node_handle a,bool*visitedNode,std::set<int> &result);
      void uniqueAboveNodesforp(node_handle p,std::set<int> &result );
      int RemoveDuplicateSet(int lvl, std::set<int>levels,std::map<int,int>map,dd_edge &e,std::set<int>&RNA,std::set<int>&RNB);
};

#endif
