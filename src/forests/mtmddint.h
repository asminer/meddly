
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

#ifndef MTMDDINT_H
#define MTMDDINT_H

#include "mtmdd.h"

namespace MEDDLY {
  class mt_mdd_int; 
};

// ******************************************************************

#ifdef NEW_MT

/** 
    Forest for multi-terminal, mdd, integer range.
*/
class MEDDLY::mt_mdd_int : public mtmdd_forest<int_terminal> {
  public:

    mt_mdd_int(int dsl, domain *d, const policies &p);
    ~mt_mdd_int();

    virtual void createEdge(int val, dd_edge &e);
    virtual void createEdge(int** vlist, int* terms, int N, dd_edge &e);
    virtual void createEdgeForVar(int vh, bool vp, const int* terms, dd_edge& a);
    virtual void evaluate(const dd_edge &f, const int* vlist, int &term)
      const;

  protected:
    virtual const char* codeChars() const;
};


#else

/** 
    Forest for multi-terminal, mdd, integer range.
*/
class MEDDLY::mt_mdd_int : public mtmdd_forest {
  public:

    mt_mdd_int(int dsl, domain *d, const policies &p);
    ~mt_mdd_int();

    virtual void createEdge(int val, dd_edge &e);
    virtual void createEdge(int** vlist, int* terms, int N, dd_edge &e);
    virtual void evaluate(const dd_edge &f, const int* vlist, int &term)
      const;

    virtual void showTerminal(FILE* s, node_handle tnode) const;
    virtual void writeTerminal(FILE* s, node_handle tnode) const;
    virtual node_handle readTerminal(FILE* s);
  protected:
    virtual const char* codeChars() const;
};

#endif

#endif

