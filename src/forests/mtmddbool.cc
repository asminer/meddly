
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


#include "mtmddbool.h"

#ifdef NEW_MT

MEDDLY::mt_mdd_bool::mt_mdd_bool(int dsl, domain *d, const policies &p)
: mtmdd_forest<int_terminal>(dsl, d, BOOLEAN, p)
{ 
  initializeForest();
}

MEDDLY::mt_mdd_bool::~mt_mdd_bool()
{ }

void MEDDLY::mt_mdd_bool::createEdge(bool term, dd_edge& e)
{
  createEdgeTempl(term, e);
}

void MEDDLY::mt_mdd_bool::createEdge(int** vlist, int N, dd_edge &e)
{
  unionOp = getOperation(UNION, this, this, this);
  enlargeVariables(vlist, N, false);
  e.set(createEdgeRT(getDomain()->getNumVariables(), vlist, (bool*) 0, N), 0);
}

void MEDDLY::mt_mdd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  createEdgeForVarTempl(vh, vp, terms, a);
}

void MEDDLY::mt_mdd_bool
::evaluate(const dd_edge &f, const int* vlist, bool &term) const
{
  evaluateTempl(f, vlist, term);
}



#else


MEDDLY::mt_mdd_bool::mt_mdd_bool(int dsl, domain *d, const policies &p)
: MEDDLY::mtmdd_forest(dsl, d, false, BOOLEAN, MULTI_TERMINAL, p)
{ 
  initializeForest();
}


MEDDLY::mt_mdd_bool::~mt_mdd_bool()
{ }


void MEDDLY::mt_mdd_bool::createEdge(bool term, dd_edge& e)
{
  if (e.getForest() != this) throw error(error::INVALID_OPERATION);
  createEdgeHelper(getTerminalNode(term), e);
}


void MEDDLY::mt_mdd_bool::createEdge(int** vlist, int N, dd_edge &e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);
  createEdgeInternal(vlist, (bool*)0, N, e);
}


void MEDDLY::mt_mdd_bool::evaluate(const dd_edge &f, const int* vlist,
    bool &term) const
{
  term = getBoolean(getTerminalNodeForEdge(f.getNode(), vlist));
}

void MEDDLY::mt_mdd_bool::showTerminal(FILE* s, int tnode) const
{
  fprintf(s, "%c", tnode ? 'T' : 'F'); 
}

void MEDDLY::mt_mdd_bool::writeTerminal(FILE* s, int tnode) const
{
  th_fprintf(s, "%c", tnode ? 'T' : 'F'); 
}

MEDDLY::node_handle MEDDLY::mt_mdd_bool::readTerminal(FILE* s)
{
  stripWS(s);
  char c = fgetc(s);
  if ('T' == c) return -1;
  if ('F' == c) return  0;
  throw error(error::INVALID_FILE);
}

#endif

const char* MEDDLY::mt_mdd_bool::codeChars() const
{
  return "dd_tvb";
}
