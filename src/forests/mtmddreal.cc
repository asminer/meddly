
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


#include "mtmddreal.h"

#ifdef NEW_MT

MEDDLY::mt_mdd_real::mt_mdd_real(int dsl, domain *d, const policies &p)
: mtmdd_forest<real_terminal>(dsl, d, REAL, p)
{ 
  initializeForest();
}

MEDDLY::mt_mdd_real::~mt_mdd_real()
{ }

void MEDDLY::mt_mdd_real::createEdge(float term, dd_edge& e)
{
  createEdgeTempl(term, e);
}

void MEDDLY::mt_mdd_real::createEdge(int** vlist, float* terms, int N, dd_edge &e)
{
  unionOp = getOperation(PLUS, this, this, this);
  enlargeVariables(vlist, N, false);
  e.set(createEdgeRT(getDomain()->getNumVariables(), vlist, terms, N), 0);
}

void MEDDLY::mt_mdd_real::
createEdgeForVar(int vh, bool vp, const float* terms, dd_edge& a)
{
  createEdgeForVarTempl(vh, vp, terms, a);
}


void MEDDLY::mt_mdd_real
::evaluate(const dd_edge &f, const int* vlist, float &term) const
{
  evaluateTempl(f, vlist, term);
}



#else

MEDDLY::mt_mdd_real::mt_mdd_real(int dsl, domain *d, const policies &p)
: MEDDLY::mtmdd_forest(dsl, d, false, REAL, MULTI_TERMINAL, p)
{ 
  initializeForest();
}


MEDDLY::mt_mdd_real::~mt_mdd_real()
{ }

void MEDDLY::mt_mdd_real::createEdge(float term, dd_edge& e)
{
  if (e.getForest() != this) throw error(error::INVALID_OPERATION);
  createEdgeHelper(getTerminalNode(term), e);
}

void MEDDLY::mt_mdd_real::createEdge(int** vlist, float* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || terms == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);

  createEdgeInternal(vlist, terms, N, e);
}

void MEDDLY::mt_mdd_real::evaluate(const dd_edge &f,
    const int* vlist, float &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  term = getReal(getTerminalNodeForEdge(f.getNode(), vlist));
}

void MEDDLY::mt_mdd_real::showTerminal(FILE* s, int tnode) const
{
  fprintf(s, "t%f", getReal(tnode)); 
}

void MEDDLY::mt_mdd_real::writeTerminal(FILE* s, int tnode) const
{
  th_fprintf(s, "t%8e", getReal(tnode)); 
}

MEDDLY::node_handle MEDDLY::mt_mdd_real::readTerminal(FILE* s)
{
  stripWS(s);
  char c = fgetc(s);
  if ('t' == c) {
    float T;
    if (1==fscanf(s, "%8e", &T)) {
      return getTerminalNode(T); 
    }
  }
  throw error(error::INVALID_FILE);
}

#endif

const char* MEDDLY::mt_mdd_real::codeChars() const
{
  return "dd_tvr";
}

