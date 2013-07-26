
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


#include "mtmddint.h"


MEDDLY::mt_mdd_int::mt_mdd_int(int dsl, domain *d, const policies &p)
: MEDDLY::mtmdd_forest(dsl, d, false, INTEGER, MULTI_TERMINAL, p)
{ 
  initializeForest();
}


MEDDLY::mt_mdd_int::~mt_mdd_int()
{ }

void MEDDLY::mt_mdd_int::createEdge(int term, dd_edge& e)
{
  if (e.getForest() != this) throw error(error::INVALID_OPERATION);
  createEdgeHelper(getTerminalNode(term), e);
}

void MEDDLY::mt_mdd_int::createEdge(const int* const* vlist,
    const int* terms, int N, dd_edge& e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || terms == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);

  createEdgeInternal(vlist, terms, N, e);
}

void MEDDLY::mt_mdd_int::evaluate(const dd_edge &f,
    const int* vlist, int &term) const
{
  // assumption: vlist does not contain any special values (-1, -2, etc).
  // vlist contains a single element.
  term = getInteger(getTerminalNodeForEdge(f.getNode(), vlist));
}

void MEDDLY::mt_mdd_int::showTerminal(FILE* s, int tnode) const
{
  fprintf(s, "t%d", getInteger(tnode)); 
}

void MEDDLY::mt_mdd_int::writeTerminal(FILE* s, int tnode) const
{
  th_fprintf(s, "t%d", getInteger(tnode)); 
}

MEDDLY::node_handle MEDDLY::mt_mdd_int::readTerminal(FILE* s)
{
  stripWS(s);
  char c = fgetc(s);
  if ('t' == c) {
    int N;
    if (1==fscanf(s, "%d", &N)) {
      if (isValidTerminalValue(N)) {
        return getTerminalNode(N);
      }
    }
  }
  throw error(error::INVALID_FILE);
}

const char* MEDDLY::mt_mdd_int::codeChars() const
{
  return "dd_tvi";
}

