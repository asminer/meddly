
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


#include "mtmxdbool.h"


MEDDLY::mt_mxd_bool::mt_mxd_bool(int dsl, domain *d, const policies &p)
: MEDDLY::mtmxd_forest(dsl, d, true, BOOLEAN, MULTI_TERMINAL, p)
{ 
  initializeForest();
}


MEDDLY::mt_mxd_bool::~mt_mxd_bool()
{ }


void MEDDLY::mt_mxd_bool::createEdge(bool val, dd_edge &e)
{
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);

  int node = createEdgeTo(getTerminalNode(val));
  e.set(node, 0);
}


void MEDDLY::mt_mxd_bool::createEdge(const int* const* vlist,
    const int* const* vplist, int N, dd_edge &e)
{
  if (e.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0 || N <= 0) 
    throw error(error::INVALID_VARIABLE);
  createEdgeInternal(vlist, vplist, (bool*)0, N, e);
}


void MEDDLY::mt_mxd_bool::evaluate(const dd_edge& f, const int* vlist,
    const int* vplist, bool &term) const
{
  MEDDLY_DCASSERT(getRangeType() == forest::BOOLEAN);
  if (f.getForest() != this) 
    throw error(error::INVALID_OPERATION);
  if (vlist == 0 || vplist == 0) 
    throw error(error::INVALID_VARIABLE);

  // assumption: vlist and vplist do not contain any special values
  // (-1, -2, etc). vlist and vplist contains a single element.
  term = getBoolean(getTerminalNodeForEdge(f.getNode(), vlist, vplist));
}

void MEDDLY::mt_mxd_bool::showTerminal(FILE* s, int tnode) const
{
  fprintf(s, "%c", tnode ? 'T' : 'F'); 
}

void MEDDLY::mt_mxd_bool::writeTerminal(FILE* s, int tnode) const
{
  th_fprintf(s, "%c", tnode ? 'T' : 'F'); 
}

MEDDLY::node_handle MEDDLY::mt_mxd_bool::readTerminal(FILE* s)
{
  stripWS(s);
  char c = fgetc(s);
  if ('T' == c) return -1;
  if ('F' == c) return  0;
  throw error(error::INVALID_FILE);
}

const char* MEDDLY::mt_mxd_bool::codeChars() const
{
  return "dd_txb";
}

