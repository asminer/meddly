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

#include "relation_node.h"
#include "forest.h"

// TBD: move this!
#define NOT_KNOWN -2

// ******************************************************************
// *                                                                *
// *                      relation_node methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::relation_node::relation_node(unsigned long sign, forest* f, int lvl, node_handle d, long enable_val, long fire_val, long inh_val) :
f(f)
{
  signature  = sign;
  level = lvl;
  down = d;
  enable = enable_val;
  fire = fire_val;
  inhibit = inh_val;
  piece_size = 0;
  token_update = NULL;
}

MEDDLY::relation_node::~relation_node()
{
}


long MEDDLY::relation_node::nextOf(long i)
{
  //to be defined for the example you use & comment this definition
  throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

bool
MEDDLY::relation_node::equals(const relation_node* n) const
{
  if((signature == n->getSignature()) && (level == n->getLevel()) && (down == n->getDown()) && (fire == n->getFire()) && (enable == n->getEnable()))
    return true;
  else
    return false;
}

void
MEDDLY::relation_node::expandTokenUpdate(long i)
{
  if(getPieceSize()==0)
  {
    token_update = (long*)malloc(1*sizeof(long));
    piece_size = 1;
    token_update[0]=NOT_KNOWN;
  }
  if(i>0)
  {
    token_update = (long*)realloc(token_update,size_t(i+1)*sizeof(long));
    for(int j = piece_size;j<=i;j++)
      token_update[j]=NOT_KNOWN;
    piece_size = i+1;
  }
}

void
MEDDLY::relation_node::setTokenUpdateAtIndex(long i,long val)
{
  MEDDLY_DCASSERT(i<getPieceSize());
  token_update[i] = val;
}
