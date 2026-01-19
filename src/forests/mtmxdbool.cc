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
#include "forest_levels.h"

MEDDLY::mt_mxd_bool::mt_mxd_bool(domain *d, const policies &p)
: mtmxd_forest(d, range_type::BOOLEAN, p)
{
}

MEDDLY::mt_mxd_bool::~mt_mxd_bool()
{ }

#ifdef ALLOW_DEPRECATED_0_18_0

void MEDDLY::mt_mxd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  // createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
  createEdgeForVarTempl<bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true, __FILE__, __LINE__);
#endif
}

#endif


bool
MEDDLY::mt_mxd_bool::checkTerminalMinterm(node_handle a,  int* from,  int* to, int levelRead, node_handle &result)
{
    terminal one_terminal(true);

   if (a <= 0 && levelRead == 0 )
    {
        result = one_terminal.getHandle();
        //printf("\n I am TERMINAL\n");
        return true;
    }

  //initial forest is empty and I am adding a single term
  #if 1
  if( getNodeLevel(a) == 0 && levelRead > 0 )
    {

        node_handle below = one_terminal.getHandle();
        for(int i = 1; i <= levelRead; i++)
          {
              int level_from = i;
              int level_to = -i;
              unsigned level_from_size = getLevelSize(level_from);
              unsigned level_to_size =  getLevelSize(level_to);
              int b_from = from[i];
              int b_to = to[i];

              //MEDDLY_DCASSERT(!isExtensibleLevel(level));

              // the level not in this subevent
              // so skip
              if((b_from == DONT_CARE) && (b_to == DONT_CARE))
                 continue;


             // the level not in this event and subevent
             // so skip
             if((b_from == DONT_CARE) && (b_to == DONT_CHANGE))
                continue;

             // else build the prime & unprime nodes
             unpacked_node *UP_node = unpacked_node::newWritable(this,
                     level_from, level_from_size, FULL_ONLY);
             unpacked_node *P_node = unpacked_node::newWritable(this,
                     level_to, level_to_size, FULL_ONLY);

              if(b_to >= 0){
              for(int j = 0; j < level_to_size; j++ )
                {
                  P_node->setFull(j, (j!=b_to) ? 0 : linkNode(below));
                }

                for(int j = 0; j < level_from_size; j++ )
                {
                  UP_node->setFull(j,
                          (j!=b_from) ? 0 : createReducedNode(j, P_node)
                  );
                }
                below = createReducedNode(-1, UP_node);
               } else {

                for(int j = 0; j < level_from_size; j++ )
                {
                  UP_node->setFull(j, (j!=b_from) ? 0 : below);
                }
                below = createReducedNode(-1, UP_node);
               }
          }

                result = below;
                return true;

    }
  #endif

    return false;
}


MEDDLY::node_handle
MEDDLY::mt_mxd_bool::unionOneMinterm(node_handle a,  int* from,  int* to, int level)
{
  node_handle result = 0;
  if(checkTerminalMinterm(a, from, to, level, result))
    return result;

  const int aLevel = getNodeLevel(a);
  const int bLevel = level;
  int resultLevel = ABS(MXD_levels::topLevel(aLevel,bLevel));
  int dwnLevel = MXD_levels::downLevel(resultLevel);

  if(aLevel < bLevel) {
    return unionOneMinterm(a, from, to, bLevel-1);
  }
  else {
    unpacked_node *A = unpacked_node::newFromNode(this, a, FULL_ONLY);
    unpacked_node *C = unpacked_node::newWritable(this, resultLevel, FULL_ONLY);
    for( int i = 0; i < getLevelSize(resultLevel); i++) {

      if(i == from[bLevel])
      {
        C->setFull(i, unionOneMinterm_r(i, dwnLevel, A->down(i), from, to));
      }
      else
      {
        C->setFull(i, linkNode(A->down(i)));
      }
    }

    result = createReducedNode(-1, C);
  }

  return result;
}


MEDDLY::node_handle
MEDDLY::mt_mxd_bool::unionOneMinterm_r(int in, int k, node_handle a,  int* from,  int* to)
{
  node_handle result = 0;
  const int aLevel = getNodeLevel(a);
  const int bLevel = -k;
  unsigned resultSize = unsigned(getLevelSize(k));

  unpacked_node* C = unpacked_node::newWritable(this, k, resultSize, FULL_ONLY);

  terminal one_terminal(true);

  //if a is terminal
  if(aLevel == 0)
  {
    if(a == one_terminal.getHandle())
    {
      return a;
    }
    else {
      if(to[bLevel] == DONT_CARE)
      {
        return one_terminal.getHandle();
      }
      else {
        for (unsigned j=0; j<resultSize; j++) {
            if(j == to[bLevel]) {
                C->setFull(j, unionOneMinterm(a, from, to, MXD_levels::downLevel(k)));
            } else {
                C->setFull(j, 0);
            }
        }
        result = createReducedNode(-1,C);
        return result;
      }
    }
  }


  int b_to  = to[bLevel];
  //if a in not at this level
  if( this->isFullyReduced() && (aLevel != k) )
  {
    if(b_to == DONT_CARE) return unionOneMinterm(a, from, to, MXD_levels::downLevel(k));
    else if(b_to == DONT_CHANGE)
    {
      C->setFull(in, unionOneMinterm(a, from, to, MXD_levels::downLevel(k)));
      return createReducedNode(-1,C);
    }
  } else if( this->isIdentityReduced() && (b_to == DONT_CHANGE) && (aLevel != k) )
  {
    return unionOneMinterm(a, from, to, MXD_levels::downLevel(k));
  }


  //unpacked_node* C = unpacked_node::newFull(this, k, resultSize);
  // if a is at this level
  unpacked_node *A =
  (aLevel == k)
  ? unpacked_node::newFromNode(this, a, FULL_ONLY)
  : this->isFullyReduced()
  ? unpacked_node::newRedundant(this, k, a, FULL_ONLY)
  : unpacked_node::newIdentity(this, k, in, a, FULL_ONLY)
  ;

  // Do computation


  for (unsigned j=0; j<resultSize; j++) {
    if(j==b_to) {
        C->setFull(j, unionOneMinterm(A->down(j), from, to, MXD_levels::downLevel(k)));
    } else {
        C->setFull(j, linkNode(A->down(j)));
    }
  }

  // cleanup
  unpacked_node::Recycle(A);

  // reduce and return result

  result = createReducedNode(-1, C);
  return result;
}
