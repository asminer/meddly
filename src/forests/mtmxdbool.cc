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

MEDDLY::mt_mxd_bool::mt_mxd_bool(domain *d, const policies &p, int* level_reduction_rule, bool tv)
: mtmxd_forest(d, range_type::BOOLEAN, p, level_reduction_rule)
{
    terminal t(tv);
    setTransparentEdge(t.getHandle());
    // setTransparentEdge(bool_Tencoder::value2handle(tv));
    initializeStorage();
}

MEDDLY::mt_mxd_bool::~mt_mxd_bool()
{ }

void MEDDLY::mt_mxd_bool::createEdge(bool term, dd_edge& e)
{
  // createEdgeTempl<bool_Tencoder, bool>(term, e);
  createEdgeTempl<bool>(term, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_bool
::createEdge(const int* const* vlist, const int* const* vplist, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(UNION, this, this, this);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);
  enlargeVariables(vplist, N, true);

  int num_vars=getNumVariables();

  // Create vlist and vplist following the mapping between variable and level
  int** ordered_vlist=static_cast<int**>(malloc(N*sizeof(int*)+(num_vars+1)*N*sizeof(int)));
  if(ordered_vlist==0){
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }
  int** ordered_vplist=static_cast<int**>(malloc(N*sizeof(int*)+(num_vars+1)*N*sizeof(int)));
  if(ordered_vplist==0){
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }

  ordered_vlist[0]=reinterpret_cast<int*>(&ordered_vlist[N]);
  ordered_vplist[0]=reinterpret_cast<int*>(&ordered_vplist[N]);
  for(int i=1; i<N; i++) {
	  ordered_vlist[i]=(ordered_vlist[i-1]+num_vars+1);
	  ordered_vplist[i]=(ordered_vplist[i-1]+num_vars+1);
  }
  for(int i=0; i<=num_vars; i++) {
	  int level=getLevelByVar(i);
	  for(int j=0; j<N; j++) {
		  ordered_vlist[j][level]=vlist[j][i];
		  ordered_vplist[j][level]=vplist[j][i];
	  }
  }

  // mtmxd_edgemaker<bool_Tencoder, bool>
  mtmxd_edgemaker<bool>
  EM(this, ordered_vlist, ordered_vplist, 0, order, N, num_vars, unionOp);

  e.set(EM.createEdge());

  free(ordered_vlist);
  free(ordered_vplist);

#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  // createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
  createEdgeForVarTempl<bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_bool::evaluate(const dd_edge &f, const int* vlist,
  const int* vplist, bool &term) const
{
    terminal t(terminal_type::BOOLEAN, evaluateRaw(f, vlist, vplist));
    term = t.getBoolean();
}

void MEDDLY::mt_mxd_bool::showEdge(output &s, const edge_value &ev,
        node_handle d) const
{
    if (d>0) {
        s.put('#');
        s.put(d);
    } else {
        terminal t(terminal_type::BOOLEAN, d);
        s.put( t.getBoolean() ? 'T' : 'F' );
    }
}


#ifdef ALLOW_DEPRECATED_0_17_3

const char* MEDDLY::mt_mxd_bool::codeChars() const
{
  return "dd_txb";
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
             unpacked_node *UP_node = unpacked_node::newFull(this, level_from, level_from_size);
             unpacked_node *P_node = unpacked_node::newFull(this, level_to, level_to_size);

              if(b_to >= 0){
              for(int j = 0; j < level_to_size; j++ )
                {
                  if (j!=b_to) P_node->d_ref(j) = 0;
                  else
                    {
                      P_node->d_ref(j) = linkNode(below);
                    }
                }

                for(int j = 0; j < level_from_size; j++ )
                {
                 if(j!=b_from) UP_node->d_ref(j) = 0;
                 else
                  {
                    UP_node->d_ref(j) = createReducedNode(j, P_node);
                  }
                }
                below = createReducedNode(-1, UP_node);
               } else {

                for(int j = 0; j < level_from_size; j++ )
                {
                 if(j!=b_from) UP_node->d_ref(j) = 0;
                 else
                  {
                    UP_node->d_ref(j) = below;
                  }
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
  int resultLevel = ABS(topLevel(aLevel,bLevel));
  int dwnLevel = downLevel(resultLevel);

  if(aLevel < bLevel) {
    return unionOneMinterm(a, from, to, bLevel-1);
  }
  else {
    unpacked_node *A = newUnpacked(a, FULL_ONLY);
    unpacked_node *C = unpacked_node::newFull(this, resultLevel, getLevelSize(resultLevel));
    for( int i = 0; i < getLevelSize(resultLevel); i++) {

      if(i == from[bLevel])
      {
        C->d_ref(i) = unionOneMinterm_r(i, dwnLevel, A->d(i), from, to);
      }
      else
        C->d_ref(i) = linkNode(A->d(i));
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

  unpacked_node* C = unpacked_node::newFull(this, k, resultSize);

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
        if(j == to[bLevel])
          C->d_ref(to[bLevel]) = unionOneMinterm(a, from, to, downLevel(k));
        else
          C->d_ref(j) = 0;
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
    if(b_to == DONT_CARE) return unionOneMinterm(a, from, to, downLevel(k));
    else if(b_to == DONT_CHANGE)
    {
      C->d_ref(in) = unionOneMinterm(a, from, to, downLevel(k));
      return createReducedNode(-1,C);
    }
  } else if( this->isIdentityReduced() && (b_to == DONT_CHANGE) && (aLevel != k) )
  {
    return unionOneMinterm(a, from, to, downLevel(k));
  }


  //unpacked_node* C = unpacked_node::newFull(this, k, resultSize);
  // if a is at this level
  unpacked_node *A =
  (aLevel == k)
  ? newUnpacked(a, FULL_ONLY)
  : this->isFullyReduced()
  ? unpacked_node::newRedundant(this, k, a, true)  // has to be a FR mxd.... make sure it is
  : unpacked_node::newIdentity(this, k, in, a, true)
  ;

  // Do computation


  for (unsigned j=0; j<resultSize; j++) {
    if(j==b_to)
      {
        C->d_ref(j) = unionOneMinterm(A->d(j), from, to, downLevel(k));
      } else
        C->d_ref(j) = linkNode(A->d(j));
  }

  // cleanup
  unpacked_node::recycle(A);

  // reduce and return result

  result = createReducedNode(-1, C);
  return result;
}
