
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


#include "../defines.h"
#include "mtmxdbool.h"

MEDDLY::mt_mxd_bool::mt_mxd_bool(unsigned dsl, domain *d, const policies &p, int* level_reduction_rule, bool tv)
: mtmxd_forest(dsl, d, BOOLEAN, p, level_reduction_rule)
{ 
  initializeForest();

  transparent=bool_Tencoder::value2handle(tv);
}

MEDDLY::mt_mxd_bool::~mt_mxd_bool()
{ }

void MEDDLY::mt_mxd_bool::createEdge(bool term, dd_edge& e)
{
  createEdgeTempl<bool_Tencoder, bool>(term, e);
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

  mtmxd_edgemaker<bool_Tencoder, bool>
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
  createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mxd_bool::evaluate(const dd_edge &f, const int* vlist, 
  const int* vplist, bool &term) const
{
  term = bool_Tencoder::handle2value(evaluateRaw(f, vlist, vplist));
}

void MEDDLY::mt_mxd_bool::showTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::show(s, tnode);
}

void MEDDLY::mt_mxd_bool::writeTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::write(s, tnode);
}

MEDDLY::node_handle MEDDLY::mt_mxd_bool::readTerminal(input &s)
{
  return bool_Tencoder::read(s);
}

const char* MEDDLY::mt_mxd_bool::codeChars() const
{
  return "dd_txb";
}

bool
MEDDLY::mt_mxd_bool::checkTerminalMinterm(node_handle a, std::vector<std::vector<int>> terms, node_handle &result)
{
  
   //printf("\n I am checking TERMINAL\n");
  if (a <= 0 && (terms[0].size() == 0 && terms[1].size() == 0) )
    {
    result = handleForValue(true);
    //printf("\n I am TERMINAL\n");
    return true;
    }
  
  //initial forest is empty and I am adding a single term
  if( getNodeLevel(a) == 0 && terms[0].size() > 0 && terms[1].size() > 0 ) 
    {
        //printf("\n I am SEMITERMINAL\n");
        int num_levels = terms[1].size();
        node_handle below = handleForValue(true);
        for(int i = 0; i < num_levels; i++)
          {
              int level_from = i + 1;
              int level_to = -level_from;
              unsigned level_from_size = getLevelSize(level_from);
              unsigned level_to_size =  getLevelSize(level_to);
              int b_from = terms[0].at(i);
              int b_to = terms[1].at(i);
          
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
                
              
              for(int j = 0; j < level_to_size; j++ )
                {
                  if(j!=b_to) P_node->d_ref(j) = 0;
                  else
                    {
                      //terms[1].erase(terms[1].begin()+level_from - 1);
                      P_node->d_ref(j) = linkNode(below);
                    }
                }
                
                for(int j = 0; j < level_from_size; j++ )
                {
                if(j!=b_from) UP_node->d_ref(j) = 0;
                else
                  {
                    //terms[0].erase(terms[0].begin()+level_from - 1);
                    UP_node->d_ref(j) = createReducedNode(j, P_node);
                  }
                }
                below = createReducedNode(-1, UP_node);
                
          }
                
                result = below;
                return true;
    
    }
                
    
    return false;
}
                
                
   
MEDDLY::node_handle 
MEDDLY::mt_mxd_bool::unionOneMinterm(node_handle a, std::vector<std::vector<int>> terms) 
{
  node_handle result = 0;
  if(checkTerminalMinterm(a, terms, result))
    return result;
  
  const int aLevel = getNodeLevel(a);
  const int bLevel = terms[0].size();
  int resultLevel = ABS(topLevel(aLevel,bLevel));
  int dwnLevel = downLevel(resultLevel);
  
  MEDDLY_DCASSERT(bLevel == resultLevel);
  
  // Existing mdd has FF kind of nodes at resultLevel
  if(aLevel < resultLevel)
   { 
    if( ( terms[0].at(bLevel-1) == DONT_CARE && ( terms[1].at(bLevel-1) == DONT_CHANGE || terms[1].at(bLevel-1) == DONT_CARE )) )
    { 
      // Minterm to be added has FF kind of nodes at resultLevel
      terms[0].erase(terms[0].begin() + bLevel - 1);
      terms[1].erase(terms[1].begin() + bLevel - 1);
      return unionOneMinterm(a, terms); 
    } else { 
      // Minterm to be added has non-skippable nodes at resultLevel
      unsigned b_from = terms[0].at(bLevel-1);
      unpacked_node *C = unpacked_node::newSparse(this, resultLevel, getLevelSize(resultLevel));
      
      terms[0].erase(terms[0].begin() + bLevel - 1);
      C->d_ref(b_from) = unionOneMinterm_r(b_from, dwnLevel, a, terms); 
      result = createReducedNode(-1, C);
     
      return result;
    }
  }
  
  if(aLevel == resultLevel)
    {
    
    if( ( terms[0].at(bLevel-1) == DONT_CARE && ( terms[1].at(bLevel-1) == DONT_CHANGE || terms[1].at(bLevel-1) == DONT_CARE )) )
      { 
        terms[0].erase(terms[0].begin() + bLevel - 1);
        unpacked_node *A = unpacked_node::newFromNode(this, a, true);
        unpacked_node *C = unpacked_node::newFull(this, resultLevel, getLevelSize(resultLevel));
        for(int i = 0; i < getLevelSize(resultLevel); i++ )
          {
            C->d_ref(i) = unionOneMinterm_r(i, dwnLevel, A->d(i), terms);
          }
         result = createReducedNode(-1, C);
        } else {
        unsigned b_from = terms[0].at(bLevel-1);
        unpacked_node *A = unpacked_node::newFromNode(this, a, true);
        unpacked_node *C = unpacked_node::newFull(this, resultLevel, getLevelSize(resultLevel));
    
          for(int i = 0; i < getLevelSize(resultLevel); i++ )
            {
              if(i==b_from) terms[0].erase(terms[0].begin() + bLevel - 1);
               C->d_ref(i) = i == b_from ? unionOneMinterm_r(i, dwnLevel, A->d(i), terms) : linkNode(A->d(i)); 
             }
          
            result = createReducedNode(-1, C);
        }
    }
  
  return result;
}
                

MEDDLY::node_handle 
MEDDLY::mt_mxd_bool::unionOneMinterm_r(int in, int k, node_handle a, std::vector<std::vector<int>> terms) 
{
  node_handle result = 0;
  const int aLevel = getNodeLevel(a);
  const int bLevel = terms[1].size();
  unsigned resultSize = unsigned(getLevelSize(k));

  unpacked_node* C = unpacked_node::newFull(this, k, resultSize);
 
  unpacked_node *A =
  (aLevel == k) 
  ? unpacked_node::newFromNode(this, a, true)
  : this->isFullyReduced()       
  ? unpacked_node::newRedundant(this, k, a, true)  // has to be a FR mxd.... make sure it is
  : unpacked_node::newIdentity(this, k, in, a, true)
  ;
  
  // Do computation
  int b_to  = terms[1].at(bLevel-1);
  
  if(b_to == DONT_CARE || b_to == DONT_CHANGE)
    {
      terms[1].erase(terms[1].begin()+bLevel-1);
      for (unsigned j=0; j<resultSize; j++) {
           C->d_ref(j) = unionOneMinterm(A->d(j), terms);
        }
    }
  else {
      for (unsigned j=0; j<resultSize; j++) {
        if(j==b_to)
          {
            terms[1].erase(terms[1].begin()+bLevel-1);
            C->d_ref(j) = unionOneMinterm(A->d(j), terms);
          } else
            C->d_ref(j) = linkNode(A->d(j));
        
      }
  }
  
  // cleanup
  unpacked_node::recycle(A);
  
  // reduce and return result
  result = createReducedNode(-1, C);
  return result;
}
