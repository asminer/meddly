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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "../forests/mt.h"
#include "convert_bdd.h"
#include<set>

namespace MEDDLY {
  class convert_mdd;
  class convert_opname;
};


// ******************************************************************
// *                                                                *
// *                     convert_mdd  class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::convert_mdd : public unary_operation {
public:
  convert_mdd(const unary_opname* cnvt, expert_forest* arg, expert_forest* res);
  
  virtual void computeDDEdge(const dd_edge& a, dd_edge& b, bool userFlag);
  
protected:
  node_handle compute_r(node_handle a);
  
  inline compute_table::entry_key* 
  findResult(node_handle a, node_handle &b) 
  {
  compute_table::entry_key* CTsrch = CT0->useEntryKey(etype[0], 0);
  MEDDLY_DCASSERT(CTsrch);
  CTsrch->writeN(a);
  CT0->find(CTsrch, CTresult[0]);
  if (!CTresult[0]) return CTsrch;
  b = resF->linkNode(CTresult[0].readN());
  CT0->recycle(CTsrch);
  return 0;
  }
  inline node_handle saveResult(compute_table::entry_key* Key, 
                                node_handle a, node_handle b) 
  {
  CTresult[0].reset();
  CTresult[0].writeN(b);
  CT0->addEntry(Key, CTresult[0]);
  return b;
  }
};

MEDDLY::convert_mdd
::convert_mdd(const unary_opname* cnvt, expert_forest* arg, expert_forest* res)
: unary_operation(cnvt, 1, arg, res)
{
  compute_table::entry_type* et = new compute_table::entry_type(cnvt->getName(), "N:N");
  et->setForestForSlot(0, arg);
  et->setForestForSlot(2, res);
  registerEntryType(0, et);
  buildCTs();
}

void MEDDLY::convert_mdd::computeDDEdge(const dd_edge& a, dd_edge& b, bool userFlag) 
{
  node_handle result = compute_r(a.getNode());
  const int num_levels = resF->getDomain()->getNumVariables();
  if (userFlag && result != resF->getTransparentNode() && resF->isQuasiReduced() &&
      resF->getNodeLevel(result) < num_levels) {
    node_handle temp = ((mt_forest*)resF)->makeNodeAtLevel(num_levels, result);
    resF->unlinkNode(result);
    result = temp;
  }
  b.set(result);
}

MEDDLY::node_handle MEDDLY::convert_mdd::compute_r(node_handle m_nh)
{
  
  printf("\n Here to find the thing for %d", m_nh);
  
  // Check compute table
  node_handle b_nh;
  compute_table::entry_key* Key = findResult(m_nh, b_nh);
  if (0==Key)
    {
       printf("\n Found bnh = %d for mnh %d",b_nh, m_nh);
      return b_nh;
    }
  
  const int node_level = argF->getNodeLevel(m_nh);
  unsigned rSize = resF->getLevelSize(node_level);
  unpacked_node *B = unpacked_node::newSparse(resF, node_level, node_level==1?1:2); // just a boolean node;
 
  // base case
  printf("\n Lvl is %d",node_level);
  if(node_level == 1)
    {
      unpacked_node *M = unpacked_node::newFromNode(argF, m_nh, true);
      int nnz = 0;
      for(int i = 0; i < M->getSize(); i++)
        {
        if(M->d(i)!=expert_forest::bool_Tencoder::handle2value(false)) nnz++;
        }
      
      printf("\n The node %d has %d non-zero children",m_nh,nnz);
      if(nnz==1)
        {
        if(M->d(0)!=expert_forest::bool_Tencoder::handle2value(false)) // single child with index 0 is in the set
          {
          printf("\n zero child is true");
          B->i_ref(0) = 0;
          B->d_ref(0) = -1;
          } else // single child with index that is NOT 0 is in the set
            {
            printf("\n non-zero child is true");
            B->i_ref(0) = 1; // if index!=0, then it is 1
            B->d_ref(0) = -1;
            }
        } else if (nnz>1){
          printf("\n non-zero child is true");
          B->i_ref(0) = 1; // if index>=0 is in the set, then it is 1
          B->d_ref(0) = -1;
        } 
      printf("\nAbout to createNode at level %d\n", node_level);
      b_nh = resF->createReducedNode(-1, B);
      printf("\nMDD %d = BDD %d", m_nh, b_nh);
    
      // Add to compute table
      return  saveResult(Key, m_nh, b_nh);
    }
  
  //recursive case
  unpacked_node *M = unpacked_node::newFromNode(argF, m_nh, true);
  int process_zero;
  std::set<node_handle> nz_children_M;
  printf("\n The node size of %d is %d", m_nh, M->getSize());
  printf("\n Converting 0th child of %d handle [%d] to bdd type", m_nh, M->d(0));
  
  if(M->d(0)!=expert_forest::bool_Tencoder::handle2value(false)) 
    process_zero  = 1;
  else 
    process_zero = 0;
  
  for(int i = 1; i < M->getSize(); i++)
    if(M->d(i)!=expert_forest::bool_Tencoder::handle2value(false)) 
          nz_children_M.insert(M->d(i));
  
  printf("\n Zero child exist? %d. Non-zero child exist %d", process_zero,nz_children_M.size() );
  
  dd_edge second(argF), all(argF);
  node_handle union_result_bnh = 0;
  
  if(nz_children_M.size()>0)
    {
      auto it = nz_children_M.begin();
      all.set(*it);
      it++;
      for(;it!=nz_children_M.end();it++)
       { 
         second.set(*it);
         apply(MEDDLY::UNION,all,second,all);
       }
      printf("\n At Level Unioned for mnh %d",all.getNode());
      union_result_bnh = compute_r(all.getNode());
    }
  
  printf("\n Level %d Union bdd handle = %d for mnh %d",node_level, union_result_bnh,all.getNode());

  if(process_zero == 1)
    {
      node_handle zero_nh =  compute_r(M->d(0));
      if((union_result_bnh != 0) && (zero_nh == union_result_bnh))
        {
          B->i_ref(0) = 1; 
          B->d_ref(0) = union_result_bnh;
        }
      else {
        printf("\n Adding zero child %d", M->d(0)); 
        B->i_ref(0) = 0; 
        B->d_ref(0) = compute_r(M->d(0));
    
        if(nz_children_M.size() > 0)
        {
          printf("\n Adding one child ---> ");
          B->i_ref(1) = 1; 
          B->d_ref(1) = union_result_bnh;
          printf("--> %d",B->d(1));
        }
      }
    
    printf("\n Shrink to do");
    if(nz_children_M.size() == 0) 
      B->shrinkSparse(1);
    else  
      {
        if(zero_nh == union_result_bnh)  B->shrinkSparse(1);
        else B->shrinkSparse(2);
      }
    } else {
      printf("\n Adding zero one index child of %d--->", m_nh);
      B->i_ref(0) = 1; 
      //B->d_ref(0) = 0;
      B->d_ref(0) = union_result_bnh;
      printf("--> %d",B->d(0));
      B->shrinkSparse(1);
    }
  
  // cleanup and Reduce
  unpacked_node::recycle(M);
  b_nh = resF->createReducedNode(-1, B);
  
  // Add to compute table
  return saveResult(Key, m_nh, b_nh);
}

// ******************************************************************
// *                                                                *
// *                       convert_opname class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::convert_opname : public unary_opname {
public:
  convert_opname();
  virtual unary_operation*
  buildOperation(expert_forest* ar, expert_forest* res) const;
};

MEDDLY::convert_opname::convert_opname()
: unary_opname("Convert_to_bdd")
{
}

MEDDLY::unary_operation*
MEDDLY::convert_opname
::buildOperation(expert_forest* arg, expert_forest* res) const
{
 /* if (0==arg) return 0;
  
  if (arg->getDomain() != res->getDomain())
    throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
  
  if (arg->getRangeType() != forest::BOOLEAN ||
      arg->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      res->getRangeType() != forest::BOOLEAN ||
      res->getEdgeLabeling() != forest::MULTI_TERMINAL ||
      arg->isForRelations() || 
      res->isForRelations() 
      ) throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
  */
 
 return new convert_mdd(this,  arg,  res);
}
// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname* MEDDLY::initializeConvertBdd()
{
  return new convert_opname;
}
