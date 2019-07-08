
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

MEDDLY::mt_mdd_bool::mt_mdd_bool(int dsl, domain *d, const policies &p, int* level_reduction_rule, bool tv)
: mtmdd_forest(dsl, d, BOOLEAN, p, level_reduction_rule)
{ 
  initializeForest();

  transparent=bool_Tencoder::value2handle(tv);
}

MEDDLY::mt_mdd_bool::~mt_mdd_bool()
{ }

void MEDDLY::mt_mdd_bool::createEdge(bool term, dd_edge& e)
{
  createEdgeTempl<bool_Tencoder, bool>(term, e);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

int MEDDLY::mt_mdd_bool::storedValue(int value) const{
	if (value==Omega){
		return 0;
	}else{
		return value+1;
	}
}

int MEDDLY::mt_mdd_bool::returnStored(int value) const{
	if (value==Omega){
		return Omega;
	}else{
		return value-1;
	}
}
void MEDDLY::mt_mdd_bool::createEdge(const int* const* vlist, int N, dd_edge &e)
{
  binary_operation* unionOp = getOperation(UNION, this, this, this);
  enlargeStatics(N);
  enlargeVariables(vlist, N, false);
  
  int num_vars=getNumVariables();

  // Create vlist following the mapping between variable and level
  int** ordered_vlist=static_cast<int**>(malloc(N*sizeof(int*)+(num_vars+1)*N*sizeof(int)));
  if(ordered_vlist==0){
	  throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
  }

  ordered_vlist[0]=reinterpret_cast<int*>(&ordered_vlist[N]);
  for(int i=1; i<N; i++) {
	  ordered_vlist[i]=(ordered_vlist[i-1]+num_vars+1);
  }
  for(int i=0; i<=num_vars; i++) {
	  int level=getLevelByVar(i);
	  for(int j=0; j<N; j++) {
		  //considering +infty if PInfty is true.
		  if (e.getHasPInfty()){
		  ordered_vlist[j][level]=storedValue(vlist[j][i]);
		  }else
		  {
			  ordered_vlist[j][level]=vlist[j][i];
		  }
	  }
  }

  mtmdd_edgemaker<bool_Tencoder, bool>
  EM(this, ordered_vlist, 0, order, N, getDomain()->getNumVariables(), unionOp);

  e.set(EM.createEdge());

  free(ordered_vlist);

#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool::
createEdgeForVar(int vh, bool vp, const bool* terms, dd_edge& a)
{
  createEdgeForVarTempl<bool_Tencoder, bool>(vh, vp, terms, a);
#ifdef DEVELOPMENT_CODE
  validateIncounts(true);
#endif
}

void MEDDLY::mt_mdd_bool
::evaluate(const dd_edge &f, const int* vlist, bool &term) const
{
	if (f.getHasPInfty()) {
		const int NumberOfLevel = f.getLevel();
		//int arraySize = sizeof(vlist) / sizeof(int);
		int* updatedVlist = new int[NumberOfLevel + 1];
		for (int ind = 1; ind <= NumberOfLevel; ind++) {
//			const int value = vlist[ind];
			updatedVlist[ind] = storedValue(vlist[ind]);
		}
		term = bool_Tencoder::handle2value(evaluateRaw(f, updatedVlist));
		free(updatedVlist);
	} else {
		term = bool_Tencoder::handle2value(evaluateRaw(f, vlist));
	}
}

void MEDDLY::mt_mdd_bool::isMarkingCovered(const dd_edge &f, const int* vlist,
		bool &term) const {
	node_handle p = f.getNode();
//	printf("XXXXF%d\n",f.getHasPInfty());
	if (!f.getHasPInfty()) {
		term = evaluateRawIsMarkingCovered(f, p, vlist);
	} else {
		const int NumberOfLevel = f.getLevel();
		//int arraySize = sizeof(vlist) / sizeof(int);
		int* updatedVlist = new int[NumberOfLevel + 1];
		for (int ind = 1; ind <= NumberOfLevel; ind++) {
			const int value = vlist[ind];
			updatedVlist[ind] = storedValue(value);
		}
//		printf("CALL WITH INFTY\n");
		term = evaluateRawIsMarkingCoveredWithInfty(f, p, updatedVlist);
		free(updatedVlist);
	}
}

void MEDDLY::mt_mdd_bool::isMarkingCovered(const dd_edge &f, const int* vlist,
		bool &term, std::map<node_handle, bool>& map) const {
	node_handle p = f.getNode();
//	printf("XXXXF%d\n",f.getHasPInfty());
	if (!f.getHasPInfty()) {
		term = evaluateRawIsMarkingCovered(f, p, vlist);
	} else {
		const int NumberOfLevel = f.getLevel();
		//int arraySize = sizeof(vlist) / sizeof(int);
		int* updatedVlist = new int[NumberOfLevel + 1];
		for (int ind = 1; ind <= NumberOfLevel; ind++) {
			const int value = vlist[ind];
			updatedVlist[ind] = storedValue(value);
		}
//		printf("CALL WITH INFTY\n");
		term = evaluateRawIsMarkingCoveredWithInfty(f, p, updatedVlist, map);
		free(updatedVlist);
	}
}
void MEDDLY::mt_mdd_bool::firstMarkingCovers(const dd_edge &f, const int* vlist,
		bool &term, int* rlist) const {
	node_handle p = f.getNode();
//	printf("XXXXF%d\n",f.getHasPInfty());
//	if (!f.getHasPInfty()){
//	term = evaluateRawIsMarkingCovered(f, p, vlist);
//	}
//	else{
	int rlistsize = getNodeLevel(p);
	int* resultlist = new int[rlistsize];
//		printf("evaluateRawFirstMarkingCoveredWithInfty WITH INFTY\n");
//	const int NumberOfLevel=f.getLevel();
	//int arraySize = sizeof(vlist) / sizeof(int);
	int* updatedVlist = new int[rlistsize + 1];
	for (int ind = 1; ind <= rlistsize; ind++) {
		const int value = vlist[ind];
		updatedVlist[ind] = storedValue(value);
	}
	term = evaluateRawFirstMarkingCoversWithInfty(f, p, updatedVlist,
			resultlist);
	free(updatedVlist);
//		printf("RESULT is %d\n",term);
	if (term) {
		for (int i = rlistsize; i >= 0; i--)
			if (i != 0)
				rlist[i] = returnStored(resultlist[i]);
			else
				rlist[i] = 0;
//			printf("%d ",rlist[i]);
	}
	free(resultlist);
//	}
}

void MEDDLY::mt_mdd_bool::firstMarkingCovers(const dd_edge &f, const int* vlist,
		bool &term, int* rlist, std::map<node_handle, int>& map) const {
	node_handle p = f.getNode();
//	printf("XXXXF%d\n",f.getHasPInfty());
//	if (!f.getHasPInfty()){
//	term = evaluateRawIsMarkingCovered(f, p, vlist);
//	}
//	else{
	int rlistsize = getNodeLevel(p);
	int* resultlist = new int[rlistsize];
//		printf("evaluateRawFirstMarkingCoveredWithInfty WITH INFTY\n");
//	const int NumberOfLevel=f.getLevel();
	//int arraySize = sizeof(vlist) / sizeof(int);
	int* updatedVlist = new int[rlistsize + 1];
	for (int ind = 1; ind <= rlistsize; ind++) {
		const int value = vlist[ind];
		updatedVlist[ind] = storedValue(value);
	}
	term = evaluateRawFirstMarkingCoversWithInfty(f, p, updatedVlist,
			resultlist, map);
	free(updatedVlist);
//		printf("RESULT is %d\n",term);
	if (term) {
		for (int i = rlistsize; i >= 0; i--)
			if (i != 0)
				rlist[i] = returnStored(resultlist[i]);
			else
				rlist[i] = 0;
//			printf("%d ",rlist[i]);
	}
	free(resultlist);
//	}
}

void MEDDLY::mt_mdd_bool::allMarkingsCover(const dd_edge &f, int* vlist,
		int* r,std::list<int*> R,bool &term, const dd_edge &w,
		std::map<node_handle, int>&map)const{

}
void MEDDLY::mt_mdd_bool::updateMarkedNode(const dd_edge &f,const int level,
		std::map<node_handle, int>& map) const {
	node_handle* p = new node_handle(f.getNode());
	node_handle* list = markNodesInSubgraph(p, 1, true);
	if (0 == list)
		return;
//	for (int k = 1; k <= getNumVariables(); ){
	for (long i = 0; list[i]; i++) {
		if (getNodeLevel(list[i]) >= level) {
			std::map<node_handle, int>::iterator it = map.find(list[i]);
			if (it != map.end()&& it->second==unpacked_node::markedWith::NC) {
				map[list[i]] = unpacked_node::markedWith::NV;
			}
		}
	}
//	k++;
//	}
	free(list);
	free(p);

}

void MEDDLY::mt_mdd_bool::showTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::show(s, tnode);
}

void MEDDLY::mt_mdd_bool::writeTerminal(output &s, node_handle tnode) const
{
  bool_Tencoder::write(s, tnode);
}

MEDDLY::node_handle MEDDLY::mt_mdd_bool::readTerminal(input &s)
{
  return bool_Tencoder::read(s);
}

const char* MEDDLY::mt_mdd_bool::codeChars() const
{
  return "dd_tvb";
}
