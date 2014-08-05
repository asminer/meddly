
// $Id$

/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "meddly.h"
#include "simple_model.h"

#include "kan_rs1.h"
#include "kan_rs2.h"
#include "kan_rs3.h"

using namespace std;

const char* kanban[] = {
  "X-+..............",  // Tin1
  "X.-+.............",  // Tr1
  "X.+-.............",  // Tb1
  "X.-.+............",  // Tg1
  "X.....-+.........",  // Tr2
  "X.....+-.........",  // Tb2
  "X.....-.+........",  // Tg2
  "X+..--+..-+......",  // Ts1_23
  "X.........-+.....",  // Tr3
  "X.........+-.....",  // Tb3
  "X.........-.+....",  // Tg3
  "X....+..-+..--+..",  // Ts23_4
  "X.............-+.",  // Tr4
  "X.............+-.",  // Tb4
  "X............+..-",  // Tout4
  "X.............-.+"   // Tg4
};

// 160 states for N=1
long expected[] = { 
  1, 160, 4600, 58400, 454475, 2546432, 11261376, 
  41644800, 133865325, 384392800, 1005927208 
};

using namespace MEDDLY;

dd_edge buildReachset(domain* d, int N)
{
  // Build initial state
  int* initial = new int[17];
  for (int i=16; i; i--) initial[i] = 0;
  initial[1] = initial[5] = initial[9] = initial[13] = N;
  forest* mdd = d->createForest(0, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge init_state(mdd);
  mdd->createEdge(&initial, 1, init_state);
  delete[] initial;

  // Build next-state function
  forest* mxd = d->createForest(1, forest::BOOLEAN, forest::MULTI_TERMINAL);
  dd_edge nsf(mxd);
  buildNextStateFunction(kanban, 16, mxd, nsf); 

  dd_edge reachable(mdd);
  apply(REACHABLE_STATES_DFS, init_state, nsf, reachable);
  
  return reachable;
}

bool matches(const char* mark, const int* minterm, int np)
{
  for (int i=0; i<np; i++) 
    if (mark[i]-48 != minterm[i]) return false;
  return true;
}

typedef struct CompLt {
	bool operator() (const int* x, const int* y){
		for(int i=15; i>=0 ;i--) {
			if(x[i]<y[i]) {
				return true;
			}
			if(x[i]>y[i]) {
				return false;
			}
		}
		return false;
	}
}CompLt;

bool checkRS(int N, const char* rs[])
{
  int sizes[16];

  for (int i=15; i>=0; i--) sizes[i] = N+1;
  domain* d = createDomainBottomUp(sizes, 16);

  dd_edge reachable = buildReachset(d, N);

  // enumerate states
  long c = 0;
  vector<int*> assignments;
  for (enumerator i(reachable); i; ++i, c++) {
    if (c>=expected[N]) {
      c++;
      break;
    }

    int* assignment = static_cast<int*>(malloc(16*sizeof(int)));
    if(assignment==0) {
    	throw error(error::INSUFFICIENT_MEMORY);
    }
    memcpy(assignment, i.getAssignments()+1, 16*sizeof(int));
    assignments.push_back(assignment);
  }

  if(c==expected[N]) {
	  sort(assignments.begin(), assignments.end(), CompLt());
	  for(int i=0; i<assignments.size(); i++) {
		  if (!matches(rs[i], assignments[i], 16)) {
		       fprintf(stderr, "Marking %ld mismatched\n", i);

		       for(int k=0; k<16; k++) {
		     	  printf("%d ",rs[i][k]-48);
		       }
		       printf("\n");
		       const int* minterm = assignments[i];
		       for(int k=0; k<16; k++) {
		     	  printf("%d ", minterm[k]);
		       }
		       printf("\n");

		       break;
		  }

		  free(assignments[i]);
	  }
  }

  destroyDomain(d);
  return c==expected[N];
}

void showRS(int N)
{
  int sizes[16];

  for (int i=15; i>=0; i--) sizes[i] = N+1;
  domain* d = createDomainBottomUp(sizes, 16);

  dd_edge reachable = buildReachset(d, N);

  // enumerate states
  long c = 0;
  for (enumerator i(reachable); i; ++i) {
    const int* minterm = i.getAssignments();
    printf("%5ld: %d", c, minterm[1]);
    for (int l=2; l<=16; l++) printf(", %d", minterm[l]);
    printf("\n");
    c++;
  }

  destroyDomain(d);
}

int main()
{
  MEDDLY::initialize();

  printf("Checking Kanban reachability set, N=1\n");
  if (!checkRS(1, kanban_rs1)) return 1;
  printf("\t%ld markings checked out\n", expected[1]);

  printf("Checking Kanban reachability set, N=2\n");
  if (!checkRS(2, kanban_rs2)) return 1;
  printf("\t%ld markings checked out\n", expected[2]);

  printf("Checking Kanban reachability set, N=3\n");
  if (!checkRS(3, kanban_rs3)) return 1;
  printf("\t%ld markings checked out\n", expected[3]);

  MEDDLY::cleanup();
  printf("Done\n");
  return 0;
}

