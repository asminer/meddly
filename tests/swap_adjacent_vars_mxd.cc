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

// $Id: test_union_mdd.cc 495 2014-03-07 01:13:30Z asminer $
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

#include <iostream>
#include <vector>
#include <algorithm>

#include "meddly.h"
#include "timer.h"

using namespace MEDDLY;

//#define VERBOSE

// Use MTMXD with real values for terminals.
#define USE_REALS 0

// Define element_type based on USE_REALS.
#if USE_REALS
  typedef float element_type;
#else
  typedef int element_type;
#endif

void printUsage(FILE *outputStream) {
	fprintf(outputStream,
			"Usage: test_union_mdd <#Variables> <VariableBound> <#Elements>\n");
}

void printState(const dd_edge& state) {
	printf("State:\n");
//	state.show(stdout, 2);
}

void printStat(const dd_edge& state) {
	printf("Peak Nodes in MDD: %ld\n", state.getForest()->getPeakNumNodes());
	printf("Current Nodes in MDD: %ld\n",
			state.getForest()->getCurrentNumNodes());
}

typedef struct CompLt {
	bool operator()(const std::vector<int>& x, const std::vector<int>& y) {
		MEDDLY_DCASSERT(x.size()==y.size());

		for (int i = x.size() - 1; i >= 0; i--) {
			if (x[i] < y[i]) {
				return true;
			}
			if (x[i] > y[i]) {
				return false;
			}
		}
		return false;
	}
} CompLt;

void printAssignments(const std::vector<std::vector<int> >& assign) {
	for (std::vector<std::vector<int> >::const_iterator iter = assign.begin();
			iter != assign.end(); iter++) {
		for (std::vector<int>::const_iterator siter = iter->begin();
				siter < iter->end(); siter++) {
			std::cout << *siter << " ";
		}
		std::cout << std::endl;
	}
}

void saveAssignments(dd_edge& state, std::vector<std::vector<int> >& assign) {
	assign.clear();

	int nVariables = state.getForest()->getDomain()->getNumVariables();
	for (enumerator iter(state); iter; ++iter) {
		assign.push_back(
				std::vector<int>(iter.getAssignments() + 1,
						iter.getAssignments() + 1 + nVariables));
	}
	std::sort(assign.begin(), assign.end(), CompLt());
}

void verifyAssignments(const std::vector<std::vector<int> >& assign1,
		const std::vector<std::vector<int> >& assign2) {
	if (assign1.size() != assign2.size()) {
		std::cerr << "Error in Size" << std::endl;
		exit(1);
	}

	for (int i = 0; i < assign1.size(); i++) {
		for (int j = 0; j < assign1[i].size(); j++) {
			if (assign1[i][j] != assign2[i][j]) {
				std::cerr << "Error in Value" << std::endl;
				exit(1);
			}
		}
	}
}

void printInCount(forest* f) {
	expert_forest* ef = static_cast<expert_forest*>(f);
	for (int i = 1; i < ef->getPeakNumNodes(); i++) {
		if (ef->isValidNodeIndex(i)
				&& (ef->isTerminalNode(i) || ef->getNode(i).offset > 0)) {
			std::cout << "Node " << i << ": " << ef->getInCount(i) << std::endl;
		}
	}
}

int main(int argc, char *argv[]) {
	unsigned int seed = time(NULL);
	srandom(seed);

	std::cout<<"Seed: "<<seed<<std::endl;

	// Initialize number of variables, their bounds and the number of elements
	// to create
	int nVariables = 50;
	int variableBound = 4;
	int nElements = 1000;

	std::cout << "#variables: " << nVariables << ", variable bound: "
			<< variableBound << ", #elements: " << nElements << std::endl;

	// create the elements randomly
	int** from = (int **) malloc(nElements * sizeof(int *));
	int** to = (int **) malloc(nElements * sizeof(int *));
	element_type* terms =
			(element_type *) malloc(nElements * sizeof(element_type));

	for (int i = 0; i < nElements; ++i)
	{
		from[i] = (int *) malloc((nVariables + 1) * sizeof(int));
		to[i] = (int *) malloc((nVariables + 1) * sizeof(int));
		from[i][0] = 0;
		to[i][0] = 0;
		for (int j = nVariables; j > 0; --j)
		{
			from[i][j] = int(float(variableBound) * random() / (RAND_MAX + 1.0));
			assert(from[i][j] >= 0 && from[i][j] < variableBound);
			to[i][j] = int(float(variableBound) * random() / (RAND_MAX + 1.0));
			assert(to[i][j] >= 0 && to[i][j] < variableBound);
		}
		terms[i] =
				element_type(int(float(variableBound) * random() / (RAND_MAX + 1.0)));
	}

	// Initialize the variable bounds array to provide to the domain
	int* bounds = (int *) malloc(nVariables * sizeof(int));
	assert(bounds != 0);
	for (int i = 0; i < nVariables; ++i) {
		bounds[i] = variableBound;
	}

	settings s;
	initialize(s);

	// Create a domain
	domain *d = createDomainBottomUp(bounds, nVariables);
	assert(d != 0);

	// Create a MTMXD forest in this domain
	forest::policies p(true);
#if USE_REALS
	forest* mtmxd =
			d->createForest(true, forest::REAL, forest::MULTI_TERMINAL, p);
#else
	forest* mtmxd =
			d->createForest(true, forest::INTEGER, forest::MULTI_TERMINAL, p);
#endif
	assert(mtmxd != 0);

	dd_edge state(mtmxd);

	timer start;

	printf("Started... ");

	// Use Meddly's batch addition to combine all elements in one step.
	mtmxd->createEdge(from, to, nElements, state);

//  printState(state);

	std::vector<std::vector<int> > assign1;
	saveAssignments(state, assign1);

	std::vector<std::vector<int> > assign2;

	printf("\n==================== Bubble Reordering ====================\n");
	int count = 1;
	start.note_time();
	for (int var = 1; var <= nVariables; var++) {
		for (int level = 1; level < nVariables; level++) {
			printf("%d. Swap variables at level %d and %d\n", count, level,
					level + 1);

			static_cast<expert_domain*>(d)->swapAdjacentVariables(level);

//			printState(state);
			saveAssignments(state, assign2);
			verifyAssignments(assign1, assign2);

			count++;
		}
	}
	start.note_time();
	printf("Time interval: %.4e seconds\n", start.get_last_interval()/1000000.0);

	printf("\n==================== Reverse Reordering ====================\n");
	count = 1;
	start.note_time();
	for (int level = 1; level < nVariables; level++) {
		printf("%d. Swap variables at level %d and %d\n", count, level,
				level + 1);
		static_cast<expert_domain*>(d)->swapAdjacentVariables(level);

//		printState(state);
		saveAssignments(state, assign2);
		verifyAssignments(assign1, assign2);

		count++;
	}
	for (int level = nVariables - 1; level >= 1; level--) {
		printf("%d. Swap variables at level %d and %d\n", count, level,
				level + 1);
		static_cast<expert_domain*>(d)->swapAdjacentVariables(level);

//		printState(state);
		saveAssignments(state, assign2);
		verifyAssignments(assign1, assign2);

		count++;
	}
	start.note_time();
	printf("Time interval: %.4e seconds\n", start.get_last_interval()/1000000.0);

	// Cleanup; in this case simply delete the domain
	destroyDomain(d);
	MEDDLY::cleanup();

	free(bounds);
	for (int i = 0; i < nElements; ++i) {
		free(elements[i]);
	}
	free(elements);

	return 0;
}

