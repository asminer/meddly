
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

#include <iostream>
#include <fstream>

#include "meddly.h"
#include "timer.h"
#include "reorder.h"

using namespace MEDDLY;
using namespace std;

void printUsage()
{
	cout<<"Usage: test_reorder <ElementFile> <OrderFile> <ReorderingHeuristic>"<<endl;
}

void printStat(const dd_edge& state)
{
	printf("Current Nodes in MDD: %ld\n", state.getForest()->getCurrentNumNodes());
}

void readElementFile(char* elementFileName, int**& elements, int& numVariable, int& numElement)
{
	assert(elements==NULL);

	ifstream in(elementFileName);
	const int MAX_LINE_LENGTH=1000;
	char lineBuffer[MAX_LINE_LENGTH];
	char literalBuffer[10];
	int x=0;
	int y=1;
	while(in.getline(lineBuffer, MAX_LINE_LENGTH)){
		if(lineBuffer[0]=='#'){
			if(lineBuffer[1]=='v'){
				// Number of variables
				sscanf(&lineBuffer[3], "%d", &numVariable);
			}
			else if(lineBuffer[1]=='e'){
				// Number of elements
				sscanf(&lineBuffer[3], "%d", &numElement);
				elements=new int*[numElement];
			}
			else if(lineBuffer[1]=='s'){
				// Seed
			}
			else {
				exit(0);
			}
		}
		else{
			elements[x]=new int[numVariable+1];
			y=1;

			char* linePtr=lineBuffer;
			char* literalPtr=literalBuffer;
			while(*linePtr!='\0'){
				literalPtr = literalBuffer;
				while(*linePtr && isspace(*linePtr)){
					linePtr++;
				}
				while(*linePtr && !isspace(*linePtr)){
					*(literalPtr++) = *(linePtr++); // copy Literal
				}
				*literalPtr = '\0'; // terminate Literal

				if(strlen(literalBuffer)>0){
					elements[x][y++] = atoi(literalBuffer);
				}
			}

			assert(y==numVariable+1);
			x++;
		}
	}
	assert(x==numElement);
}

int main(int argc, char *argv[])
{
	if(argc!=4){
		printUsage();
		exit(0);
	}

	char elementFileName[200];
	sscanf(argv[1], "%s", &elementFileName);
	char orderFileName[200];
	sscanf(argv[2], "%s", &orderFileName);
	char heuristic[30];
	sscanf(argv[3], "%s", &heuristic);

	int numVariable=0;
	int numElement=0;
	int** elements=NULL;
	int* order=NULL;

	readElementFile(elementFileName, elements, numVariable, numElement);
	readOrderFile(orderFileName, order, numVariable);

	int* bounds = new int[numVariable];
	for (int i = 0; i < numVariable; ++i) {
		bounds[i] = 2;
	}

	settings s;
	initialize(s);

	// Create a domain
	domain *d = createDomainBottomUp(bounds, numVariable);

	// Create an MDD forest in this domain (to store states)
	forest::policies p(false);
	if(strcmp(heuristic, "LI")==0) {
		p.setLowestInversion();
	}
	else if(strcmp(heuristic, "HI")==0) {
		p.setHighestInversion();
	}
	else if(strcmp(heuristic, "BD")==0) {
		p.setBubbleDown();
	}
	else if(strcmp(heuristic, "BU")==0) {
		p.setBubbleUp();
	}
	else if(strcmp(heuristic, "LC")==0) {
		p.setLowestCost();
	}

	forest* f = d->createForest(false, forest::BOOLEAN,
			forest::MULTI_TERMINAL, p);

	dd_edge state(f);

	cout<<"Started..."<<endl;

	f->createEdge(elements, numElement, state);

	static_cast<expert_forest*>(f)->removeAllComputeTableEntries();
	printStat(state);

	timer runTimer;
	runTimer.note_time();

	static_cast<expert_domain*>(d)->reorderVariables(order);

	runTimer.note_time();
	cout<<"Time: "<<runTimer.get_last_interval()/1000000.0<<" s"<<endl;
	printStat(state);

	destroyDomain(d);
	MEDDLY::cleanup();

	delete[] bounds;
	delete[] order;
	for (int i = 0; i < numElement; i++) {
		delete[] elements[i];
	}
	delete[] elements;

	return 0;
}
