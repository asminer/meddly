
// $Id: rubiks_cube.cc 653 2016-02-17 00:00:51Z cjiang1209 $

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



/*
  State space generation for a 3x3 Rubik's Cube
  -- Junaid Babar

	The model has 6 square faces with 9 squares (of the same color) in each face
	for a total of 54 squares.
  
	These 54 squares can be grouped together as some of them move together.

	Type 1: Single square (6 such squares)  ,  6 * 1 =  6
	Type 2: Two squares (12 such L-corners) , 12 * 2 = 24
	Type 3: Three squares (8 such corners)  ,  8 * 3 = 24
	
	The model thus has 26 components and 26 component locations.
	Each component location is represented as a MDD level. The size of the MDD
	level depends on the number of components that can fit into that location.
	For example, Type 2 component locations are of size 12.
  
	Level 0        : Terminals
	Level 01 - 12  : Type 2 component locations (size 12)
	Level 13 - 20  : Type 3 component locations (size 8)

  The 6 Type 1 component locations (size 6) are not represented since they
  never move.	Previously at levels 21 - 26.

	Up:     In order (going right starting from left-upper corner)
          of components (type:id:var),
          (3:4:17, 2:5:6, 3:5:18, 2:6:7, 3:1:14, 2:0:1, 3:0:13, 2:4:5)
          Note: (1:0) also belongs to this row but since it stays in the
          same location when this row is moved left or right, it is ignored.
	Down:   (3:3:16, 2:2:3, 3:2:15, 2:8:9, 3:6:19, 2:9:10, 3:7:20, 2:10:11) (1:5 ignored)
	Left:   (3:4:17, 2:4:5, 3:0:13, 2:3:4, 3:3:16, 2:10:11, 3:7:20, 2:11:12) (1:4 ignored)
	Right:  (3:1:14, 2:6:7, 3:5:18, 2:7:8, 3:6:19, 2:8:9, 3:2:15, 2:1:2) (1:2 ignored)
	Front:  (3:0:13, 2:0:1, 3:1:14, 2:1:2, 3:2:15, 2:2:3, 3:3:16, 2:3:4) (1:1 ignored)
	Back:   (3:5:18, 2:5:6, 3:4:17, 2:11:12, 3:7:20, 2:9:10, 3:6:19, 2:7:8) (1:3 ignored)

	Initially components are placed in components locations that match their
	Ids.
*/

#define SATURATE 1
#define NSF_CARDINALITY 1

#define REORDER_CUBE 0
#define FACED_ORDERING 0

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cassert>

#include "meddly.h"
#include "meddly_expert.h"
#include "timer.h"

using namespace std;
using namespace MEDDLY;

enum class Face {F, B, L, R, U, D};
enum class Direction {CW, CCW, FLIP};

// level *variables = NULL;
int* sizes = nullptr;
int** initst = NULL;

static int num_components = 20;

// Number of variables of each type
const int type1 = 6;
const int type2 = 12;
const int type3 = 8;

// Number of levels
const int num_levels = type2 + type3;

// Domain handle
domain *d;

// Forest storing the next state function
// forest_hndl relation;
//forest* relations[3];
vector<forest*> relations;

// Forest storing the set of states
// forest_hndl states;
forest* states;

static int comp_type[] =
  { 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 };
static int comp3_map[] = { 13, 14, 15, 16, 17, 18, 19, 20 };
static int comp2_map[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

// Specify index of each variable
// Variable 0 is not used
static int level2vars[3][21] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

static int num_phases = 3;

struct Action
{
  Face face;
  Direction direction;
};

struct Component
{
  int type;
  int id;
};

static Component level2components[3][21] = {
  {
    // Not used
    { 0, 0 },
    // Front
    { 3, 0 }, { 3, 1 }, { 3, 2 }, { 3, 3 },
    { 2, 0 }, { 2, 1 }, { 2, 2 }, { 2, 3 },
    // Back
    { 3, 5 }, { 3, 4 }, { 3, 7 }, { 3, 6 },
    { 2, 5 }, { 2, 11 }, { 2, 9 }, { 2, 7 },
    // Not affected
    { 2, 4 }, { 2, 10 }, { 2, 8 }, { 2, 6 }
  },
  {
    // Not used
    { 0, 0 },
    // Left
    { 3, 4 }, { 3, 0 }, { 3, 3 }, { 3, 7 },
    { 2, 4 }, { 2, 3 }, { 2, 10 }, { 2, 11 },
    // Right
    { 3, 1 }, { 3, 5 }, { 3, 6 }, { 3, 2 },
    { 2, 6 }, { 2, 7 }, { 2, 8 }, { 2, 1 },
    // Not affected
    { 2, 5 }, { 2, 9 }, { 2, 2 }, { 2, 0 }
  },
  {
    // Not used
    { 0, 0 },
    // Up
    { 3, 4 }, { 3, 5 }, { 3, 1 }, { 3, 0 },
    { 2, 5 }, { 2, 6 }, { 2, 0 }, { 2, 4 },
    // Down
    { 3, 6 }, { 3, 7 }, { 3, 3 }, { 3, 2 },
    { 2, 9 }, { 2, 10 }, { 2, 2 }, { 2, 8 },
    // Not affected
    { 2, 11 }, { 2, 3 }, { 2, 1 }, { 2, 7 }
  }
};

//static vector<vector<Action>> actions = {
//  { { Face::F, Direction::CW }, { Face::B, Direction::CW } },
//  { { Face::L, Direction::CW }, { Face::R, Direction::CW } },
//  { { Face::U, Direction::CW }, { Face::D, Direction::CW } }
//};

static vector<vector<Action>> actions = {
  { { Face::F, Direction::CW }, { Face::B, Direction::CW } },
  { { Face::L, Direction::CW }, { Face::R, Direction::CW } },
  { { Face::U, Direction::CW }, { Face::D, Direction::CW } }
};

//static vector<vector<Action>> actions ={
//  { { Face::F, Direction::CW }, { Face::F, Direction::CCW }, { Face::F, Direction::FLIP },
//    { Face::B, Direction::CW }, { Face::B, Direction::CCW }, { Face::B, Direction::FLIP } },
//  { { Face::L, Direction::CW }, { Face::L, Direction::CCW }, { Face::L, Direction::FLIP },
//    { Face::R, Direction::CW }, { Face::R, Direction::CCW }, { Face::R, Direction::FLIP } },
//  { { Face::U, Direction::CW }, { Face::U, Direction::CCW }, { Face::U, Direction::FLIP },
//    { Face::D, Direction::CW }, { Face::D, Direction::CCW }, { Face::D, Direction::FLIP } }
//};

//void SetUpArrays()
//{
//  // Set up comp_type arrays based on order[]
//  for (int i = 1; i < num_levels + 1; i++) {
//    comp_type[i] = order[i] > 11? 3: 2;
//    if (comp_type[i] == 3) {
//      comp3_map[order[i] - 12] = i;
//    } else {
//      comp2_map[order[i]] = i;
//    }
//  }
//
//  for (int i = 0; i < 12; i++) {
//    fprintf(stderr, "%d, ", comp2_map[i]);
//  }
//  fprintf(stderr, "\n");
//  for (int i = 0; i < 8; i++) {
//    fprintf(stderr, "%d, ", comp3_map[i]);
//  }
//  fprintf(stderr, "\n");
//}

int get_component_type (int comp_level)
{
  assert(comp_level > 0 && comp_level <= num_levels);
  return comp_type[comp_level];
}

int get_component_size (int c_type)
{
  assert(c_type == 2 || c_type == 3);
  return (c_type == 3)? type3: type2;
}

int getLevelSize (int comp_level)
{
  return get_component_size(get_component_type(comp_level));
}

int get_component_level (int c_type, int index)
{
  assert(c_type == 2 || c_type == 3);
  assert((c_type == 3 && index >= 0 && index < 8) ||
      (c_type == 2 && index >= 0 && index < 12));

  return (c_type == 3)? comp3_map[index]: comp2_map[index];
}


void Init()
{
  assert(num_levels == (type2 + type3));

  // node size for each level
  sizes = (int *) malloc(num_levels * sizeof(int));
  assert(sizes != NULL);
  
  assert(num_levels == 20);
  for (int i = 0; i < num_levels; i++) {
    sizes[i] = getLevelSize(i+1);
//    fprintf(stderr, "sizes[%d] = %d\n", i, sizes[i]);
  }
//  fflush(stderr);

  // sets of states
  initst = (int**) malloc(1 * sizeof(int*));
  assert(initst != NULL);
  initst[0] = (int *) malloc((num_levels + 1) * sizeof(int));
  assert(initst[0] != NULL);
  // all start at state 0
  memset(initst[0], 0, (num_levels + 1) * sizeof(int));
}


void CheckVars(domain *d)
{
  // Make sure that the level handles are the same as the level heights.
  // If this is not the case, the program needs to be modified to be
  // more flexible.
  int nVars = d->getNumVariables();
  int currHeight = nVars;
  int topVar = d->getNumVariables();
  while (topVar > 0)
  {
    if (topVar != currHeight) {
      // complain
      exit(1);
    }
    topVar--;
    currHeight--;
  }
}

void PrintIntArray(int p[], const int p_size)
{
  printf("[");
  for (int i=0; i<p_size; i++) printf("%d ", p[i]);
  printf("]\n");
}

void show_num_nodes(const dd_edge& e)
{
  static_cast<expert_forest*>(e.getForest())->removeAllComputeTableEntries();

  printf("#Nodes: %d\n", e.getForest()->getCurrentNumNodes());
  FILE_output out(stdout);
  e.show(out, 2);
}

dd_edge BuildMoveHelper(
  forest* relation,
  int type3_a,
  int type2_a,
  int type3_b,
  int type2_b,
  int type3_c,
  int type2_c,
  int type3_d,
  int type2_d
  )
{
  assert(relation->isForRelations());

  // transform to variables
  int a2 = get_component_level(2, type2_a);
  int b2 = get_component_level(2, type2_b);
  int c2 = get_component_level(2, type2_c);
  int d2 = get_component_level(2, type2_d);
  int a3 = get_component_level(3, type3_a);
  int b3 = get_component_level(3, type3_b);
  int c3 = get_component_level(3, type3_c);
  int d3 = get_component_level(3, type3_d);

//  fprintf(stderr, "type2_a, a2 = %d, %d\n", type2_a, a2);
//  fprintf(stderr, "type2_b, b2 = %d, %d\n", type2_b, b2);
//  fprintf(stderr, "type2_c, c2 = %d, %d\n", type2_c, c2);
//  fprintf(stderr, "type2_d, d2 = %d, %d\n", type2_d, d2);
//  fprintf(stderr, "type3_a, a3 = %d, %d\n", type3_a, a3);
//  fprintf(stderr, "type3_b, b3 = %d, %d\n", type3_b, b3);
//  fprintf(stderr, "type3_c, c3 = %d, %d\n", type3_c, c3);
//  fprintf(stderr, "type3_d, d3 = %d, %d\n", type3_d, d3);

  // face is ordered like this:
  // type3, type2, type3, type2, type3, type2, type3, type2

  const int sz = num_levels + 1;

  // Adding all the elements in one go
  //
  // 4 type 2 additions = 4 * 12 = 48
  // 4 type 3 additions = 4 * 8 = 32
  // total additions = 4 * 12 + 4 * 8 = 4 * 20 = 80
  //
  int nElements = 4 * type2 + 4 * type3;
  assert(nElements == 80);
  int** from = new int*[nElements];
  int** to = new int*[nElements];
  for (int i = 0; i < nElements; i++)
  {
    // allocate elements
    from[i] = new int[sz];
    to[i] = new int [sz];

    // initialize elements
    from[i][0] = 0;
    to[i][0] = 0;

    std::fill_n(&from[i][1], sz - 1, DONT_CARE);
    std::fill_n(&to[i][1], sz - 1, DONT_CHANGE);

    to[i][a3] = DONT_CARE;
    to[i][b3] = DONT_CARE;
    to[i][c3] = DONT_CARE;
    to[i][d3] = DONT_CARE;
    to[i][a2] = DONT_CARE;
    to[i][b2] = DONT_CARE;
    to[i][c2] = DONT_CARE;
    to[i][d2] = DONT_CARE;
  }

  int currElement = 0;

  // a3' <- d3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][d3] = i;
    to[currElement][a3] = i;
    currElement++;
  }

  // b3' <- a3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][a3] = i;
    to[currElement][b3] = i;
    currElement++;
  }

  // c3' <- b3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][b3] = i;
    to[currElement][c3] = i;
    currElement++;
  }

  // d3' <- c3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][c3] = i;
    to[currElement][d3] = i;
    currElement++;
  }


  // a2' <- d2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][d2] = i;
    to[currElement][a2] = i;
    currElement++;
  }

  // b2' <- a2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][a2] = i;
    to[currElement][b2] = i;
    currElement++;
  }

  // c2' <- b2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][b2] = i;
    to[currElement][c2] = i;
    currElement++;
  }

  // d2' <- c2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][c2] = i;
    to[currElement][d2] = i;
    currElement++;
  }

  // compute result = union elements to create an mxd for each component
  // and then intersect the mxds.

  dd_edge result(relation);
  dd_edge temp(relation);
  int offset = 0;

  // a3' <- d3

  relation->createEdge(from + offset, to + offset, type3, result);
  offset += type3;

  // b3' <- a3
  relation->createEdge(from + offset, to + offset, type3, temp);
  result *= temp;
  offset += type3;

  // c3' <- b3
  relation->createEdge(from + offset, to + offset, type3, temp);
  result *= temp;
  offset += type3;

  // d3' <- c3
  relation->createEdge(from + offset, to + offset, type3, temp);
  result *= temp;
  offset += type3;

  // a2' <- d2
  relation->createEdge(from + offset, to + offset, type2, temp);
  result *= temp;
  offset += type2;

  // b2' <- a2
  relation->createEdge(from + offset, to + offset, type2, temp);
  result *= temp;
  offset += type2;

  // c2' <- b2
  relation->createEdge(from + offset, to + offset, type2, temp);
  result *= temp;
  offset += type2;

  // d2' <- c2
  relation->createEdge(from + offset, to + offset, type2, temp);
  result *= temp;
  offset += type2;

  assert(offset == nElements);

  // delete arrays
  for (int i = 0; i < nElements; i++)
  {
    delete[] from[i];
    delete[] to[i];
  }
  delete[] from;
  delete[] to;
  
  return result;
}


// Modified
dd_edge BuildFlipMoveHelper(
  forest* relation,
  int type3_a,
  int type2_a,
  int type3_b,
  int type2_b,
  int type3_c,
  int type2_c,
  int type3_d,
  int type2_d
  )
{
  assert(relation->isForRelations());

  // transform to levels
  int a2 = get_component_level(2, type2_a);
  int b2 = get_component_level(2, type2_b);
  int c2 = get_component_level(2, type2_c);
  int d2 = get_component_level(2, type2_d);
  int a3 = get_component_level(3, type3_a);
  int b3 = get_component_level(3, type3_b);
  int c3 = get_component_level(3, type3_c);
  int d3 = get_component_level(3, type3_d);

//  fprintf(stderr, "type2_a, a2 = %d, %d\n", type2_a, a2);
//  fprintf(stderr, "type2_b, b2 = %d, %d\n", type2_b, b2);
//  fprintf(stderr, "type2_c, c2 = %d, %d\n", type2_c, c2);
//  fprintf(stderr, "type2_d, d2 = %d, %d\n", type2_d, d2);
//  fprintf(stderr, "type3_a, a3 = %d, %d\n", type3_a, a3);
//  fprintf(stderr, "type3_b, b3 = %d, %d\n", type3_b, b3);
//  fprintf(stderr, "type3_c, c3 = %d, %d\n", type3_c, c3);
//  fprintf(stderr, "type3_d, d3 = %d, %d\n", type3_d, d3);

  // face is ordered like this:
  // type3, type2, type3, type2, type3, type2, type3, type2

  const int sz = num_levels + 1;

  // Adding all the elements in one go
  //
  // 4 type 2 additions = 4 * 12 = 48
  // 4 type 3 additions = 4 * 8 = 32
  // total additions = 4 * 12 + 4 * 8 = 4 * 20 = 80
  //
  int nElements = 4 * type2 + 4 * type3;
  assert(nElements == 80);
  int** from = (int**) malloc(nElements * sizeof(int *));
  int** to = (int**) malloc(nElements * sizeof(int *));
  for (int i = 0; i < nElements; i++)
  {
    // allocate elements
    from[i] = (int*) malloc(sz * sizeof(int));
    to[i] = (int*) malloc(sz * sizeof(int));

    // initialize elements
    from[i][0] = 0;
    to[i][0] = 0;

    std::fill_n(&from[i][1], sz - 1, DONT_CARE);
    std::fill_n(&to[i][1], sz - 1, DONT_CHANGE);

    to[i][a3] = DONT_CARE;
    to[i][b3] = DONT_CARE;
    to[i][c3] = DONT_CARE;
    to[i][d3] = DONT_CARE;
    to[i][a2] = DONT_CARE;
    to[i][b2] = DONT_CARE;
    to[i][c2] = DONT_CARE;
    to[i][d2] = DONT_CARE;
  }


  int currElement = 0;

  // a3' <- c3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][c3] = i;
    to[currElement][a3] = i;
    currElement++;
  }

  // b3' <- d3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][d3] = i;
    to[currElement][b3] = i;
    currElement++;
  }

  // c3' <- a3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][a3] = i;
    to[currElement][c3] = i;
    currElement++;
  }

  // d3' <- b3

  for (int i = 0; i < type3; i++)
  {
    from[currElement][b3] = i;
    to[currElement][d3] = i;
    currElement++;
  }


  // a2' <- c2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][c2] = i;
    to[currElement][a2] = i;
    currElement++;
  }

  // b2' <- d2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][d2] = i;
    to[currElement][b2] = i;
    currElement++;
  }

  // c2' <- a2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][a2] = i;
    to[currElement][c2] = i;
    currElement++;
  }

  // d2' <- b2

  for (int i = 0; i < type2; i++)
  {
    from[currElement][b2] = i;
    to[currElement][d2] = i;
    currElement++;
  }


  // compute result = union elements to create an mxd for each component
  // and then intersect the mxds.

  dd_edge result(relation);
  int offset = 0;

  // a3' <- c3
  {
    dd_edge temp(relation);
    relation->createEdge(from + offset, to + offset, type3, temp);
    result = temp;
    offset += type3;
  }

  // b3' <- d3
  // c3' <- a3
  // d3' <- b3
  for (int i = 0; i < 3; i++)
  {
    dd_edge temp(relation);
    relation->createEdge(from + offset, to + offset, type3, temp);
    result *= temp;
    offset += type3;
  }

  // a2' <- c2
  // b2' <- d2
  // c2' <- a2
  // d2' <- b2
  for (int i = 0; i < 4; i++)
  {
    dd_edge temp(relation);
    relation->createEdge(from + offset, to + offset, type2, temp);
    result *= temp;
    offset += type2;
  }

  assert(offset == nElements);

  // delete arrays
  for (int i = 0; i < nElements; i++)
  {
    free(from[i]);
    free(to[i]);
  }
  free(from);
  free(to);
  
  return result;
}


dd_edge BuildMove(forest* relation, Face f, Direction d) {
  assert(relation->isForRelations());

  dd_edge result(relation);
  switch (f) {
    case Face::U:
      if (d == Direction::CW) {
        result = BuildMoveHelper(relation, 4, 5, 5, 6, 1, 0, 0, 4);
      } else if (d == Direction::CCW) {
        result = BuildMoveHelper(relation, 4, 4, 0, 0, 1, 6, 5, 5);
      } else {
        result = BuildFlipMoveHelper(relation, 4, 5, 5, 6, 1, 0, 0, 4);
      }
      break;
    case Face::D:
      if (d == Direction::CW) {
        result =  BuildMoveHelper(relation, 3, 2, 2, 8, 6, 9, 7, 10);
      } else if (d == Direction::CCW) {
        result =  BuildMoveHelper(relation, 3, 10, 7, 9, 6, 8, 2, 2);
      } else {
        result =  BuildFlipMoveHelper(relation, 3, 2, 2, 8, 6, 9, 7, 10);
      }
      break;
    case Face::L:
      if (d == Direction::CW) {
        result =  BuildMoveHelper(relation, 4, 4, 0, 3, 3, 10, 7, 11);
      } else if (d == Direction::CCW) {
        result =  BuildMoveHelper(relation, 4, 11, 7, 10, 3, 3, 0, 4);
      } else {
        result =  BuildFlipMoveHelper(relation, 4, 4, 0, 3, 3, 10, 7, 11);
      }
      break;
    case Face::R:
      if (d == Direction::CW) {
        result =  BuildMoveHelper(relation, 1, 6, 5, 7, 6, 8, 2, 1);
      } else if (d == Direction::CCW) {
        result =  BuildMoveHelper(relation, 1, 1, 2, 8, 6, 7, 5, 6);
      } else {
        result =  BuildFlipMoveHelper(relation, 1, 6, 5, 7, 6, 8, 2, 1);
      }
      break;
    case Face::F:
      if (d == Direction::CW) {
        result = BuildMoveHelper(relation, 0, 0, 1, 1, 2, 2, 3, 3);
      } else if (d == Direction::CCW) {
        result = BuildMoveHelper(relation, 0, 3, 3, 2, 2, 1, 1, 0);
      } else {
        result = BuildFlipMoveHelper(relation, 0, 0, 1, 1, 2, 2, 3, 3);
      }
      break;
    case Face::B:
      if (d == Direction::CW) {
        result =  BuildMoveHelper(relation, 5, 5, 4, 11, 7, 9, 6, 7);
      } else if (d == Direction::CCW) {
        result =  BuildMoveHelper(relation, 5, 7, 6, 9, 7, 11, 4, 5);
      } else {
        result =  BuildFlipMoveHelper(relation, 5, 5, 4, 11, 7, 9, 6, 7);
      }
      break;
  }
  return result;
}


const char* face_to_string(Face f){
  switch(f) {
    case Face::F: return "Front";
    case Face::B: return "Back";
    case Face::L: return "Left";
    case Face::R: return "Right";
    case Face::U: return "Up";
    case Face::D: return "Down";
    default: return "Invalid Face";
  }
}


void usage() {
  fprintf(stderr, "Usage: rubik_cube [-bfs|-dfs|-l<key>|-p]\n");
  fprintf(stderr, "-bfs   : use traditional algorithm to compute reachable states\n");
  fprintf(stderr, "-msat  : use saturation with monolithic relation compute reachable states\n");
  fprintf(stderr, "-esat  : use saturation with event-wise relation compute reachable states\n");
  fprintf(stderr, "-ksat  : use saturation with level-wise relation compute reachable states\n");
  fprintf(stderr, "-kspsat: same as ksat, with additional pre-processing of events to improve saturation\n");
  fprintf(stderr, "-l<key>: key can be any combination of\n");
  fprintf(stderr, "         A: Front face clock-wise rotation,\n");
  fprintf(stderr, "         a: Front face counter clock-wise rotation,\n");
  fprintf(stderr, "         1: Front face flip,\n");
  fprintf(stderr, "         Back face (B, b, 2),\n");
  fprintf(stderr, "         Left face (C, c, 3),\n");
  fprintf(stderr, "         Right face (D, d, 4),\n");
  fprintf(stderr, "         Up face (E, e, 5),\n");
  fprintf(stderr, "         Down face (F, f, 6),\n");
  fprintf(stderr, "\n");
}


class moves {
  public:

    // ?CW, ?CCW, ?F stand for clock-wise, counter clock-wise and flip resp.

    bool FCW;
    bool FCCW;
    bool FF;
    bool BCW;
    bool BCCW;
    bool BF;
    bool LCW;
    bool LCCW;
    bool LF;
    bool RCW;
    bool RCCW;
    bool RF;
    bool UCW;
    bool UCCW;
    bool UF;
    bool DCW;
    bool DCCW;
    bool DF;

    moves()
      : FCW(false), FCCW(false), FF(false), BCW(false), BCCW(false),
      BF(false), LCW(false), LCCW(false), LF(false), RCW(false),
      RCCW(false), RF(false), UCW(false), UCCW(false), UF(false),
      DCW(false), DCCW(false), DF(false) {}

    int countEnabledMoves() const {
      int count = 0;
      if (FCW) count++;
      if (BCW) count++;
      if (UCW) count++;
      if (DCW) count++;
      if (LCW) count++;
      if (RCW) count++;
      if (FCCW) count++;
      if (BCCW) count++;
      if (UCCW) count++;
      if (DCCW) count++;
      if (LCCW) count++;
      if (RCCW) count++;
      if (FF) count++;
      if (BF) count++;
      if (UF) count++;
      if (DF) count++;
      if (LF) count++;
      if (RF) count++;
      return count;
    }
};

void show_order(int* level2var)
{
  for (int i = 1; i <= num_components; i++) {
    printf("%d, ", level2var[i]);
  }
  printf("\n");
}

void generate_variable_orders()
{
  for (int i = 0; i< num_phases; i++) {
    for (int j = 1 ; j <= num_components; j++) {
      Component& comp = level2components[i][j];
      level2vars[i][j] = get_component_level(comp.type, comp.id);
    }
  }

//  // F and B
//  orders[1][get_component_level(3, 0)] = 1;
//  orders[1][get_component_level(3, 1)] = 2;
//  orders[1][get_component_level(3, 2)] = 3;
//  orders[1][get_component_level(3, 3)] = 4;
//  orders[1][get_component_level(2, 0)] = 5;
//  orders[1][get_component_level(2, 1)] = 6;
//  orders[1][get_component_level(2, 2)] = 7;
//  orders[1][get_component_level(2, 3)] = 8;
//
//  orders[1][get_component_level(3, 5)] = 9;
//  orders[1][get_component_level(3, 4)] = 10;
//  orders[1][get_component_level(3, 7)] = 11;
//  orders[1][get_component_level(3, 6)] = 12;
//  orders[1][get_component_level(2, 5)] = 13;
//  orders[1][get_component_level(2, 11)] = 14;
//  orders[1][get_component_level(2, 9)] = 15;
//  orders[1][get_component_level(2, 7)] = 16;
//
//  orders[1][get_component_level(2, 4)] = 17;
//  orders[1][get_component_level(2, 10)] = 18;
//  orders[1][get_component_level(2, 8)] = 19;
//  orders[1][get_component_level(2, 6)] = 20;
//
//  // L and R
//  orders[2][get_component_level(3, 4)] = 1;
//  orders[2][get_component_level(3, 0)] = 2;
//  orders[2][get_component_level(3, 3)] = 3;
//  orders[2][get_component_level(3, 7)] = 4;
//  orders[2][get_component_level(2, 4)] = 5;
//  orders[2][get_component_level(2, 3)] = 6;
//  orders[2][get_component_level(2, 10)] = 7;
//  orders[2][get_component_level(2, 11)] = 8;
//
//  orders[2][get_component_level(3, 1)] = 9;
//  orders[2][get_component_level(3, 5)] = 10;
//  orders[2][get_component_level(3, 6)] = 11;
//  orders[2][get_component_level(3, 2)] = 12;
//  orders[2][get_component_level(2, 6)] = 13;
//  orders[2][get_component_level(2, 7)] = 14;
//  orders[2][get_component_level(2, 8)] = 15;
//  orders[2][get_component_level(2, 1)] = 16;
//
//  orders[2][get_component_level(2, 5)] = 17;
//  orders[2][get_component_level(2, 9)] = 18;
//  orders[2][get_component_level(2, 2)] = 19;
//  orders[2][get_component_level(2, 0)] = 20;
//
//  // U and D
//  orders[0][get_component_level(3, 4)] = 1;
//  orders[0][get_component_level(3, 5)] = 2;
//  orders[0][get_component_level(3, 1)] = 3;
//  orders[0][get_component_level(3, 0)] = 4;
//  orders[0][get_component_level(2, 5)] = 5;
//  orders[0][get_component_level(2, 6)] = 6;
//  orders[0][get_component_level(2, 0)] = 7;
//  orders[0][get_component_level(2, 4)] = 8;
//
//  orders[0][get_component_level(3, 3)] = 9;
//  orders[0][get_component_level(3, 2)] = 10;
//  orders[0][get_component_level(3, 6)] = 11;
//  orders[0][get_component_level(3, 7)] = 12;
//  orders[0][get_component_level(2, 2)] = 13;
//  orders[0][get_component_level(2, 8)] = 14;
//  orders[0][get_component_level(2, 9)] = 15;
//  orders[0][get_component_level(2, 10)] = 16;
//
//  orders[0][get_component_level(2, 11)] = 17;
//  orders[0][get_component_level(2, 3)] = 18;
//  orders[0][get_component_level(2, 1)] = 19;
//  orders[0][get_component_level(2, 7)] = 20;

  for (int i = 0; i < 3; i++) {
    show_order(level2vars[i]);
  }
}

static int phase = 0;

void execute_phase(const dd_edge& initial, const dd_edge& nsf, dd_edge& result)
{
  fprintf(stdout, "-------------------------------------------------------------------\n");
  fprintf(stdout, "Phase %d:\n", phase);
  phase++;

  expert_forest* relation = static_cast<expert_forest*>(nsf.getForest());
  int* level2var = new int[relation->getNumVariables() + 1];
  relation->getVariableOrder(level2var);

  fprintf(stdout, "Reordering...\n");

  timer start;

  start.note_time();
  static_cast<expert_forest*>(initial.getForest())->reorderVariables(level2var);
  delete[] level2var;
  start.note_time();

  fprintf(stdout, "Time: %.4e seconds\n", start.get_last_interval()/1000000.0);

  fprintf(stdout, "Computing the reachable states...\n");
  show_num_nodes(initial);

  start.note_time();
  apply(REACHABLE_STATES_DFS, initial, nsf, result);
  start.note_time();

//  show_num_nodes(result);
  fprintf(stdout, "Time: %.4e seconds\n", start.get_last_interval()/1000000.0);
}

int phased_saturate()
{
  fprintf(stdout, "Constructing the relations...\n");

  generate_variable_orders();
  // Reorder the relations
  static_cast<expert_forest*>(relations[0])->reorderVariables(&level2vars[0][0]);
  static_cast<expert_forest*>(relations[1])->reorderVariables(&level2vars[1][0]);
  static_cast<expert_forest*>(relations[2])->reorderVariables(&level2vars[2][0]);

//  {
//    dd_edge e(relations[2]);
//    for (const auto& action : actions[2]) {
//      e += BuildMove(relations[2], action.face, action.direction);
//    }
//    static_cast<expert_forest*>(relations[2])->removeAllComputeTableEntries();
//  }

  vector<dd_edge> nsfs;
  for (int i = 0; i < num_phases; i++) {
//  for (int i = num_phases - 1; i >= 0; i--) {
    dd_edge e(relations[i]);
    for (const auto& action : actions[i]) {
      e += BuildMove(relations[i], action.face, action.direction);
    }
    nsfs.push_back(e);

    show_num_nodes(e);
  }

//  show_num_nodes(nsfs[0]);
  static_cast<expert_forest*>(relations[0])->dynamicReorderVariables(16, 1);
////  show_num_nodes(nsfs[0]);
//  static_cast<expert_forest*>(relations[0])->dynamicReorderVariables(8, 5);
////  show_num_nodes(nsfs[0]);
//  static_cast<expert_forest*>(relations[0])->dynamicReorderVariables(12, 9);
////  show_num_nodes(nsfs[0]);
//  static_cast<expert_forest*>(relations[0])->dynamicReorderVariables(16, 13);
  show_num_nodes(nsfs[0]);

//  show_num_nodes(nsfs[1]);
//  static_cast<expert_forest*>(relations[1])->dynamicReorderVariables(4, 1);
//  static_cast<expert_forest*>(relations[1])->dynamicReorderVariables(8, 5);
//  static_cast<expert_forest*>(relations[1])->dynamicReorderVariables(12, 9);
//  static_cast<expert_forest*>(relations[1])->dynamicReorderVariables(16, 13);
//  show_num_nodes(nsfs[1]);

//  show_num_nodes(nsfs[2]);
//  static_cast<expert_forest*>(relations[2])->dynamicReorderVariables(4, 1);
//  static_cast<expert_forest*>(relations[2])->dynamicReorderVariables(8, 5);
//  static_cast<expert_forest*>(relations[2])->dynamicReorderVariables(12, 9);
//  static_cast<expert_forest*>(relations[2])->dynamicReorderVariables(16, 13);
//  show_num_nodes(nsfs[2]);

  // Build the initial states
  assert(states);
  dd_edge initial(states);
  states->createEdge(initst, 1, initial);

  timer start;
  start.note_time();

  // Fixed point
  bool fp = false;
  dd_edge result = initial;
  while (!fp) {
    fp = true;
    for (int i = 0; i < num_phases; i++) {
      execute_phase(initial, nsfs[i], result);
      if (result != initial) {
        fp = false;
        initial = result;
      }
      show_num_nodes(result);
    }
  }

  start.note_time();
  fprintf(stdout, " done!\n");
  fflush(stdout);
  fprintf(stdout, "Time for constructing reachability set: %.4e seconds\n",
      start.get_last_interval()/1000000.0);
  fprintf(stdout, "# of reachable states: %1.6e\n",
      initial.getCardinality());
  fflush(stdout);

  return 0;
}

int main(int argc, char *argv[])
{
  Init();

  // Initialize MEDDLY
  MEDDLY::settings s;
  s.ctSettings.style = MonolithicChainedHash;
  s.ctSettings.maxSize = 16 * 16777216;
  // s.ctSettings.staleRemoval =
  //   MEDDLY::settings::computeTableSettings::Lazy;
  // s.ctSettings.staleRemoval =
  //   MEDDLY::settings::computeTableSettings::Moderate;
  // s.ctSettings.staleRemoval =
  //   MEDDLY::settings::computeTableSettings::Aggressive;

  MEDDLY::initialize(s);

  // Set up the state variables, as described earlier
  d = createDomainBottomUp(sizes, num_levels);
  if (d == nullptr) {
    fprintf(stderr, "Couldn't create domain\n");
    return 1;
  }
  CheckVars(d);

  // Create forests
  forest::policies p(false);
  p.setSinkDown();
  states = d->createForest(false, forest::BOOLEAN, forest::MULTI_TERMINAL, p);
  relations.push_back(d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL));
  relations.push_back(d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL));
  relations.push_back(d->createForest(true, forest::BOOLEAN, forest::MULTI_TERMINAL));

  if (states == nullptr) {
    fprintf(stderr, "Couldn't create forest of states\n");
    return 1;
  } else {
    fprintf(stderr, "Created forest of states\n");
  }

  for (const auto& relation : relations) {
    if (relation == nullptr) {
      fprintf(stderr, "Couldn't create forest of relations\n");
      return 1;
    }
  }

  fprintf(stderr, "Created forest of relations\n");

  // Build set of initial states

  for (int j = 0; j < type2; j++) {
    initst[0][get_component_level(2, j)] = j;
  }
  for (int j = 0; j < type3; j++) {
    initst[0][get_component_level(3, j)] = j;
  }
  initst[0][0] = 0;

  phased_saturate();

  destroyDomain(d);
  cleanup();
  fprintf(stderr, "\n\nDONE\n");
  return 0;
}


