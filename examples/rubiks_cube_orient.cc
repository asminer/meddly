
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

	Up:     In order (going right starting from front face left-upper corner)
          of components (Type:id),
          (3:0, 2:0, 3:1, 2:1, 3:2, 2:2, 3:3, 2:3)
          Note: (1:0) also belongs to this row but since it stays in the
          same location when this row is moved left or right, it is ignored.
	Down:   (3:4, 2:8, 3:5, 2:9, 3:6, 2:10, 3:7, 2:11) (1:5 ignored)
	Left:   (3:0, 2:4, 3:4, 2:11, 3:7, 2:7, 3:3, 2:3) (1:4 ignored)
	Right:  (3:1, 2:5, 3:5, 2:9, 3:6, 2:6, 3:2, 2:1) (1:2 ignored)
	Front:  (3:0, 2:0, 3:1, 2:5, 3:5, 2:8, 3:4, 2:4) (1:1 ignored)
	Back:   (3:3, 2:2, 3:2, 2:6, 3:6, 2:10, 3:7, 2:7) (1:3 ignored)

	Initially components are placed in components locations that match their
	Ids.
*/


#define SATURATE 1
#define NO_CHOICE 1
#define COUNT_TIME 1
#define COUNT_STATES 1

#define REORDER_CUBE 0
#define FACED_ORDERING 0

#define FROM_TO 0

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include <domain.h>
#include <forest.h>
#include <dd_edge.h>
#include <ophandle.h>

typedef enum {F, B, L, R, U, D} face;
typedef enum {CW, CCW} direction;

level *variables = NULL;
int *sizes = NULL;
int *initst = NULL;

// Domain handle
domain *d;

// Forest storing the next state function
forest_hndl relation;

// Forest storing the set of states
forest_hndl states;

// Edge representing the initial set of states
dd_edge *initial;

// Edge representing the next state function
dd_edge *nsf;

// Number of variables of each type
const int type1 = 6;
const int type2 = 12;
const int type3 = 8;

// Number of levels
const int num_levels = type2 + type3;


int get_component_type (int comp_level)
{
  assert(comp_level > 0 && comp_level <= num_levels);

#if FACED_ORDERING
  // type3 = 1, 3, 5, 7, 10, 12, 15, 18
  // type2 = the rest; except 0
  static int comp_type[] =
      {0, 3, 2, 3, 2, 3, 2, 3, 2, 2, 3, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2};
#else
  
#if REORDER_CUBE
  static int comp_type[] =
      {0, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
#else
  static int comp_type[] =
      {0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3};
#endif

#endif
  return comp_type[comp_level];
}

int get_component_size (int comp_type)
{
  assert(comp_type > 0 && comp_type < 4);  // 1..3 depending on # of facets
  assert(comp_type != 1);
  return (comp_type == 3)? type3 * 3: type2 * 2;
}

int getLevelSize (int comp_level)
{
  return get_component_size(get_component_type(comp_level));
}

int get_component_level (int comp_type, int index)
{
  assert(comp_type > 0 && comp_type < 4);  // 1..3 depending on # of facets
  assert(comp_type != 1);
  assert((comp_type == 3 && index >= 0 && index <= 7) ||
      (comp_type == 2 && index >= 0 && index <= 11));

#if FACED_ORDERING
  static int comp3_map[] = {1, 3, 5, 7, 10, 12, 15, 18};
  static int comp2_map[] = {2, 4, 6, 8, 9, 11, 13, 14, 16, 17, 19, 20};
#else

#if REORDER_CUBE
  static int comp3_map[] = {1, 2, 3, 4, 5, 6, 7, 8};
  static int comp2_map[] = {9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
#else
  static int comp3_map[] = {13, 14, 15, 16, 17, 18, 19, 20};
  static int comp2_map[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
#endif

#endif

  return (comp_type == 3)? comp3_map[index]: comp2_map[index];
}

void Init()
{
  assert(num_levels == (type2 + type3));

  // store for level handles
  variables = (level *) malloc((num_levels + 1) * sizeof(level));
  assert(variables != NULL);
  memset(variables, 0, (num_levels + 1) * sizeof(level));
  
  // node size for each level
  sizes = (int *) malloc((num_levels + 1) * sizeof(int));
  assert(sizes != NULL);
  
  assert(num_levels == 20);
  sizes[0] = 0;
  for (int i = 1; i < (num_levels + 1); i++) {
    sizes[i] = getLevelSize(i);
    fprintf(stderr, "sizes[%d] = %d\n", i, sizes[i]);
  }
  fflush(stderr);

  // sets of states
  initst = (int *) malloc((num_levels + 1) * sizeof(int));
  assert(initst != NULL);
  // all start at state 0
  memset(initst, 0, (num_levels + 1) * sizeof(int));
}

void CheckVars()
{
  // Sanity check
  for (int i = num_levels; i > 0; i--) 
    if ((variables[i] > num_levels) || (variables[i] < 1)) {
      fprintf(stderr, "Level handle for variable %d is %d, ",
          i, variables[i]);
      fprintf(stderr, "outside of expected range\n");
      exit(1);
    }
}

void SetIntArray(int p[], const int p_size, const int c)
{
  for (int i=0; i<p_size; i++) p[i] = c;
}

int plus1_mod2[] = {1, 0};
int plus1_mod3[] = {1, 2, 0};
int plus2_mod3[] = {2, 0, 1};

// type3: along Z axis
// type2: along Y axis; then along Z axis
//
// type3: along Z axis
// type3:
//   F : no change
//   B : no change
//   L : +1 (3, 4), +2 (0, 7)
//   R : +1 (1, 6), +2 (2, 5)
//   U : +1 (0, 5), +2 (1, 4)
//   D : +1 (2, 7), +2 (3, 6)
//
// type2: along Y axis; then along Z axis
// type2:
//   F : +1
//   B : +1
//   L : no change
//   R : no change
//   U : no change
//   D : no change
int chor_type2 (face f, int from, int curr_or) {
  return (f == F || f == B)? plus1_mod2[curr_or]: curr_or;
}

// type3: along Z axis
// type3:
//   F : no change
//   B : no change
//   L : +1 (3, 4), +2 (0, 7)
//   R : +1 (1, 6), +2 (2, 5)
//   U : +1 (0, 5), +2 (1, 4)
//   D : +1 (2, 7), +2 (3, 6)
//
int chor_type3 (face f, int from, int curr_or) {
  switch (f) {
    case L:
      //   L : +1 (3, 4), +2 (0, 7)
      return (from == 3 || from == 4)?
        plus1_mod3[curr_or]: plus2_mod3[curr_or];
    case R:
      //   R : +1 (1, 6), +2 (2, 5)
      return (from == 1 || from == 6)?
        plus1_mod3[curr_or]: plus2_mod3[curr_or];
    case U:
      //   U : +1 (0, 5), +2 (1, 4)
      return (from == 0 || from == 5)?
        plus1_mod3[curr_or]: plus2_mod3[curr_or];
    case D:
      //   D : +1 (2, 7), +2 (3, 6)
      return (from == 2 || from == 7)?
        plus1_mod3[curr_or]: plus2_mod3[curr_or];
    case F:
    case B:
      //   F : no change
      //   B : no change
    default:
      return curr_or;
  }
}

dd_edge* DoMoveHelper(
  face f,
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
  const int sz = num_levels + 1;
  int from[sz];  // num_levels is a constant
  int to[sz];

  // face is ordered like this:
  // type3, type2, type3, type2, type3, type2, type3, type2

  // transform to levels
  int a2 = get_component_level(2, type2_a);
  int b2 = get_component_level(2, type2_b);
  int c2 = get_component_level(2, type2_c);
  int d2 = get_component_level(2, type2_d);
  int a3 = get_component_level(3, type3_a);
  int b3 = get_component_level(3, type3_b);
  int c3 = get_component_level(3, type3_c);
  int d3 = get_component_level(3, type3_d);

  fprintf(stderr, "type2_a, a2 = %d, %d\n", type2_a, a2);
  fprintf(stderr, "type2_b, b2 = %d, %d\n", type2_b, b2);
  fprintf(stderr, "type2_c, c2 = %d, %d\n", type2_c, c2);
  fprintf(stderr, "type2_d, d2 = %d, %d\n", type2_d, d2);
  fprintf(stderr, "type3_a, a3 = %d, %d\n", type3_a, a3);
  fprintf(stderr, "type3_b, b3 = %d, %d\n", type3_b, b3);
  fprintf(stderr, "type3_c, c3 = %d, %d\n", type3_c, c3);
  fprintf(stderr, "type3_d, d3 = %d, %d\n", type3_d, d3);
  fflush(stderr);

  // create node at level 13
  dd_tempedge *temp13 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type3; i++) {
    for (int j=0; j<3; j++) {
      from[variables[d3]] = i + (j * type3);
      to[variables[a3]] = i + (chor_type3(f, type3_d, j) * type3);
#if FROM_TO
      AddMatrixElement(temp13, from, to, sz, true);
#else
      AddMatrixElement(temp13, to, from, sz, true);
#endif
    }
  }

  dd_edge *e13 = CreateEdge(temp13);

  // create node at level 14
  dd_tempedge *temp14 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type3; i++) {
    for (int j=0; j<3; j++) {
      from[variables[a3]] = i + (j * type3);
      to[variables[b3]] = i + (chor_type3(f, type3_a, j) * type3);
#if FROM_TO
      AddMatrixElement(temp14, from, to, sz, true);
#else
      AddMatrixElement(temp14, to, from, sz, true);
#endif
    }
  }

  dd_edge *e14 = CreateEdge(temp14);

  // create node at level 15
  dd_tempedge *temp15 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type3; i++) {
    for (int j=0; j<3; j++) {
      from[variables[b3]] = i + (j * type3);
      to[variables[c3]] = i + (chor_type3(f, type3_b, j) * type3);
#if FROM_TO
      AddMatrixElement(temp15, from, to, sz, true);
#else
      AddMatrixElement(temp15, to, from, sz, true);
#endif
    }
  }

  dd_edge *e15 = CreateEdge(temp15);

  // create node at level 16
  dd_tempedge *temp16 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type3; i++) {
    for (int j=0; j<3; j++) {
      from[variables[c3]] = i + (j * type3);
      to[variables[d3]] = i + (chor_type3(f, type3_c, j) * type3);
#if FROM_TO
      AddMatrixElement(temp16, from, to, sz, true);
#else
      AddMatrixElement(temp16, to, from, sz, true);
#endif
    }
  }

  dd_edge *e16 = CreateEdge(temp16);
  
  // create node at level 1
  dd_tempedge *temp1 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type2; i++) {
    for (int j=0; j<2; j++) {
      from[variables[d2]] = i + (j * type2);
      to[variables[a2]] = i + (chor_type2(f, type2_d, j) * type2);
#if FROM_TO
      AddMatrixElement(temp1, from, to, sz, true);
#else
      AddMatrixElement(temp1, to, from, sz, true);
#endif
    }
  }

  dd_edge *e1 = CreateEdge(temp1);

  // create node at level 2
  dd_tempedge *temp2 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type2; i++) {
    for (int j=0; j<2; j++) {
      from[variables[a2]] = i + (j * type2);
      to[variables[b2]] = i + (chor_type2(f, type2_a, j) * type2);
#if FROM_TO
      AddMatrixElement(temp2, from, to, sz, true);
#else
      AddMatrixElement(temp2, to, from, sz, true);
#endif
    }
  }

  dd_edge *e2 = CreateEdge(temp2);

  // create node at level 3
  dd_tempedge *temp3 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type2; i++) {
    for (int j=0; j<2; j++) {
      from[variables[b2]] = i + (j * type2);
      to[variables[c2]] = i + (chor_type2(f, type2_b, j) * type2);
#if FROM_TO
      AddMatrixElement(temp3, from, to, sz, true);
#else
      AddMatrixElement(temp3, to, from, sz, true);
#endif
    }
  }

  dd_edge *e3 = CreateEdge(temp3);
  
  // create node at level 4
  dd_tempedge *temp4 = CreateTempEdge(relation, NULL);

  // Set all levels (except term) to don't care
  SetIntArray(from + 1, sz - 1, -2);
  SetIntArray(to + 1, sz - 1, -2);
  from[variables[a3]] = -1; to[variables[a3]] = -1;
  from[variables[b3]] = -1; to[variables[b3]] = -1;
  from[variables[c3]] = -1; to[variables[c3]] = -1;
  from[variables[d3]] = -1; to[variables[d3]] = -1;
  from[variables[a2]] = -1; to[variables[a2]] = -1;
  from[variables[b2]] = -1; to[variables[b2]] = -1;
  from[variables[c2]] = -1; to[variables[c2]] = -1;
  from[variables[d2]] = -1; to[variables[d2]] = -1;

  for (int i=0; i<type2; i++) {
    for (int j=0; j<2; j++) {
      from[variables[c2]] = i + (j * type2);
      to[variables[d2]] = i + (chor_type2(f, type2_c, j) * type2);
#if FROM_TO
      AddMatrixElement(temp4, from, to, sz, true);
#else
      AddMatrixElement(temp4, to, from, sz, true);
#endif
    }
  }

  dd_edge *e4 = CreateEdge(temp4);
  // ShowDDEdge(stderr, e4);

  dd_edge *result = NULL;

#if 0
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, e1, e4, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e3, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e2, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e13, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e16, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e15, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e14, result));
#else    
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, e1, e2, result));
#if 0
  fprintf(stderr, "e1: from [%d] to [%d]\n", variables[d2], variables[a2]);
  ShowDDEdge(stderr, e1);
  fprintf(stderr, "e2: from [%d] to [%d]\n", variables[a2], variables[b2]);
  ShowDDEdge(stderr, e2);
  fprintf(stderr, "e1 union e2:\n");
  ShowDDEdge(stderr, result);
  fprintf(stderr, "\n");
  exit(0);
#endif
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e3, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e4, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e13, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e14, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e15, result));
  assert(SUCCESS == ApplyBinary(OP_INTERSECTION, result, e16, result));
#if 0
  ReleaseEdge(e1);
  ReleaseEdge(e2);
  ReleaseEdge(e3);
  ReleaseEdge(e4);
  ReleaseEdge(e13);
  ReleaseEdge(e14);
  ReleaseEdge(e15);
  ReleaseEdge(e16);
#endif
#endif
  return result;
}

# if 0
dd_edge* DoMoveHelper(int a3, int b3, int c3, int d3,
    int a1, int b1, int c1, int d1);
dd_edge* DoMove(face f, direction d);

dd_edge* Front(direction dir);
dd_edge* Back(direction dir);
dd_edge* Left(direction dir);
dd_edge* Right(direction dir);
dd_edge* Up(direction dir);
dd_edge* Down(direction dir);
#endif

dd_edge* DoMove(face f, direction d) {
  dd_edge* result = NULL;
  switch (f) {
    case U:
      if (d == CW) {
        // result = DoMoveHelper(16, 15, 14, 13, 4, 3, 2, 1);
        result = DoMoveHelper(f, 4, 5, 5, 6, 1, 0, 0, 4);
      } else {
        assert(d == CCW);
        // result = DoMoveHelper(13, 14, 15, 16, 1, 2, 3, 4);
        result = DoMoveHelper(f, 0, 4, 1, 0, 5, 6, 4, 5);
      }
      break;
    case D:
      if (d == CW) {
        // result =  DoMoveHelper(20, 19, 18, 17, 12, 11, 10, 9);
        result =  DoMoveHelper(f, 3, 2, 2, 8, 6, 9, 7, 10);
      } else {
        assert(d == CCW);
        // result =  DoMoveHelper(17, 18, 19, 20, 9, 10, 11, 12);
        result =  DoMoveHelper(f, 7, 10, 6, 9, 2, 8, 3, 2);
      }
      break;
    case L:
      if (d == CW) {
        // result =  DoMoveHelper(13, 17, 20, 16, 5, 12, 8, 4);
        result =  DoMoveHelper(f, 0, 3, 3, 10, 7, 11, 4, 4);
      } else {
        assert(d == CCW);
        // result =  DoMoveHelper(16, 20, 17, 13, 4, 8, 12, 5);
        result =  DoMoveHelper(f, 4, 4, 7, 11, 3, 10, 0, 3);
      }
      break;
    case R:
      if (d == CW) {
        // result =  DoMoveHelper(14, 15, 19, 18, 2, 7, 10, 6);
        result =  DoMoveHelper(f, 1, 1, 5, 6, 6, 7, 2, 8);
      } else {
        assert(d == CCW);
        // result =  DoMoveHelper(18, 19, 15, 14, 6, 10, 7, 2);
        result =  DoMoveHelper(f, 2, 8, 6, 7, 5, 6, 1, 1);
      }
      break;
    case F:
      if (d == CW) {
        // result =  DoMoveHelper(13, 14, 18, 17, 1, 6, 9, 5);
        result = DoMoveHelper(f, 0, 0, 1, 1, 2, 2, 3, 3);
      } else {
        assert(d == CCW);
        // result =  DoMoveHelper(17, 18, 14, 13, 5, 9, 6, 1);
        result = DoMoveHelper(f, 3, 3, 2, 2, 1, 1, 0, 0);
      }
      break;
    case B:
      if (d == CW) {
        // result =  DoMoveHelper(15, 16, 20, 19, 3, 8, 11, 7);
        result =  DoMoveHelper(f, 5, 5, 4, 11, 7, 9, 6, 7);
      } else {
        assert(d == CCW);
        // result =  DoMoveHelper(19, 20, 16, 15, 7, 11, 8, 3);
        result =  DoMoveHelper(f, 6, 7, 7, 9, 4, 11, 5, 5);
      }
      break;
  }
  return result;
}

dd_edge* Front(direction dir) { return DoMove(F, dir); }
dd_edge* Back(direction dir) { return DoMove(B, dir); }
dd_edge* Left(direction dir) { return DoMove(L, dir); }
dd_edge* Right(direction dir) { return DoMove(R, dir); }
dd_edge* Up(direction dir) { return DoMove(U, dir); }
dd_edge* Down(direction dir) { return DoMove(D, dir); }

const char* face_to_string(face f){
  switch(f) {
    case F: return "Front";
    case B: return "Back";
    case L: return "Left";
    case R: return "Right";
    case U: return "Up";
    case D: return "Down";
    default: return "Invalid Face";
  }
}


void usage() {
  fprintf(stderr, "Usage: rubik_cube [-m <MB>|-dfs|-l<int>|-p|-pgif]\n");
  fprintf(stderr, "-m<MB> : sets the total memory to be used (in MB)\n");
  fprintf(stderr, "-dfs   : use depth-first algorithm to compute reachable states\n");
  fprintf(stderr, "-l<int>: combo to use {2:FU, 3:FUR, 4:FURL, 5:FURLB, 6:FURLBD}\n");
  fprintf(stderr, "-p     : prints initial states, nsf, reachable states on stderr\n");
  fprintf(stderr, "-pgif  : writes GIFs for initial states, nsf, reachable states\n");
  fprintf(stderr, "\n");
  fflush(stderr);
}

int main(int argc, char *argv[])
{
  bool make_gifs = false;
  bool pretty_print = false;
  bool no_choice = false;
  bool enable_F = false;
  bool enable_U = false;
  bool enable_R = false;
  bool enable_L = false;
  bool enable_B = false;
  bool enable_D = false;
  if (argc > 1) {
    assert(argc <= 5);
    for (int i=1; i<argc; i++) {
      char *cmd = argv[i];
      if (strncmp(cmd, "-pgif", 6) == 0) make_gifs = true;
      else if (strncmp(cmd, "-p", 3) == 0) pretty_print = true;
      else if (strncmp(cmd, "-dfs", 5) == 0) no_choice = true;
      else if (strncmp(cmd, "-l", 2) == 0) {
        for (unsigned j = 2; j < strlen(cmd); j++) {
          switch (cmd[j]) {
            case 'F':
                enable_F = true;
                break;
            case 'U':
                enable_U = true;
                break;
            case 'R':
                enable_R = true;
                break;
            case 'L':
                enable_L = true;
                break;
            case 'B':
                enable_B = true;
                break;
            case 'D':
                enable_D = true;
                break;
          }
        }
      } else if (strncmp(cmd, "-m", 2) == 0) {
        int mem_total = strtol(&cmd[2], NULL, 10);
        if (mem_total < 1 || mem_total > 100*1024) { // 10 GB!
          usage();
          exit(1);
        }
        // set up memory available
        InitMemoryManager(mem_total*1024*1024);
      }
      else {
        usage();
        exit(1);
      }
    }
  }

  // set up arrays based on number of levels
  Init();

  // Initialize MEDDLY
  initialize();

  // Set up the state variables, as described earlier
  d = CreateDomain(num_levels, sizes, variables);
  if (NULL == d) {
    fprintf(stderr, "Couldn't create domain\n");
    return 1;
  }
  CheckVars();

  // Create forests
#if 1
  // states = CreateForest(d, MDD, false, forest::QUASI_REDUCED, FULL_OR_SPARSE_STORAGE);
  states = CreateForest(d, MDD, false, forest::FULLY_REDUCED, FULL_OR_SPARSE_STORAGE);
  relation = CreateForest(d, MXD, false, forest::IDENTITY_REDUCED, FULL_OR_SPARSE_STORAGE);
  // relation = CreateForest(d, MXD, false, forest::QUASI_REDUCED, FULL_OR_SPARSE_STORAGE);
#else
  states = CreateForest(d, MDD, false, forest::FULLY_REDUCED, FULL_STORAGE);
  relation = CreateForest(d, MXD, false, forest::IDENTITY_REDUCED, FULL_STORAGE);
#endif

  if (INVALID_FOREST == states) {
    fprintf(stderr, "Couldn't create forest of states\n");
    return 1;
  } else {
    fprintf(stderr, "Created forest of states\n");
  }
  if (INVALID_FOREST == relation) {
    fprintf(stderr, "Couldn't create forest of relations\n");
    return 1;
  } else {
    fprintf(stderr, "Created forest of relations\n");
  }

  // Build set of initial states
  for (int j = 0; j < type2; j++) {
    initst[variables[get_component_level(2, j)]] = j;
  }
  for (int j = 0; j < type3; j++) {
    initst[variables[get_component_level(3, j)]] = j;
  }

  int i = num_levels;

#if FACED_ORDERING
  i = 0;
#else
#if REORDER_CUBE
  for (int j = type2; j > 0; i--) { assert(initst[variables[i]] == --j); }
  for (int j = type3; j > 0; i--) { assert(initst[variables[i]] == --j); }
#else
  for (int j = type3; j > 0; i--) { assert(initst[variables[i]] == --j); }
  for (int j = type2; j > 0; i--) { assert(initst[variables[i]] == --j); }
#endif
  assert (i == 0);
#endif

  initst[0] = 0;
  initial = CreateVectorElement(states, initst, num_levels + 1, true); 

  dd_edge *curr = NULL;

  if (NULL==initial) {
    fprintf(stderr, "Couldn't create set of initial states\n");
    return 1;
  } else {
    fprintf(stderr, "Created set of initial states\n");
#if 0
    fprintf(stderr, "Initial state: ");
    ShowDDEdge(stderr, initial);
    exit(0);
#endif
  }

  // Build transitions
  const int num_moves = 6;
  dd_edge* move[num_moves][2];

#if !(NO_CHOICE)

  move[F][CW] = Front(CW);
#if 0
  fprintf(stderr, "Overall transition relation: ");
  ShowDDEdge(stderr, move[F][CW]);
  exit(0);
#endif

  move[B][CW] = Back(CW);
  move[L][CW] = Left(CW);
  move[R][CW] = Right(CW);
  move[U][CW] = Up(CW);
  move[D][CW] = Down(CW);
#if 0
  move[F][CCW] = Front(CCW);
  move[B][CCW] = Back(CCW);
  move[L][CCW] = Left(CCW);
  move[R][CCW] = Right(CCW);
  move[U][CCW] = Up(CCW);
  move[D][CCW] = Down(CCW);
#else
  move[F][CCW] = move[F][CW];
  move[B][CCW] = move[B][CW];
  move[L][CCW] = move[L][CW];
  move[R][CCW] = move[R][CW];
  move[U][CCW] = move[U][CW];
  move[D][CCW] = move[D][CW];
#endif

#endif
  
  char dummy;

  // Build overall next-state function
  if(no_choice) {
    nsf = NULL;
    
    if (enable_F) {
      move[F][CW] = Front(CW);
      if (nsf == NULL) {
        nsf = move[F][CW];
      } else {
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, nsf, move[F][CW], nsf));
      }
      fprintf(stderr, "Union f-cw: ");
    }

    if (enable_U) {
      move[U][CW] = Up(CW);
      if (nsf == NULL) {
        nsf = move[U][CW];
      } else {
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, nsf, move[U][CW], nsf));
      }
      fprintf(stderr, "Union u-cw: ");
    }

    if (enable_R) {
      move[R][CW] = Right(CW);
      if (nsf == NULL) {
        nsf = move[R][CW];
      } else {
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, nsf, move[R][CW], nsf));
      }
      fprintf(stderr, "Union r-cw: ");
    }

    if (enable_L) {
      move[L][CW] = Left(CW);
      if (nsf == NULL) {
        nsf = move[L][CW];
      } else {
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, nsf, move[L][CW], nsf));
      }
      fprintf(stderr, "Union l-cw: ");
    }

    if (enable_B) {
      move[B][CW] = Back(CW);
      if (nsf == NULL) {
        nsf = move[B][CW];
      } else {
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, nsf, move[B][CW], nsf));
      }
      fprintf(stderr, "Union b-cw: ");
    }

    if (enable_D) {
      move[D][CW] = Down(CW);
      if (nsf == NULL) {
        nsf = move[D][CW];
      } else {
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, nsf, move[D][CW], nsf));
      }
      fprintf(stderr, "Union d-cw: ");
    }

    if (NULL == nsf) {
      fprintf(stderr, "Couldn't create next-state function\n");
      return 1;
    } else {
      fprintf(stderr, "Created next-state function\n");
    }
  }

#ifdef INTERMEDIATE_PRINT
  fprintf(stderr, "Initial states: ");
  ShowDDEdge(stderr, initial);
  fprintf(stderr, "\nTransition relation nodes\n");
  for (int i=0; i<num_moves; i++) {
    fprintf(stderr, "Move %d, CW: ", i);
    ShowDDEdge(stderr, move[i][0]);
    fprintf(stderr, "\n"); 
    fprintf(stderr, "Move %d, CCW: ", i);
    ShowDDEdge(stderr, move[i][1]);
    fprintf(stderr, "\n"); 
  }
  fprintf(stderr, "Overall transition relation: ");
  ShowDDEdge(stderr, nsf);
  fprintf(stderr, "\n"); 
#endif

  double card = 0;

#ifdef SATURATE

#ifdef COUNT_TIME
  struct rusage start, stop;
  assert(getrusage(RUSAGE_SELF, &start) == 0);
#endif

  if (no_choice) {
#if 0
    Saturate(initial, nsf, curr);
#else
#if 0
    vector<dd_edge *> xd;
    // FURLBD
    xd.push_back(move[F][CW]);
    xd.push_back(move[U][CW]);
    /*
       xd.push_back(move[R][CW]);
       xd.push_back(move[L][CW]);
       xd.push_back(move[B][CW]);
       xd.push_back(move[D][CW]);
     */
    assert(SUCCESS == Saturate(initial, &xd, curr));
#else
    fflush(stderr);
    vector<dd_edge *> *xd = NULL;
    assert(SUCCESS == SplitMxd(nsf, xd));
    assert(SUCCESS == Saturate(initial, xd, curr));
#endif
#endif
  } else {
    int choice = 0;
    curr = initial;
    card = Cardinality(curr);
    fprintf(stderr, "# of reachable states: %1.6e\n", card);
    while (true) {
      fprintf(stderr, "\n\nSaturate using...\n");
      for (i = 0; i < num_moves; i++) {
        fprintf(stderr, "\t%d. %s\n", i, face_to_string((face)i));
      }
      fprintf(stderr, "\t6. Print reachable states\n");
      fprintf(stderr, "\t7. Stop\n");
      fprintf(stderr, "Enter Choice (0-7): ");
      fscanf(stdin, "%d%c", &choice, &dummy);
      if (choice >= 0 && choice < num_moves) {
#if 0
        Saturate(curr, temp_nsf[choice], curr);
#else
#if 1
        Saturate(curr, move[choice][0], curr);
#else
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, move[F][CW], move[U][CW], nsf));
#if 1
        assert(SUCCESS == 
            ApplyBinary(OP_UNION, nsf, move[R][CW], nsf));
#endif
        Saturate(initial, nsf, curr);
        break;
#endif
#endif
        card = Cardinality(curr);
        fprintf(stderr, "# of reachable states: %1.6e\n", card);
      } else if (choice == 6) {
        fprintf(stderr, "Reachable states: ");
        ShowDDEdge(stderr, curr);
        fprintf(stderr, "\n");
      } else if (choice == 7) {
        fprintf(stderr, "Stopping... \n");
        break;
      }
    }
  }

#ifdef COUNT_TIME
  assert(getrusage(RUSAGE_SELF, &stop) == 0);

  // stop.ru_utime - start.ru_utime
  suseconds_t u_sat_time =
    (stop.ru_utime.tv_sec * 1000000 + stop.ru_utime.tv_usec) -
    (start.ru_utime.tv_sec * 1000000 + start.ru_utime.tv_usec);
  suseconds_t s_sat_time =
    (stop.ru_stime.tv_sec * 1000000 + stop.ru_stime.tv_usec) -
    (start.ru_stime.tv_sec * 1000000 + start.ru_stime.tv_usec);

  fprintf(stderr, "\nTime for constructing initial states and nsf:\n");
  fprintf(stderr, "  %ld.%ld sec user, %ld.%ld sec system\n",
      (long int)start.ru_utime.tv_sec, (long int)start.ru_utime.tv_usec,
      (long int)start.ru_stime.tv_sec, (long int)start.ru_stime.tv_usec);

  fprintf(stderr, "\nTime for constructing reachability set:\n");
  fprintf(stderr, "  %06f sec user, %06f system\n",
      u_sat_time/1000000.0, s_sat_time/1000000.0);
#endif  // COUNT_TIME

#endif  // SATURATE

#ifdef COUNT_STATES
  card = Cardinality(curr);
  fprintf(stderr, "\n# of reachable states: %1.12e\n", card);
  fflush(stderr);
  // fprintf(stderr, "\n# of reachable states: %s\n", ll_to_pa(card));
#endif

#ifdef MEM_USAGE
  fprintf(stderr, "\nPeak memory usage: %d\n", stop.ru_maxrss);
  fprintf(stderr, "\nIntegral shared memory usage: %d\n", stop.ru_ixrss);
#endif

  const char *fn[] = {"reachable", "initial", "nsf", "gif"};
  if (make_gifs) {
    CreateDDEdgePic(fn[1], fn[3], initial);
    printf("Wrote initial states to %s.%s\n", fn[1], fn[3]);
    CreateDDEdgePic(fn[2], fn[3], nsf);
    printf("Wrote next-state function to %s.%s\n", fn[2], fn[3]);
    CreateDDEdgePic(fn[0], fn[3], curr);
    printf("Wrote reachable states to %s.%s\n", fn[0], fn[3]);
  } else if (pretty_print) {
    fprintf(stderr, "\nInitial States: ");
    ShowDDEdge(stderr, initial);
    fprintf(stderr, "\n");
    fprintf(stderr, "\nNext-State Function: ");
    ShowDDEdge(stderr, nsf);
    fprintf(stderr, "\n");
    fprintf(stderr, "\nReachable States: ");
    ShowDDEdge(stderr, curr);
    fprintf(stderr, "\n");
  }

#if 0
  ReleaseEdge(move[F][CW]);
  ReleaseEdge(move[U][CW]);
  ReleaseEdge(move[R][CW]);
  ReleaseEdge(move[L][CW]);
  ReleaseEdge(move[B][CW]);
  ReleaseEdge(move[D][CW]);
  ReleaseEdge(move[F][CCW]);
  ReleaseEdge(move[B][CCW]);
  ReleaseEdge(move[L][CCW]);
  ReleaseEdge(move[R][CCW]);
  ReleaseEdge(move[U][CCW]);
  ReleaseEdge(move[D][CCW]);
#endif

  DestroyForest(states);
  if (INVALID_FOREST != states) {
    fprintf(stderr, "Couldn't destroy forest of states\n");
    return 1;
  } else {
    fprintf(stderr, "Destroyed forest of states\n");
  }

  DestroyForest(relation);
  if (INVALID_FOREST != relation) {
    fprintf(stderr, "Couldn't destroy forest of relations\n");
    return 1;
  } else {
    fprintf(stderr, "Destroyed forest of relations\n");
  }

  DestroyDomain(d);
  cleanup();
  fprintf(stderr, "\n\nDONE\n");
  return 0;
}


