

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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../defines.h"
#include "init_builtin.h"

#include "copy.h"
#include "cardinality.h"
#include "complement.h"
#include "maxmin_range.h"
#include "mdd2index.h"

#include "union.h"
#include "intersection.h"
#include "difference.h"
#include "cross.h"

#include "maxmin.h"
#include "plus.h"
#include "minus.h"
#include "multiply.h"
#include "divide.h"

#include "comp_eq.h"
#include "comp_ne.h"
#include "comp_lt.h"
#include "comp_le.h"
#include "comp_gt.h"
#include "comp_ge.h"

#include "prepostimage.h"
#include "reach_bfs.h"
#include "reach_dfs.h"

#include "vect_matr.h"

#include "mm_mult.h"

template <class T>
inline void initP(const T* &global, T* &local, T* init)
{
  global = (local = init);
}

void MEDDLY::builtin_initializer::init(const settings &s)
{
  initP(MEDDLY::COPY,                 COPY,       initializeCopy(s)         );
  initP(MEDDLY::CARDINALITY,          CARD,       initializeCardinality(s)  );
  initP(MEDDLY::COMPLEMENT,           COMPL,      initializeComplement(s)   );
  initP(MEDDLY::MAX_RANGE,            MAXRANGE,   initializeMaxRange(s)     );
  initP(MEDDLY::MIN_RANGE,            MINRANGE,   initializeMaxRange(s)     );
  initP(MEDDLY::CONVERT_TO_INDEX_SET, MDD2INDEX,  initializeMDD2INDEX(s)    );

  initP(MEDDLY::UNION,                UNION,      initializeUnion(s)        );
  initP(MEDDLY::INTERSECTION,         INTERSECT,  initializeIntersection(s) );
  initP(MEDDLY::DIFFERENCE,           DIFFERENCE, initializeDifference(s)   );
  initP(MEDDLY::CROSS,                CROSS,      initializeCross(s)        );

  initP(MEDDLY::MAXIMUM,              MAX,        initializeMaximum(s)      );
  initP(MEDDLY::MINIMUM,              MIN,        initializeMinimum(s)      );
  initP(MEDDLY::PLUS,                 PLUS,       initializePlus(s)         );
  initP(MEDDLY::MINUS,                MINUS,      initializeMinus(s)        );
  initP(MEDDLY::MULTIPLY,             MULTIPLY,   initializeMultiply(s)     );
  initP(MEDDLY::DIVIDE,               DIVIDE,     initializeMultiply(s)     );

  initP(MEDDLY::EQUAL,                EQ,           initializeEQ(s)         );
  initP(MEDDLY::NOT_EQUAL,            NE,           initializeNE(s)         );
  initP(MEDDLY::LESS_THAN,            LT,           initializeLT(s)         );
  initP(MEDDLY::LESS_THAN_EQUAL,      LE,           initializeLE(s)         );
  initP(MEDDLY::GREATER_THAN,         GT,           initializeGT(s)         );
  initP(MEDDLY::GREATER_THAN_EQUAL,   GE,           initializeGE(s)         );

  initP(MEDDLY::PRE_IMAGE,            PRE_IMAGE,    initializePreImage(s)   );
  initP(MEDDLY::POST_IMAGE,           POST_IMAGE,   initializePostImage(s)  );
  initP(MEDDLY::REACHABLE_STATES_DFS, FORWARD_DFS,  initializeForwardDFS(s) );
  initP(MEDDLY::REACHABLE_STATES_BFS, FORWARD_BFS,  initializeForwardBFS(s) );
  initP(MEDDLY::REVERSE_REACHABLE_DFS,BACKWARD_DFS, initializeBackwardDFS(s));
  initP(MEDDLY::REVERSE_REACHABLE_BFS,BACKWARD_BFS, initializeBackwardBFS(s));

  initP(MEDDLY::VM_MULTIPLY,          VM_MULTIPLY,  initializeVMmult(s)     );
  initP(MEDDLY::MV_MULTIPLY,          MV_MULTIPLY,  initializeMVmult(s)     );

  initP(MEDDLY::EXPLVECT_MATR_MULT, EXPLVECT_MATR_MULT, initExplVectorMatrixMult(s));
  initP(MEDDLY::MATR_EXPLVECT_MULT, MATR_EXPLVECT_MULT, initMatrixExplVectorMult(s));

  initP(MEDDLY::MM_MULTIPLY,          MM_MULTIPLY,  initializeMMMultiply(s) );
}

template <class T>
inline void cleanPair(T *local, const T* &global)
{
  delete local;
  global = 0;
}

void MEDDLY::builtin_initializer::cleanup()
{
  cleanPair(COPY,           MEDDLY::COPY);
  cleanPair(CARD,           MEDDLY::CARDINALITY);
  cleanPair(COMPL,          MEDDLY::COMPLEMENT);
  cleanPair(MAXRANGE,       MEDDLY::MAX_RANGE);
  cleanPair(MINRANGE,       MEDDLY::MIN_RANGE);
  cleanPair(MDD2INDEX,      MEDDLY::CONVERT_TO_INDEX_SET);

  cleanPair(UNION,          MEDDLY::UNION);
  cleanPair(INTERSECT,      MEDDLY::INTERSECTION);
  cleanPair(DIFFERENCE,     MEDDLY::DIFFERENCE);
  cleanPair(CROSS,          MEDDLY::CROSS);

  cleanPair(MAX,            MEDDLY::MAXIMUM);
  cleanPair(MIN,            MEDDLY::MINIMUM);
  cleanPair(PLUS,           MEDDLY::PLUS);
  cleanPair(MINUS,          MEDDLY::MINUS);
  cleanPair(MULTIPLY,       MEDDLY::MULTIPLY);
  cleanPair(DIVIDE,         MEDDLY::DIVIDE);

  cleanPair(EQ,             MEDDLY::EQUAL);
  cleanPair(NE,             MEDDLY::NOT_EQUAL);
  cleanPair(LT,             MEDDLY::LESS_THAN);
  cleanPair(LE,             MEDDLY::LESS_THAN_EQUAL);
  cleanPair(GT,             MEDDLY::GREATER_THAN);
  cleanPair(GE,             MEDDLY::GREATER_THAN_EQUAL);

  cleanPair(PRE_IMAGE,      MEDDLY::PRE_IMAGE);
  cleanPair(POST_IMAGE,     MEDDLY::POST_IMAGE);
  cleanPair(FORWARD_DFS,    MEDDLY::REACHABLE_STATES_DFS);
  cleanPair(FORWARD_BFS,    MEDDLY::REACHABLE_STATES_BFS);
  cleanPair(BACKWARD_DFS,   MEDDLY::REVERSE_REACHABLE_DFS);
  cleanPair(BACKWARD_BFS,   MEDDLY::REVERSE_REACHABLE_BFS);

  cleanPair(EXPLVECT_MATR_MULT, MEDDLY::EXPLVECT_MATR_MULT);
  cleanPair(MATR_EXPLVECT_MULT, MEDDLY::MATR_EXPLVECT_MULT);

  cleanPair(MM_MULTIPLY,    MEDDLY::MM_MULTIPLY);
}

