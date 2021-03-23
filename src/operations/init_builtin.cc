
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
#include "cycle.h"
#include "select.h"
#include "incoming_edge_count.h"
#include "above_count.h"
#include "below_count.h"
#include "highest_unique.h"
#include "unique_count.h"

#include "union.h"
#include "intersection.h"
#include "difference.h"
#include "cross.h"

#include "maxmin.h"
#include "plus.h"
#include "minus.h"
#include "multiply.h"
#include "divide.h"
#include "modulo.h"

#include "comp_eq.h"
#include "comp_ne.h"
#include "comp_lt.h"
#include "comp_le.h"
#include "comp_gt.h"
#include "comp_ge.h"

#include "prepostplus.h"

#include "prepostimage.h"
#include "reach_bfs.h"
#include "reach_dfs.h"
#include "sat_pregen.h"
#include "sat_otf.h"
#include "sat_impl.h"
#include "sat_hyb.h"

#include "vect_matr.h"

#include "mm_mult.h"

#include "constrained.h"
#include "transitive_closure.h"

#include "mpz_object.h"

namespace MEDDLY {

  // unary operation "codes"

  const unary_opname* COPY = 0;
  const unary_opname* CARDINALITY = 0;
  const unary_opname* COMPLEMENT = 0;
  const unary_opname* MAX_RANGE = 0;
  const unary_opname* MIN_RANGE = 0;
  const unary_opname* CONVERT_TO_INDEX_SET = 0;
  const unary_opname* CYCLE = 0;
  const unary_opname* SELECT = 0;
  const unary_opname* IEC = 0;
  const unary_opname* AC=0;
  const unary_opname* BC=0;
  const unary_opname* HU=0;
  const unary_opname* UC=0;

  // binary operation "codes"

  const binary_opname* UNION = 0;
  const binary_opname* INTERSECTION = 0;
  const binary_opname* DIFFERENCE = 0;
  const binary_opname* CROSS = 0;

  const binary_opname* MINIMUM = 0;
  const binary_opname* MAXIMUM = 0;
  const binary_opname* PLUS = 0;
  const binary_opname* MINUS = 0;
  const binary_opname* MULTIPLY = 0;
  const binary_opname* DIVIDE = 0;
  const binary_opname* MODULO = 0;

  const binary_opname* EQUAL = 0;
  const binary_opname* NOT_EQUAL = 0;
  const binary_opname* LESS_THAN = 0;
  const binary_opname* LESS_THAN_EQUAL = 0;
  const binary_opname* GREATER_THAN = 0;
  const binary_opname* GREATER_THAN_EQUAL = 0;

  const binary_opname* PRE_PLUS = 0;
  const binary_opname* POST_PLUS = 0;

  const binary_opname* PRE_IMAGE = 0;
  const binary_opname* POST_IMAGE = 0;
  const binary_opname* TC_POST_IMAGE = 0;
  const binary_opname* REACHABLE_STATES_DFS = 0;
  const binary_opname* REACHABLE_STATES_BFS = 0;
  const binary_opname* REVERSE_REACHABLE_DFS = 0;
  const binary_opname* REVERSE_REACHABLE_BFS = 0;
  const binary_opname* REACHABLE_STATES_BFS_UA = 0;
  const binary_opname* REACHABLE_STATES_BFS_HUA = 0;

  const binary_opname* VM_MULTIPLY = 0;
  const binary_opname* MV_MULTIPLY = 0;

  const binary_opname* MM_MULTIPLY = 0;

  // numerical operation "codes"

  const numerical_opname* EXPLVECT_MATR_MULT = 0;
  const numerical_opname* MATR_EXPLVECT_MULT = 0;

  // saturation operation "codes"

  const satpregen_opname* SATURATION_FORWARD = 0;
  const satpregen_opname* SATURATION_BACKWARD = 0;
  const satotf_opname* SATURATION_OTF_FORWARD = 0;
  const satimpl_opname* SATURATION_IMPL_FORWARD = 0;
  const sathyb_opname* SATURATION_HYB_FORWARD = 0;

  // minimum witness operation "codes"
  const constrained_opname* CONSTRAINED_BACKWARD_BFS = 0;
  const constrained_opname* CONSTRAINED_FORWARD_DFS = 0;
  const constrained_opname* CONSTRAINED_BACKWARD_DFS = 0;
  const constrained_opname* TRANSITIVE_CLOSURE_DFS = 0;
};



MEDDLY::builtin_initializer::builtin_initializer(initializer_list *p)
 : initializer_list(p)
{
}



template <class T>
inline void initP(const T* &global, T* &local, T* init)
{
  global = (local = init);
}


void MEDDLY::builtin_initializer::setup()
{
  initP(MEDDLY::COPY,                 COPY,       initializeCopy()          );
  initP(MEDDLY::CARDINALITY,          CARD,       initializeCardinality()   );
  initP(MEDDLY::COMPLEMENT,           COMPL,      initializeComplement()    );
  initP(MEDDLY::MAX_RANGE,            MAXRANGE,   initializeMaxRange()      );
  initP(MEDDLY::MIN_RANGE,            MINRANGE,   initializeMaxRange()      );
  initP(MEDDLY::CONVERT_TO_INDEX_SET, MDD2INDEX,  initializeMDD2INDEX()     );
  initP(MEDDLY::CYCLE,                CYCLE,      initializeCycle()         );
  initP(MEDDLY::SELECT,               SELECT,     initializeSelect()        );
  initP(MEDDLY::IEC,                  IEC,        initializeIncomingEdgeCount());
  initP(MEDDLY::AC,                   AC,         initializeAboveCount());
  initP(MEDDLY::BC,                   BC,         initializeBelowCount());
  initP(MEDDLY::HU,                   HU,         initializeHighestUnique());
  initP(MEDDLY::UC,                   UC,         initializeUniqueCount());

  initP(MEDDLY::UNION,                UNION,      initializeUnion()         );
  initP(MEDDLY::INTERSECTION,         INTERSECT,  initializeIntersection()  );
  initP(MEDDLY::DIFFERENCE,           DIFFERENCE, initializeDifference()    );
  initP(MEDDLY::CROSS,                CROSS,      initializeCross()         );

  initP(MEDDLY::MAXIMUM,              MAX,        initializeMaximum()       );
  initP(MEDDLY::MINIMUM,              MIN,        initializeMinimum()       );
  initP(MEDDLY::PLUS,                 PLUS,       initializePlus()          );
  initP(MEDDLY::MINUS,                MINUS,      initializeMinus()         );
  initP(MEDDLY::MULTIPLY,             MULTIPLY,   initializeMultiply()      );
  initP(MEDDLY::DIVIDE,               DIVIDE,     initializeDivide()        );
  initP(MEDDLY::MODULO,               MODULO,     initializeModulo()        );

  initP(MEDDLY::EQUAL,                EQ,           initializeEQ()          );
  initP(MEDDLY::NOT_EQUAL,            NE,           initializeNE()          );
  initP(MEDDLY::LESS_THAN,            LT,           initializeLT()          );
  initP(MEDDLY::LESS_THAN_EQUAL,      LE,           initializeLE()          );
  initP(MEDDLY::GREATER_THAN,         GT,           initializeGT()          );
  initP(MEDDLY::GREATER_THAN_EQUAL,   GE,           initializeGE()          );

  initP(MEDDLY::PRE_PLUS,             PRE_PLUS,     initializePrePlus()     );
  initP(MEDDLY::POST_PLUS,            POST_PLUS,    initializePostPlus()    );

  initP(MEDDLY::PRE_IMAGE,            PRE_IMAGE,    initializePreImage()    );
  initP(MEDDLY::POST_IMAGE,           POST_IMAGE,   initializePostImage()   );
  initP(MEDDLY::TC_POST_IMAGE,        TC_POST_IMAGE,initializeTCPostImage() );
  initP(MEDDLY::REACHABLE_STATES_DFS, FORWARD_DFS,  initializeForwardDFS()  );
  initP(MEDDLY::REACHABLE_STATES_BFS, FORWARD_BFS,  initializeForwardBFS()  );
  initP(MEDDLY::REVERSE_REACHABLE_DFS,BACKWARD_DFS, initializeBackwardDFS() );
  initP(MEDDLY::REVERSE_REACHABLE_BFS,BACKWARD_BFS, initializeBackwardBFS() );
  initP(MEDDLY::REACHABLE_STATES_BFS_UA, FORWARD_BFS_UA, initializeForwardBFSUA());
  initP(MEDDLY::REACHABLE_STATES_BFS_HUA, FORWARD_BFS_HUA, initializeForwardBFSHUA());

  initP(MEDDLY::VM_MULTIPLY,          VM_MULTIPLY,  initializeVMmult()      );
  initP(MEDDLY::MV_MULTIPLY,          MV_MULTIPLY,  initializeMVmult()      );

  initP(MEDDLY::MM_MULTIPLY,          MM_MULTIPLY,  initializeMMMultiply()  );

  initP(MEDDLY::EXPLVECT_MATR_MULT, EXPLVECT_MATR_MULT, initExplVectorMatrixMult()  );
  initP(MEDDLY::MATR_EXPLVECT_MULT, MATR_EXPLVECT_MULT, initMatrixExplVectorMult()  );

  initP(MEDDLY::SATURATION_FORWARD,   SATURATION_FORWARD,   initSaturationForward()   );
  initP(MEDDLY::SATURATION_BACKWARD,  SATURATION_BACKWARD,  initSaturationBackward()  );
  initP(MEDDLY::SATURATION_OTF_FORWARD,   SATURATION_OTF_FORWARD,   initOtfSaturationForward()  );
  initP(MEDDLY::SATURATION_IMPL_FORWARD, SATURATION_IMPL_FORWARD, initImplSaturationForward()  );
  initP(MEDDLY::SATURATION_HYB_FORWARD, SATURATION_HYB_FORWARD, initHybSaturationForward()  );
  initP(MEDDLY::CONSTRAINED_BACKWARD_BFS,   CONSTRAINED_BACKWARD_BFS,   initConstrainedBFSBackward()  );
  initP(MEDDLY::CONSTRAINED_FORWARD_DFS,   CONSTRAINED_FORWARD_DFS,   initConstrainedDFSForward()  );
  initP(MEDDLY::CONSTRAINED_BACKWARD_DFS,   CONSTRAINED_BACKWARD_DFS,   initConstrainedDFSBackward()  );
  initP(MEDDLY::TRANSITIVE_CLOSURE_DFS,   TRANSITIVE_CLOSURE_DFS,   initTransitiveClosureDFS()  );

#ifdef HAVE_LIBGMP
  mpz_object::initBuffer();
#endif
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
  cleanPair(IEC,            MEDDLY::IEC);
  cleanPair(AC,             MEDDLY::AC);
  cleanPair(BC,             MEDDLY::BC);
  cleanPair(HU,             MEDDLY::HU);
  cleanPair(UC,             MEDDLY::UC);
  cleanPair(COMPL,          MEDDLY::COMPLEMENT);
  cleanPair(MAXRANGE,       MEDDLY::MAX_RANGE);
  cleanPair(MINRANGE,       MEDDLY::MIN_RANGE);
  cleanPair(MDD2INDEX,      MEDDLY::CONVERT_TO_INDEX_SET);
  cleanPair(CYCLE,          MEDDLY::CYCLE);

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
  cleanPair(MODULO,         MEDDLY::MODULO);

  cleanPair(EQ,             MEDDLY::EQUAL);
  cleanPair(NE,             MEDDLY::NOT_EQUAL);
  cleanPair(LT,             MEDDLY::LESS_THAN);
  cleanPair(LE,             MEDDLY::LESS_THAN_EQUAL);
  cleanPair(GT,             MEDDLY::GREATER_THAN);
  cleanPair(GE,             MEDDLY::GREATER_THAN_EQUAL);

  cleanPair(PRE_PLUS,       MEDDLY::PRE_PLUS);
  cleanPair(POST_PLUS,      MEDDLY::POST_PLUS);

  cleanPair(PRE_IMAGE,      MEDDLY::PRE_IMAGE);
  cleanPair(POST_IMAGE,     MEDDLY::POST_IMAGE);
  cleanPair(TC_POST_IMAGE,  MEDDLY::TC_POST_IMAGE);
  cleanPair(FORWARD_DFS,    MEDDLY::REACHABLE_STATES_DFS);
  cleanPair(FORWARD_BFS,    MEDDLY::REACHABLE_STATES_BFS);
  cleanPair(BACKWARD_DFS,   MEDDLY::REVERSE_REACHABLE_DFS);
  cleanPair(BACKWARD_BFS,   MEDDLY::REVERSE_REACHABLE_BFS);
  cleanPair(FORWARD_BFS_UA, MEDDLY::REACHABLE_STATES_BFS_UA);
  cleanPair(FORWARD_BFS_HUA, MEDDLY::REACHABLE_STATES_BFS_HUA);

  cleanPair(SATURATION_BACKWARD,      MEDDLY::SATURATION_BACKWARD );
  cleanPair(SATURATION_FORWARD,       MEDDLY::SATURATION_FORWARD  );
  cleanPair(SATURATION_OTF_FORWARD,   MEDDLY::SATURATION_OTF_FORWARD  );
  cleanPair(SATURATION_IMPL_FORWARD,   MEDDLY::SATURATION_IMPL_FORWARD  );
  cleanPair(SATURATION_HYB_FORWARD,   MEDDLY::SATURATION_HYB_FORWARD  );

  cleanPair(EXPLVECT_MATR_MULT, MEDDLY::EXPLVECT_MATR_MULT);
  cleanPair(MATR_EXPLVECT_MULT, MEDDLY::MATR_EXPLVECT_MULT);

  cleanPair(MM_MULTIPLY,    MEDDLY::MM_MULTIPLY);


#ifdef HAVE_LIBGMP
  mpz_object::clearBuffer();
#endif
}

