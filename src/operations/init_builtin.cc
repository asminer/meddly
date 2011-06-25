

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
#include "mdd2evplus.h"

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

void MEDDLY::builtin_initializer::init(const settings &s)
{
  MEDDLY::COPY                  = (COPY       = initializeCopy(s)         );
  MEDDLY::CARDINALITY           = (CARD       = initializeCardinality(s)  );
  MEDDLY::COMPLEMENT            = (COMPL      = initializeComplement(s)   );
  MEDDLY::MAX_RANGE             = (MAXRANGE   = initializeMaxRange(s)     );
  MEDDLY::MIN_RANGE             = (MINRANGE   = initializeMaxRange(s)     );
  MEDDLY::CONVERT_TO_INDEX_SET  = (MDD2EVPLUS = initializeMDD2EVPLUS(s)   );

  MEDDLY::UNION                 = (UNION      = initializeUnion(s)        );
  MEDDLY::INTERSECTION          = (INTERSECT  = initializeIntersection(s) );
  MEDDLY::DIFFERENCE            = (DIFFERENCE = initializeDifference(s)   );
  MEDDLY::CROSS                 = (CROSS      = initializeCross(s)        );

  MEDDLY::MAXIMUM               = (MAX        = initializeMaximum(s)      );
  MEDDLY::MINIMUM               = (MIN        = initializeMinimum(s)      );
  MEDDLY::PLUS                  = (PLUS       = initializePlus(s)         );
  MEDDLY::MINUS                 = (MINUS      = initializeMinus(s)        );
  MEDDLY::MULTIPLY              = (MULTIPLY   = initializeMultiply(s)     );
  MEDDLY::DIVIDE                = (DIVIDE     = initializeMultiply(s)     );

  MEDDLY::EQUAL                 = (EQ         = initializeEQ(s)           );
  MEDDLY::NOT_EQUAL             = (NE         = initializeNE(s)           );
  MEDDLY::LESS_THAN             = (LT         = initializeLT(s)           );
  MEDDLY::LESS_THAN_EQUAL       = (LE         = initializeLE(s)           );
  MEDDLY::GREATER_THAN          = (GT         = initializeGT(s)           );
  MEDDLY::GREATER_THAN_EQUAL    = (GE         = initializeGE(s)           );

  MEDDLY::PRE_IMAGE             = (PRE_IMAGE  = initializePreImage(s)     );
  MEDDLY::POST_IMAGE            = (POST_IMAGE = initializePreImage(s)     );
}

template <class T>
inline void cleanPair(T *local, const T* &global)
{
  delete local;
  global = 0;
}

void MEDDLY::builtin_initializer::cleanup()
{
  cleanPair(COPY,         MEDDLY::COPY);
  cleanPair(CARD,         MEDDLY::CARDINALITY);
  cleanPair(COMPL,        MEDDLY::COMPLEMENT);
  cleanPair(MAXRANGE,     MEDDLY::MAX_RANGE);
  cleanPair(MINRANGE,     MEDDLY::MIN_RANGE);
  cleanPair(MDD2EVPLUS,   MEDDLY::CONVERT_TO_INDEX_SET);

  cleanPair(UNION,        MEDDLY::UNION);
  cleanPair(INTERSECT,    MEDDLY::INTERSECTION);
  cleanPair(DIFFERENCE,   MEDDLY::DIFFERENCE);
  cleanPair(CROSS,        MEDDLY::CROSS);

  cleanPair(MAX,          MEDDLY::MAXIMUM);
  cleanPair(MIN,          MEDDLY::MINIMUM);
  cleanPair(PLUS,         MEDDLY::PLUS);
  cleanPair(MINUS,        MEDDLY::MINUS);
  cleanPair(MULTIPLY,     MEDDLY::MULTIPLY);
  cleanPair(DIVIDE,       MEDDLY::DIVIDE);

  cleanPair(EQ,           MEDDLY::EQUAL);
  cleanPair(NE,           MEDDLY::NOT_EQUAL);
  cleanPair(LT,           MEDDLY::LESS_THAN);
  cleanPair(LE,           MEDDLY::LESS_THAN_EQUAL);
  cleanPair(GT,           MEDDLY::GREATER_THAN);
  cleanPair(GE,           MEDDLY::GREATER_THAN_EQUAL);

  cleanPair(PRE_IMAGE,    MEDDLY::PRE_IMAGE);
  cleanPair(POST_IMAGE,   MEDDLY::POST_IMAGE);
}

