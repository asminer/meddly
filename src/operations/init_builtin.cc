

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

#include "cardinality.h"
#include "complement.h"
#include "maxmin_range.h"
#include "mdd2evplus.h"

void MEDDLY::builtin_initializer::init(const settings &s)
{
  MEDDLY::CARDINALITY           = (CARD       = initializeCardinality(s));
  MEDDLY::COMPLEMENT            = (COMPL      = initializeComplement(s) );
  MEDDLY::MAX_RANGE             = (MAXRANGE   = initializeMaxRange(s)   );
  MEDDLY::MIN_RANGE             = (MINRANGE   = initializeMaxRange(s)   );
  MEDDLY::CONVERT_TO_INDEX_SET  = (MDD2EVPLUS = initializeMDD2EVPLUS(s) );
}

template <class T>
inline void cleanPair(T *local, const T* &global)
{
  delete local;
  global = 0;
}

void MEDDLY::builtin_initializer::cleanup()
{
  cleanPair(CARD,         MEDDLY::CARDINALITY);
  cleanPair(COMPL,        MEDDLY::COMPLEMENT);
  cleanPair(MAXRANGE,     MEDDLY::MAX_RANGE);
  cleanPair(MINRANGE,     MEDDLY::MIN_RANGE);
  cleanPair(MDD2EVPLUS,   MEDDLY::CONVERT_TO_INDEX_SET);
}
