
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

#ifndef MEDDLY_INIT_BUILTIN
#define MEDDLY_INIT_BUILTIN

#include "../opname.h"
#include "../opname.h"
#include "../opname_numer.h"
#include "../opname_satur.h"
#include "../initializer.h"

namespace MEDDLY {
  class builtin_initializer;
};

class MEDDLY::builtin_initializer : public initializer_list {

    /*
  numerical_opname* EXPLVECT_MATR_MULT;
  numerical_opname* MATR_EXPLVECT_MULT;
  */

  satpregen_opname* SATURATION_FORWARD;
  satpregen_opname* SATURATION_BACKWARD;
  satotf_opname* SATURATION_OTF_FORWARD;
  satimpl_opname* SATURATION_IMPL_FORWARD;
  sathyb_opname* SATURATION_HYB_FORWARD;

  constrained_opname* CONSTRAINED_BACKWARD_BFS;
  constrained_opname* CONSTRAINED_FORWARD_DFS;
  constrained_opname* CONSTRAINED_BACKWARD_DFS;
  constrained_opname* TRANSITIVE_CLOSURE_DFS;

public:
  builtin_initializer(initializer_list *p);
protected:
  virtual void setup();
  virtual void cleanup();
};

#endif // #include guard
