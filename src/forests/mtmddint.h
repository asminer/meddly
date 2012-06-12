
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

#ifndef MTMDDINT_H
#define MTMDDINT_H

#include "mtmdd.h"

namespace MEDDLY {
  class mt_mdd_int; 
};

// ******************************************************************

/** 
    Forest for multi-terminal, mdd, integer range.
*/
class MEDDLY::mt_mdd_int : public mtmdd_forest {
  public:

    mt_mdd_int(int dsl, domain *d);
    ~mt_mdd_int();
};

#endif

