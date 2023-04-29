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

#ifdef _MEDDLY_NOINST_

#include "old_meddly.h"
#include "old_meddly.hh"

#include "error.h"
#include "io.h"
#include "memstats.h"
#include "variable.h"
#include "varorder.h"
#include "domain.h"
#include "forest.h"
#include "relation_node.h"
#include "memory.h"

#else

#include "meddly/old_meddly.h"
#include "meddly/old_meddly.hh"

#include "meddly/error.h"
#include "meddly/io.h"
#include "meddly/memstats.h"
#include "meddly/variable.h"
#include "meddly/varorder.h"
#include "meddly/domain.h"
#include "meddly/forest.h"
#include "meddly/relation_node.h"
#include "meddly/memory.h"

#endif
