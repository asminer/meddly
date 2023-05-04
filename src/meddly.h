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

#define PATH(X) X

#else

#define PATH(X) "meddly/" ## X

#endif


#include PATH("old_meddly.h")

#include PATH("error.h")
#include PATH("io.h")
#include PATH("memstats.h")
#include PATH("variable.h")
#include PATH("varorder.h")
#include PATH("domain.h")
#include PATH("forest.h")
#include PATH("relation_node.h")
#include PATH("memory.h")

#include PATH("ops_builtin.h")
#include PATH("opname.h")
#include PATH("opname_numer.h")
#include PATH("opname_satur.h")

#include PATH("oper.h")
#include PATH("oper_unary.h")
#include PATH("oper_binary.h")
#include PATH("oper_special.h")

#include PATH("ct_initializer.h")

#include PATH("global_rebuilder.h")

