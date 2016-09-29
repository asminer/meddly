
// $Id: init_builtin.cc 700 2016-07-07 21:06:50Z asminer $

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
#include "init_managers.h"

// include specific manager headers here

/*
namespace MEDDLY {
  // Global memory manager constants here
};
*/




MEDDLY::memman_initializer::memman_initializer(initializer_list *p)
 : initializer_list(p)
{
  // TBD - zero out any memory_manager_factory members here
}

void MEDDLY::memman_initializer::setup()
{
  memory_manager::resetGlobalStats();
  
  // TBD - initialize Global and member factories here
}

void MEDDLY::memman_initializer::cleanup()
{
  // TBD - destroy and zero member factories
  // TBD - zero globals
}

