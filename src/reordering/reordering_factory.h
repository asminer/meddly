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

#ifndef MEDDLY_REORDERING_FACTORY_H
#define MEDDLY_REORDERING_FACTORY_H

#include <memory>

#include "reordering_base.h"

#include "lowest_inversion_reordering.h"
#include "highest_inversion_reordering.h"
#include "sink_down_reordering.h"
#include "bring_up_reordering.h"
#include "lowest_cost_reordering.h"
#include "lowest_memory_reordering.h"
#include "random_reordering.h"
#include "larc_reordering.h"

namespace MEDDLY {

class reordering_factory
{
private:
  reordering_factory();

public:
  static std::unique_ptr<reordering_base> create(policies::reordering_type r);
};

inline std::unique_ptr<reordering_base> reordering_factory::create(policies::reordering_type r)
{
  switch(r) {
  case policies::reordering_type::LOWEST_INVERSION:
    return std::unique_ptr<reordering_base>(new lowest_inversion_reordering());
  case policies::reordering_type::HIGHEST_INVERSION:
    return std::unique_ptr<reordering_base>(new highest_inversion_reordering());
  case policies::reordering_type::SINK_DOWN:
    return std::unique_ptr<reordering_base>(new sink_down_reordering());
  case policies::reordering_type::BRING_UP:
    return std::unique_ptr<reordering_base>(new bring_up_reordering());
  case policies::reordering_type::LOWEST_COST:
    return std::unique_ptr<reordering_base>(new lowest_cost_reordering());
  case policies::reordering_type::LOWEST_MEMORY:
    return std::unique_ptr<reordering_base>(new lowest_memory_reordering());
  case policies::reordering_type::RANDOM:
    return std::unique_ptr<reordering_base>(new random_reordering());
  case policies::reordering_type::LARC:
    return std::unique_ptr<reordering_base>(new larc_reordering());
  default:
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
  }
}

}

#endif
