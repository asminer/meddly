// $Id$

#ifndef REORDERING_FACTORY_H
#define REORDERING_FACTORY_H

#include <memory>

#include "reordering_type.h"
#include "reordering_base.h"

#include "lowest_inversion_reordering.h"
#include "highest_inversion_reordering.h"
#include "sink_down_reordering.h"
#include "bubble_up_reordering.h"
#include "lowest_cost_reordering.h"
#include "lowest_memory_reordering.h"
#include "random_reordering.h"
#include "LARC_reordering.h"

namespace MEDDLY {

class reordering_factory
{
private:
  reordering_factory();

public:
  static std::unique_ptr<reordering_base> create(reordering_type r);
};

inline std::unique_ptr<reordering_base> reordering_factory::create(reordering_type r)
{
  switch(r) {
  case reordering_type::LOWEST_INVERSION:
    return std::unique_ptr<reordering_base>(new lowest_inversion_reordering());
  case reordering_type::HIGHEST_INVERSION:
    return std::unique_ptr<reordering_base>(new highest_inversion_reordering());
  case reordering_type::SINK_DOWN:
    return std::unique_ptr<reordering_base>(new sink_down_reordering());
  case reordering_type::BUBBLE_UP:
    return std::unique_ptr<reordering_base>(new bubble_up_reordering());
  case reordering_type::LOWEST_COST:
    return std::unique_ptr<reordering_base>(new lowest_cost_reordering());
  case reordering_type::LOWEST_MEMORY:
    return std::unique_ptr<reordering_base>(new lowest_memory_reordering());
  case reordering_type::RANDOM:
    return std::unique_ptr<reordering_base>(new random_reordering());
  case reordering_type::LARC:
    return std::unique_ptr<reordering_base>(new LARC_reordering());
  default:
    throw error(error::INVALID_ARGUMENT);
  }
}

}

#endif
