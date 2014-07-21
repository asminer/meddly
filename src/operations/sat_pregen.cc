
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

#include "../defines.h"
#include "sat_pregen.h"
#include <typeinfo> // for "bad_cast" exception

namespace MEDDLY {
  class fb_saturation_opname;
};

// ******************************************************************
// *                                                                *
// *                   fb_saturation_opname class                   *
// *                                                                *
// ******************************************************************

class MEDDLY::fb_saturation_opname : public satpregen_opname {
    bool forward;
  public:
    fb_saturation_opname(bool fwd);
    virtual specialized_operation* buildOperation(arguments* a) const;
};

MEDDLY::fb_saturation_opname::fb_saturation_opname(bool fwd)
 : satpregen_opname(fwd ? "SaturationFwd" : "SaturationBack")
{
  forward = fwd;
}

MEDDLY::specialized_operation*
MEDDLY::fb_saturation_opname::buildOperation(arguments* a) const
{
  pregen_relation* rel = dynamic_cast<pregen_relation*>(a);
  if (0==rel) throw error(error::INVALID_ARGUMENT);

  //
  // No sanity checks needed here; we did them already when constructing a.
  //

  if (forward) {

    // 
    // TBD. 
    // Build appropriate forward operation here, transferring
    // rel to the operation.
    // Note that rel specifies the forests for the input and output sets
    // (as well as the relations).
    //

    throw error(error::NOT_IMPLEMENTED);
  } else {

    // 
    // TBD. 
    // Build appropriate backward operation here, transferring
    // rel to the operation.
    // Note that rel specifies the forests for the input and output sets
    // (as well as the relations).
    //

    throw error(error::NOT_IMPLEMENTED);
  }
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::satpregen_opname* MEDDLY::initSaturationForward(const settings &s)
{
  return new fb_saturation_opname(true);
}

MEDDLY::satpregen_opname* MEDDLY::initSaturationBackward(const settings &s)
{
  return new fb_saturation_opname(false);
}

