
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

#include "defines.h"

/*
  The error class is implemented here.

  For now, this is overkill, because it is a very simple class.
  But someday we may want a class hierarchy...

*/

MEDDLY::error::error(MEDDLY::error::code c)
{ 
  errcode = c; 
  fname = 0;
  lineno = 0;
}

MEDDLY::error::error(MEDDLY::error::code c, const char* fn, int ln)
{ 
  errcode = c; 
  fname = fn;
  lineno = ln;
}

const char* MEDDLY::error::getName() const 
{
  switch (errcode) {
      case  MEDDLY::error::UNINITIALIZED:        return "Uninitialized";
      case  MEDDLY::error::ALREADY_INITIALIZED:  return "Already initialized";
      case  MEDDLY::error::NOT_IMPLEMENTED:      return "Not implemented";
      case  MEDDLY::error::INSUFFICIENT_MEMORY:  return "Insufficient memory";
      case  MEDDLY::error::INVALID_OPERATION:    return "Invalid operation";
      case  MEDDLY::error::INVALID_VARIABLE:     return "Invalid variable";
      case  MEDDLY::error::INVALID_LEVEL:        return "Invalid level";
      case  MEDDLY::error::INVALID_BOUND:        return "Invalid bound";
      case  MEDDLY::error::DOMAIN_NOT_EMPTY:     return "Domain not empty";
      case  MEDDLY::error::UNKNOWN_OPERATION:    return "Unknown operation";
      case  MEDDLY::error::DOMAIN_MISMATCH:      return "Domain mismatch";
      case  MEDDLY::error::FOREST_MISMATCH:      return "Forest mismatch";
      case  MEDDLY::error::TYPE_MISMATCH:        return "Type mismatch";
      case  MEDDLY::error::WRONG_NUMBER:         return "Wrong number";
      case  MEDDLY::error::VALUE_OVERFLOW:       return "Overflow";
      case  MEDDLY::error::DIVIDE_BY_ZERO:       return "Divide by zero";
      case  MEDDLY::error::INVALID_POLICY:       return "Invalid policy";
      case  MEDDLY::error::INVALID_ASSIGNMENT:   return "Invalid assignment";
      case  MEDDLY::error::INVALID_ARGUMENT:     return "Invalid argument";
      case  MEDDLY::error::INVALID_FILE:         return "Invalid file";
      case  MEDDLY::error::COULDNT_WRITE:        return "Couldn't write to file";
      case  MEDDLY::error::COULDNT_READ:         return "Couldn't read from file";
      case  MEDDLY::error::MISCELLANEOUS:        return "Miscellaneous";
      default:                           return "Unknown error";
  }
}
