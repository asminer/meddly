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
#include "error.h"

MEDDLY::error::error(MEDDLY::error::code c, const char* fn, unsigned ln)
{
    errcode = c;
    fname = fn;
    lineno = ln;
}

const char* MEDDLY::error::getName() const
{
    switch (errcode) {
        case  error::UNINITIALIZED:         return "Uninitialized";
        case  error::ALREADY_INITIALIZED:   return "Already initialized";
        case  error::NOT_IMPLEMENTED:       return "Not implemented";
        case  error::INSUFFICIENT_MEMORY:   return "Insufficient memory";
        case  error::INVALID_OPERATION:     return "Invalid operation";
        case  error::INVALID_VARIABLE:      return "Invalid variable";
        case  error::INVALID_LEVEL:         return "Invalid level";
        case  error::INVALID_BOUND:         return "Invalid bound";
        case  error::INVALID_ITERATOR:      return "Invalid iterator";
        case  error::DOMAIN_NOT_EMPTY:      return "Domain not empty";
        case  error::UNKNOWN_OPERATION:     return "Unknown operation";
        case  error::DOMAIN_MISMATCH:       return "Domain mismatch";
        case  error::FOREST_MISMATCH:       return "Forest mismatch";
        case  error::TYPE_MISMATCH:         return "Type mismatch";
        case  error::WRONG_NUMBER:          return "Wrong number";
        case  error::VALUE_OVERFLOW:        return "Overflow";
        case  error::DIVIDE_BY_ZERO:        return "Divide by zero";
        case  error::SUBTRACT_INFINITY:     return "Subtract infinity";
        case  error::INVALID_POLICY:        return "Invalid policy";
        case  error::INVALID_ASSIGNMENT:    return "Invalid assignment";
        case  error::INVALID_ARGUMENT:      return "Invalid argument";
        case  error::INVALID_FILE:          return "Invalid file";
        case  error::COULDNT_WRITE:         return "Couldn't write to file";
        case  error::COULDNT_READ:          return "Couldn't read from file";
        case  error::MISCELLANEOUS:         return "Miscellaneous";
        default:                            return "Unknown error";
    }
}
