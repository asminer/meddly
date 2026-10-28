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

#ifndef MEDDLY_LOG_JSON_H
#define MEDDLY_LOG_JSON_H

#include "logger.h"
#include "io.h"

namespace MEDDLY {

  class json_logger;

};


// ******************************************************************
// *                                                                *
// *                        json_logger class                       *
// *                                                                *
// ******************************************************************

/** For logging stats in json format.
*/

class MEDDLY::json_logger : public MEDDLY::logger {
    output &out;
  public:
    json_logger(output &s);
    virtual ~json_logger();

    virtual void addComment(const char* comment);
    virtual void newPhase(const forest* f, const char* comment);
    virtual void logForestInfo(const forest* f, const char* name);
    virtual void addToActiveNodeCount(const forest* f, int level, long delta);
};

#endif // #include guard
