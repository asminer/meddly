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

#ifndef MEDDLY_LOG_SIMPLE_H
#define MEDDLY_LOG_SIMPLE_H

#include "logger.h"
#include "io.h"

namespace MEDDLY {

  class simple_logger;

};


// ******************************************************************
// *                                                                *
// *                       simple_logger class                      *
// *                                                                *
// ******************************************************************

/** For logging stats in a flexible and compact text-based format.
    No timestamps.
*/

class MEDDLY::simple_logger : public MEDDLY::logger {
    output &out;
    long** active_delta;
    int* left;
    int* right;
    unsigned batch_forests;
    int aggregate;
    int ucount;
  public:
    simple_logger(output &s, int agg=16);
    virtual ~simple_logger();

    virtual void addComment(const char* comment);
    virtual void newPhase(const forest* f, const char* comment);
    virtual void logForestInfo(const forest* f, const char* name);
    virtual void addToActiveNodeCount(const forest* f, int level, long delta);

  protected:
    inline long* activeArray(unsigned fid) const {
      return (left[fid]<0) ? active_delta[fid] - left[fid] : active_delta[fid];
    }

    void flushLog();
};

#endif // #include guard
