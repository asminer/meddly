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
#include "log_json.h"
#include "forest.h"

#define BATCH

// ******************************************************************
// *                                                                *
// *                       json_logger methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::json_logger::json_logger(output &s)
  : out(s)
{
}

MEDDLY::json_logger::~json_logger()
{
}

void MEDDLY::json_logger::addComment(const char*)
{
  // Completely ignored
}

void MEDDLY::json_logger::newPhase(const forest*, const char*)
{
  // Completely ignored
}

void MEDDLY::json_logger::logForestInfo(const forest* f, const char* name)
{
  if (!f) return;

  fixLogger();

  int L = f->getNumVariables();
  int smallest = f->isForRelations() ? -L : 1;

  out << "{ \"forest_id\":" << f->FID() << ", ";
  if (name) {
    out << "\"name\":\"" << name << "\", ";
  }
  out << "\"left\":" << smallest << ", ";
  out << "\"right\":" << L;
  if (recordingNodeCounts()) {
    long* raw_active;
    long* active;
    if (f->isForRelations()) {
      raw_active = new long[2*L+1];
      active = raw_active+L;
    } else {
      raw_active = new long[L+1];
      active = raw_active;
    }
    f->countNodesByLevel(active);

    out << ", \"an\":[";
    for (int l=smallest; l<=L; l++) {
      out << active[l];
      if (l<L) out << ", ";
    }
    out << "]";
  }
  out << " }\n";
  out.flush();
}

void MEDDLY::json_logger::addToActiveNodeCount(const forest* f, int level, long delta)
{
  if (0==f) return;
  out << "{ \"f\":" << f->FID() << ", \"l\":" << level << ", \"anc\":" << delta << " }\n";
  out.flush();
}

