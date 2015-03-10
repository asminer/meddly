
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


#include "defines.h"
#include "loggers.h"


// ******************************************************************
// *                                                                *
// *                       json_logger methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::json_logger::json_logger(std::ostream &s)
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

void MEDDLY::json_logger::logForestInfo(const forest* f, const char* name)
{
  const expert_forest* ef = dynamic_cast<const expert_forest*>(f);
  MEDDLY_DCASSERT(ef);
  if (0==ef) return;

  fixLogger();

  int L = ef->getNumVariables();
  int smallest = ef->isForRelations() ? -L : 1;

  out << "{ \"forest_id\":" << ef->FID() << ", ";
  if (name) {
    out << "\"name\":\"" << name << "\", ";
  }
  out << "\"left\":" << smallest << ", ";
  out << "\"right\":" << L;
  if (recordingNodeCounts()) {
    long* raw_active;
    long* active;
    if (ef->isForRelations()) {
      raw_active = new long[2*L+1];
      active = raw_active+L;
    } else {
      raw_active = new long[L+1];
      active = raw_active;
    }
    ef->countNodesByLevel(active);

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

// ******************************************************************
// *                                                                *
// *                      simple_logger methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::simple_logger::simple_logger(std::ostream &s)
  : out(s)
{
  out << "T simple\n";
}

MEDDLY::simple_logger::~simple_logger()
{
}

void MEDDLY::simple_logger::addComment(const char* str)
{
  //
  // TBD - fancify - scan str, and after every newline,
  // start the next line with #
  //
  if (str) {
    out << "# " << str << "\n";
  }
}

void MEDDLY::simple_logger::logForestInfo(const forest* f, const char* name)
{
  const expert_forest* ef = dynamic_cast<const expert_forest*>(f);
  MEDDLY_DCASSERT(ef);
  if (0==ef) return;

  fixLogger();

  int L = ef->getNumVariables();
  int smallest = ef->isForRelations() ? -L : 1;

  out << "F " << ef->FID();
  if (name) {
    out << " \"" << name << "\"";
  }
  out << " " << smallest; // Left value
  out << " " << L;        // Right value

  if (recordingNodeCounts()) {
    long* raw_active;
    long* active;
    if (ef->isForRelations()) {
      raw_active = new long[2*L+1];
      active = raw_active+L;
    } else {
      raw_active = new long[L+1];
      active = raw_active;
    }
    ef->countNodesByLevel(active);

    out << " [";
    for (int l=smallest; l<=L; l++) {
      out << active[l];
      if (l<L) out << ", ";
    }
    out << "]";
  }
  out << "\n";
  out.flush();
}

void MEDDLY::simple_logger
::addToActiveNodeCount(const forest* f, int level, long delta)
{
  if (0==f) return;
  out << "a " << f->FID() << " " << level << " " << delta << "\n";
  out.flush();
}

