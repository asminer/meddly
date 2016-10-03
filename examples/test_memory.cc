
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

#include <iostream>
#include "meddly_expert.h"

const int MAX_ID = 32;

int usage(const char* what)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) {
      name = ptr+1;
    }
  }

  cerr << "\nUsage: " << name << " [options] script script ...\n\n";
  cerr << "Run testing scripts.  If none specified, read scripts from standard input.\n";
  cerr << "Possible options:\n";
  cerr << "\t-mm=MANAGER:     set the memory manager to the specified style.\n";
  cerr << "\t                 Possible managers:\n";
  cerr << "\t                     TBD\n";
  cerr << "\nScripts:\n";
  cerr << "Scripts are free-form text files, with the following statements:\n";
  cerr << "    # for comments (ignore until end of line)\n";
  cerr << "    + ident size;  to allocate an identifier of a given size.\n";
  cerr << "    - ident;       to free an identifier.\n";
  cerr << "    ?;             to show active identifiers and memory management information\n";
  cerr << "Allocating an identifier makes it active, and freeing an identifier also\n";
  cerr << "deactivates it.  It is an error to allocate an already active identifier,\n";
  cerr << "or to free an identifier that is not active.  Identifiers follow the usual\n";
  cerr << "naming restrictions: can contain letters, digits, or underscores, and must\n";
  cerr << "not start with a digit.  Only the first " << MAX_ID << " characters of an identifier\n";
  cerr << "are used.\n\n";

  return 0;
}

