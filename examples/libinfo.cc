
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

/*
    This program simply displays info about MEDDLY itself.
*/

#include <iostream>
#include <stdlib.h>
#include "../src/meddly.h"

int main(int argc, const char** argv)
{
  using namespace std;
  using namespace MEDDLY;
  int i;
  if (1==argc) {
    cout << getLibraryInfo(0) << "\n";
  } else {
    for (i=1; i<argc; i++) {
      long N = atol(argv[i]);
      const char* info = getLibraryInfo(N);
      if (info) {
        cout << info << "\n";
      } else {
        cout << "(null string)\n";
      }
    }
  }
  return 0;
}
