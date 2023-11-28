/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#ifndef REPORTING_H
#define REPORTING_H

#include <fstream>

class timer;

struct two_strings {
    const char* first;
    const char* second;

    two_strings(const char* f, const char* s) {
        first = f;
        second = s;
    }
};


///
/// Looks through command-line switches for
///     -r  filename
/// and sets the report stream accordingly.
/// Will not touch any other arguments.
///
///     @param      argc    Arg count (copy from main)
///     @param      argv    Arg values (copy from main)
///
///     @throws     A two_strings object on error.
///
void setReport(int argc, const char** argv);


///
/// Start reporting, if the report stream is set up.
///
///     @param  T       Timer, with stop time noted already.
///     @param  exe     Executable name, possibly with an extra extension,
///                     so you can pass __FILE__ as the argument
///
///     @return true, if the report stream is set up.
///
bool startReport(const timer &T, const char* exe);


extern std::ofstream report;

#endif
