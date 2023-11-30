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

#include "timer.h"
#include "reporting.h"

std::ofstream report;

void setReport(int argc, const char** argv)
{
    const char* fname = nullptr;
    // Process command-line args
    for (int i=1; i<argc; i++) {

        if (0==strcmp("-r", argv[i])) {
            ++i;
            fname = argv[i];
            continue;
        }

        //
        // Ignore any other args
    }

    if (fname) {
        report.open(fname);
        if (!report) {
            throw two_strings("Couldn't write to file ", fname);
        }
    }
}

bool startReport(const timer &T, const char* exe)
{
    if (report) {
        report << T.get_last_seconds() << "  ";
        for (unsigned i=0; exe[i]; i++) {
            if ('.' == exe[i]) break;
            report << exe[i];
        }
        report << "  ";
        return true;
    } else {
        return false;
    }
}

void show_sec(std::ostream &s, const timer &T, int a, int b)
{
    long sec, usec;
    T.get_last(sec, usec);
    s << std::setw(a) << sec;
    if (b) {
        s << '.';
        long d = 100000;
        for (int i=0; i<b; i++) {
            s << (usec / d);
            usec %= d;
            d /= 10;
        }
    }
    s << " seconds\n";
}



