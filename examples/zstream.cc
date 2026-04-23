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

#include "zstream.h"

// **********************************************************************

/**
    Finds the extension for a file, and returns an appropriate character:

        'x' : for .xz files
        'g' : for .gz files
        ' ' : for all others

 */
char getExtension(const char* fname)
{
    if (!fname) return ' ';

    const char* ext = nullptr;
    for (const char* ptr = fname; *ptr; ++ptr) {
        if ('.' == *ptr) ext = ptr+1;
    }

    // No extension? bail
    if (!ext) return ' ';

    // Extension length != 2? Bail
    if (0==ext[0] || 0==ext[1] || ext[2]) return ' ';

    // Extension should end in z
    if ('z' != ext[1]) return ' ';

    switch (ext[0]) {
        case 'g':
        case 'x':
                    return ext[0];

        default:
                    return ' ';
    }
}

// **********************************************************************
// *                                                                    *
// *                      compressed_input methods                      *
// *                                                                    *
// **********************************************************************

compressed_input::compressed_input(const char* fname)
{
    char type = getExtension(fname);

    if (' ' == type) {
        ispipe = false;
        fin = fopen(fname, "r");
    } else {
        ispipe = true;

        char cmd[1024]; // lazy

        snprintf(cmd, 1024, "%s %s",
                ('x' == type) ? "xzcat" : "gzcat",
                fname
        );

        fin = popen(cmd, "r");
    }

    instrm.setFILE(fin);
}

compressed_input::~compressed_input()
{
    if (fin) {
        if (ispipe) {
            pclose(fin);
        } else {
            fclose(fin);
        }
        fin = nullptr;
    }
}

// **********************************************************************
// *                                                                    *
// *                     compressed_output  methods                     *
// *                                                                    *
// **********************************************************************

compressed_output::compressed_output(const char* fname)
{
    char type = getExtension(fname);

    if (' ' == type) {
        ispipe = false;
        fout = fopen(fname, "w");
    } else {
        ispipe = true;

        char cmd[1024]; // lazy

        snprintf(cmd, 1024, "%s > %s",
                ('x' == type) ? "xz" : "gzip",
                fname
        );

        fout = popen(cmd, "w");
    }

    outstrm.setFILE(fout);
}

compressed_output::~compressed_output()
{
    if (fout) {
        if (ispipe) {
            pclose(fout);
        } else {
            fclose(fout);
        }
        fout = nullptr;
    }
}

