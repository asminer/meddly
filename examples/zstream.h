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

#ifndef ZSTREAM_H
#define ZSTREAM_H

#include "../src/io.h"


/**
    Decompress and read an input file.
    The type of compression is determined by the file name extension.
    Supported compression types:
        .gz :   uses gzcat
        .xz :   uses xzcat
 */
class compressed_input {
    public:
        compressed_input(const char* fname);
        ~compressed_input();

        inline operator bool() const {
            return fin;
        }

        inline MEDDLY::input& instream() {
            return instrm;
        }
    private:
        MEDDLY::FILE_input instrm;
        FILE* fin;
        bool  ispipe;
};


/**
    Compress and write an output file.
    The type of compression is determined by the file name extension.
    Supported compression types:
        .gz :   uses gzcat
        .xz :   uses xzcat
 */
class compressed_output {
    public:
        compressed_output(const char* fname);
        ~compressed_output();

        inline operator bool() const {
            return fout;
        }

        inline MEDDLY::output& outstream() {
            return outstrm;
        }

    private:
        MEDDLY::FILE_output outstrm;
        FILE* fout;
        bool  ispipe;
};


#endif
