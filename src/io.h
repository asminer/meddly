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

#ifndef MEDDLY_IO_H
#define MEDDLY_IO_H

#include "defines.h"

#ifndef _MEDDLY_WITHOUT_CSTDIO_
#include <cstdio>
#endif

#ifndef _MEDDLY_WITHOUT_IOSTREAM_
#include <iostream>
#endif

namespace MEDDLY {
    class input;
    class FILE_input;
    class istream_input;
    class output;
    class FILE_output;
    class ostream_output;
};

//
// We are using our own I/O classes because we want to support both
// C++ streams and C FILE*s, and this seems to be the cleanest way
// to do that.  Also it allows users to (easily?) define their own
// I/O mechanisms.
//
//

// ******************************************************************
// *                                                                *
// *                          input  class                          *
// *                                                                *
// ******************************************************************

/**
    Input class.
    Abstract base class for file input.
*/
class MEDDLY::input {
    public:
        input();
        virtual ~input();

        /**
            Return true if and only if the input stream has hit EOF.
        */
        virtual bool eof() const = 0;

        /**
            Read exactly one character from the input stream.
                @return   The character consumed, or EOF.
                @throws   an appropriate error
        */
        virtual int get_char() = 0;

        /**
            Put the last consumed character back on the stream.
            Does nothing if no character was consumed.
                @param  x   Character to put back.
                            x will be the last consumed character.
                @throws     an appropriate error
        */
        virtual void unget(char x) = 0;

        /**
            Read an integer in signed decimal notation.
            TBD - integer is terminated by whitespace, or any character?
                @return   The integer consumed
                @throws   an appropriate error
        */
        virtual long get_integer() = 0;

        /**
            Read a floating point value in signed decimal notation.
            TBD - value is terminated by whitespace, or any character?
                @return   The value consumed
                @throws   an appropriate error
        */
        virtual double get_real() = 0;

        /**
            Read raw bytes into a memory location.
                @param  bytes   Number of bytes requested
                @param  buffer  Pointer to store the bytes

                @return Number of bytes actually read
        */
        virtual size_t read(size_t bytes, unsigned char* buffer) = 0;


        /*
            Handy input stream manipulation
        */

        /**
            Consume whitespace (if any) from the input stream,
            including comments of the form #...\n
        */
        void stripWS();

        /**
            Consume a keyword from the input stream.
            If the keyword does not match the input stream,
            we throw an INVALID_FILE error.
        */
        void consumeKeyword(const char* keyword);

};  // end of input class

// ******************************************************************
// *                                                                *
// *                        FILE_input class                        *
// *                                                                *
// ******************************************************************

#ifndef _MEDDLY_WITHOUT_CSTDIO_

/** FILE_Input class.
    Use for C-style FILE* input.
*/
class MEDDLY::FILE_input : public MEDDLY::input {
    public:
        FILE_input(FILE* _inf = nullptr);
        virtual ~FILE_input();

        inline void setFILE(FILE* _inf) {
            inf = _inf;
        }
        inline operator bool() const {
            return inf;
        }

        virtual bool eof() const;
        virtual int get_char();
        virtual void unget(char);
        virtual long get_integer();
        virtual double get_real();
        virtual size_t read(size_t bytes, unsigned char* buffer);

    private:
        FILE* inf;

};  // end of FILE_input class

#endif // ifndef _MEDDLY_WITHOUT_CSTDIO_

// ******************************************************************
// *                                                                *
// *                      istream_input  class                      *
// *                                                                *
// ******************************************************************

#ifndef _MEDDLY_WITHOUT_IOSTREAM_

/** istream_Input class.
    Use for C++-style istream input.
*/
class MEDDLY::istream_input : public MEDDLY::input {
    public:
        istream_input(std::istream &in);
        virtual ~istream_input();

        virtual bool eof() const;
        virtual int get_char();
        virtual void unget(char);
        virtual long get_integer();
        virtual double get_real();
        virtual size_t read(size_t bytes, unsigned char* buffer);

    private:
        std::istream &in;

};  // end of istream_input class

#endif // ifndef _MEDDLY_WITHOUT_IOSTREAM_

// ******************************************************************
// *                                                                *
// *                          output class                          *
// *                                                                *
// ******************************************************************

/** Output class.
    Abstract base class.
*/
class MEDDLY::output {
        unsigned tabstops;
    public:
        output();
        virtual ~output();

        /// Increase the indentation level by one.
        inline void indent_more() {
            ++tabstops;
        }

        /// Decrease the indentation level by one.
        inline void indent_less() {
            if (tabstops) --tabstops;
        }

        /// Get the indentation level.
        inline unsigned indentation() const {
            return tabstops;
        }

        /// Set the indentation level.
        inline void indentation(unsigned ts) {
            tabstops = ts;
        }

        /**
            Write exactly one character to the output stream.
                @param  x   Character to write
                @throws     An appropriate error
        */
        virtual void put(char x) = 0;

        /**
            Write a string to the output stream.
                @param  x   String to write
                @param  w   Width for formatting
                @throws     An appropriate error
        */
        virtual void put(std::string s, int w=0) = 0;

        /**
            Write a signed, decimal integer to the output stream.
                @param  x   Integer to write
                @param  w   Width for formatting
                @throws     An appropriate error
        */
        virtual void put(long x, int w=0) = 0;

        inline void put(int x, int w=0) {
            put(long(x), w);
        }

        /**
            Write an unsigned, decimal integer to the output stream.
                @param  x   Integer to write
                @param  w   Width for formatting
                @throws     An appropriate error
        */
        virtual void put(unsigned long x, int w=0) = 0;

        inline void put(unsigned x, int w=0) {
            put((unsigned long) x, w);
        }

        /**
            Write hex digits to the output stream.
                @param  x   Value to write
                @param  w   Width for formatting
                @throws     An appropriate error
        */
        virtual void put_hex(unsigned long x, int w=0) = 0;

        /**
            Write a signed, decimal, floating-point value
            to the output stream.
                @param  x   Value to write
                @param  w   Width
                @param  p   Precision
                @param  f   Format, either 'e', 'f', or 'g'
                            (for the style of printf)
                @throws     An appropriate error
        */
        virtual void put(double x, int w=0, int p=6, char f='g') = 0;

        /**
            Write raw bytes from a memory location.
                @param  bytes   Number of bytes in the buffer
                @param  buffer  Pointer to memory location

                @return Number of bytes actually written
        */
        virtual size_t write(size_t bytes, const unsigned char* buffer) = 0;

        /**
            Flush the output stream.
        */
        virtual void flush() = 0;

        /*
            Handy output stream manipulation
        */

        /**
            Write number of bytes, with units.
                @param  m       Number of bytes
                @param  human   If false, units will be "bytes".
                                If true, scale units so that output
                                is betwen one and 1000.
        */
        void put_mem(size_t m, bool human);

};  // end of output class


// ******************************************************************
// *                                                                *
// *                       FILE_output  class                       *
// *                                                                *
// ******************************************************************

#ifndef _MEDDLY_WITHOUT_CSTDIO_

/** FILE_output class.
    Use for C-style FILE* input.
*/
class MEDDLY::FILE_output : public MEDDLY::output {
    public:
        FILE_output(FILE* outf = nullptr);
        virtual ~FILE_output();

        inline void setFILE(FILE* _outf) {
            outf = _outf;
        }

        inline operator bool() const {
            return outf;
        }

        virtual void put(char x);
        virtual void put(std::string, int w);
        virtual void put(long x, int w);
        virtual void put(unsigned long x, int w);
        virtual void put_hex(unsigned long x, int w);
        virtual void put(double x, int w, int p, char f);
        virtual size_t write(size_t bytes, const unsigned char* buffer);
        virtual void flush();

    private:
        FILE* outf;

};  // end of FILE_output class

#endif // ifndef _MEDDLY_WITHOUT_CSTDIO_

// ******************************************************************
// *                                                                *
// *                      ostream_output class                      *
// *                                                                *
// ******************************************************************

#ifndef _MEDDLY_WITHOUT_IOSTREAM_

/** ostream_output class.
    Use for C++-style ostream input.
*/
class MEDDLY::ostream_output : public MEDDLY::output {
    public:
        ostream_output(std::ostream &out);
        virtual ~ostream_output();

        virtual void put(char x);
        virtual void put(std::string, int w);
        virtual void put(long x, int w);
        virtual void put(unsigned long x, int w);
        virtual void put_hex(unsigned long x, int w);
        virtual void put(double x, int w, int p, char f);
        virtual size_t write(size_t bytes, const unsigned char* buffer);
        virtual void flush();

    private:
        std::ostream &out;

};  // end of ostream_output class

#endif // ifndef _MEDDLY_WITHOUT_IOSTREAM_

#endif
