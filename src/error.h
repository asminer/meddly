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

#ifndef MEDDLY_ERROR_H
#define MEDDLY_ERROR_H

namespace MEDDLY {
    class error;
};

/**
    Class for errors thrown by MEDDLY.

 */
class MEDDLY::error {
    public:
        /// Error codes.
        enum code {
            /// The library was not initialized.
            UNINITIALIZED,
            /// The library was already initialized.
            ALREADY_INITIALIZED,
            /// The requested operation is not yet implemented!
            NOT_IMPLEMENTED,
            /// An operation failed due to lack of memory.
            INSUFFICIENT_MEMORY,
            /// An operation is not supported for the given forest.
            INVALID_OPERATION,
            /// A provided variable is erroneous.
            INVALID_VARIABLE,
            /// A provided level number is erroneous.
            INVALID_LEVEL,
            /// A provided variable bound is out of range.
            INVALID_BOUND,
            /// Attempt to dereference an invalid iterator.
            INVALID_ITERATOR,
            /// We expected an empty domain, but it wasn't
            DOMAIN_NOT_EMPTY,
            /// Unknown operation (bad operation handle).
            UNKNOWN_OPERATION,
            /// Requested operation requires same domains, they weren't.
            DOMAIN_MISMATCH,
            /// Requested operation requires same forest, it wasn't.
            FOREST_MISMATCH,
            /// Requested operation unsupported for operand or result type.
            TYPE_MISMATCH,
            /// Requested operation requires different number of operands.
            WRONG_NUMBER,
            /// A result won't fit in an integer / float.
            VALUE_OVERFLOW,
            /// Integer division by 0 is invalid.
            DIVIDE_BY_ZERO,
            /// Subtracting infinity is invalid.
            SUBTRACT_INFINITY,
            /// Infinity / infinity and infinity % infinity are invalid.
            INFINITY_DIV_INFINITY,
            /// Invalid policy setting.
            INVALID_POLICY,
            /// Bad value for something.
            INVALID_ASSIGNMENT,
            /// Invalid argument (for specialized operations)
            INVALID_ARGUMENT,
            /// File format error.
            INVALID_FILE,
            /// File input error.
            COULDNT_READ,
            /// File output error.
            COULDNT_WRITE,
            /// Miscellaneous error
            MISCELLANEOUS
        };
    public:
        error(code c, const char* fn, unsigned ln);

        /// Return a human-readable error message
        const char* getName() const;

        inline operator code() const        { return errcode; }
        inline code getCode() const         { return errcode; }
        inline const char* getFile() const  { return fname; }
        inline unsigned getLine() const     { return lineno; }

    private:
        code errcode;
        const char* fname;
        unsigned lineno;
};


#endif
