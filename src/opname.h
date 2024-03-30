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

#ifndef MEDDLY_OPNAME_H
#define MEDDLY_OPNAME_H

#include "defines.h"

#include "oper_unary.h"
#include "oper_binary.h"

namespace MEDDLY {
    class opname;
    // class unary_opname;
    // class binary_opname;
    class specialized_opname;

    // typedef unary_opname* (*unary_handle)();
    // typedef binary_opname* (*binary_handle)();

    class operation;
    class unary_operation;
    class binary_operation;
    class specialized_operation;

    class ct_object;
    class dd_edge;
    class forest;

    class initializer_list;
};

// ******************************************************************
// *                                                                *
// *                          opname class                          *
// *                                                                *
// ******************************************************************

/**
    Class (hierarchy) for names of operations.
    These are used as factories, to build specific operations
    (tied to forests) to be applied to decision diagrams
*/
class MEDDLY::opname {
    public:
        opname(const char* n);
        virtual ~opname();

        inline const char* getName() const { return name; }

    private:
        const char* name;
};

// ******************************************************************
// *                                                                *
// *                       unary_opname class                       *
// *                                                                *
// ******************************************************************

/// Unary operation names.
/*
class MEDDLY::unary_opname : public opname {
    public:
        unary_opname(const char* n);
        virtual ~unary_opname();

        unary_operation* getOperation(const dd_edge &arg, const dd_edge &res);
        unary_operation* getOperation(const dd_edge &arg, opnd_type res);

        inline void removeOperationFromCache(unary_operation* uop) {
            cache.remove(uop);
        }

    protected:
        virtual unary_operation*
            buildOperation(unary_list& c, forest* arg, forest* res);

        virtual unary_operation*
            buildOperation(unary_list& c, forest* arg, opnd_type res);

    private:
        unary_list cache;
};
*/

// ******************************************************************
// *                                                                *
// *                      binary_opname  class                      *
// *                                                                *
// ******************************************************************

/// Binary operation names.
/*
class MEDDLY::binary_opname : public opname {
    public:
        binary_opname(const char* n);
        virtual ~binary_opname();

        binary_operation* getOperation(const dd_edge &arg1,
                const dd_edge &arg2, const dd_edge &res);

        inline void removeOperationFromCache(binary_operation* bop) {
            cache.remove(bop);
        }

    protected:
        virtual binary_operation* buildOperation(binary_list &c,
            forest* arg1, forest* arg2, forest* res) = 0;

    private:
        binary_list cache;
};
*/

// ******************************************************************
// *                                                                *
// *                    specialized_opname class                    *
// *                                                                *
// ******************************************************************

/// Specialized operation names.
class MEDDLY::specialized_opname : public opname {
    public:
        /**
            Abstract base class for arguments to buildOperation().
            Derived operation names must provide derived classes for arguments.
        */
        class arguments {
            public:
                arguments();
                virtual ~arguments();

                /**
                    Specify if arguments should be destroyed or not.
                    If yes (the default), the operation will destroy
                    the arguments once they are no longer needed.
                */
                inline void setAutoDestroy(bool destroy) {
                    destroyWhenDone = destroy;
                }
                inline bool autoDestroy() const {
                    return destroyWhenDone;
                }

            private:
                bool destroyWhenDone;
        };

    public:
        specialized_opname(const char* n);
        virtual ~specialized_opname();

        /** Note - unlike the more general binary and unary ops,
            a specialized operation might be crafted to the specific
            arguments (passed as an abstract class).
            Examples are:
                -   operations that will be called several times with
                    several of the arguments unchanged, and we want
                    to do some preprocessing on the unchanged arguments.

                -   operations with very bizarre and/or user-defined
                    parameters (like on-the-fly saturation).

            @param  a   Arguments.  Will be destroyed when we are finished,
                        if autoDestroy() is set for the arguments.
        */
        virtual specialized_operation* buildOperation(arguments* a) = 0;
};


#endif // #include guard
