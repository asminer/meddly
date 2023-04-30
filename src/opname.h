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

#include "old_meddly.h"

namespace MEDDLY {
    class opname;
    class unary_opname;
    class binary_opname;
    class specialized_opname;

    class specialized_operation;

    class initializer_list;
    class expert_forest;

    void initialize(initializer_list *);
    void cleanup();
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
    const char* name;
    int index;
    static int next_index;

    friend void MEDDLY::initialize(initializer_list *);
    friend void MEDDLY::cleanup();
  public:
    opname(const char* n);
    virtual ~opname();

    int getIndex() const;
    const char* getName() const;
};

// ******************************************************************
// *                                                                *
// *                       unary_opname class                       *
// *                                                                *
// ******************************************************************

/// Unary operation names.
class MEDDLY::unary_opname : public opname {
public:
  unary_opname(const char* n);
  virtual ~unary_opname();

  virtual unary_operation*
  buildOperation(expert_forest* arg, expert_forest* res) const;

  virtual unary_operation*
  buildOperation(expert_forest* arg, opnd_type res) const;
};

// ******************************************************************
// *                                                                *
// *                      binary_opname  class                      *
// *                                                                *
// ******************************************************************

/// Binary operation names.
class MEDDLY::binary_opname : public opname {
  public:
    binary_opname(const char* n);
    virtual ~binary_opname();

    virtual binary_operation* buildOperation(expert_forest* arg1,
      expert_forest* arg2, expert_forest* res) const = 0;
};


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
        void setAutoDestroy(bool destroy);
        bool autoDestroy() const;

      private:
        bool destroyWhenDone;
    };

    specialized_opname(const char* n);
    virtual ~specialized_opname();

    /** Note - unlike the more general binary and unary ops,
        a specialized operation might be crafted to the specific
        arguments (passed as an abstract class).
        Examples are:
          - operations that will be called several times with
            several of the arguments unchanged, and we want
            to do some preprocessing on the unchanged arguments.

          - operations with very bizarre and/or user-defined
            parameters (like on-the-fly saturation).

        @param  a   Arguments.  Will be destroyed when we are finished,
                    if autoDestroy() is set for the arguments.
    */
    virtual specialized_operation* buildOperation(arguments* a) const = 0;
};

// ******************************************************************
// *                                                                *
// *                     inlined opname methods                     *
// *                                                                *
// ******************************************************************

inline int
MEDDLY::opname::getIndex() const
{
  return index;
}

inline const char*
MEDDLY::opname::getName() const
{
  return name;
}

// ******************************************************************
// *                                                                *
// *               inlined specialized_opname methods               *
// *                                                                *
// ******************************************************************

inline void
MEDDLY::specialized_opname::arguments::setAutoDestroy(bool destroy)
{
  destroyWhenDone = destroy;
}

inline bool
MEDDLY::specialized_opname::arguments::autoDestroy() const
{
  return destroyWhenDone;
}


#endif // #include guard
