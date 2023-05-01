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

#include "defines.h"
#include "opname.h"
#include "error.h"

#include "initializer.h"

// These are needed for builtin opname initialization

//
// Unary operations
//
#include "operations/copy.h"
#include "operations/cardinality.h"
#include "operations/complement.h"
#include "operations/maxmin_range.h"
#include "operations/mdd2index.h"
#include "operations/cycle.h"
#include "operations/select.h"

// ******************************************************************
// *                                                                *
// *                       opname_init  class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::opname_init : public initializer_list {
    public:
        opname_init(initializer_list* p);
    protected:
        virtual void setup();
        virtual void cleanup();
    private:
        inline void mydelete(unary_opname* &p) {
            delete p;
            p = nullptr;
        }
};

MEDDLY::opname_init::opname_init(initializer_list* p)
    : initializer_list(p)
{
    opname::next_index = 0;

    opname::_COPY       = nullptr;
    opname::_CARD       = nullptr;
    opname::_COMPL      = nullptr;
    opname::_MAXRANGE   = nullptr;
    opname::_MINRANGE   = nullptr;
    opname::_MDD2INDEX  = nullptr;
    opname::_CYCLE      = nullptr;
    opname::_SELECT     = nullptr;
}

void MEDDLY::opname_init::setup()
{
    //
    // Unary ops
    //
    opname::_COPY       = initializeCopy();
    opname::_CARD       = initializeCardinality();
    opname::_COMPL      = initializeComplement();
    opname::_MAXRANGE   = initializeMaxRange();
    opname::_MINRANGE   = initializeMaxRange();
    opname::_MDD2INDEX  = initializeMDD2INDEX();
    opname::_CYCLE      = initializeCycle();
    opname::_SELECT     = initializeSelect();

}

void MEDDLY::opname_init::cleanup()
{
    mydelete(opname::_COPY);
    mydelete(opname::_CARD);
    mydelete(opname::_COMPL);
    mydelete(opname::_MAXRANGE);
    mydelete(opname::_MINRANGE);
    mydelete(opname::_MDD2INDEX);
    mydelete(opname::_CYCLE);
    mydelete(opname::_SELECT);

    opname::next_index = 0;
}

// ******************************************************************
// *                                                                *
// *                         opname methods                         *
// *                                                                *
// ******************************************************************

size_t MEDDLY::opname::next_index;

MEDDLY::unary_opname* MEDDLY::opname::_COPY;
MEDDLY::unary_opname* MEDDLY::opname::_CARD;
MEDDLY::unary_opname* MEDDLY::opname::_COMPL;
MEDDLY::unary_opname* MEDDLY::opname::_MAXRANGE;
MEDDLY::unary_opname* MEDDLY::opname::_MINRANGE;
MEDDLY::unary_opname* MEDDLY::opname::_MDD2INDEX;
MEDDLY::unary_opname* MEDDLY::opname::_CYCLE;
MEDDLY::unary_opname* MEDDLY::opname::_SELECT;

MEDDLY::opname::opname(const char* n)
{
    name = n;
    index = next_index;
    next_index++;
}

MEDDLY::opname::~opname()
{
    // library must be closing
}

MEDDLY::initializer_list*
MEDDLY::opname::makeInitializer(initializer_list* prev)
{
    return new opname_init(prev);
}

// ******************************************************************
// *                                                                *
// *                      unary_opname methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::unary_opname::unary_opname(const char* n) : opname(n)
{
}

MEDDLY::unary_opname::~unary_opname()
{
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(expert_forest* ar, expert_forest* rs) const
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}

MEDDLY::unary_operation*
MEDDLY::unary_opname::buildOperation(expert_forest* ar, opnd_type res) const
{
    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
}


// ******************************************************************
// *                                                                *
// *                     binary_opname  methods                     *
// *                                                                *
// ******************************************************************

MEDDLY::binary_opname::binary_opname(const char* n) : opname(n)
{
}

MEDDLY::binary_opname::~binary_opname()
{
}

// ******************************************************************
// *                                                                *
// *                   specialized_opname methods                   *
// *                                                                *
// ******************************************************************

MEDDLY::specialized_opname::arguments::arguments()
{
    setAutoDestroy(true);
}

MEDDLY::specialized_opname::arguments::~arguments()
{
}


MEDDLY::specialized_opname::specialized_opname(const char* n) : opname(n)
{
}

MEDDLY::specialized_opname::~specialized_opname()
{
}

// ******************************************************************
// *                                                                *
// *                           front  end                           *
// *                                                                *
// ******************************************************************

const MEDDLY::unary_opname* MEDDLY::COPY()
{
    return opname::COPY();
}

const MEDDLY::unary_opname* MEDDLY::CARDINALITY()
{
    return opname::CARDINALITY();
}

const MEDDLY::unary_opname* MEDDLY::COMPLEMENT()
{
    return opname::COMPLEMENT();
}

const MEDDLY::unary_opname* MEDDLY::MAX_RANGE()
{
    return opname::MAX_RANGE();
}

const MEDDLY::unary_opname* MEDDLY::MIN_RANGE()
{
    return opname::MIN_RANGE();
}

const MEDDLY::unary_opname* MEDDLY::CONVERT_TO_INDEX_SET()
{
    return opname::CONVERT_TO_INDEX_SET();
}

const MEDDLY::unary_opname* MEDDLY::CYCLE()
{
    return opname::CYCLE();
}

const MEDDLY::unary_opname* MEDDLY::SELECT()
{
    return opname::SELECT();
}

