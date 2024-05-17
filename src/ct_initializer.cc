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

#include "error.h"
#include "ct_initializer.h"
#include "compute_table.h"
#include "ct_entry_type.h"

#include "storage/ct_styles.h"

// ******************************************************************
// *                                                                *
// *                     ct_initializer  methods                    *
// *                                                                *
// ******************************************************************

MEDDLY::ct_settings MEDDLY::ct_initializer::the_settings;
const MEDDLY::compute_table_style* MEDDLY::ct_initializer::ct_factory;
MEDDLY::compute_table_style* MEDDLY::ct_initializer::builtin_ct_factory;

MEDDLY::ct_initializer::ct_initializer(initializer_list* prev)
    : initializer_list(prev)
{
    ct_factory = nullptr;
    builtin_ct_factory = nullptr;

    setBuiltinStyle(MonolithicUnchainedHash);
    setMaxSize(16777216);
    setStaleRemoval(staleRemovalOption::Moderate);
    // setCompression(compressionOption::None);
    setCompression(compressionOption::TypeBased);

    //
    // Set to null for now.
    // Set to the proper manager in setup()
    //
    setMemoryManager(nullptr);
}

MEDDLY::ct_initializer::~ct_initializer()
{
    delete builtin_ct_factory;
    builtin_ct_factory = nullptr;
}

void MEDDLY::ct_initializer::setup()
{
    if (!ct_factory) {
        throw error(error::INVALID_ASSIGNMENT, __FILE__, __LINE__);
    }

    MEDDLY_DCASSERT(FREELISTS);
    setMemoryManager(FREELISTS);

    compute_table::initStatics(ct_factory, the_settings);
    /*

    if (ct_factory->usesMonolithic()) {
        compute_table::initStatics( ct_factory->create(the_settings) );
    //     operation::Monolithic_CT = ct_factory->create(the_settings);
    } else {
        compute_table::initStatics(nullptr);
    }

    */
}

void MEDDLY::ct_initializer::cleanup()
{
    // delete operation::Monolithic_CT;
    // operation::Monolithic_CT = nullptr;

    compute_table::doneStatics();
}

void MEDDLY::ct_initializer::setMemoryManager(const memory_manager_style* mms)
{
    the_settings.MMS = mms;
}


void MEDDLY::ct_initializer::setStaleRemoval(staleRemovalOption sro)
{
    the_settings.staleRemoval = sro;
}

void MEDDLY::ct_initializer::setMaxSize(unsigned long ms)
{
    the_settings.maxSize = ms;
    the_settings.allowHugeTables = (ms > 1000000000);
}

void MEDDLY::ct_initializer::setBuiltinStyle(builtinCTstyle cts)
{
    delete builtin_ct_factory;
    builtin_ct_factory = nullptr;
    switch (cts) {
        case MonolithicUnchainedHash:
            builtin_ct_factory = new monolithic_unchained_style;
            break;

        case MonolithicChainedHash:
            builtin_ct_factory = new monolithic_chained_style;
            break;

        case OperationUnchainedHash:
            builtin_ct_factory = new operation_unchained_style;
            break;

        case OperationChainedHash:
            builtin_ct_factory = new operation_chained_style;
            break;
    }

    ct_factory = builtin_ct_factory;
}

void MEDDLY::ct_initializer::setUserStyle(const compute_table_style* cts)
{
    delete builtin_ct_factory;
    builtin_ct_factory = nullptr;
    ct_factory = cts;
}

void MEDDLY::ct_initializer::setCompression(compressionOption co)
{
    the_settings.compression = co;
}

void MEDDLY::ct_initializer::setHugeTables(bool on)
{
    the_settings.allowHugeTables = on;
}

MEDDLY::compute_table* MEDDLY::ct_initializer::createForOp(const ct_entry_type* et)
{
    MEDDLY_DCASSERT(et);
    if (ct_factory) {
        return ct_factory->create(the_settings, et->getID());
    } else {
        return nullptr;
    }
}


