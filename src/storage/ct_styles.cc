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


#include "ct_styles.h"


// #define DEBUG_SLOW

// #define DEBUG_CT_SEARCHES
// #define DEBUG_CT_SLOTS

// #define DEBUG_REHASH
// #define DEBUG_TABLE2LIST
// #define DEBUG_LIST2TABLE
// #define DEBUG_CTALLOC

// #define DEBUG_ISDEAD
// #define DEBUG_ISSTALE
// #define DEBUG_REMOVESTALES
// #define DEBUG_REMOVESTALES_DETAILS

// #define DEBUG_CT_SCAN
// #define DEBUG_CT

// #define DEBUG_VALIDATE_COUNTS

// #define INTEGRATED_MEMMAN

#ifndef INTEGRATED_MEMMAN
#include "../memory.h"
#endif

#include "ct_typebased.h"
#include "ct_none.h"


#include <climits>

#include "../hash_stream.h"
#include "../ct_initializer.h"

/*
#include "../node_storage.h"
#include "../compute_table.h"
#include "../ct_entry_key.h"
#include "../ct_entry_result.h"
#include "../ct_vector.h"
*/

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***                                                                     ***
// ***                    New templatized compute table                    ***
// ***                                                                     ***
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

namespace MEDDLY {
    template <bool MONOLITHIC, bool CHAINED, bool INTSLOTS>
    class ct_tmpl : public compute_table {
        public:
            ct_tmpl(const ct_settings &s, operation* op, unsigned slot);
            virtual ~ct_tmpl();

        private:
            /// Hash table
            std::vector <unsigned long> table;

            /// When to next expand the table
            unsigned long tableExpand;

            /// When to next shrink the table
            unsigned long tableShrink;

    };
};

// **********************************************************************
// *                                                                    *
// *                        ct_template  methods                        *
// *                                                                    *
// **********************************************************************

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************
// ***                                                                     ***
// ***                             Front  ends                             ***
// ***                                                                     ***
// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// **********************************************************************
// *                                                                    *
// *                  monolithic_chained_style methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::monolithic_chained_style::monolithic_chained_style()
    : compute_table_style(true)
{
}

MEDDLY::compute_table*
MEDDLY::monolithic_chained_style::create(const ct_settings &s) const
{
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<true, true>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_typebased<true, true>(s, 0, 0);
        default:
                return 0;
    }
}

// **********************************************************************
// *                                                                    *
// *                 monolithic_unchained_style methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::monolithic_unchained_style::monolithic_unchained_style()
    : compute_table_style(true)
{
}

MEDDLY::compute_table*
MEDDLY::monolithic_unchained_style::create(const ct_settings &s) const
{
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<true, false>(s, 0, 0);
        case compressionOption::TypeBased:
                return new ct_typebased<true, false>(s, 0, 0);
        default:
                return 0;
    }
}

// **********************************************************************
// *                                                                    *
// *                  operation_chained_style  methods                  *
// *                                                                    *
// **********************************************************************

MEDDLY::operation_chained_style::operation_chained_style()
    : compute_table_style(false)
{
}

MEDDLY::compute_table*
MEDDLY::operation_chained_style::create(const ct_settings &s, operation* op, unsigned slot) const
{
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<false, true>(s, op, slot);
        case compressionOption::TypeBased:
                return new ct_typebased<false, true>(s, op, slot);
        default:
                return 0;
    }
}

// **********************************************************************
// *                                                                    *
// *                 operation_unchained_style  methods                 *
// *                                                                    *
// **********************************************************************


MEDDLY::operation_unchained_style::operation_unchained_style()
    : compute_table_style(false)
{
}

MEDDLY::compute_table*
MEDDLY::operation_unchained_style::create(const ct_settings &s, operation* op, unsigned slot) const
{
    switch (s.compression) {
        case compressionOption::None:
                return new ct_none<false, false>(s, op, slot);
        case compressionOption::TypeBased:
                return new ct_typebased<false, false>(s, op, slot);
        default:
                return 0;
    }
}

