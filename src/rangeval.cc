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

#include "meddly/rangeval.h"
#include "meddly/io.h"

void MEDDLY::rangeval::write(output &s) const
{
    switch (the_type) {
        case range_type::BOOLEAN:
            MEDDLY_DCASSERT(range_special::NORMAL == s_value);
            s.put("b ");
            s.put( l_value ? 'T' : 'F' );
            break;

        case range_type::INTEGER:
            s.put("i ");
            switch (s_value) {
                case range_special::NORMAL:
                    s.put(l_value);
                    break;

                case range_special::PLUS_INFINITY:
                    s.put("oo");
                    break;

                default:
                    FAIL(__FILE__, __LINE__, "Unknown special value");
            }
            break;

        case range_type::REAL:
            MEDDLY_DCASSERT(range_special::NORMAL == s_value);
            s.put("r ");
            s.put(d_value, 0, 10, 'e');
            break;

        default:
            FAIL(__FILE__, __LINE__, "Unknown range type");
    }
    s.put(' ');
}

void MEDDLY::rangeval::read(input &s)
{
    s.stripWS();
    char T = s.get_char();
    char C;

    switch (T) {
        case 'b':
            // Boolean.
            s.stripWS();
            C = s.get_char();
            if ('T' == C) {
                setBoolean(true);
                return;
            }
            if ('F' == C) {
                setBoolean(false);
                return;
            }
            throw error(error::INVALID_FILE, __FILE__, __LINE__);

        case 'i':
            // Integer.
            // Check for infinity.
            s.stripWS();
            if ('o' == s.next_char()) {
                s.consumeKeyword("oo");
                setSpecial(
                    range_special::PLUS_INFINITY,
                    range_type::INTEGER
                );
                return;
            }
            setInteger(s.get_integer());
            return;

        case 'r':
            setReal(s.get_real());
            return;

        default:
            throw error(error::INVALID_FILE, __FILE__, __LINE__);
    }
}

