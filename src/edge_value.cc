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



#include "edge_value.h"
#include "io.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                       edge_value methods                       *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::edge_value::edge_value()
{
    set();
}

MEDDLY::edge_value::edge_value(int v)
{
    set(v);
}

MEDDLY::edge_value::edge_value(long v)
{
    set(v);
}

MEDDLY::edge_value::edge_value(float v)
{
    set(v);
}

MEDDLY::edge_value::edge_value(double v)
{
    set(v);
}

MEDDLY::edge_value::edge_value(edge_type newtype, const edge_value &v)
{
    switch (v.getType()) {
        case edge_type::VOID:
            if (edge_type::VOID == newtype) {
                set();
                return;
            }
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

        case edge_type::INT:
            setTempl(newtype, v.getInt());
            return;

        case edge_type::LONG:
            setTempl(newtype, v.getLong());
            return;

        case edge_type::FLOAT:
            setTempl(newtype, v.getFloat());
            return;

        case edge_type::DOUBLE:
            setTempl(newtype, v.getDouble());
            return;

        default:
            throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
}


void MEDDLY::edge_value::read(input &s)
{
#ifdef DEBUG_READ_DD
    std::cerr << "\tin edge_value::read\n";
#endif
    s.stripWS();
    // Get the type
    //      v: void
    //      i: integer
    //      l: long
    //      f: float
    //      d: double
    int t = s.get_char();
    s.stripWS();

    switch (t) {
        case 'v':
            set();
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got void\n";
#endif
            break;

        case 'i':
            set(int(s.get_integer()));
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got int " << ev_int << "\n";
#endif
            break;

        case 'l':
            set(long(s.get_integer()));
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got long " << ev_long << "\n";
#endif
            break;

        case 'f':
            set(float(s.get_real()));
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got float " << ev_float << "\n";
#endif
            break;

        case 'd':
            set(double(s.get_real()));
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got double " << ev_double << "\n";
#endif
            break;

        default:
            throw error(error::INVALID_FILE, __FILE__, __LINE__);

    }
#ifdef DEBUG_READ_DD
    std::cerr << "\tdone edge_value::read\n";
#endif

}

void MEDDLY::edge_value::write(output &s) const
{
    switch (mytype) {
        case edge_type::VOID:
            s.put("v ");
            return;

        case edge_type::INT:
            s.put("i ");
            s.put(ev_int);
            s.put(' ');
            return;

        case edge_type::LONG:
            s.put("l ");
            s.put(ev_long);
            s.put(' ');
            return;

        case edge_type::FLOAT:
            s.put("f ");
            s.put(ev_float);
            s.put(' ');
            return;

        case edge_type::DOUBLE:
            s.put("d ");
            s.put(ev_double);
            s.put(' ');
            return;

        default:
            MEDDLY_DCASSERT(false);
    }
}

void MEDDLY::edge_value::show(output &s) const
{
    switch (mytype) {
        case edge_type::VOID:
            return;

        case edge_type::INT:
            s.put(ev_int);
            return;

        case edge_type::LONG:
            s.put(ev_long);
            return;

        case edge_type::FLOAT:
            s.put(ev_float);
            return;

        case edge_type::DOUBLE:
            s.put(ev_double);
            return;

        default:
            MEDDLY_DCASSERT(false);
    }
}

