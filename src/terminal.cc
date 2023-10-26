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



#include "terminal.h"
#include "io.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        terminal methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************

MEDDLY::terminal::terminal()
{
    setOmega();
}

MEDDLY::terminal::terminal(bool t)
{
    setBoolean(t);
}

MEDDLY::terminal::terminal(int t)
{
    setInteger(t);
}

MEDDLY::terminal::terminal(long t)
{
    // TBD: check for overflow?
    setInteger(t);
}

MEDDLY::terminal::terminal(float t)
{
    setReal(t);
}



void MEDDLY::terminal::read(input &s)
{
#ifdef DEBUG_READ_DD
    std::cerr << "\tin terminal::read\n";
#endif
    s.stripWS();
    // Get the type
    //      (never use: n for non-terminal node handle)
    //      w: omega
    //      b: boolean
    //      i: integer
    //      r: real
    int t = s.get_char();
    int c;

#ifdef DEBUG_READ_DD
    if (t>31) {
        std::cerr << "\t  got char '" << char(t) << "'\n";
    } else {
        std::cerr << "\t  got char ascii " << t << "\n";
    }
#endif

    s.stripWS();

    switch (t) {
        case 'w':
            setOmega(s.get_integer());
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got omega " << t_omega << "\n";
#endif
            break;

        case 'b':
            c = s.get_char();
            if ('T' == c) {
                setBoolean(true);
#ifdef DEBUG_READ_DD
                std::cerr << "\t  got true\n";
#endif
                break;
            }
            if ('F' == c) {
                setBoolean(false);
#ifdef DEBUG_READ_DD
                std::cerr << "\t  got false\n";
#endif
                break;
            }
            throw error(error::INVALID_FILE, __FILE__, __LINE__);

        case 'i':
            setInteger(s.get_integer());
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got integer " << t_integer << "\n";
#endif
            break;

        case 'r':
            setReal(s.get_real());
#ifdef DEBUG_READ_DD
            std::cerr << "\t  got real " << t_real << "\n";
#endif
            break;

        default:
            throw error(error::INVALID_FILE, __FILE__, __LINE__);
    }
#ifdef DEBUG_READ_DD
    std::cerr << "\tdone terminal::read\n";
#endif
}

void MEDDLY::terminal::write(output &s) const
{
    switch (mytype) {
        case terminal_type::OMEGA:
            s.put("w ");
            s.put(t_omega);
            break;

        case terminal_type::BOOLEAN:
            s.put("b ");
            s.put( t_boolean ? 'T' : 'F' );
            break;

        case terminal_type::INTEGER:
            s.put("i ");
            s.put(t_integer);
            break;

        case terminal_type::REAL:
            s.put("r ");
            s.put(t_real, 0, 10, 'e');
            break;

        default:
            MEDDLY_DCASSERT(false);
    }
    s.put(' ');
}

