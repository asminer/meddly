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



#include "meddly/terminal.h"
#include "meddly/io.h"

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        terminal methods                        *
// *                                                                *
// *                                                                *
// ******************************************************************



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
            FAIL(__FILE__, __LINE__, "Unknown terminal type");
    }
    s.put(' ');
}

void MEDDLY::terminal::verify(input &s) const
{
    s.stripWS();
    int t = s.get_char();
    if (EOF == t) {
        throw error(error::INVALID_FILE, __FILE__, __LINE__);
    }
    node_handle fomega;
    int fbool;
    long finteger;
    float freal;

    switch (mytype) {
        case terminal_type::OMEGA:
            if (t != 'w') {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  type 'w'\n";
                std::cerr << "\tfile says type '" << char(t) << "'\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            fomega = s.get_integer();
            if (fomega != t_omega) {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  omega: " << t_omega << "\n";
                std::cerr << "\tfile says omega: " << fomega << "\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            break;

        case terminal_type::BOOLEAN:
            if (t != 'b') {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  type 'b'\n";
                std::cerr << "\tfile says type '" << char(t) << "'\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            s.stripWS();
            fbool = s.get_char();
            if (EOF == fbool) {
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            if (t_boolean != ('T' == fbool)) {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  bool: " << (t_boolean ? 'T' : 'F') << "\n";
                std::cerr << "\tfile says bool: " << char(fbool) << "\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            break;

        case terminal_type::INTEGER:
            if (t != 'i') {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  type 'i'\n";
                std::cerr << "\tfile says type '" << char(t) << "'\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            finteger = s.get_integer();
            if (t_integer != finteger) {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  integer: " << t_integer << "\n";
                std::cerr << "\tfile says integer: " << finteger << "\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            break;

        case terminal_type::REAL:
            if (t != 'r') {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  type 'r'\n";
                std::cerr << "\tfile says type '" << char(t) << "'\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            freal = s.get_real();
            // probably need to do within epsilon comparison
            if (t_real != freal) {
#ifdef VERBOSE_VERIFY_FAIL
                std::cerr << "\nVerify error in terminal\n";
                std::cerr << "\texpected  real: " << t_real << "\n";
                std::cerr << "\tfile says real: " << freal << "\n";
#endif
                throw error(error::INVALID_FILE, __FILE__, __LINE__);
            }
            break;

        default:
            FAIL(__FILE__, __LINE__, "Unknown terminal type");
    }
}


void MEDDLY::terminal::show(output &s) const
{
    switch (mytype) {
        case terminal_type::OMEGA:
            s.put("w ");
            s.put(t_omega);
            break;

        case terminal_type::BOOLEAN:
            s.put( t_boolean ? 'T' : 'F' );
            break;

        case terminal_type::INTEGER:
            s.put('t');
            s.put(t_integer);
            break;

        case terminal_type::REAL:
            s.put('t');
            s.put(t_real, 0, 10, 'e');
            break;

        default:
            FAIL(__FILE__, __LINE__, "Unknown terminal type");
    }
}

