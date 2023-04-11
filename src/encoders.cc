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

#include "encoders.h"

// ******************************************************************
// *                     bool_Tencoder  methods                     *
// ******************************************************************

void MEDDLY::bool_Tencoder::show(output &s, node_handle h)
{
  s.put(handle2value(h) ? 'T' : 'F');
}

void MEDDLY::bool_Tencoder::write(output &s, node_handle h)
{
  s.put(handle2value(h) ? 'T' : 'F');
}

MEDDLY::node_handle MEDDLY::bool_Tencoder::read(input &s)
{
  s.stripWS();
  int c = s.get_char();
  if ('T' == c) return value2handle(true);
  if ('F' == c) return value2handle(false);
  throw error(error::INVALID_FILE, __FILE__, __LINE__);
}

// ******************************************************************
// *                      int_Tencoder methods                      *
// ******************************************************************


void MEDDLY::int_Tencoder::show(output &s, node_handle h)
{
  s << "t" << handle2value(h);
}

void MEDDLY::int_Tencoder::write(output &s, node_handle h)
{
  s << "t" << handle2value(h);
}

MEDDLY::node_handle MEDDLY::int_Tencoder::read(input &s)
{
  s.stripWS();
  int c = s.get_char();
  if ('t' != c) throw error(error::INVALID_FILE, __FILE__, __LINE__);
  return value2handle(s.get_integer());
}

// ******************************************************************
// *                     float_Tencoder methods                     *
// ******************************************************************

void MEDDLY::float_Tencoder::show(output &s, node_handle h)
{
  s << "t" << handle2value(h);
}

void MEDDLY::float_Tencoder::write(output &s, node_handle h)
{
  s.put('t');
  s.put(handle2value(h), 8, 8, 'e');
}

MEDDLY::node_handle MEDDLY::float_Tencoder::read(input &s)
{
  s.stripWS();
  int c = s.get_char();
  if ('t' != c) throw error(error::INVALID_FILE, __FILE__, __LINE__);
  return value2handle(s.get_real());
}
