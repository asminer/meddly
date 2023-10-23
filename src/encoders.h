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

#ifndef MEDDLY_ENCODER_H
#define MEDDLY_ENCODER_H

#include "error.h"
#include "io.h"
#include "edge_value.h"

#include <cstring>

/*
    Classes to encode and decode terminal values and edge values.
    As classes, because we can use them in template functions.
*/

namespace MEDDLY {

    class bool_Tencoder;
    class int_Tencoder;
    class float_Tencoder;

    template <typename T>
    class EVencoder;
};

    /** Encoding for booleans into (terminal) node handles */
    class MEDDLY::bool_Tencoder {
        // -1 true
        //  0 false
        public:
            static node_handle value2handle(bool v);
            static bool handle2value(node_handle h);
            static void show(output &s, node_handle h);
            static void write(output &s, node_handle h);
            static node_handle read(input &s);
    };

    /** Encoding for integers into (terminal) node handles */
    class MEDDLY::int_Tencoder {
        public:
            static node_handle value2handle(int v);
            static int handle2value(node_handle h);
            static void show(output &s, node_handle h);
            static void write(output &s, node_handle h);
            static node_handle read(input &s);
    };

    /** Encoding for floats into (terminal) node handles */
    class MEDDLY::float_Tencoder {
            union intfloat {
                float real;
                int integer;
            };
        public:
            static node_handle value2handle(float v);
            static float handle2value(node_handle h);
            static void show(output &s, node_handle h);
            static void write(output &s, node_handle h);
            static node_handle read(input &s);
    };

    /** Encoding types into edge values */
    template<typename T>
    class MEDDLY::EVencoder {
        public:
            static size_t edgeBytes();
            static void writeValue(edge_value &ptr, T val);
            static void readValue(const edge_value &ptr, T &val);
            static void show(output &s, const edge_value &ptr);
            static void write(output &s, const edge_value &ptr);
            static void read(input &s, edge_value &ptr);
    };

inline MEDDLY::node_handle
MEDDLY::bool_Tencoder::value2handle(bool v)
{
  return v ? -1 : 0;
}

inline bool
MEDDLY::bool_Tencoder::handle2value(MEDDLY::node_handle h)
{
  if (-1 == h)
    return true;
  if (0 == h)
    return false;
  throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
}

inline MEDDLY::node_handle
MEDDLY::int_Tencoder::value2handle(int v)
{
  MEDDLY_DCASSERT(4 == sizeof(MEDDLY::node_handle));
  if (v < -1073741824 || v > 1073741823) {
    // Can't fit in 31 bits (signed)
    throw error(error::VALUE_OVERFLOW, __FILE__, __LINE__);
  }
  if (v)
    v |= -2147483648; // sets the sign bit
  return v;
}

inline int
MEDDLY::int_Tencoder::handle2value(MEDDLY::node_handle h)
{
  // << 1 kills the sign bit
  // >> 1 puts us back, and extends the (new) sign bit
  return (h << 1) >> 1;
}

inline MEDDLY::node_handle
MEDDLY::float_Tencoder::value2handle(float v)
{
  MEDDLY_DCASSERT(4 == sizeof(MEDDLY::node_handle));
  MEDDLY_DCASSERT(sizeof(float) <= sizeof(MEDDLY::node_handle));
  if (0.0 == v)
    return 0;
  intfloat x;
  x.real = v;
  // strip lsb in fraction, and add sign bit
  return node_handle( (unsigned(x.integer) >> 1) | 0x80000000 );
}

inline float
MEDDLY::float_Tencoder::handle2value(MEDDLY::node_handle h)
{
  MEDDLY_DCASSERT(4 == sizeof(MEDDLY::node_handle));
  MEDDLY_DCASSERT(sizeof(float) <= sizeof(MEDDLY::node_handle));
  if (0 == h)
    return 0.0;
  intfloat x;
  x.integer = (h << 1); // remove sign bit
  return x.real;
}

template<typename T>
inline size_t
MEDDLY::EVencoder<T>::edgeBytes()
{
  return sizeof(T);
}
template<typename T>
inline void
MEDDLY::EVencoder<T>::writeValue(edge_value &ptr, T val)
{
    ptr.set(val);
}
template<typename T>
inline void
MEDDLY::EVencoder<T>::readValue(const edge_value &ptr, T &val)
{
    ptr.get(val);
}

template<typename T>
inline void MEDDLY::EVencoder<T>::show(output &s, const edge_value &ptr)
{
  T val;
  readValue(ptr, val);
  s << val;
}

namespace MEDDLY {

template<>
inline void EVencoder<int>::write(output &s, const edge_value &ptr)
{
  int val;
  readValue(ptr, val);
  s.put(val);
}

template<>
inline void EVencoder<long>::write(output &s, const edge_value &ptr)
{
  long val;
  readValue(ptr, val);
  s.put(val);
}

template<>
inline void EVencoder<float>::write(output &s, const edge_value &ptr)
{
  float val;
  readValue(ptr, val);
  s.put(val, 8, 8, 'e');
}

template<>
inline void EVencoder<double>::write(output &s, const edge_value &ptr)
{
  double val;
  readValue(ptr, val);
  s.put(val, 8, 8, 'e');
}

template<>
inline void EVencoder<int>::read(input &s, edge_value &ptr)
{
  writeValue(ptr, int(s.get_integer()));
}

template<>
inline void EVencoder<long>::read(input &s, edge_value &ptr)
{
  writeValue(ptr, long(s.get_integer()));
}

template<>
inline void EVencoder<float>::read(input &s, edge_value &ptr)
{
  writeValue(ptr, float(s.get_real()));
}

template<>
inline void EVencoder<double>::read(input &s, edge_value &ptr)
{
  writeValue(ptr, double(s.get_real()));
}

}


#endif
