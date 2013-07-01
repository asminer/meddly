
// $Id$

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


#ifndef BYTEPACK_H
#define BYTEPACK_STORAGE_H

// ******************************************************************
// *                                                                *
// *                                                                *
// *     Utilities to compact integers  into fewer bits by size     *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                                                                *
// *                  Determine space requirements                  *
// *                                                                *
// ******************************************************************

template <class INT>
inline int bytesRequired4(INT a)
{
  // proceed with byte examinations
  static const unsigned long byte2 = 0xff00;
  static const unsigned long byte3 = byte2<<8;
  static const unsigned long byte4 = byte3<<8;
  if (a & (byte4 | byte3)) {
    if (a & byte4) {
      return 4;
    } else {
      return 3;
    }
  } else {
    if (a & byte2) {
      return 2;
    } else {
      return 1;
    }
  }
}

template <class INT>
inline int bytesRequired8(INT a)
{
  // proceed with byte examinations
  static const unsigned long byte2 = 0xff00;
  static const unsigned long byte3 = byte2<<8;
  static const unsigned long byte4 = byte3<<8;
  static const unsigned long byte5 = byte4<<8;
  static const unsigned long byte6 = byte5<<8;
  static const unsigned long byte7 = byte6<<8;
  static const unsigned long byte8 = byte7<<8;
  if (a & (byte8 | byte7 | byte6 | byte5)) {
    if (a & (byte8 | byte7)) {
      if (a & byte8) {
        return 8;
      } else {
        return 7;
      }
    } else {
      if (a & byte6) {
        return 6;
      } else {
        return 5;
      }
    }
  } else {
    if (a & (byte4 | byte3)) {
      if (a & byte4) {
        return 4;
      } else {
        return 3;
      }
    } else {
      if (a & byte2) {
        return 2;
      } else {
        return 1;
      }
    }
  }
}

template <int N, class INT>
inline int bytesRequired(INT a);

template <>
inline int bytesRequired<4>(int a)
{
  return bytesRequired4(a);
}

template <>
inline int bytesRequired<4>(long a)
{
  return bytesRequired4(a);
}

template <>
inline int bytesRequired<8>(long a)
{
  return bytesRequired8(a);
}


template <class INT>
inline void stripDownEncodingForSizing(INT &a)
{
  if (a<0) {
    // terminal value
    a <<= 1;
    if (a<0) {
      a <<= 1;
      a = -a;
    } else {
      a <<= 1;
    }
  } else {
    // nonterminal node
    a <<= 1;
  }
}

template <class INT>
inline void stripSignedEncodingForSizing(INT &a)
{
  a = (a<0) ? -(a<<1) : a<<1; 
}

/**
      Determine the bytes required for downward pointer \a a.

      Downward pointers follow the rules:
        if positive, it points to a node
        if negative, it encodes a terminal node
        
      so the implementation takes that encoding into account.
*/
inline int bytesRequiredForDown(int a)
{
  stripDownEncodingForSizing(a);
  return bytesRequired<sizeof(int)>(a);
}

/**
      Determine the bytes required for downward pointer \a a.

      Downward pointers follow the rules:
        if positive, it points to a node
        if negative, it encodes a terminal node
        
      so the implementation takes that encoding into account.
*/
inline int bytesRequiredForDown(long a)
{
  stripDownEncodingForSizing(a);
  return bytesRequired<sizeof(long)>(a);
}

/**
      Determine the bytes required for a signed value \a a.
      We are effectively assuming 2s complement storage.
*/
inline int bytesRequiredForSigned(int a)
{
  stripSignedEncodingForSizing(a);
  return bytesRequired<sizeof(int)>(a);
}

/**
      Determine the bytes required for a signed value \a a.
      We are effectively assuming 2s complement storage.
*/
inline int bytesRequiredForSigned(long a)
{
  stripSignedEncodingForSizing(a);
  return bytesRequired<sizeof(long)>(a);
}


// ******************************************************************
// *                                                                *
// *       Convert from natural  to smaller space requirement       *
// *                                                                *
// ******************************************************************

template <int bytes, class INT>
inline void rawToData(INT a, unsigned char* b)
{
  //
  // The compiler is smart enough to optimize out these if's
  //
  b[0] = a & 0xff;
  if (bytes > 1)  { a >>= 8;  b[1] = a & 0xff;  }
  if (bytes > 2)  { a >>= 8;  b[2] = a & 0xff;  }
  if (bytes > 3)  { a >>= 8;  b[3] = a & 0xff;  }
  if (bytes > 4)  { a >>= 8;  b[4] = a & 0xff;  }
  if (bytes > 5)  { a >>= 8;  b[5] = a & 0xff;  }
  if (bytes > 6)  { a >>= 8;  b[6] = a & 0xff;  }
  if (bytes > 7)  { a >>= 8;  b[7] = a & 0xff;  }
}

template <int bytes, class INT>
inline void signedToData(INT L, unsigned char* d)
{
  rawToData<bytes>(L, d);
}

template <int bytes, class INT>
inline void downToData(INT  P, unsigned char* d)
{
  // positive P: as usual.
  if (P >= 0) {
    rawToData<bytes>(P, d);
    return;
  }

  // negative P: this is a terminal pointer.
  //              next msb set - terminal value is negative.
  //
  //  No conversion necessary because msb propogates when we shift
  static const unsigned long nmsb = (0x40L) << ((sizeof(INT)-1)*8);
  if (P & nmsb) {
    rawToData<bytes>(P, d);
    return;
  }

  // negative P: this is a terminal pointer.
  //              next msb clr - terminal value is positive.
  //              
  //  The thing to do here is deal with the msb manually:
  //  clear msb, encode, set msb.
  static const unsigned long msboff = ~ ((0x80L) << ((sizeof(INT)-1)*8));
  rawToData<bytes>(P & msboff, d);
  d[bytes-1] |= 0x80;
}


// ******************************************************************
// *                                                                *
// *       Convert from smaller  to natural space requirement       *
// *                                                                *
// ******************************************************************

template <int bytes, class INT>
inline void dataToRaw(const unsigned char* b, INT &a)
{
  //
  // The compiler is smart enough to optimize out these if's
  //
  if (bytes > 7)  { a |= b[7];  a <<= 8;  }
  if (bytes > 6)  { a |= b[6];  a <<= 8;  }
  if (bytes > 5)  { a |= b[5];  a <<= 8;  }
  if (bytes > 4)  { a |= b[4];  a <<= 8;  }
  if (bytes > 3)  { a |= b[3];  a <<= 8;  }
  if (bytes > 2)  { a |= b[2];  a <<= 8;  }
  if (bytes > 1)  { a |= b[1];  a <<= 8;  }
  a |= b[0];
}

template <int bytes, class INT>
inline void dataToSigned(const unsigned char* d, INT& L)
{
  // deal with negatives properly
  if (d[bytes-1] & 0x80) {
    L = (~0L) << 8;
  } else {
    L = 0;
  }
  dataToRaw<bytes>(d, L);
}


template <int bytes, class INT>
inline void dataToDown(const unsigned char* d, INT& P)
{
  // Is this a terminal value?
  if (d[bytes-1] & 0x80) {
    // YES.
    // Is this a negative terminal value?
    if (d[bytes-1] & 0x40) {
      // YES.
      // Easy case: same as ordinary negatives.
      P = (~0L) << 8;
      dataToRaw<bytes>(d, P);
      return;
    }
    // NO.
    // Positive terminal value.
    P = 0;
    dataToRaw<bytes>(d, P);
    if (bytes != 8) {
      // Move MSB
      static const unsigned long bmsboff = ~ ((0x80L) << ((bytes-1)*8));
      static const unsigned long msbon = (0x80L) << ((sizeof(INT)-1)*8);
      P = (P & bmsboff) | msbon;
    }
    return;
  }

  // non-terminal value: as usual
  P = 0;
  dataToRaw<bytes>(d, P);
}

// ******************************************************************
// *                                                                *
// *                           Old  stuff                           *
// *                                                                *
// ******************************************************************

template <int bytes>
inline void longToData(long L, unsigned char* d)
{
  rawToData<bytes>(L, d);
}

template <int bytes>
inline void dataToLong(const unsigned char* d, long& L)
{
  // deal with negatives properly
  if (d[bytes-1] & 0x80) {
    L = (~0L) << 8;
  } else {
    L = 0;
  }
  dataToRaw<bytes>(d, L);
}

template <int bytes>
inline void intToData(int L, unsigned char* d)
{
  rawToData<bytes>(L, d);
}

template <int bytes>
inline void dataToInt(const unsigned char* d, int& L)
{
  // deal with negatives properly
  if (d[bytes-1] & 0x80) {
    L = (~0) << 8;
  } else {
    L = 0;
  }
  dataToRaw<bytes>(d, L);
}

template <class INT>
inline void bytesRequiredForVal(INT a, int& bytes)
{
  if (sizeof(INT) == bytes) return;
  long la = (a<0) ? -(a<<1) : a<<1; 
  static unsigned long msbyte = (0xffL) << (sizeof(INT)-1)*8;
  unsigned long mask = msbyte;
  for (int b=sizeof(INT); b>bytes; b--) {
    if (la & mask) {
      bytes = b;
      return;
    }
    mask = (mask >> 8) & ~msbyte;
  }
}

inline void bytesRequiredForDown(long a, int& bytes)
{
  if (a<0) {
    // terminal value
    a <<= 1;
    if (a<0) {
      a <<= 1;
      a = -a;
    } else {
      a <<= 1;
    }
  } else {
    a <<= 1;
  }
  static unsigned long msbyte = (0xffL) << (sizeof(long)-1)*8;
  unsigned long mask = msbyte;
  for (int b=sizeof(long); b>bytes; b--) {
    if (a & mask) {
      bytes = b;
      return;
    }
    mask = (mask >> 8) & ~msbyte;
  }
}

#endif  // include guard

