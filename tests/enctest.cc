
// $Id$

/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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

#include <stdio.h>
#include <stdlib.h>

#include "../src/timer.h"

// super hack

#define TEST_ENCODING
#include "../src/node_wrappers.cc"

#define FULL_TEST

// now, test those functions :^)

unsigned char buffer[10];

void encodeInt(int a, int bytes)
{
  switch (bytes) {
    case 1:   intToData<1>(a, buffer);  return;
    case 2:   intToData<2>(a, buffer);  return;
    case 3:   intToData<3>(a, buffer);  return;
    case 4:   intToData<4>(a, buffer);  return;
    case 5:   intToData<5>(a, buffer);  return;
    case 6:   intToData<6>(a, buffer);  return;
    case 7:   intToData<7>(a, buffer);  return;
    case 8:   intToData<8>(a, buffer);  return;

    default:
      fprintf(stderr, "Bad number of bytes for encoding: %d\n", bytes);
      exit(1);
  }
}

void decodeInt(int &a, int bytes)
{
  switch (bytes) {
    case 1:   dataToInt<1>(buffer, a);  return;
    case 2:   dataToInt<2>(buffer, a);  return;
    case 3:   dataToInt<3>(buffer, a);  return;
    case 4:   dataToInt<4>(buffer, a);  return;
    case 5:   dataToInt<5>(buffer, a);  return;
    case 6:   dataToInt<6>(buffer, a);  return;
    case 7:   dataToInt<7>(buffer, a);  return;
    case 8:   dataToInt<8>(buffer, a);  return;

    default:
      fprintf(stderr, "Bad number of bytes for encoding: %d\n", bytes);
      exit(1);
  }
}

void encodeLong(long a, int bytes)
{
  switch (bytes) {
    case 1:   longToData<1>(a, buffer);  return;
    case 2:   longToData<2>(a, buffer);  return;
    case 3:   longToData<3>(a, buffer);  return;
    case 4:   longToData<4>(a, buffer);  return;
    case 5:   longToData<5>(a, buffer);  return;
    case 6:   longToData<6>(a, buffer);  return;
    case 7:   longToData<7>(a, buffer);  return;
    case 8:   longToData<8>(a, buffer);  return;

    default:
      fprintf(stderr, "Bad number of bytes for encoding: %d\n", bytes);
      exit(1);
  }
}

void decodeLong(long &a, int bytes)
{
  switch (bytes) {
    case 1:   dataToLong<1>(buffer, a);  return;
    case 2:   dataToLong<2>(buffer, a);  return;
    case 3:   dataToLong<3>(buffer, a);  return;
    case 4:   dataToLong<4>(buffer, a);  return;
    case 5:   dataToLong<5>(buffer, a);  return;
    case 6:   dataToLong<6>(buffer, a);  return;
    case 7:   dataToLong<7>(buffer, a);  return;
    case 8:   dataToLong<8>(buffer, a);  return;

    default:
      fprintf(stderr, "Bad number of bytes for encoding: %d\n", bytes);
      exit(1);
  }
}

void encodeDown(long a, int bytes)
{
  switch (bytes) {
    case 1:   downToData<1>(a, buffer);  return;
    case 2:   downToData<2>(a, buffer);  return;
    case 3:   downToData<3>(a, buffer);  return;
    case 4:   downToData<4>(a, buffer);  return;
    case 5:   downToData<5>(a, buffer);  return;
    case 6:   downToData<6>(a, buffer);  return;
    case 7:   downToData<7>(a, buffer);  return;
    case 8:   downToData<8>(a, buffer);  return;

    default:
      fprintf(stderr, "Bad number of bytes for encoding: %d\n", bytes);
      exit(1);
  }
}

void decodeDown(long &a, int bytes)
{
  switch (bytes) {
    case 1:   dataToDown<1>(buffer, a);  return;
    case 2:   dataToDown<2>(buffer, a);  return;
    case 3:   dataToDown<3>(buffer, a);  return;
    case 4:   dataToDown<4>(buffer, a);  return;
    case 5:   dataToDown<5>(buffer, a);  return;
    case 6:   dataToDown<6>(buffer, a);  return;
    case 7:   dataToDown<7>(buffer, a);  return;
    case 8:   dataToDown<8>(buffer, a);  return;

    default:
      fprintf(stderr, "Bad number of bytes for encoding: %d\n", bytes);
      exit(1);
  }
}

inline void fprinthex(FILE* s, int i)
{
  for (i--; i>=0; i--) fprintf(s, "%02x ", buffer[i]);
}

void testInt(int a)
{
  int br = 1;
  bytesRequiredForVal(a, br);
  
  for (int b=sizeof(int); b>=br; b--) {
    int x;
    encodeInt(a, b);
    decodeInt(x, b);
    if (a != x) {
      fprintf(stdout, "%d -> ", a);
      fprinthex(stdout, b);
      fprintf(stdout, " -> %d\n", x);
      exit(2);
    }
  }
}

void testLong(long a)
{
  int br = 1;
  bytesRequiredForVal(a, br);
  
  for (int b=sizeof(long); b>=br; b--) {
    long x;
    encodeLong(a, b);
    decodeLong(x, b);
    if (a != x) {
      fprintf(stdout, "%ld -> ", a);
      fprinthex(stdout, b);
      fprintf(stdout, " -> %ld\n", x);
      exit(2);
    }
  }
}

void testDown(long a)
{
  int br = 1;
  bytesRequiredForDown(a, br);
  
  for (int b=sizeof(long); b>=br; b--) {
    long x;
    encodeDown(a, b);
    decodeDown(x, b);
    if (a != x) {
      fprintf(stdout, "%ld -> ", a);
      fprinthex(stdout, b);
      fprintf(stdout, " -> %ld\n", x);
      exit(2);
    }
  }
}

int main()
{
  printf("Testing integer encodings...\n");
  testInt(0);
  int check = 1024;
  int count = check;
  timer foo;
#ifdef FULL_TEST
  for (int i=1; i!=2147483647; i++) {
#else
  for (int i=1; i<1048576; i++) {
#endif
    testInt(i);
    testInt(-i);
    count--;
    if (count) continue;
    printf("I %d\n", i);
    if (foo.get_last_interval() < 5000000) {
      if (check < 1073741824) check *= 2;
    }
    count = check;
    foo.note_time();
  }
  printf("Testing long encodings (with shifted ints)...\n");
  testLong(0);
  check = 1024;
  count = check;
  foo.note_time();
#ifdef FULL_TEST
  for (int i=1; i!=2147483647; i++) {
#else
  for (int i=1; i<1048576; i++) {
#endif
    long L = i;
    while (L << 1) {
      testLong(L);
      testLong(-L);
      L <<= 1;
    }
    count--;
    if (count) continue;
    printf("L %d\n", i);
    if (foo.get_last_interval() < 5000000) {
      if (check < 1073741824) check *= 2;
    }
    count = check;
    foo.note_time();
  }
  printf("Testing down encodings (with shifted ints)...\n");
  testDown(0);
  check = 1024;
  count = check;
  foo.note_time();
#ifdef FULL_TEST
  for (int i=1; i!=2147483647; i++) {
#else
  for (int i=1; i<1048576; i++) {
#endif
    long L = i;
    while (L << 2) {
      testDown(L);
      testDown(-L);
      static const long msb = (0x80L) << ((sizeof(long)-1)*8);
      testDown(L | msb);
      L <<= 1;
    }
    count--;
    if (count) continue;
    printf("D %d\n", i);
    if (foo.get_last_interval() < 5000000) {
      if (check < 1073741824) check *= 2;
    }
    count = check;
    foo.note_time();
  }
  printf("Done.\n");
}
