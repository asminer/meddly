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

#ifndef MEDDLY_HASH_STREAM_H
#define MEDDLY_HASH_STREAM_H

#include "error.h"

namespace MEDDLY {
  class hash_stream;
};

// #define DEBUG_HASH

/**
    Class to hash a stream of unsigned integers.
    Length of the stream can be unknown.
    Based on Bob Jenkin's Hash; see
    http://burtleburtle.net/bob/hash/doobs.html
*/
class MEDDLY::hash_stream {
    unsigned z[3];
    int slot;
  public:
    hash_stream() { }
  protected:
    inline static unsigned rot(unsigned x, int k) {
      return (((x)<<(k)) | ((x)>>(32-(k))));
    }
    inline static void mix(unsigned &a, unsigned &b, unsigned &c) {
        a -= c;  a ^= rot(c, 4);  c += b;
        b -= a;  b ^= rot(a, 6);  a += c;
        c -= b;  c ^= rot(b, 8);  b += a;
        a -= c;  a ^= rot(c,16);  c += b;
        b -= a;  b ^= rot(a,19);  a += c;
        c -= b;  c ^= rot(b, 4);  b += a;
    }
    inline static void final_mix(unsigned &a, unsigned &b, unsigned &c) {
        c ^= b; c -= rot(b,14);
        a ^= c; a -= rot(c,11);
        b ^= a; b -= rot(a,25);
        c ^= b; c -= rot(b,16);
        a ^= c; a -= rot(c,4);
        b ^= a; b -= rot(a,14);
        c ^= b; c -= rot(b,24);
    }
    inline void mix()   { mix(z[2], z[1], z[0]); }
    inline void final_mix() { final_mix(z[2], z[1], z[0]); }
  public:
    inline void start(unsigned init) {
#ifdef DEBUG_HASH
        printf("hash_stream::start %u\n", init);
#endif
        z[2] = init;
        z[1] = 0;
        z[0] = 0xdeadbeef;
        slot = 2;
    }
    inline void start() {
#ifdef DEBUG_HASH
        printf("hash_stream::start\n");
#endif
        z[2] = 0;
        z[1] = 0;
        z[0] = 0xdeadbeef;
        slot = 3;
    }
    inline unsigned finish() {
        final_mix();
#ifdef DEBUG_HASH
        printf("hash_stream::finish: %u\n", z[0]);
#endif
        return z[0];
    }
    inline unsigned long finish64() {
        final_mix();
#ifdef DEBUG_HASH
        printf("hash_stream::finish: %u,%u\n", z[0], z[1]);
#endif
        unsigned long foo = z[0];
        foo <<= 32;
        foo |= z[1];
        return foo;
    }
    inline void push(unsigned v) {
#ifdef DEBUG_HASH
        printf("    push %u\n", v);
#endif
        if (slot) {
          slot--;
          z[slot] += v;
        } else {
          mix();
          z[2] += v;
          slot = 2;
        }
    }
    inline void push(unsigned v1, unsigned v2) {
#ifdef DEBUG_HASH
        printf("    push %u, %u\n", v1, v2);
#endif
        switch (slot) {
            case 0:
                mix();
                z[2] += v1;
                z[1] += v2;
                slot = 1;
                return;

            case 1:
                z[0] += v1;
                mix();
                z[2] += v2;
                slot = 2;
                return;

            case 2:
                z[1] += v1;
                z[0] += v2;
                slot = 0;
                return;

            case 3:
                z[2] += v1;
                z[1] += v2;
                slot = 1;
                return;

            default: throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
        };
    }
    inline void push(unsigned v1, unsigned v2, unsigned v3) {
#ifdef DEBUG_HASH
        printf("    push %u, %u, %u\n", v1, v2, v3);
#endif
        switch (slot) {
            case 0:
                mix();
                z[2] += v1;
                z[1] += v2;
                z[0] += v3;
                return;

            case 1:
                z[0] += v1;
                mix();
                z[2] += v2;
                z[1] += v3;
                return;

            case 2:
                z[1] += v1;
                z[0] += v2;
                mix();
                z[2] += v3;
                return;

            case 3:
                z[2] += v1;
                z[1] += v2;
                z[0] += v3;
                return;

            default: throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
        };
    }
    inline void push(const void* data, std::size_t bytes) {
      const unsigned* hack = (const unsigned*) data;
      std::size_t num_unsigneds = (bytes / sizeof(unsigned));
      const unsigned* hackend = hack + num_unsigneds;
      for (; hack < hackend; hack++) {
        push(*hack);
      }
      //
      // Leftover bytes that don't fill an unsigned.
      //
      const unsigned char* lastfew = (const unsigned char*) hackend;
      unsigned leftover;
      switch (bytes % sizeof(unsigned)) {
        case 0:
            return;

        case 1:
            push(lastfew[0]);
            return;

        case 2:
            leftover = lastfew[0];
            leftover <<= 8;
            leftover |= lastfew[1];
            push(leftover);
            return;

        case 3:
            leftover = lastfew[0];
            leftover <<= 8;
            leftover |= lastfew[1];
            leftover <<= 8;
            leftover |= lastfew[2];
            push(leftover);
            return;

        // how large can an unsigned be?
        // just to be sure:

        case 4:
            leftover = lastfew[0];
            leftover <<= 8;
            leftover |= lastfew[1];
            leftover <<= 8;
            leftover |= lastfew[2];
            leftover <<= 8;
            leftover |= lastfew[3];
            push(leftover);
            return;

        case 5:
            leftover = lastfew[0];
            leftover <<= 8;
            leftover |= lastfew[1];
            leftover <<= 8;
            leftover |= lastfew[2];
            leftover <<= 8;
            leftover |= lastfew[3];
            leftover <<= 8;
            leftover |= lastfew[4];
            push(leftover);
            return;

        case 6:
            leftover = lastfew[0];
            leftover <<= 8;
            leftover |= lastfew[1];
            leftover <<= 8;
            leftover |= lastfew[2];
            leftover <<= 8;
            leftover |= lastfew[3];
            leftover <<= 8;
            leftover |= lastfew[4];
            leftover <<= 8;
            leftover |= lastfew[5];
            push(leftover);
            return;

        case 7:
            leftover = lastfew[0];
            leftover <<= 8;
            leftover |= lastfew[1];
            leftover <<= 8;
            leftover |= lastfew[2];
            leftover <<= 8;
            leftover |= lastfew[3];
            leftover <<= 8;
            leftover |= lastfew[4];
            leftover <<= 8;
            leftover |= lastfew[5];
            leftover <<= 8;
            leftover |= lastfew[6];
            push(leftover);
            return;

        default:
            throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
      }
    }

    //
    // Hash an array of unsigneds
    //
    static inline unsigned raw_hash(const unsigned* k, int len) {
        unsigned a, b, c;
    //    a = b = c = 0xdeadbeef;
        a = b = 0;
        c = 0xdeadbeef;

        // handle most of the key
        while (len > 3)
        {
          a += *k++;
          b += *k++;
          c += *k++;
          mix(a,b,c);
          len -= 3;
        }

        // handle the last 3 uint32_t's
        switch(len)
        {
          // all the case statements fall through
          case 3: c += k[2];
          case 2: b += k[1];
          case 1: a += k[0];
          case 0: // nothing left to add (shouldn't get to this case)
                  final_mix(a,b,c);
                  break;
        }

        return c;
    }

    //
    // Hash an array of unsigneds, return an unsigned long
    //
    static inline unsigned long raw_hash64(const unsigned* k, int len) {
        unsigned a, b, c;
    //    a = b = c = 0xdeadbeef;
        a = b = 0;
        c = 0xdeadbeef;

        // handle most of the key
        while (len > 3)
        {
          a += *k++;
          b += *k++;
          c += *k++;
          mix(a,b,c);
          len -= 3;
        }

        // handle the last 3 uint32_t's
        switch(len)
        {
          // all the case statements fall through
          case 3: c += k[2];
          case 2: b += k[1];
          case 1: a += k[0];
          case 0: // nothing left to add (shouldn't get to this case)
                  final_mix(a,b,c);
                  break;
        }

        unsigned long answer = b;
        answer <<= 32;
        answer |= c;
        return answer;
    }
};

#endif // #include guard
