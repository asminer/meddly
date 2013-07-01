
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
#include <assert.h>

#define MEDDLY_DCASSERT(X) assert(X)

#include "../src/storage/heap.h"

// ----------------------------------------------------------------------

class myHeapMan {
    struct node {
      unsigned long ID;
      long left;
      long right;

      inline void reset(unsigned long id, long l, long r) {
        ID = id;
        left = l;
        right = r;
      };
    };

    node* nodePile;
    int size;
  public:
    myHeapMan(int sz);
    ~myHeapMan();

    inline long left(long h)  const { return nodePile[h].left; };
    inline long right(long h) const { return nodePile[h].right; };
    inline unsigned long ID(long h) const { return nodePile[h].ID; };
    inline long key(long h) const { return h; }

    inline long& left(long h)  { return nodePile[h].left; };
    inline long& right(long h) { return nodePile[h].right; };
    inline unsigned long& ID(long h) { return nodePile[h].ID; };

    inline void reset(long h, unsigned long id, long l, long r) {
      nodePile[h].reset(id, l, r); 
    }

};

myHeapMan::myHeapMan(int sz)
{
  size = sz;
  nodePile = new node[sz];
}

myHeapMan::~myHeapMan()
{
  delete[] nodePile;
}

// ----------------------------------------------------------------------

int main()
{
  using namespace MEDDLY;
  const int MAX = 1000;
  myHeapMan hm(MAX);

  myheap<myHeapMan> H(hm);

  printf("Enter values to add to / remove from the heap (0 removes smallest), one per line\n");
  printf("Legal values are between 1 and %d\n", MAX-1);

  for (;;) {
    int N;
    if (0==scanf("%d", &N)) break;
    if (feof(stdin)) break;
    if (N>=MAX) {
      printf("Value %d too large, ignoring\n", N);
      continue;
    }

    if (0==N) {
      printf("Removing smallest\n");
      long which = H.removeByID(1);
      printf("Removed %ld\n", which);
      hm.reset(which, 0, 0, 0);
    } else if (hm.ID(N)) {
      if (hm.ID(N)) {
        printf("Removing %d\n", N);
        long which = H.removeByID(hm.ID(N));
        printf("Removed %ld\n", which);
        hm.reset(which, 0, 0, 0);
      }
    } else {
      printf("Adding %d\n", N);
      H.addToHeap(N);
    }
    printf("Heap:");
    H.dumpHeap(stdout);
  }

  printf("\nDone\n");
  return 0;
}
