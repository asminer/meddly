
/* $Id$ */

/*
    termview: terminal viewing utility for Meddly trace data.
    Copyright (C) 2015, Iowa State University Research Foundation, Inc.

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by 
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along 
    with this software.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "system.h"
#include "parse.h"

#define DEBUG_PARSER

int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) name = ptr+1;
  }
  printf("\nUsage: %s [options]\n\n", name);
  printf("Read trace data from standard input, and display to the terminal.\n");
  printf("If conflicting options are passed, then the last one takes precedence.\n\n");
  printf("Options:\n");
  printf("    -f  file:  read from the given file, instead of standard input.\n");
  printf("\n");
  return 1;
}

int main(int argc, const char** argv)
{
  const char* pathname = 0;

  /* 
    Parse arguments 
  */
  for (int i=1; i<argc; i++) {
    if ('-' != argv[i][0]) {
      return usage(argv[0]);
    }

    /* Future expansion - long options */
    if ('-' == argv[i][1]) {
      
      return usage(argv[0]);
    }
    
    /* Still going?  Short options only. */
    if (argv[i][2]) {
      return usage(argv[0]);
    }

    switch (argv[i][1]) {

      case 'f':
        pathname = argv[i+1];
        i++;
        continue;

      default:
        return usage(argv[0]);
    }


  } /* For i (looping over arguments) */



  /*
    Initialize parser
  */
  FILE* inf = stdin;
  if (pathname) {
    inf = fopen(pathname, "r");
    if (0==inf) {
      printf("\nError, couldn't open file `%s', giving up\n\n", pathname);
      return 2;
    }
  }

  /*
    Structures for grabbing parse data
  */
  const int plength = 256;
  char pbuffer[plength];
  update_t* alist = 0;
  forest_t F;
  initialize(&F);

  /*
    For now - just a loop to parse file to test the parser
  */
  for (int line=1;;line++) {
    int c = fgetc(inf);
    int t;
    if (EOF == c) break;
    if ('\n' == c) continue;  /* In case of empty lines */

    switch (c) {
      case 'T':
#ifdef DEBUG_PARSER
        printf("Parsing T record\n");
#endif
        t = parse_T(inf);
        if (t<0) {
          printf("Parse error line %d\n", line);
          return 2;
        }
        if (t!=1) {
          printf("Unknown file type\n");
          return 3;
        }
        continue;

      case 'F':
#ifdef DEBUG_PARSER
        printf("Parsing F record\n");
#endif
        parse_F(inf, &F);
#ifdef DEBUG_PARSER
        printf("Got F structure:\n");
        printf("  fid %d\n", F.fid);
        if (F.name) printf("  name `%s'\n", F.name);
        printf("  left %d\n", F.left);
        printf("  right %d\n", F.right);
        if (F.counts) {
          printf("  counts [%d", F.counts[F.left]);
          for (int k=F.left+1; k<=F.right; k++) {
            printf(", %d", F.counts[k]);
          }
          printf("]\n");
        }
#endif
        continue;

      
      case 'a': /* TBD */
#ifdef DEBUG_PARSER
        printf("Parsing a record\n");
#endif
        alist = parse_a(inf);
#ifdef DEBUG_PARSER
        printf("Got a structure:\n");
        for (update_t* curr = alist; curr; curr=curr->next) {
          printf("\tfid: %2d level: %3d delta: %3d    ", 
                  curr->fid, curr->level, curr->delta);
          if (curr->next) printf("->");
          printf("\n");
        }
#endif
        kill_update(alist);
        continue;

      case 'p': /* TBD */
#ifdef DEBUG_PARSER
        printf("Parsing p record\n");
#endif
        parse_p(inf, pbuffer, plength);
#ifdef DEBUG_PARSER
        printf("Got p string: `%s'\n", pbuffer);
#endif
        continue;
      
      default:
#ifdef DEBUG_PARSER
        if ('#' == c) {
          printf("Comment, ignoring line\n");
        } else {
          printf("Unknown record `%c', ignoring line\n", c);
        }
#endif
        ignoreLine(inf);
    }
  }

  return 0;
}
