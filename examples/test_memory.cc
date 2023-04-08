
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

#define _MEDDLY_NOINST_
#include <iostream>
#include <fstream>
#include <map>  // for symbol table
#include "../src/meddly.h"
#include "../src/meddly_expert.h"

// #define DEBUG_PARSER

const int MAX_ID = 32;
const char* filename;
int lineno;

struct handle_size {
  unsigned long handle;
  size_t size;

  handle_size(unsigned long h, size_t s) {
    handle = h;
    size = s;
  }
};

void showFileLine(std::ostream &out)
{
  if (filename)   out << filename;
  else            out << "standard input";
  out << " near line " << lineno;
}

void runError(const char* ident, const char* what)
{
  using namespace std;
  cerr << "Runtime error in ";
  showFileLine(cerr);
  cerr << ":\n\tIdentifier `" << ident << "' " << what << "\n";
}

void parseError(char c, const char* what)
{
  using namespace std;
  cerr << "Parse error in ";
  showFileLine(cerr);
  cerr << ":\n\t";
  if (c) cerr << "`" << c << "' ";
  cerr << what << "\n";
  exit(1);
}

char skipWhite(std::istream &in)
{
  bool comment = false;
  for (;;) {
    char c = in.get();
    if (!in) return ' ';
    switch (c) {
      case '\n':
            comment = false;
            lineno++;
            continue;

      case ' ': continue;
      case '\t': continue;
      case '\r': continue;

      case '#':
            comment = true;
            continue;

      default:
            if (comment) continue;
            return c;
    };
  }
}

void readIdent(std::istream &in, char* buffer, int buflen)
{
  buflen--;
  buffer[0] = skipWhite(in);
  if (!in) {
    parseError(0, "identifier expected");
  }
  if ( ('_' != buffer[0]) &&
       ('a' > buffer[0] || 'z' < buffer[0]) &&
       ('A' > buffer[0] || 'Z' < buffer[0]) ) {
    parseError(0, "identifier expected");
  }
  int i = 1;
  for (;;) {
    char c = in.get();
    if (!in) break;
    if (i<buflen) {
      buffer[i] = c;
      i++;
    }
    if ('_' == c) continue;
    if (('0' <= c) && ('9' >= c)) continue;
    if (('a' <= c) && ('z' >= c)) continue;
    if (('A' <= c) && ('Z' >= c)) continue;
    in.unget();
    break;
  }
  i--;
  buffer[i] = 0;
#ifdef DEBUG_PARSER
  std::cerr << "Consumed identifier `" << buffer << "'\n";
#endif
}

size_t readNumber(std::istream &in)
{
  size_t num = 0;
  char c = skipWhite(in);
  if (!in) {
    parseError(0, "number expected");
  }
  if (('0' > c) || ('9' < c)) {
    parseError(0, "number expected");
  }
  num = size_t(c - '0');
  for (;;) {
    c = in.get();
    if (!in) break;
    if (('0' > c) || ('9' < c)) {
      in.unget();
      break;
    }
    num *= 10;
    num += size_t(c - '0');
  }
#ifdef DEBUG_PARSER
  std::cerr << "Consumed number " << num << "\n";
#endif
  return num;
}

void matchChar(std::istream &in, char key)
{
  char c = skipWhite(in);
  if (c!=key) {
    parseError(key, "expected");
  }
}

int usage(const char* who)
{
  /* Strip leading directory, if any: */
  const char* name = who;
  for (const char* ptr=who; *ptr; ptr++) {
    if ('/' == *ptr) {
      name = ptr+1;
    }
  }

  using namespace std;

  cerr << "\nUsage: " << name << " [options] script script ...\n\n";
  cerr << "Run testing scripts.  If none specified, read scripts from standard input.\n";
  cerr << "Possible options:\n";
  cerr << "\t--              Remaining arguments are filenames.\n\n";
  cerr << "\t-g N            set the granularity of the memory manager to N.\n";
  cerr << "\t                Valid sizes are from 1 to 16; default is 4.\n\n";
  cerr << "\t-m MANAGER:     set the memory manager to the specified style.\n";
  cerr << "\t                 Possible managers:\n";
  cerr << "\t                     NONE  (default; test input file only)\n";
  cerr << "\t                     ORIGINAL_GRID\n";

  // TBD - granularity

  cerr << "\nScripts:\n";
  cerr << "Scripts are free-form text files, with the following statements:\n";
  cerr << "    # for comments (ignore until end of line)\n";
  cerr << "    + ident size;  to allocate an identifier of a given size.\n";
  cerr << "    - ident;       to free an identifier.\n";
  cerr << "    ?;             to show active identifiers and memory management information\n";
  cerr << "Allocating an identifier makes it active, and freeing an identifier also\n";
  cerr << "deactivates it.  It is an error to allocate an already active identifier,\n";
  cerr << "or to free an identifier that is not active.  Identifiers follow the usual\n";
  cerr << "naming restrictions: can contain letters, digits, or underscores, and must\n";
  cerr << "not start with a digit.  Only the first " << MAX_ID << " characters of an identifier\n";
  cerr << "are used.\n\n";

  return 0;
}

int main(int argc, const char** argv)
{
  MEDDLY::initialize();

  char buffer[MAX_ID+1];
  const char* program = argv[0];

  const MEDDLY::memory_manager_style* mst = 0;  // set a default
  unsigned char granularity = 4;  // make this an option
  unsigned char minsize = 0;      // make this an option

  for (; argc>1; argv++, argc--) {

    //
    // Process --
    //
    if (strcmp("--", argv[1]) == 0) {
      argv++;
      argc--;
      break;
    }

    //
    // Process -g
    //
    if (strcmp("-g", argv[1]) == 0) {
      argv++;
      argc--;

      if (!argv[1]) {
        std::cerr << "Missing argument to -g\n";
        return usage(program)+1;
      }
      int G = atoi(argv[1]);
      if (G<1 || G>16) {
        std::cerr << "Value of " << G << " for -g is out of range 1..16, ignoring\n";
      } else {
        granularity = G;
      }

      continue;
    }

    //
    // Process -m
    //
    if (strcmp("-m", argv[1]) == 0) {
      argv++;
      argc--;

      if (0==argv[1]) {
        std::cerr << "Missing argument to -m\n";
        return usage(program)+1;
      }

      if (0==strcmp("NONE", argv[1])) {
        mst = 0;
        continue;
      }
      // DON'T continue for these,
      // so we can check for errors
      if (0==strcmp("ORIGINAL_GRID", argv[1])) {
        mst = MEDDLY::ORIGINAL_GRID;
      }
      // else if ... else if ...
      else {
        std::cerr << "Unknown memory manager " << argv[1] << "\n";
        return usage(program)+1;
      }

      if (0==mst) {
        std::cerr << "Hmm, memory manager style " << argv[1] << " does not seem to be\n";
        std::cerr << "    available; falling back to NONE\n";
      }
      continue;
    }

    //
    // Is this an attempt at a switch?
    //
    if ('-' == argv[1][0]) {
      std::cerr << "Unknown switch: " << argv[1] << "\n";
      return usage(program)+1;
    }

    //
    // Done with switches
    //
    break;

  } // for switches

  //
  // Build memory manager
  //

  MEDDLY::memory_manager* Mmm = 0;

  MEDDLY::memstats stats;

  if (mst) {
    std::cout << "Using memory manager of type " << mst->getName() << "\n";
    Mmm = mst->initManager(granularity, minsize, stats);
    if (0==Mmm) {
      std::cout << "  error, couldn't create memory manager with\n";
      std::cout << "      granularity = " << int(granularity) << "\n";
      std::cout << "      minsize = " << int(minsize) << "\n";
      return 0;
    }
  }

  //
  // We have the memory manager; start parsing!
  //



  std::istream* infile = 0;
  std::map <std::string, handle_size> symbols;

  lineno = 1;
  if (argc < 2) {
    filename = 0;
#ifdef DEBUG_PARSER
    std::cerr << "Reading from standard input\n";
#endif
    infile = &std::cin;
  } else {
    filename = argv[1];
    argv++;
    argc--;
#ifdef DEBUG_PARSER
    std::cerr << "Reading from " << filename << "\n";
#endif
    infile = new std::ifstream(filename);
    if (! *infile) {
      std::cerr << "Couldn't open file `" << filename << "'\n";
    }
  }

  for (;;) {
    char c = skipWhite(*infile);
    size_t slots;
    if (' ' == c) {
      if (filename) delete infile;
      // open next file, if any
      if (argc > 1) {
        filename = argv[1];
        lineno = 1;
        argv++;
        argc--;
#ifdef DEBUG_PARSER
        std::cerr << "Reading from " << filename << "\n";
#endif
        infile = new std::ifstream(filename);
        if (! *infile) {
          std::cerr << "Couldn't open file `" << filename << "'\n";
        }
        continue;
      } else {
        break;
      }
    }

    switch (c) {
      case '+':
#ifdef DEBUG_PARSER
          std::cerr << "Parsing + statement\n";
#endif
          readIdent(*infile, buffer, MAX_ID);
          slots = readNumber(*infile);
          matchChar(*infile, ';');
          if (symbols.count(buffer)) {
            runError(buffer, " is already active; ignoring");
          } else {
            unsigned long handle = 0;
            if (Mmm) {
              handle = Mmm->requestChunk(slots);
              if (0==handle) {
                runError(buffer, "could not be allocated");
              }
            }
            symbols.insert(std::pair<std::string, handle_size>
              (buffer, handle_size(handle, slots))
            );
          }
          continue;

      case '-':
#ifdef DEBUG_PARSER
          std::cerr << "Parsing - statement\n";
#endif
          readIdent(*infile, buffer, MAX_ID);
          matchChar(*infile, ';');
          if (!symbols.count(buffer)) {
            runError(buffer, " is not active; ignoring");
          } else {
            std::map<std::string, handle_size>::iterator pos = symbols.find(buffer);
            if (pos == symbols.end()) {
              std::cerr << "Internal error - can't find item " << buffer << "\n";
              exit(2);
            }
            if (Mmm) {
              Mmm->recycleChunk(pos->second.handle, pos->second.size);
            }
            symbols.erase(pos);
          }
          continue;

      case '?':
#ifdef DEBUG_PARSER
          std::cerr << "Parsing ? statement\n";
#endif
          matchChar(*infile, ';');

          std::cout << "Current active identifiers, as of ";
          showFileLine(std::cout);
          std::cout << ":\n";
          for (std::map<std::string, handle_size>::iterator pos= symbols.begin();
               pos != symbols.end();
               ++pos)
          {
            std::cout << "\t" << pos->first;
            std::cout << "    size: " << pos->second.size;
            std::cout << "    handle: " << pos->second.handle << "\n";
          }

          if (Mmm) {
            std::cout << "Memory manager internals:\n";
            MEDDLY::ostream_output out(std::cout);
            Mmm->dumpInternal(out);
          }
          continue;

      default:
          parseError(c, "is not a valid statement");

    }
  }
#ifdef DEBUG_PARSER
  std::cerr << "End of input\n";
#endif

  MEDDLY::cleanup();
}

