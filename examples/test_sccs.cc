
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

#include "../src/meddly.h"
#include "../src/operations/sccgraph.h"

#include <iostream>

class my_edge_label : public MEDDLY::sccgraph::edge_label {
  public:
    my_edge_label(char L);
    virtual ~my_edge_label();

    virtual void append_and_recycle(edge_label* L);
    virtual void show(MEDDLY::output &out);

  private:
    std::string the_label;
};

//
//

my_edge_label::my_edge_label(char L)
{
  the_label = std::string(1, L);
}

my_edge_label::~my_edge_label()
{
}

void my_edge_label::append_and_recycle(edge_label* L)
{
  my_edge_label* ML = dynamic_cast<my_edge_label*> (L);
  the_label += ", ";
  the_label += ML->the_label;
  delete L;
}

void my_edge_label::show(MEDDLY::output &out)
{
  out.put(the_label.c_str());
}

//
//

int run_test(std::istream &in, const char* answer)
{
  MEDDLY::sccgraph G;

  char key;
  unsigned from, to;
  char label;
  MEDDLY::ostream_output out(std::cout);

  const bool interactive = !answer;

  for (;;) {
    if (interactive) {
      std::cout << "Options:\n\n  V: add vertices\n  E: add edge  \n  D: display graph\n  U: update sccs\n  Q: quit\n\n";
    }
    in >> key;

    switch (key) {

      case 'V':
      case 'v':
                if (interactive) std::cout << "Enter vertex number N; will add vertices 0..N as needed\n";
                in >> to;
                std::cout << "add_vertex(" << to << ")\n";
                G.add_vertex(to);
                continue;

      case 'E':
      case 'e':
                if (interactive) std::cout << "Enter and edge as SRC DEST LABEL(character)\n";
                in >> from >> to >> label;
                std::cout << "add_edge(" << from << ", " << to << ", " << label << ")\n";
                G.add_edge(from, to, new my_edge_label(label));
                continue;

      case 'D':
      case 'd':
                G.dumpGraph(out);
                continue;

      case 'U':
      case 'u':
                std::cout << "updating sccs\n";
                G.update_SCCs();
                continue;

      case 'Q':
      case 'q':
                break;

      default:
                if (!interactive) return 1;
                std::cout << "Not a valid option\n";
                continue;
    }

    break;
  };

  //
  // If interactive - dump SCCs.
  // Otherwise, we need to pass in expected SCCs and check them.
  //
  if (interactive) {
    for (unsigned s = 0; s<G.get_num_SCCs(); s++) {
      const unsigned* vertex_list = G.get_SCC_vertices(s);
      for (unsigned i = 0; i<G.get_SCC_size(s); i++) {
        if (i) std::cout << ", ";
        std::cout << vertex_list[i];
      }
      std::cout << ";\n";
    }
  }

  return 0;
}

int main()
{
  using namespace std;

  run_test(std::cin, 0);

  cout << "Done\n";

  return 0;
}

