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

#ifndef MEDDLY_IO_DOT_H
#define MEDDLY_IO_DOT_H

#include <fstream>
#include <string>

namespace MEDDLY {
    class forest;
    class expert_forest;
    class node_marker;
    class dd_edge;

    class dot_maker;
};

/**
    Helper object for creating graphs from MDDs.

    TBD:
        we could have methods for setting various graph attributes
        different from the defaults
        (size, edge colors, BDD-style, etc.)
        that would be called between construction and startGraph().
 */
class MEDDLY::dot_maker {

    public:
        /**
            Constructor.
                @param  F           Forest
                @param  basename    Base name of the file;
                                    will add .dot for the dot input file.
         */
        dot_maker(forest* F, const char* basename);

        /// Destructor for cleanup.
        ~dot_maker();

        /**
            Add (possibly labeled) root edges to the .dot file.
            This should be done first.
            Can be called several times.
         */
        void addRootEdge(const dd_edge &E);

        /**
            Finish and close the .dot file.
        */
        void doneGraph();

        /**
            Run dot and produce an output with given extension.
            The type name and extension must match.
            Can be called several times, but only after doneGraph().
                @param  ext     Extension, e.g., "pdf" or "png"
        */
        void runDot(const char* ext);

    private:
        std::string basename;
        std::ofstream outfile;

        std::vector <dd_edge> roots;

        expert_forest* For;
        node_marker* nm;
};


#endif // include guard
