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

#ifndef MEDDLY_IO_MDDS_H
#define MEDDLY_IO_MDDS_H

#include "defines.h"

namespace MEDDLY {
    class output;
    class input;
    class forest;
    class dd_edge;

    class mdd_writer;
    class mdd_reader;
};


/**     Object for writing MDDs to an exchange file.
*/
class MEDDLY::mdd_writer {
    public:
        mdd_writer(output &s, const forest* F);
        ~mdd_writer();
        inline void writeRootEdge(const dd_edge &E) {
            roots.push_back(E);
        }
        void finish();

    private:
        output &out;
        std::vector <dd_edge> roots;
        const forest* For;
        bool finished;
};

/**     Object for reading MDDs from an exchange file.
*/
class MEDDLY::mdd_reader {
    public:
        mdd_reader(input &s, forest* F);
        void readRootEdge(dd_edge &E);
        inline unsigned numRoots() const { return roots.size(); }
    private:
        std::vector <dd_edge> roots;
        unsigned rptr;
};

#endif
