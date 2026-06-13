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

#include <iostream>
#include <fstream>

#include "../src/meddly.h"
#include "zstream.h"

// **********************************************************************
// print usage instructions
// **********************************************************************

int usage(const char* exe)
{
    using namespace std;
    const char* base = exe;
    for (; *exe; ++exe)
    {
        if ('/' == *exe)
        {
            base = exe+1;
        }
    }

    cerr << "\n";
    cerr << "Usage: " << base << " file file ...\n";
    cerr << "\n";
    cerr << "Read in DDs from each input file, into a single forest. Files with\n"
         << "extensions .gz or .xz are decompressed using gzcat and xzcat, respectively.\n";
    cerr << "Displays information about each DD, and the overall forest.\n";
    cerr << "\n";

    return 1;
}


// **********************************************************************
// main
// **********************************************************************

int main(int argc, const char** argv)
{
    using namespace std;
    using namespace MEDDLY;

    //
    // For now: no switches
    //
    if (argc < 2) {
        return usage(argv[0]);
    }

    try {
        MEDDLY::initialize();
        std::vector <dd_edge> all_roots;
        domain* d = nullptr;
        forest* F = nullptr;
        ostream_output mout(cout);

        for (int i=1; i<argc; ++i) {
            mout << "\n======================================================================\n";
            mout.indent_more();
            mout << argv[i] << "\n";
            compressed_input reader(argv[i]);
            if (!reader) {
                mout << "couldn't read, skipping\n\n";
                continue;
            }
            mdd_reader* mddr = nullptr;
            if (d) {
                d->verify(reader.instream());
                mddr = new mdd_reader(reader.instream(), F, &mout);
            } else {
                d = domain::create(reader.instream());
                mddr = new mdd_reader(reader.instream(), d, &mout);
                F = mddr->getForest();
            }
            mout << mddr->getFileNodes() << " nodes, " << mddr->numRoots() << " roots\n\n";

            for (unsigned r=0; r<mddr->numRoots(); r++) {
                mout << "Root #" << r << ": "
                     << mddr->getRoot(r).getNodeCount() << " nodes "
                     << mddr->getRoot(r).getEdgeCount() << " edges\n";

                all_roots.push_back(mddr->getRoot(r));
            }

            delete mddr;

            mout.indent_less();
        }

        mout << "\n======================================================================\n\n";

        mout << "Forest stats:\n";
        F->reportStats(mout, "    ",
            HUMAN_READABLE_MEMORY | BASIC_STATS | EXTRA_STATS
        );


        MEDDLY::cleanup();
        mout << "Done!\n";
        return 0;
    }
    catch (MEDDLY::error e) {
        cerr    << "\nCaught meddly error '" << e.getName()
                << "'\n    thrown in " << e.getFile()
                << " line " << e.getLine() << "\n";
        return 2;
    }
    catch (const char* e) {
        cerr    << "\nCaught our own error: " << e << "\n";
        return 4;
    }
    cerr << "\nSome other error?\n";
    return 6;
}
