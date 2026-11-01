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

#include "meddly/meddly.h"
#include "utils/zstream.h"

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
    cerr << "Usage: " << base << " [options] file file ...\n";
    cerr << "\n";
    cerr << "Read in DDs from each input file, into a single forest. Files with\n"
         << "extensions .gz or .xz are decompressed using gzcat and xzcat, respectively.\n";
    cerr << "Displays information about each DD, and the overall forest.\n";
    cerr << "\n";
    cerr << "Switches:\n";
    cerr << "\n";
    cerr << "    -disp: Display the entire forest after it is built.\n";
    cerr << "\n";
    cerr << "    -dot file: Produce a dot file for the DDs.\n";
    cerr << "\n";
    return 1;
}

// **********************************************************************
// strip extension
// **********************************************************************
char* strip_extension(const char* inpath)
{
    unsigned len = strlen(inpath);
    if (len>4) {
        const char* end = inpath+len-4;
        if (0==strcmp(".dot", end))
            len -= 4;
    }
    char* buf = new char[len+1];
    strncpy(buf, inpath, len);
    return buf;
}

// **********************************************************************
// main
// **********************************************************************

int main(int argc, const char** argv)
{
    using namespace std;
    using namespace MEDDLY;

    bool disp = false;
    char* dotfile = nullptr;

    try {
        MEDDLY::initialize();
        std::vector <dd_edge> all_roots;
        domain* d = nullptr;
        forest* F = nullptr;
        ostream_output mout(cout);

        for (int i=1; i<argc; ++i) {
            if (0==strcmp("-disp", argv[i])) {
                disp = true;
                continue;
            }
            if (0==strcmp("-dot", argv[i])) {
                i++;
                dotfile = strip_extension(argv[i]);
                continue;
            }

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
                const dd_edge& re = mddr->getRoot(r);
                mout << "Root #" << r;
                if (disp) {
                    mout << " -> " << re << ":\n    ";
                } else {
                    mout << ": ";
                }
                mout << re.getNodeCount() << " nodes "
                     << re.getEdgeCount() << " edges\n";

                all_roots.push_back(mddr->getRoot(r));
            }

            delete mddr;

            mout.indent_less();
        } // for i

        if (!F) {
            return usage(argv[0]);
        }

        mout << "\n======================================================================\n\n";

        if (disp) {
            F->dump(mout, 0);
            mout << "\n======================================================================\n\n";
        }

        mout << "Forest stats:\n";
        F->reportStats(mout, "    ",
                HUMAN_READABLE_MEMORY | BASIC_STATS | EXTRA_STATS
        );

        if (dotfile) {
            mout << "Writing to " << dotfile << ".dot\n";
            dot_maker dotwrite(F, dotfile); // tbd convert to base
            for (unsigned i=0; i<all_roots.size(); i++) {
                dotwrite.addRootEdge(all_roots[i]);
            }
            dotwrite.doneGraph();
            delete[] dotfile;
            dotfile = nullptr;
        }



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
