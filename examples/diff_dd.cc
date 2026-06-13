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
    cerr << "Usage: " << base << " file1 file2\n";
    cerr << "\n";
    cerr << "Read in and compare two DDs from files. Files with extensions .gz or .xz\n";
    cerr << "are decompressed using gzcat and xzcat, respectively.\n";
    cerr << "The DDs must be the same type.\n";
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
    if (argc != 3) {
        return usage(argv[0]);
    }

    cout << "Comparing files " << argv[1] << " and " << argv[2] << "\n";

    try {
        MEDDLY::initialize();

        compressed_input read1(argv[1]);
        if (!read1) {
            cerr << "Couldn't open file: '" << argv[1] << "'\n";
            return 1;
        }

        compressed_input read2(argv[2]);
        if (!read2) {
            cerr << "Couldn't open file: '" << argv[2] << "'\n";
            return 1;
        }

        cout << "Reading " << argv[1] << "...\n";
        domain *d = domain::create(read1.instream());
        mdd_reader mdd1(read1.instream(), d);
        forest* F = mdd1.getForest();
        cout << "... done, " << mdd1.getFileNodes() << " nodes, "
             << mdd1.numRoots() << " roots\n";


        cout << "Reading " << argv[2] << "...\n";
        d->verify(read2.instream());
        mdd_reader mdd2(read2.instream(), F);
        cout << "... done, " << mdd2.getFileNodes() << " nodes, "
             << mdd2.numRoots() << " roots\n";

        for (unsigned i=0; ; i++) {
            if (i >= mdd1.numRoots()) break;
            if (i >= mdd2.numRoots()) break;

            cout << "Roots #" << i << ": ";

            if (mdd1.getRoot(i) == mdd2.getRoot(i)) {
                cout << "MATCH\n";
            } else {
                cout << "differ\n";
            }

        }

        MEDDLY::cleanup();
        cout << "Done!\n";
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
