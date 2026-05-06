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
#include <string.h>

#include "../src/meddly.h"
#include "zstream.h"


void myConv(const MEDDLY::rangeval &x, MEDDLY::rangeval &y)
{
    using namespace MEDDLY;
    if (x) {
        y.setInteger(0);
    } else {
        y.setSpecial(range_special::PLUS_INFINITY, range_type::INTEGER);
    }
}

// **********************************************************************
// get distance from a pathname.
// If the pathname does not match the required pattern, returns negative.
// **********************************************************************

int extract_distance(const char* path)
{
    int dots[3];
    dots[0] = -1;
    dots[1] = -1;
    dots[2] = -1;
    for (int i=0; path[i]; i++) {
        if ('.' == path[i]) {
            dots[0] = dots[1];
            dots[1] = dots[2];
            dots[2] = i;
        }
    }
    if (dots[1]<0) return -1;

    if (0==strcmp(".txt", path+dots[2])) {
        return atol(path+dots[1]+1);
    }

    if (dots[0]<0) return -2;
    if (0==strcmp(".txt.gz", path+dots[1])) {
        return atol(path+dots[0]+1);
    }
    if (0==strcmp(".txt.xz", path+dots[1])) {
        return atol(path+dots[0]+1);
    }
    return -3;
}


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
    cerr << "Usage: " << base << " -o outfile -i file file ...\n";
    cerr << "\n";
    cerr << "Read in a series of Boolean DDs, one from each input file, and construct\n"
         << "a single distance EV+MDD from them. The input files must have pathname\n"
         << "that ends in .distance.txt, possibly followed by .gz or .xz if\n"
         << "compressed, where distance is an integer. Interpretation is, that the\n"
         << "file contains states reachable in at most 'distance' steps.\n";
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
    // Check command line
    //
    if (argc < 5) {
        return usage(argv[0]);
    }
    if (strcmp("-o", argv[1])) {
        cerr << "First switch should be -o\n";
        return usage(argv[0]);
    }
    if (strcmp("-i", argv[3])) {
        cerr << "Second switch should be -i\n";
        return usage(argv[0]);
    }

    user_unary_factory MyConversion("MyConv", myConv);

    cout << "Output file " << argv[2] << "\n";

    try {
        MEDDLY::initialize();
        std::vector <dd_edge> all_roots;
        domain* d = nullptr;
        forest* inF = nullptr;
        forest* outF = nullptr;
        dd_edge outdd;

        for (int i=4; i<argc; i++) {
            long dist = extract_distance(argv[i]);
            if (dist<0) {
                cout << "Input file " << argv[i] << ": cannot determine distance\n";
                return 1;
            }
            cout << "Input file " << argv[i] << ": distance " << dist << "\n";

            compressed_input reader(argv[i]);
            if (!reader) {
                cout << "    couldn't read\n";
                throw "Bad input file";
            }

            mdd_reader* mddr = nullptr;
            if (d) {
                d->verify(reader.instream());
                mddr = new mdd_reader(reader.instream(), inF);
            } else {
                d = domain::create(reader.instream());
                mddr = new mdd_reader(reader.instream(), d);
                inF = mddr->getForest();
                outF = forest::create(d, SET, range_type::INTEGER,
                            edge_labeling::EVPLUS
                       );

                outdd.attach(outF);
                rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);
                outF->createConstant(infty, outdd);
            }
            cout << "    " << mddr->getFileNodes() << " nodes\n";

            dd_edge curr(outF);
            apply(MyConversion, mddr->getRoot(0), curr);
            curr.setEdgeValue(dist);

            /*
            ostream_output mout(cout);

            cout << "MDD we read in:\n";
            mddr->getRoot(0).showGraph(mout);

            cout << "EV+MDD converted:\n";
            curr.showGraph(mout);
            */

            delete mddr;

            apply(MINIMUM, outdd, curr, outdd);

        } // for i

        cout << "Writing to " << argv[2] << "\n";
        compressed_output writer(argv[2]);
        if (!writer) {
            cout << "    couldn't write\n";
            throw "bad writer";
        }
        d->write(writer.outstream());
        mdd_writer mywriter(writer.outstream(), outF);
        mywriter.writeRootEdge(outdd);
        mywriter.finish();

        //
        // Reporting
        //
        ostream_output meddlyout(cout);
        meddlyout << "\n======================================================================\n\n";
        meddlyout << "MDD stats:\n";
        inF->reportStats(meddlyout, "    ",
                HUMAN_READABLE_MEMORY | BASIC_STATS
        );
        meddlyout << "EV+MDD stats:\n";
        outF->reportStats(meddlyout, "    ",
                HUMAN_READABLE_MEMORY | BASIC_STATS
        );

        compute_table::showAll(meddlyout, 2);


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
