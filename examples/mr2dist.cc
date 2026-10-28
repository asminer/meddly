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
#include <iomanip>
#include <fstream>
#include <string.h>
#include <cassert>

#include "meddly/meddly.h"
#include "utils/zstream.h"


void myConv(const MEDDLY::rangeval &x, MEDDLY::rangeval &y)
{
    using namespace MEDDLY;
    if (x) {
        y.setInteger(0);
    } else {
        y.setSpecial(range_special::PLUS_INFINITY, range_type::INTEGER);
    }
}

void show_node_count(std::ostream &sout, long x, const char* desc)
{
    sout << "        ";
    if (x<0) {
        sout << "-";
        x = -x;
    } else {
        sout << " ";
    }
    long grouping = 1000000000000L;
    char filler = ' ';
    while (grouping) {
        if (filler != ' ') {
            sout << ",";
        } else {
            sout << " ";
        }
        const long d = x / grouping;
        if (d) {
            if (d < 100) {
                sout << filler;
            }
            if (d < 10) {
                sout << filler;
            }
            sout << d;
            x %= grouping;
            filler = '0';
        } else {
            sout << filler << filler << filler;
        }
        grouping /= 1000;
    }
    sout << " nodes " << desc << "\n";
    sout.flush();
}

void show_forest_memory(std::ostream &sout, long mem, const char* desc)
{
    sout << "        ";
    sout << "                ";
    unsigned units = 0;
    while (mem >= 1024) {
        mem /= 1024;
        ++units;
    }
    const char* uletters = " KMGTPEZYabcdefgh";
    sout << std::setw(5) << mem << " " << uletters[units] << "bytes " << desc << "\n";
    sout.flush();
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
    cerr << "Usage: " << base << " -o outfile [-d distfile] -i file file ...\n";
    cerr << "\n";
    cerr << "Read in a series of Boolean DDs, one from each input file, and construct\n"
         << "a single distance EV+MDD from them. The input files must have pathname\n"
         << "that ends in .distance.txt, possibly followed by .gz or .xz if\n"
         << "compressed, where distance is an integer. Interpretation is, that the\n"
         << "file contains states reachable in at most 'distance' steps.\n";
    cerr << "\n";
    cerr << "Switches:\n";
    cerr << "    -d dfile:  Start with the given distance function. If omitted,\n"
         << "               start with infinity.\n";
    cerr << "\n";
    cerr << "    -i:        Everything from here on is an input file.\n";
    cerr << "\n";
    cerr << "    -o ofile:  Specify the output file.\n";
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

    const char* outfile = nullptr;
    const char* distfile = nullptr;
    std::vector <const char*> inputs;

    //
    // Check command line
    //
    char status = ' ';
    for (int i=1; i<argc; i++) {
        if (0==strcmp("-o", argv[i])) {
            if (outfile) {
                cerr << "Error: two -o switches?\n";
                return usage(argv[0]);
            }
            outfile = argv[i+1];
            i++;
            continue;
        }
        if (0==strcmp("-d", argv[i])) {
            if (distfile) {
                cerr << "Error: two -d switches?\n";
                return usage(argv[0]);
            }
            distfile = argv[i+1];
            i++;
            continue;
        }
        if (0==strcmp("-i", argv[i])) {
            status = 'i';
            continue;
        }
        if (' ' == status) {
            cerr << "Expected switch -d, -i, or -o; got " << argv[i] << "\n";
            return usage(argv[0]);
        }
        inputs.push_back(argv[i]);
    }

    if (!outfile) {
        cerr << "Unspecified output file (missing -o outfile)\n";
        return usage(argv[0]);
    }

    user_unary_factory MyConversion("MyConv", myConv);

    cout << "Output file " << outfile << "\n";

    try {
        MEDDLY::initialize();
        std::vector <dd_edge> all_roots;
        domain* d = nullptr;
        forest* inF = nullptr;
        forest* outF = nullptr;
        mdd_reader* mddr = nullptr;
        dd_edge outdd;

        if (distfile) {
            cout << "Distance file " << distfile << "\n";
            compressed_input dreader(distfile);
            if (!dreader) {
                cout << "    couldn't read distance file\n";
                throw "Bad distance file";
            }
            d = domain::create(dreader.instream());
            mddr = new mdd_reader(dreader.instream(), d);
            outF = mddr->getForest();
            outdd.attach(outF);
            outdd = mddr->getRoot(0);
            cout << "    " << mddr->getFileNodes() << " nodes\n";
            cout.flush();
            delete mddr;
            mddr = nullptr;
        }

        for (unsigned i=0; i<inputs.size(); i++) {
            long dist = extract_distance(inputs[i]);
            if (dist<0) {
                cout << "Input file " << inputs[i] << ": cannot determine distance\n";
                return 1;
            }
            cout << "Input file " << inputs[i] << ": distance " << dist << "\n";
            cout.flush();

            compressed_input reader(inputs[i]);
            if (!reader) {
                cout << "    couldn't read\n";
                throw "Bad input file";
            }

            if (inF) {
                d->verify(reader.instream());
                mddr = new mdd_reader(reader.instream(), inF);
            } else {
                if (d) {
                    d->verify(reader.instream());
                    assert(outF);
                } else {
                    d = domain::create(reader.instream());
                    outF = forest::create(d, SET, range_type::INTEGER,
                            edge_labeling::EVPLUS
                    );
                    outdd.attach(outF);
                    rangeval infty(range_special::PLUS_INFINITY,
                            range_type::INTEGER
                    );
                    outF->createConstant(infty, outdd);
                }
                mddr = new mdd_reader(reader.instream(), d);
                inF = mddr->getForest();
            }

            cout << "    " << mddr->getFileNodes() << " nodes\n";
            cout.flush();

            dd_edge curr(outF);
            apply(MyConversion, mddr->getRoot(0), curr);
            curr.setEdgeValue(dist);

            cout << "    converted to single ev+mdd\n";
            cout.flush();

            delete mddr;
            mddr = nullptr;

            apply(MINIMUM, outdd, curr, outdd);
            cout << "    incorporated into distance function\n";
            show_node_count(cout, outdd.getNodeCount(), " (function)");
            if (outF) {
                show_node_count(cout, outF->getCurrentNumNodes(), " (forest)");
                show_forest_memory(cout, outF->getCurrentMemoryUsed(), "(forest)");
            }

        } // for i

        cout << "Writing to " << outfile << "\n";
        compressed_output writer(outfile);
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
