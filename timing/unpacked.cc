/*
    Meddly: Multi-terminal and Edge-valued Decision Diagram LibrarY.
    Copyright (C) 2011, Iowa State University Research Foundation, Inc.

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
#include "timer.h"
#include "reporting.h"
#include "park_random.h"

#include <iostream>
#include <iomanip>
#include <vector>

// #define REPORTING

using namespace MEDDLY;

#ifdef DEVELOPMENT_CODE
const unsigned CHANGES = 128 * 1024;
#else
const unsigned CHANGES = 8 * 1024 * 1024;
#endif
const unsigned DOTS = 16;

const unsigned USIZE = 128;

const unsigned LEVELS = 10;
const unsigned VARSIZE = 25;

void make_blank(const std::vector <forest*> &FA, node_storage_flags fs)
{
    const char* nodetype = (FULL_ONLY == fs) ? "full" : "sparse";
    std::cout   << "    Creating"
                << std::setw(11) << "blank,"
                << std::setw(7) << nodetype << " nodes ";

    unpacked_node* UA[USIZE];
    for (unsigned i=0; i<USIZE; i++) UA[i] = nullptr;
    unsigned u = 0;

    const unsigned FAszm1 = FA.size() - 1;

    timer T;
    for (unsigned d=0; d<DOTS; d++) {
        std::cout << '.';
        std::cout.flush();
        for (unsigned c=0; c<CHANGES; c++) {
            //
            // Randomly choose forest
            //
            forest* f = FA[Equilikely(0, FAszm1)];
            if (0==f) throw "null forest";

            //
            // Randomly choose level
            //
            int k = (int) Equilikely(1, LEVELS);

            //
            // recycle old in slot u,
            // and build a new blank node there
            //
            unpacked_node::Recycle(UA[u]);
            const unsigned sz = (FULL_ONLY == fs) ? VARSIZE : 2;
            UA[u] = unpacked_node::newWritable(f, k, sz, fs);

            //
            // Next slot
            //
            u = (u+1) % USIZE;
        }
    }
    std::cout << 'd';
    std::cout.flush();
    for (unsigned i=0; i<USIZE; i++) {
        unpacked_node::Recycle(UA[i]);
    }
    T.note_time();

    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "blank ", nodetype)) {
        report  << "Created " << DOTS * CHANGES
                << " blank " << nodetype << " nodes" << std::endl;
    }

}


void make_redundant(const std::vector <forest*> &FA, node_storage_flags fs)
{
    const char* nodetype = (FULL_ONLY == fs) ? "full" : "sparse";
    std::cout   << "    Creating"
                << std::setw(11) << "redundant,"
                << std::setw(7) << nodetype << " nodes ";

    unpacked_node* UA[USIZE];
    for (unsigned i=0; i<USIZE; i++) UA[i] = nullptr;
    unsigned u = 0;

    const unsigned FAszm1 = FA.size() - 1;

    timer T;
    for (unsigned d=0; d<DOTS; d++) {
        std::cout << '.';
        std::cout.flush();
        for (unsigned c=0; c<CHANGES; c++) {
            //
            // Randomly choose forest
            //
            forest* f = FA[Equilikely(0, FAszm1)];
            if (0==f) throw "null forest";

            //
            // Randomly choose level
            //
            int k = (int) Equilikely(1, LEVELS);

            //
            // recycle old in slot u,
            // and build a new redundant node there
            //
            unpacked_node::Recycle(UA[u]);
            if (f->isMultiTerminal()) {
                UA[u] = unpacked_node::newRedundant(f, k, -1, fs);
            } else if (f->isEVPlus()) {
                UA[u] = unpacked_node::newRedundant(f, k, 0L, -1, fs);
            } else {
                UA[u] = unpacked_node::newRedundant(f, k, 1.0f, -1, fs);
            }

            //
            // Next slot
            //
            u = (u+1) % USIZE;
        }
    }
    std::cout << 'd';
    std::cout.flush();
    for (unsigned i=0; i<USIZE; i++) {
        unpacked_node::Recycle(UA[i]);
    }
    T.note_time();

    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "rednt ", nodetype)) {
        report  << "Created " << DOTS * CHANGES
                << " redundant " << nodetype << " nodes" << std::endl;
    }

}


void make_identity(const std::vector <forest*> &FA, node_storage_flags fs)
{
    const char* nodetype = (FULL_ONLY == fs) ? "full" : "sparse";
    std::cout   << "    Creating"
                << std::setw(11) << "identity,"
                << std::setw(7) << nodetype << " nodes ";

    unpacked_node* UA[USIZE];
    for (unsigned i=0; i<USIZE; i++) UA[i] = nullptr;
    unsigned u = 0;

    const unsigned FAszm1 = FA.size() - 1;

    timer T;
    for (unsigned d=0; d<DOTS; d++) {
        std::cout << '.';
        std::cout.flush();
        for (unsigned c=0; c<CHANGES; c++) {
            //
            // Randomly choose forest
            //
            forest* f = FA[Equilikely(0, FAszm1)];
            if (0==f) throw "null forest";

            //
            // Randomly choose level
            //
            int k = (int) Equilikely(1, LEVELS);

            //
            // Randomly choose identity index
            //
            unsigned i = Equilikely(0, VARSIZE-1);

            //
            // recycle old in slot u,
            // and build a new redundant node there
            //
            unpacked_node::Recycle(UA[u]);
            if (f->isMultiTerminal()) {
                UA[u] = unpacked_node::newIdentity(f, k, i, -1, fs);
            } else if (f->isEVPlus()) {
                UA[u] = unpacked_node::newIdentity(f, k, i, 0L, -1, fs);
            } else {
                UA[u] = unpacked_node::newIdentity(f, k, i, 1.0f, -1, fs);
            }

            //
            // Next slot
            //
            u = (u+1) % USIZE;
        }
    }
    std::cout << 'd';
    std::cout.flush();
    for (unsigned i=0; i<USIZE; i++) {
        unpacked_node::Recycle(UA[i]);
    }
    T.note_time();

    show_sec(std::cout, T, 3, 3);

    if (startReport(T, __FILE__, "idnty ", nodetype)) {
        report  << "Created " << DOTS * CHANGES
                << " identity " << nodetype << " nodes" << std::endl;
    }

}

int main(int argc, const char** argv)
{
    try {
        setReport(argc, argv);

        MEDDLY::initialize();

        seed = 123456789;

        std::cout << "Initializing domain and forests\n";

        int bounds[LEVELS+1];
        for (unsigned i=0; i<=LEVELS; i++) {
            bounds[i] = VARSIZE;
        }
        domain *D = domain::createBottomUp(bounds, LEVELS);

        //
        // Make all kinds of relation forests.
        //

        std::vector <forest*> relF;
        relF.push_back(
            forest::create(D, true, range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL)
        );
        relF.push_back(
            forest::create(D, true, range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL)
        );
        relF.push_back(
            forest::create(D, true, range_type::REAL,
                edge_labeling::MULTI_TERMINAL)
        );
        relF.push_back(
            forest::create(D, true, range_type::INTEGER,
                edge_labeling::EVPLUS)
        );
        relF.push_back(
            forest::create(D, true, range_type::REAL,
                edge_labeling::EVTIMES)
        );

        //
        // Make all kinds of set or relation forests.
        //

        std::vector <forest*> allF;
        // Copy relation forests
        for (unsigned i=0; i<relF.size(); i++) {
            allF.push_back(relF[i]);
        }
        allF.push_back(
            forest::create(D, false, range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL)
        );
        allF.push_back(
            forest::create(D, false, range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL)
        );
        allF.push_back(
            forest::create(D, false, range_type::REAL,
                edge_labeling::MULTI_TERMINAL)
        );
        allF.push_back(
            forest::create(D, false, range_type::INTEGER,
                edge_labeling::EVPLUS)
        );

        //
        // Check forests
        //

        bool found = false;
        for (unsigned i=0; i<allF.size(); i++) {
            if (!allF[i]) {
                std::cerr << "forest " << i << " is null\n";
                found = true;
            }
        }
        if (found) throw "null forest";


        //
        // Run tests
        //

        std::cout << "Timing tests to build/recycle " << DOTS * CHANGES
                  << " nodes\n";

        make_blank(allF, SPARSE_ONLY);
        make_blank(allF, FULL_ONLY);

        make_redundant(allF, SPARSE_ONLY);
        make_redundant(allF, FULL_ONLY);

        make_identity(relF, SPARSE_ONLY);
        make_identity(relF, FULL_ONLY);

        std::cout << "Done; cleaning up\n";

        MEDDLY::cleanup();
        return 0;
    }

    catch (MEDDLY::error e) {
        std::cerr   << "\nCaught meddly error '" << e.getName()
                    << "'\n    thrown in " << e.getFile()
                    << " line " << e.getLine() << "\n";
        return 1;
    }
    catch (const char* e) {
        std::cerr << "\nCaught our own error: " << e << "\n";
        return 2;
    }
    catch (two_strings p) {
        std::cerr << "\t" << p.first << p.second << "\n";
        return 3;
    }
    std::cerr << "\nSome other error?\n";
    return 4;
}
