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

// #include <cstdlib>
#include <iostream>
// #include <iomanip>
// #include <string.h>

#include "kan_base.h"

// ============================================

long expected[] = {
    1, 160, 4600, 58400, 454475, 2546432, 11261376,
    41644800, 133865325, 384392800, 1005927208
};

const int nstart = 1;
const int nstop = 10;

// ============================================

class kanban : public kanban_base {
    public:
        kanban();

        inline void setN(int n) {
            N = n;
        }
        inline void fill_initial(MEDDLY::dd_edge &init_state)
        {
            MEDDLY::minterm init(mdd);
            build_initial(init);
            init.buildFunction(false, init_state);
        }

        long writeReachset(FILE* s);
        bool generateAndRead(FILE* s);
        bool readAndGenerate(FILE* s);

    protected:
        void make_forests();
        void done_forests();

    private:
        MEDDLY::forest* mdd;
        MEDDLY::forest* mxd;
        MEDDLY::forest* evmdd;
};

// ============================================

kanban::kanban()
{
    disable_output_options();
    mdd = nullptr;
    mxd = nullptr;
    evmdd = nullptr;
}

// ============================================

/*
    Generate transition relation, reachability set for given N,
    and write them to a file.
*/
long kanban::writeReachset(FILE* s)
{
    using namespace MEDDLY;

    make_forests();

    // Build initial state
    dd_edge init_state(mdd);
    fill_initial(init_state);

#ifdef PROGRESS
    fputc('i', stdout);
    fflush(stdout);
#endif

    // Build next-state function
    dd_edge nsf(mxd);
    build_monolithic(nsf);

#ifdef PROGRESS
    fputc('n', stdout);
    fflush(stdout);
#endif

    // Build reachable states
    dd_edge reachable(mdd);
    apply(REACHABLE_TRAD_NOFS(true), init_state, nsf, reachable);

#ifdef PROGRESS
    fputc('r', stdout);
    fflush(stdout);
#endif

    // Build index set for reachable states
    dd_edge reach_index(evmdd);
    apply(CONVERT_TO_INDEX_SET, reachable, reach_index);

#ifdef PROGRESS
    fputc('x', stdout);
    fflush(stdout);
#endif

    long c = -1;
    apply(CARDINALITY, reachable, c);

#ifdef PROGRESS
    fputc('c', stdout);
    fflush(stdout);
#endif

    FILE_output mys(s);
    get_domain()->write(mys);

#ifdef WRITE_MXD
    mdd_writer wmxd(mys, mxd);
    wmxd.writeRootEdge(nsf);
    wmxd.finish();
#endif

#ifdef WRITE_MDD
    // Write initial state + reachable states
    mdd_writer wmdd(mys, mdd);
    wmdd.writeRootEdge(init_state);
    wmdd.writeRootEdge(reachable);
    wmdd.finish();
#endif

#ifdef WRITE_EVMDD
    mdd_writer wevmdd(mys, evmdd);
    wevmdd.writeRootEdge(reach_index);
    wevmdd.finish();
#endif

    //
    // clean up forests
    //
    done_forests();
    return c;
}

// ============================================

/*
    Generate transition relation, reachability set for given N,
    then read from a file and check for equality.
*/
bool kanban::generateAndRead(FILE* s)
{
    using namespace MEDDLY;

    make_forests();

    // Build initial state
    dd_edge init_state(mdd);
    fill_initial(init_state);

    // Build next-state function
    dd_edge nsf(mxd);
    build_monolithic(nsf);

    // Build reachable states
    dd_edge reachable(mdd);
    apply(REACHABLE_TRAD_NOFS(true), init_state, nsf, reachable);

    // Build index set for reachable states
    dd_edge reach_index(evmdd);
    apply(CONVERT_TO_INDEX_SET, reachable, reach_index);

    // Now, read from the file and verify
    FILE_input mys(s);
    get_domain()->verify(mys);

    dd_edge list[2];
#ifdef WRITE_MXD
    mdd_reader R_mxd(mys, mxd);
    dd_edge r_nsf;
    R_mxd.readRootEdge(r_nsf);
    if (r_nsf != nsf) {
        printf("Failed to generate and read MXD\n");
        return false;
    }
#endif

#ifdef WRITE_MDD
    mdd_reader R_mdd(mys, mdd);
    dd_edge r_init, r_reach;
    R_mdd.readRootEdge(r_init);
    R_mdd.readRootEdge(r_reach);
    if (r_init != init_state) {
        printf("Failed to generate and read initial state\n");
        return false;
    }
    if (r_reach != reachable) {
        printf("Failed to generate and read reachable states\n");
        return false;
    }
#endif

#ifdef WRITE_EVMDD
    mdd_reader R_evmdd(mys, evmdd);
    dd_edge r_index;
    R_evmdd.readRootEdge(r_index);
    if (r_index != reach_index) {
        printf("Failed to generate and read reachable state indexes\n");
        FILE_output myout(stdout);
        myout << "reach_index: " << reach_index << "\n";
        myout << "from file  : " << r_index << "\n";
#ifdef DEBUG_FILE
        reach_index.showGraph(myout);
        r_index.showGraph(myout);
#endif
        return false;
    }
#endif

    done_forests();
    return true;
}

// ============================================

/*
    Read the transition relation, reachability set from a file,
    then generate them for a given N and check for equality.
*/
bool kanban::readAndGenerate(FILE* s)
{
    using namespace MEDDLY;

    // Build (empty) forests
    make_forests();

    dd_edge mxdfile;
    dd_edge mddfile0, mddfile1;
    dd_edge indexfile;

    // Read from the file
    FILE_input mys(s);
    get_domain()->verify(mys);

#ifdef WRITE_MXD
    mdd_reader R_mxd(mys, mxd);
    R_mxd.readRootEdge(mxdfile);
#endif
#ifdef WRITE_MDD
    mdd_reader R_mdd(mys, mdd);
    R_mdd.readRootEdge(mddfile0);
    R_mdd.readRootEdge(mddfile1);
#endif
#ifdef WRITE_EVMDD
    mdd_reader R_evmdd(mys, evmdd);
    R_evmdd.readRootEdge(indexfile);
#endif

    // Build initial state
    dd_edge init_state(mdd);
    fill_initial(init_state);

    // Build next-state function
    dd_edge nsf(mxd);
    build_monolithic(nsf);

    // Build reachable states
    dd_edge reachable(mdd);
    apply(REACHABLE_TRAD_NOFS(true), init_state, nsf, reachable);

    // Build index set for reachable states
    dd_edge reach_index(evmdd);
    apply(CONVERT_TO_INDEX_SET, reachable, reach_index);


    // Verify
#ifdef WRITE_MXD
    if (mxdfile != nsf) {
        printf("Failed to read and generate MXD\n");
        return false;
    }
#endif
#ifdef WRITE_MDD
    if (mddfile0 != init_state) {
        printf("Failed to generate and read initial state\n");
        return false;
    }
    if (mddfile1 != reachable) {
        printf("Failed to generate and read reachable states\n");
        return false;
    }
#endif
#ifdef WRITE_EVMDD
    if (indexfile != reach_index) {
        printf("Failed to generate and read reachable state indexes\n");
        return false;
    }
#endif

    done_forests();
    return true;
}

// ============================================

void kanban::make_forests()
{
    using namespace MEDDLY;
    mdd = forest::create(get_domain(), SET, range_type::BOOLEAN,
                            edge_labeling::MULTI_TERMINAL);
    mxd = forest::create(get_domain(), RELATION, range_type::BOOLEAN,
                            edge_labeling::MULTI_TERMINAL);
    evmdd = forest::create(get_domain(), SET, range_type::INTEGER,
                            edge_labeling::INDEX_SET);
}

// ============================================

void kanban::done_forests()
{
    using namespace MEDDLY;
    forest::destroy(evmdd);
    forest::destroy(mxd);
    forest::destroy(mdd);
}

// ============================================

/*
    Main program, of course
*/
int main()
{
#ifdef DEBUG_FILE
    const int N = 2;
#else
    const int N = 11;
#endif
    MEDDLY::initialize();
    kanban K;

    char filename[20];
    strcpy(filename, "kan_io.data.XXXXXX");
    mkstemp(filename);

    MEDDLY::FILE_output out(stdout);

    try {
        for (int n=1; n<N; n++) {
            out.indent_more();
            out << "N=" << n << "\n";
            out << "Building model...\n";
            K.setN(n);
            K.build_model(out);

            out << "Saving dds to file " << filename << "... ";
            FILE* s = fopen(filename, "w");
            if (!s) {
                out << "    Couldn't open for writing\n";
                return 2;
            }

            long c = K.writeReachset(s);
            out << "\n" << c << " states\n";
            if (c !=expected[n]) {
                out << "    wrong number of states!\n";
                throw 1;
            }
            fclose(s);

            // read back the file into already built forest

            out << "Generate and read... ";
            s = fopen(filename, "r");
            if (!s) {
                out << "    Couldn't open for reading\n";
                throw 2;
            }

            if (!K.generateAndRead(s)) throw 3;
            out << "    verified\n";
            fclose(s);

            // read back the file and then build up forest

            out << "Read and generate... ";
            s = fopen(filename, "r");
            if (!s) {
                out << "    Couldn't open for reading\n";
                throw 2;
            }
            if (!K.readAndGenerate(s)) throw 4;
            out << "    verified\n";
            fclose(s);

            // cleanup before next iteration

            K.destroy_model();
            remove(filename);

            out.indent_less();
            out << "\n";
        } // for n
        out << "Done\n";
        MEDDLY::cleanup();
        return 0;
    }
    catch (int e) {
        // cleanup
#ifndef DEBUG_FILE
        remove(filename);
#endif
        MEDDLY::cleanup();
        printf("\nError %d\n", e);
        return e;
    }
    catch (MEDDLY::error e) {
        // cleanup
#ifndef DEBUG_FILE
        remove(filename);
#endif
        MEDDLY::cleanup();
        printf("\nError: %s, thrown at %s line %u\n", e.getName(),
                e.getFile(), e.getLine());
        return 1;
    }
    catch (...) {
        // cleanup
#ifndef DEBUG_FILE
        remove(filename);
#endif
        MEDDLY::cleanup();
        printf("\nFailed (caught exception)\n");
        return 1;
    }
}

