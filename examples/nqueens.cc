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

/*
    Builds the set of solutions to the N-queens problem for
    user-specified N.

    In other words, finds all possible ways to put N queens onto
    an NxN chessboard so that no queen is attacking any other queen.
    Clearly, we cannot have two queens in any given row, therefore
    every such solution must have exactly one queen in each row.
    We use an MDD with N variables, with variable x_i indicating
    the column position (from 1 to N) for the queen in row i.
*/

#include <cstdio>
#include <cassert>

#include "../src/meddly.h"
#include "../timing/timer.h"
#include "../src/log_simple.h"

// #define SHOW_ALL_SOLUTIONS

using namespace MEDDLY;

int N;
long* scratch;
const char* lfile;
bool build_pdf;
bool pessimistic;
bool mark_sweep;
unsigned ct_max;

bool use_folding;


void intersect_acc(dd_edge** A, int L)
{
    if (0==A[0]) {
        // find first non-zero and swap
        for (int i=1; i<L; i++) {
            if (A[i]) {
                A[0] = A[i];
                A[i] = 0;
                break;
            }
        }
        if (0==A[0]) return;
    }
    fprintf(stderr, "\t");
    for (int i=1; i<L; i++) {
        fprintf(stderr, "%d ", L-i);
        if (A[i]) {
            apply(INTERSECTION, *A[0], *A[i], *A[0]);
            delete A[i];
            A[i] = 0;
            // operation::removeStalesFromMonolithic();
        }
    }
    fprintf(stderr, "\n");
}


void intersect_fold(dd_edge** A, int L)
{
    while (L>1) {
        fprintf(stderr, "\t%2d terms to combine ", L);
        // combine adjacent pairs
        for (int i=0; i<L; i+=2) {
            if (A[i] && A[i+1]) {
                apply(INTERSECTION, *A[i], *A[i+1], *A[i]);
                delete A[i+1];
                A[i+1] = 0;
                fprintf(stderr, ".");
            }
        } // for i
        fprintf(stderr, "\n");
        // compact
        int p = 0;
        for (int i=0; i<L; i++) {
            if (A[i]) {
                if (i>p) {
                    A[p] = A[i];
                    A[i] = 0;
                }
                p++;
            }
        } // for i
        L = p;
    } // while
}


inline void intersect(dd_edge** A, int L)
{
    if (use_folding)
        intersect_fold(A, L);
    else
        intersect_acc(A, L);
}


void createQueenNodes(forest* f, int q, dd_edge &col, dd_edge &cp, dd_edge &cm)
{
    assert(q>0);
    assert(q<=N);
    f->createEdgeForVar(q, false, col);
    for (int i=0; i<N; i++) {
        scratch[i] = i+q;
    }
    f->createEdgeForVar(q, false, scratch, cp);
    for (int i=0; i<N; i++) {
        scratch[i] = i-q;
    }
    f->createEdgeForVar(q, false, scratch, cm);
}

bool processArgs(int argc, const char** argv)
{
    lfile = 0;
    build_pdf = false;
    mark_sweep = true;
    pessimistic = true;
    ct_max = 0;
    bool setN = false;
    use_folding = false;
    for (int i=1; i<argc; i++) {
        if ('-' == argv[i][0]) {
            if (strcmp("-acc", argv[i])==0) {
                use_folding = false;
                continue;
            }
            if (strcmp("-fold", argv[i])==0) {
                use_folding = true;
                continue;
            }
            if (strcmp("-ms", argv[i])==0) {
                mark_sweep = true;
                continue;
            }
            if (strcmp("-opt", argv[i])==0) {
                mark_sweep = false;
                pessimistic = false;
                continue;
            }
            if (strcmp("-pess", argv[i])==0) {
                mark_sweep = false;
                pessimistic = true;
                continue;
            }
            if (strcmp("-ct", argv[i])==0) {
                long ctm = atol(argv[i+1]);
                i++;
                if (ctm > 0) {
                    ct_max = ctm;
                }
                continue;
            }
            if (strcmp("-l", argv[i])==0) {
                lfile = argv[i+1];
                i++;
                continue;
            }
            if (strcmp("-pdf", argv[i])==0) {
                build_pdf = true;
                continue;
            }
            return false;
        }
        if (setN) return false;
        N = atoi(argv[i]);
        setN = true;
    }

    if (!setN || N<1) return false;
    return true;
}

int usage(const char* who)
{
    /* Strip leading directory, if any: */
    const char* name = who;
    for (const char* ptr=who; *ptr; ptr++) {
        if ('/' == *ptr) name = ptr+1;
    }
    printf("Usage: %s [options] N\n\n", name);
    printf("\t       N:  board dimension\n");
    printf("\t    -acc:  Accumulate constraints in order (default)\n");
    printf("\t   -fold:  Accumulate constraints in pairs\n");
    printf("\t-ct size:  Maximum entries in compute table\n");
    printf("\t     -ms:  Mark and sweep node deletion (default)\n");
    printf("\t    -opt:  Reference counts and Optimistic node deletion\n");
    printf("\t   -pess:  Reference counts and Pessimistic node deletion\n");
    printf("\t-l lfile:  Write logging information to specified file\n");
    printf("\t    -pdf:  Write MDD of solutions to out.pdf\n\n");
    return 1;
}

int main(int argc, const char** argv)
{
    if (!processArgs(argc, argv)) return usage(argv[0]);

    try {
        initializer_list* INIT = defaultInitializerList(0);
        if (ct_max) {
            ct_initializer::setMaxSize(ct_max);
        }
        initialize(INIT);

        policies p(false);
        if (mark_sweep) {
            p.useReferenceCounts = false;
        } else if (pessimistic) {
            p.useReferenceCounts = true;
            p.setPessimistic();
        } else {
            p.useReferenceCounts = true;
            p.setOptimistic();
        }
        printf("Using %s\n", getLibraryInfo(0));

        timer watch;
        printf("%d-Queens solutions.\n", N);
        scratch = new long[N+1];

        watch.note_time();
        printf("Initializing domain and forest\n");
        const char* ndp = "unknown node deletion";
        if (!p.useReferenceCounts) {
            // TBD - can we use "never delete" with mark and sweep?
            ndp = "mark and sweep node deletion";
        } else {
            switch (p.deletion) {
                case policies::node_deletion::NEVER:
                    ndp = "never delete";
                    break;

                case policies::node_deletion::OPTIMISTIC:
                    ndp = "optimistic node deletion";
                    break;

                case policies::node_deletion::PESSIMISTIC:
                    ndp = "pessimistic node deletion";
                    break;
            }
        }
        printf("\tUsing %s policy\n", ndp);

        int* varsizes = new int[N+1];
        for (int i=0; i<N; i++) {
            varsizes[i] = N;
        }
        domain* d = domain::createBottomUp(varsizes, unsigned(N));
        assert(d);
        delete[] varsizes;
        forest* fi = forest::create(d, MEDDLY::SET, range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL, p);
        assert(fi);
        forest* fb = forest::create(d, MEDDLY::SET, range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL, p);
        assert(fb);

        FILE_output log;
        logger* LOG = 0;
        if (lfile) {
            FILE* fout = fopen(lfile, "w");
            if (!fout) {
                printf("Couldn't open %s for writing, no logging\n", lfile);
            } else {
                log.setFILE(fout);
                LOG = new simple_logger(log);
                LOG->recordNodeCounts();
                LOG->recordTimeStamps();
                char comment[80];
                snprintf(comment, 80, "Automatically generated by nqueens (N=%d)", N);
                LOG->addComment(comment);
                fb->setLogger(LOG, "forest name");
                LOG->newPhase(fb, "Initializing");
            }
        }

        printf("Building nodes for queen column and diagonals\n");

        dd_edge** col = new dd_edge*[N];
        dd_edge** dgp = new dd_edge*[N];
        dd_edge** dgm = new dd_edge*[N];
        dd_edge** constr = new dd_edge*[N+1];


        for (int i=0; i<N; i++) {
            col[i] = new dd_edge(fi);
            dgp[i] = new dd_edge(fi);
            dgm[i] = new dd_edge(fi);
            createQueenNodes(fi, i+1, *col[i], *dgp[i], *dgm[i]);
            constr[i] = new dd_edge(fb);
            fb->createConstant(true, *constr[i]);
        }
        constr[N] = 0;

        for (int i=0; i<N-1; i++) {
            char buffer[80];
            snprintf(buffer, 80, "Building queen %2d constraints\n", i+1);
            printf("%s\n", buffer);
            if (LOG) LOG->newPhase(fb, buffer);
            for (int j=N-1; j>i; j--) {
                dd_edge uniq_col(fb);
                apply(NOT_EQUAL, *col[i], *col[j], uniq_col);
                dd_edge uniq_dgp(fb);
                apply(NOT_EQUAL, *dgp[i], *dgp[j], uniq_dgp);
                dd_edge uniq_dgm(fb);
                apply(NOT_EQUAL, *dgm[i], *dgm[j], uniq_dgm);
                // build overall "not attacking each other" set...
                apply(INTERSECTION, uniq_col, uniq_dgp, uniq_col);
                apply(INTERSECTION, uniq_col, uniq_dgm, uniq_col);
                int k = uniq_col.getLevel()-1;
                if (k<0) k=0;
                assert(k<N);
                apply(INTERSECTION, *constr[k], uniq_col, *constr[k]);
            } // for j
        } // for i
        forest::destroy(fi);

        printf("Building solutions\n");
        fflush(stdout);
        if (LOG) LOG->newPhase(fb, "Building solutions");
        intersect(constr, N);
        assert(constr[0]);
        dd_edge* solutions = constr[0];
        constr[0] = 0;
        watch.note_time();

        printf("Elapsed time: %lf seconds\n", watch.get_last_seconds());

        printf("Cleanup\n");
        if (LOG) LOG->newPhase(fb, "Cardinality");
        for (int i=0; i<N; i++) {
            delete col[i];
            delete dgp[i];
            delete dgm[i];
            delete constr[i];
        }
        delete[] col;
        delete[] dgp;
        delete[] dgm;
        delete[] constr;

        printf("Done.\n\n");

        printf("Set of solutions requires %lu nodes\n", solutions->getNodeCount());
        printf("Forest stats:\n");
        FILE_output myout(stdout);
        fb->reportStats(myout, "\t",
                HUMAN_READABLE_MEMORY  |
                BASIC_STATS | EXTRA_STATS |
                STORAGE_STATS | STORAGE_DETAILED |
                HOLE_MANAGER_STATS | HOLE_MANAGER_DETAILED
                );

        long c;
        apply(CARDINALITY, *solutions, c);
        printf("\nThere are %ld solutions to the %d-queens problem\n\n", c, N);

        FILE_output out(stdout);
        dd_edge::iterator first = solutions->begin();
#ifdef SHOW_ALL_SOLUTIONS
        c = 0;
        for (; first; ++first) {
            c++;
            out << "Solution ";
            out.put(c, 6);
            out << ": ";
            (*first).show(out);
            out << "\n";
        }
#else
        // show one of the solutions
        if (first) {
            out << "One of the solutions:\n";
            for (int i=1; i<=N; i++) {
                out << "\tQueen for row " << i << " in column " << (*first).from(i)+1 << "\n";
            }
        }
#endif
        if (build_pdf) {
            solutions->setLabel("queens");
            dot_maker dm(fb, "out");
            dm.addRootEdge(*solutions);
            dm.doneGraph();
            dm.runDot("pdf");
        }
        delete solutions;
        fb->validateCacheCounts();
        compute_table::showAll(myout, 3);
        if (LOG) {
            LOG->newPhase(fb, "Cleanup");
            domain::destroy(d);
        }
        cleanup();
        delete LOG;
        return 0;
    }
    catch (MEDDLY::error e) {
        fprintf(stderr, "\nCaught meddly error '%s'\n", e.getName());
        fprintf(stderr, "    thrown in %s line %u\n", e.getFile(), e.getLine());
        return 1;
    }
    catch (const char* e) {
        fprintf(stderr, "\nCaught our own error: '%s'\n", e);
        return 2;
    }
    fprintf(stderr, "\nSome other error?\n");
    return 4;
}
