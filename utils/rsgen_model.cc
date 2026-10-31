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

#include "utils/rsgen_model.h"
#include "utils/timer.h"

// #define DEBUG_EVENTS

//
// Add commas to long integer x
//
inline void show_node_count(long x)
{
    if (x<0) {
        std::cerr << "-";
        x = -x;
    } else {
        std::cerr << " ";
    }
    long grouping = 1000000000000L;
    char filler = ' ';
    while (grouping) {
        if (filler != ' ') {
            std::cerr << ",";
        } else {
            std::cerr << " ";
        }
        const long d = x / grouping;
        if (d) {
            if (d < 100) {
                std::cerr << filler;
            }
            if (d < 10) {
                std::cerr << filler;
            }
            std::cerr << d;
            x %= grouping;
            filler = '0';
        } else {
            std::cerr << filler << filler << filler;
        }
        grouping /= 1000;
    }
}


//
// progress for BFS iterations
//

MEDDLY::forest* F_prog = nullptr;

void my_progress(unsigned iter, char st)
{
    static timer watch;
    if (' ' == st) {
        watch.note_time();
        std::cerr << "    Iter" << std::setw(5) << iter << ": ";
        return;
    }
    if (';' == st) {
        watch.note_time();
        std::cerr   << std::fixed
                    << std::setprecision(4)
                    << std::setw(12)
                    << watch.get_last_seconds()
                    << " sec";

        if (F_prog) {
            show_node_count(F_prog->getCurrentNumNodes());
            std::cerr << " nodes";

            unsigned long mem = F_prog->getCurrentMemoryUsed();
            unsigned units = 0;
            while (mem >= 1024) {
                mem /= 1024;
                ++units;
            }
            const char* uletters = " KMGTPEZYabcdefgh";
            std::cerr << std::setw(5) << mem << " " << uletters[units] << "bytes";
        }

        std::cerr << std::endl;
        return;
    }
    std::cerr << st << " ";
}

// ******************************************************************
// *                   rsgen_model::event methods                   *
// ******************************************************************

rsgen_model::event::event(const char* n)
{
    name = n;
}

rsgen_model::event::~event()
{
}

// ******************************************************************
// *               rsgen_model::simple_event  methods               *
// ******************************************************************

rsgen_model::simple_event::simple_event(const char* n,
        const std::vector <level_delta> &_upd)
        : rsgen_model::event(n), updates(_upd)
{
#ifdef DEBUG_EVENTS
    std::cerr << "Creating event " << get_name() << " [";
    for (unsigned u=0; u<updates.size(); u++) {
        if (u) std::cerr << ", ";
        std::cerr << "(" << updates[u].level << ":" << updates[u].delta << ")";
    }
    std::cerr << "]\n";
#endif
}

bool rsgen_model::simple_event::fire(const MEDDLY::minterm &curr,
        MEDDLY::minterm &next)
{
    for (unsigned i=next.getNumVars(); i; i--) {
        next.from(i) = curr.from(i);
    }
    for (unsigned u=0; u<updates.size(); u++) {
        const unsigned i = updates[u].level;
        if (updates[u].delta < 0) {
            if (curr.from(i) < -updates[u].delta) return false;
        }
        next.from(i) = curr.from(i) + updates[u].delta;
    }
    return true;
}

bool rsgen_model::simple_event::build_nsf(MEDDLY::forest* mtmxd,
        MEDDLY::dd_edge &R)
{
    using namespace MEDDLY;

    static rangeval* temp = nullptr;
    static unsigned maxBound = 1;

    const domain* d = mtmxd->getDomain();
    const unsigned nVars = d->getNumVariables();
    if (!temp) {
        // build a temporary array, of dimension largest variable
        for (unsigned i=1; i<=nVars; i++) {
            const unsigned xib = d->getVariableBound(i, false);
            if (xib > maxBound) maxBound = xib;
        }
        temp = new rangeval[maxBound];
    }

    forest* mxd = R.getForest();

    //
    // Build mask for this event.
    // Levels not in the updates vector: don't change
    // Levels in the updates vector: don't care
    //
    minterm mask(mtmxd);
    for (unsigned i=1; i<=nVars; i++) {
        mask.setVars(i, DONT_CARE, DONT_CHANGE);
    }
    for (unsigned u=0; u<updates.size(); u++) {
        mask.setVars(updates[u].level, DONT_CARE, DONT_CARE);
    }
    mask.buildFunction(false, R);

    //
    // Build a term for each updated variable,
    // AND it with the relations so far.
    //
    dd_edge xun(mtmxd);
    dd_edge xpr(mtmxd);
    dd_edge term(mxd);
    for (unsigned u=0; u<updates.size(); u++) {
        const unsigned i = updates[u].level;
        //
        // Build xi'
        //
        mtmxd->createEdgeForVar(int(i), true, xpr);

        //
        // Build xi + delta
        //
        for (unsigned v=0; v<maxBound; v++) {
            temp[v].setInteger(int(v) + updates[u].delta);
        }
        mtmxd->createEdgeForVar(int(i), false, xun);

        //
        // update for this variable is: xi' == xi + delta
        //
        apply(EQUAL, xpr, xun, term);

        //
        // AND with relation
        //
        R *= term;
    }

    return true;
}

// ******************************************************************
// *                      rsgen_model  methods                      *
// ******************************************************************

rsgen_model::rsgen_model()
{
    D = nullptr;
    Fset = nullptr;
    Fmxd = nullptr;
    Fmtmxd = nullptr;

    // default options

    method = '1';
    approx = false;
    verb = 0;
}

rsgen_model::~rsgen_model()
{
    // MEDDLY::forest::destroy(Fmtmxd);
    // MEDDLY::forest::destroy(Fmxd);
    // MEDDLY::forest::destroy(Fset);
    // MEDDLY::domain::destroy(D);
}

bool rsgen_model::process_arg(unsigned &i, const char** argv)
{
    if (nullptr == argv[i]) {
        return false;
    }
    if ('-' != argv[i][0]) {
        return false;
    }

    if ('v' == argv[i][1]) {
        // verbosity
        if ((argv[i][2] >= '0') && (argv[i][2] <= '9')) {
            verb = argv[i][2] - '0';
            i++;
            return true;
        }
        return false;
    }

    if (0==strcmp(argv[i], "-approx")) {
        approx = true;
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-exact")) {
        approx = false;
        i++;
        return true;
    }

    if (0==strcmp(argv[i], "-front")) {
        method = 'f';
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-trad")) {
        method = 't';
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-dfs")) {
        method = '1';
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-sat1")) {
        method = '1';
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-sat2")) {
        method = '2';
        i++;
        return true;
    }


    return false;
}

void rsgen_model::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options]\n\n";
    }

    out << "\t-approx: approximate cardinalities";
    if (approx) {
        out << " (default)";
    }
    out << "\n";
    out << "\t-exact : exact cardinalities";
    if (!approx) {
        out << " (default)";
    }
    out << "\n";
    out << "\n";
    out << "\t-front : use traditional iterations, with frontier\n";
    out << "\t-trad  : use traditional iterations, without frontier\n";
    out << "\t-dfs   : use default saturation (currently, -sat1)\n";
    out << "\t-sat1  : Saturation v1, new implementation\n";
    out << "\t-sat2  : Saturation v2\n";
    out << "\n";
    out << "\t-vL    : Set verbosity (0..9) for traditional (default 0).\n";
    out << "\n";
}

void rsgen_model::generate()
{
    using namespace MEDDLY;

    //
    // Policies. Might change from command line switches.
    //
    policies pmdd(false), pmxd(true);

    // TBD: policy changes?

    //
    // build forests, even if we don't need them
    //

    // TBD: distance vs reachable states
    std::cout << "Building forests\n";

    Fset = forest::create(D, SET,
                range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL,
                pmdd
    );
    Fmxd = forest::create(D, RELATION,
                range_type::BOOLEAN,
                edge_labeling::MULTI_TERMINAL,
                pmxd
    );
    Fmtmxd = forest::create(D, RELATION,
                range_type::INTEGER,
                edge_labeling::MULTI_TERMINAL,
                pmxd
    );

    //
    // build initial states
    //
    std::cout << "Building initial\n";
    minterm initial(Fset);
    build_initial(initial);
    initial.setValue(true); // unless distance...

    dd_edge init_state(Fset);
    initial.buildFunction(false, init_state);

    //
    // build transition relation
    //
    dd_edge nsf(Fmxd);
    timer watch;
    watch.note_time();

    // tbd - different things based on options
    std::cout << "Building monolithic\n";
    build_monolithic(nsf);

    forest::destroy(Fmtmxd);
    watch.note_time();
    std::cout   << "Next-state function construction took "
                << watch.get_last_seconds() << " seconds\n";
    std::cout << "MxD for relation requires "
              << nsf.getNodeCount() << " nodes" << std::endl;
    printStats("MxD", Fmxd);

    //
    // Set progress updates according to verbosity
    //
    F_prog = (verb>1) ? Fset : nullptr;

    //
    // build reachset / distance
    //
    dd_edge reachable(Fset);
    watch.note_time();
    binary_operation* rsgen = nullptr;
    switch (method) {
        case 't':
            std::cout << "Building reachability set using traditional algorithm (no frontier)" << std::endl;
            rsgen = build(REACHABLE_TRAD_NOFS(true), Fset, Fmxd, Fset);
            if (verb) rsgen->setProgressNotifier(my_progress);
            rsgen->compute(init_state, nsf, reachable);
            break;

        case 'f':
            std::cout << "Building reachability set using traditional algorithm (with frontier)" << std::endl;
            rsgen = build(REACHABLE_TRAD_FS(true), Fset, Fmxd, Fset);
            if (verb) rsgen->setProgressNotifier(my_progress);
            rsgen->compute(init_state, nsf, reachable);
            break;

        case 'd':
            std::cout << "Building reachability set using default saturation" << std::endl;
            apply(REACHABLE_SATUR(true), init_state, nsf, reachable);
            break;

        case '1':
            std::cout << "Building reachability set using saturation v1, monolithic relation" << std::endl;
            apply(REACHABLE_SATUR(true, 1), init_state, nsf, reachable);
            break;

        case '2':
            std::cout << "Building reachability set using saturation v2, monolithic relation" << std::endl;
            apply(REACHABLE_SATUR(true, 2), init_state, nsf, reachable);
            break;

            /*
            case 'e':
                cout << "Building reachability set using explicit search\n";
                cout << "Using batch size: " << batchsize << endl;
                explicitReachset(kanban, 16, mdd, init_state, reachable, batchsize);
                break;

            case 'k':
            case 's':
                cout << "Building reachability set using saturation, relation";
                if ('k'==method)  cout << " by levels" << endl;
                else              cout << " by events" << endl;
                sat = SATURATION_FORWARD(mdd, ensf, mdd);
                if (0==sat) {
                    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
                }
                sat->compute(init_state, reachable);
                break;
                */

            default:
                std::cout << "Error - unknown method '" << method << "'\n";
                exit(2);
        };
        watch.note_time();
        std::cout << "Done\n";
        std::cout   << "Reachability set construction took "
                    << watch.get_last_seconds() << " seconds" << std::endl;
        std::cout   << "Reachability set uses "
                    << reachable.getNodeCount() << " nodes" << std::endl;

        //
        // Write to terminal
        //
        // if (show_dd) {
            // cout << "Reachable states MDD:\n";
            // reachable.showGraph(meddlyout);
        // }

        //
        // Write to file
        //
        // compressed_output writer(outfile);
        // if (writer) {
            // cout << "Writing to " << outfile << "\n";
            // d->write(writer.outstream());
            // mdd_writer mywriter(writer.outstream(), mdd);
            // mywriter.writeRootEdge(reachable);
            // mywriter.finish();
        // }

        printStats("MDD", Fset);



}

void rsgen_model::build_monolithic(MEDDLY::dd_edge &R)
{
    Fmxd->createConstant(false, R);
    MEDDLY::dd_edge Re(Fmxd);
    for (unsigned e=0; e<events.size(); e++) {
        if (events[e]) {
            events[e]->build_nsf(Fmtmxd, Re);
            R += Re;
        }
    }
}


void rsgen_model::printStats(const char* who, const MEDDLY::forest* F)
{
    using namespace MEDDLY;
    if (F) {
        ostream_output meddlyout(std::cout);
        meddlyout << who << " stats:\n";
        F->reportStats(meddlyout, "\t",
            HUMAN_READABLE_MEMORY | BASIC_STATS
            | EXTRA_STATS
            | STORAGE_STATS
            // | STORAGE_DETAILED
            // | HOLE_MANAGER_STATS
            // | HOLE_MANAGER_DETAILED
        );
    }
}

