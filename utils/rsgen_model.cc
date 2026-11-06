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
#include "utils/zstream.h"

#include <iomanip>

// #define DEBUG_UPDATE_SORT
// #define DEBUG_UPDATE_STRING
// #define DEBUG_EVENTS
// #define DEBUG_EVENT_BUILD
// #define DEBUG_GENERATE
// #define SHOW_MXD_FOREST
// #define WRITE_MXD
// #define WRITE_EVENTS

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
        std::cerr << "    Iter" << std::setw(6) << iter << ": ";
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
        const std::vector <level_axb> &_upd)
        : rsgen_model::event(n), updates(_upd)
{
    //
    // put update vector in order. using something like insertion sort
    //
    for (unsigned i=1; i<updates.size(); i++) {
        for (unsigned j=i; j; ) {
            --j;
            if (updates[j].level < updates[j+1].level) break;
            updates[j].swap(updates[j+1]);
        }
    }
#ifdef DEBUG_UPDATE_SORT
    show();
#endif
}

rsgen_model::simple_event::simple_event(const char* n, const char* pmstr)
        : rsgen_model::event(n)
{
    if (pmstr) {
        MEDDLY_DCASSERT('X' == pmstr[0]);
        unsigned u = 0;
        for (unsigned i=1; pmstr[i]; i++) {
            if (('+' == pmstr[i]) || ('-' == pmstr[i])) {
                u++;
            }
        }
        if (u) {
            updates.resize(u);
            u = 0;
            for (unsigned i=1; pmstr[i]; i++) {
                if ('.' == pmstr[i]) continue;
                const int delta = ('+' == pmstr[i]) ? 1 : -1;
#ifdef DEVELOPMENT_CODE
                updates.at(u).set(i, 1, delta);
#else
                updates[u].set(i, 1, delta);
#endif
                u++;
            }
        }
    }
#ifdef DEBUG_UPDATE_STRING
    show();
#endif
}

bool rsgen_model::simple_event::fire(const MEDDLY::minterm &curr,
        MEDDLY::minterm &next)
{
    for (unsigned i=next.getNumVars(); i; i--) {
        next.setVar(i, curr.getVar(i));
    }
    for (unsigned u=0; u<updates.size(); u++) {
        const unsigned i = updates[u].level;
        next.setVar(i,
            updates[u].scale * curr.getVar(i) + updates[u].delta
        );
        if (next.getVar(i) < 0) return false;
    }
    return true;
}

#ifdef OLD_SIMPLE_MODEL

bool rsgen_model::simple_event::build_nsf(MEDDLY::forest* mtmxd,
        MEDDLY::dd_edge &R)
{
    using namespace MEDDLY;

#ifdef DEBUG_EVENTS
    ostream_output out(std::cerr);
    out << "Building nsf for event " << get_name() << " [";
    for (unsigned u=0; u<updates.size(); u++) {
        if (u) out << ", ";
        out << "(" << updates[u].level << ":" << updates[u].delta << ")";
    }
    out << "]\n";
#endif

    forest* mxd = R.getForest();

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

#ifdef DEBUG_EVENT_BUILD
    out << "Event mask:\n";
    R.showGraph(out);
#endif

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
            temp[v].setInteger(updates[u].scale * int(v) + updates[u].delta);
        }
        mtmxd->createEdgeForVar(int(i), false, temp, xun);

        //
        // update for this variable is: xi' == xi + delta
        //
        apply(EQUAL, xpr, xun, term);

#ifdef DEBUG_EVENT_BUILD
        out << "Term at level " << i << ":\n";
        term.showGraph(out);
#endif

        //
        // AND with relation
        //
        R *= term;
    }
    term.set(0);

#ifdef DEBUG_EVENTS
    out << "Overall event:\n";
    R.show(out);
#endif

    return true;
}


#else


bool rsgen_model::simple_event::build_nsf(MEDDLY::dd_edge &R)
{
    using namespace MEDDLY;

#ifdef DEBUG_EVENTS
    ostream_output out(std::cerr);
    out << "Building nsf for event " << get_name() << " [";
    for (unsigned u=0; u<updates.size(); u++) {
        if (u) out << ", ";
        out << "(" << updates[u].level << ":" << updates[u].delta << ")";
    }
    out << "]\n";
#endif

    forest* mxd = R.getForest();

    node_handle root = mxd->handleForValue(true);
    unpacked_node* Uun;
    unpacked_node* Upr;
    edge_value ev;
    const bool IR = mxd->isIdentityReduced();
    int K = 0;
    for (unsigned u=0; u<updates.size(); u++) {
        int L = int(updates[u].level);
        root = mxd->makeIdentitiesTo(root, K, L-1, ~0);
        Uun = unpacked_node::newWritable(mxd, L, FULL_ONLY);
        for (unsigned i=0; i<Uun->getSize(); i++) {
            const int j = updates[u].scale * int(i) + updates[u].delta;
            if (j < 0 || j >= int(Uun->getSize())) {
                Uun->setFull(i, 0);
            } else if (unsigned(j)==i && IR) {
                Uun->setFull(i, mxd->linkNode(root));
            } else {
                // build singleton
                Upr = unpacked_node::newWritable(mxd, -L, 1, SPARSE_ONLY);
                Upr->setSparse(0, j, mxd->linkNode(root));
                node_handle h;
                mxd->createReducedNode(Upr, ev, h);
                MEDDLY_DCASSERT(ev.isVoid());
                Uun->setFull(i, h);
            }
        }
        mxd->unlinkNode(root);
        mxd->createReducedNode(Uun, ev, root);
        MEDDLY_DCASSERT(ev.isVoid());
        K = L;
    }
    root = mxd->makeIdentitiesTo(root, K, mxd->getNumVariables(), ~0);
    R.set(root);

#ifdef DEBUG_EVENTS
    out << "Overall event:\n";
    R.show(out);
#endif

    return true;
}

#endif // OLD_SIMPLE_MODEL

void rsgen_model::simple_event::show() const
{
    std::cerr << "Event ";
    if (get_name()) std::cerr << get_name() << " ";
    std::cerr << ": [";
    for (unsigned u=0; u<updates.size(); u++) {
        if (u) std::cerr << ", ";
        std::cerr << "(" << updates[u].level << ": *" << updates[u].scale;
        if (updates[u].delta < 0) {
            std::cerr << " " << updates[u].delta << ")";
        } else {
            std::cerr << " +" << updates[u].delta << ")";
        }
    }
    std::cerr << "]\n";
}

// ******************************************************************
// *                   rsgen_model public methods                   *
// ******************************************************************

rsgen_model::rsgen_model()
{
    D = nullptr;
    Fset = nullptr;
    Fmxd = nullptr;
#ifdef OLD_SIMPLE_MODEL
    Fmtmxd = nullptr;
#endif

    // default options

    set_set_type('b');
    set_gen_method('1');
    set_garbage('o');
    forwd = true;
#ifdef HAVE_LIBGMP
    approx = false;
#else
    approx = true;
#endif
    verb = 0;
    expl_batch = 1024;

    allow_output_files = true;
    file_rel = nullptr;
    file_set = nullptr;
    check_dir = nullptr;
}

rsgen_model::~rsgen_model()
{
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

#ifdef HAVE_LIBGMP
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
#endif

    //
    // What are we generating
    //
    if (0==strcmp(argv[i], "-evdist")) {
        set_set_type('+');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-mtdist")) {
        set_set_type('i');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-reach")) {
        set_set_type('b');
        i++;
        return true;
    }

    if (0==strcmp(argv[i], "-fwd")) {
        forwd = true;
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-bwd")) {
        forwd = false;
        i++;
        return true;
    }

    //
    // Generation algorithms
    //

    if (0==strcmp(argv[i], "-exp")) {
        set_gen_method('x');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-eb")) {
        i++;
        if (argv[i]) {
            long x = atol(argv[i]);
            if (x>0) {
                expl_batch = unsigned(x);
            }
            i++;
        }
        return true;
    }
    if (0==strcmp(argv[i], "-front")) {
        set_gen_method('f');
        i++;
        check_dir = nullptr;
        return true;
    }
    if (0==strcmp(argv[i], "-trad")) {
        set_gen_method('t');
        i++;
        check_dir = nullptr;
        return true;
    }
    if (0==strcmp(argv[i], "-tradc")) {
        set_gen_method('t');
        i++;
        check_dir = argv[i];
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-dfs")) {
        set_gen_method('d');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-sat1")) {
        set_gen_method('1');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-sat2")) {
        set_gen_method('2');
        i++;
        return true;
    }

    //
    // Garbage collection
    //
    if (0==strcmp(argv[i], "-ms")) {
        set_garbage('m');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-opt")) {
        set_garbage('o');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-pess")) {
        set_garbage('p');
        i++;
        return true;
    }

    //
    // File output
    //
    if (allow_output_files) {
        if (0==strcmp(argv[i], "-wr")) {
            i++;
            file_rel = argv[i];
            i++;
            return true;
        }
        if (0==strcmp(argv[i], "-ws")) {
            i++;
            file_set = argv[i];
            i++;
            return true;
        }
    }

    return false;
}

void rsgen_model::show_usage(std::ostream &out, const char* arg0)
{
    if (arg0) {
        out << "\nUsage: " << basename(arg0) << " [options]\n";
    }
    out << "\n";
    out << "    GENERATION TYPES\n";
    out << "\n";
    out << "\t-evdist    : build the distance function using EV+MDDs.\n";
    out << "\t-mtdist    : build the distance function.using MTMDDs.\n";
    out << "\t-reach     : build the reachability set (default).\n";
    out << "\n";
    out << "\t-fwd       : use forward generation (default).\n";
    out << "\t-bwd       : use backward generation.\n";

    out << "\n";
    out << "    GENERATION METHODS\n";
    out << "\n";
    out << "\t-exp       : use explicit generation (forward only)\n";
    out << "\t-front     : use traditional iterations, with frontier\n";
    out << "\t-trad      : use traditional iterations, without frontier\n";
    out << "\t-tradc dir : -trad with checkpoints saved in directory dir\n";
    out << "\t-dfs       : use default saturation (currently, -sat1)\n";
    out << "\t-sat1      : Saturation v1, new implementation\n";
    out << "\t-sat2      : Saturation v2\n";
    out << "\n";
    out << "\t-eb N      : Set explicit batch size (default "
        << expl_batch << ")\n";
    out << "\t-vL        : Set verbosity (0..9) for traditional (default 0).\n";
    out << "\n";
    out << "    GARBAGE COLLECTION\n";
    out << "\n";
    out << "\t-ms       : mark and sweep\n";
    out << "\t-opt      : optimistic reference counts\n";
    out << "\t-pess     : pessimistic reference counts\n";
    out << "\n";
    if (allow_output_files) {
        out << "    FILE OUTPUT\n";
        out << "\n";
        out << "\t-wr file  : write the transition relation to file\n";
        out << "\t-ws file  : write the generated set to file\n";
        out << "\n";
        out << "\t            Output filenames ending in .gz or .xz will be compressed.\n";
        out << "\n";
    }
#ifdef HAVE_LIBGMP
    out << "    CARDINALITY\n";
    out << "\n";
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
#endif
}

void rsgen_model::generate(MEDDLY::output &out, bool mxdstats, bool mddstats)
{
    using namespace MEDDLY;

    //
    // Policies. Might change from command line switches.
    //
    out << "Building forests\n";
    policies pmdd(false), pmxd(true);

    switch (garbage) {
        case 'm':
            out << "    Using mark and sweep\n";
            pmdd.useReferenceCounts = false;
            pmxd.useReferenceCounts = false;
            break;

        case 'p':
            out << "    Using pessimistic node deletion\n";
            pmdd.useReferenceCounts = true;
            pmxd.useReferenceCounts = true;
            pmdd.setPessimistic();
            pmxd.setPessimistic();
            break;

        default:
            out << "    Using optimistic node deletion\n";
            pmdd.useReferenceCounts = true;
            pmxd.useReferenceCounts = true;
            pmdd.setOptimistic();
            pmxd.setOptimistic();
            break;
    }

    //
    // build forests, even if we don't need them
    //

    switch (set_type) {
        case 'b':
            Fset =  forest::create(D, SET,
                    range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL,
                    pmdd
                );
            break;

        case 'i':
            Fset =  forest::create(D, SET,
                    range_type::INTEGER,
                    edge_labeling::MULTI_TERMINAL,
                    pmdd
                );
            break;

        case '+':
            Fset =  forest::create(D, SET,
                    range_type::INTEGER,
                    edge_labeling::EVPLUS,
                    pmdd
                );
            break;

        default:
            throw "internal: unknown set_type";
    }

    if (' ' != reltype) {
        Fmxd = forest::create(D, RELATION,
                    range_type::BOOLEAN,
                    edge_labeling::MULTI_TERMINAL,
                    pmxd
        );
#ifdef OLD_SIMPLE_MODEL
        Fmtmxd = forest::create(D, RELATION,
                    range_type::INTEGER,
                    edge_labeling::MULTI_TERMINAL,
                    pmxd
        );
#endif
    }

    //
    // build initial states
    //
    out << "Building initial\n";
    minterm initial(Fset);
    build_initial(initial);
    dd_edge init_state(Fset);
    init_state.setLabel("initial");
    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);

    switch (set_type) {
        case 'i':
            initial.setValue(0);
            initial.buildFunction(-1, init_state);
            break;

        case '+':
            initial.setValue(0);
            initial.buildFunction(infty, init_state);
            break;

        default:
            initial.setValue(true);
            initial.buildFunction(false, init_state);
    }

    //
    // build transition relation
    //
    timer watch;
    dd_edge nsf(Fmxd);
    watch.note_time();
    switch (reltype) {
        case 'm':   // monolithic
            out << "Building monolithic relation\n";
            build_monolithic(nsf);
#ifdef OLD_SIMPLE_MODEL
            forest::destroy(Fmtmxd);
#endif
            break;

        case ' ':   // no relation
            break;

        default:
            out << "Unknown relation type: '" << reltype << "'\n";
            throw "internal error";
    };
    watch.note_time();
    if (Fmxd) {
        out   << "Next-state function construction took "
              << watch.get_last_seconds() << " seconds\n";
        out << "MxD for relation requires "
            << nsf.getNodeCount() << " nodes\n";
        if (mxdstats) {
            printStats(out, Fmxd);
        }

        write_to_file("relation", nsf, file_rel);
    }

    //
    // Set progress updates according to verbosity
    //
    F_prog = (verb>1) ? Fset : nullptr;

    //
    // build reachset / distance
    //
    dd_edge reachable(Fset);
    reachable.setLabel((set_type != 'b') ? "distfunc" : "reachset");
    watch.note_time();
    binary_operation* rsgen = nullptr;
    out << "Building " << get_direction() << " " << get_set_name();
    switch (gen_method) {
        case 't':
            out << " using traditional algorithm (no frontier)\n";
            if (check_dir) {
                trad_with_checkpoints(nsf, init_state, reachable);
            } else {
                rsgen = build(REACHABLE_TRAD_NOFS(forwd), Fset, Fmxd, Fset);
                if (verb) rsgen->setProgressNotifier(my_progress);
                rsgen->compute(init_state, nsf, reachable);
            }
            break;

        case 'f':
            if (set_type != 'b') {
                throw "Cannot use traditional with frontier for distance";
            }
            out   << " using traditional algorithm (with frontier)\n";
            rsgen = build(REACHABLE_TRAD_FS(forwd), Fset, Fmxd, Fset);
            if (verb) rsgen->setProgressNotifier(my_progress);
            rsgen->compute(init_state, nsf, reachable);
            break;

        case 'd':
            out << " using default saturation\n";
            apply(REACHABLE_SATUR(forwd), init_state, nsf, reachable);
            break;

        case '1':
            out << " using saturation v1\n";
            apply(REACHABLE_SATUR(forwd, 1), init_state, nsf, reachable);
            break;

        case '2':
            out << " using saturation v2\n";
            apply(REACHABLE_SATUR(forwd, 2), init_state, nsf, reachable);
            break;

        case 'x':
            out << " using explicit generation\n";
            out << "    using batch size " << expl_batch << "\n";
            gen_explicit(init_state, reachable);
            break;

            /*
            case 'k':
            case 's':
                cout << "Building reachability set using saturation, relation";
                if ('k'==gen_method)  cout << " by levels" << endl;
                else              cout << " by events" << endl;
                sat = SATURATION_FORWARD(mdd, ensf, mdd);
                if (0==sat) {
                    throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
                }
                sat->compute(init_state, reachable);
                break;
                */

        default:
            out << "Error - unknown method '" << gen_method << "'\n";
            exit(2);
    };
    watch.note_time();
    out << "Done\n"
        << get_Set_name() << " construction took "
        << watch.get_last_seconds() << " seconds\n"
        << get_Set_name() << " uses "
        << reachable.getNodeCount() << " nodes\n";
    out.flush();

    after_generation(out, reachable);

    write_to_file(get_set_name(), reachable, file_set);

    if (mddstats) {
        printStats(out, Fset);
    }
}

void rsgen_model::build_monolithic(MEDDLY::dd_edge &R)
{
#ifdef WRITE_MXD
    MEDDLY::FILE_output myout(fopen("mxd.txt", "w"));
    D->write(myout);
    MEDDLY::mdd_writer mywriter(myout, Fmxd);
#endif
    Fmxd->createConstant(false, R);
    MEDDLY::dd_edge Re(Fmxd);
    for (unsigned e=0; e<events.size(); e++) {
        if (events[e]) {
#ifdef OLD_SIMPLE_MODEL
            events[e]->build_nsf(Fmtmxd, Re);
#else
            events[e]->build_nsf(Re);
#endif
            R += Re;
#ifdef WRITE_EVENTS
            mywriter.writeRootEdge(Re);
#endif
        }
    }
    Re.set(0);
#ifdef WRITE_MXD
    mywriter.writeRootEdge(R);
    mywriter.finish();
#endif
#ifdef SHOW_MXD_FOREST
    MEDDLY::ostream_output out(std::cerr);
    out << "Monolithic relation: " << R << "\n";
    out << "Complete forest:\n";
    Fmxd->dump(out, MEDDLY::SHOW_DETAILS);
#endif
    R.setLabel("nsf-monolithic");
}

void rsgen_model::after_generation(MEDDLY::output &out,
        const MEDDLY::dd_edge &e)
{
    if ('b' == get_set_type()) {
        show_cardinality(out, e);
    }
}

void rsgen_model::destroy_model()
{
#ifdef OLD_SIMPLE_MODEL
    MEDDLY::forest::destroy(Fmtmxd);
#endif
    MEDDLY::forest::destroy(Fmxd);
    MEDDLY::forest::destroy(Fset);
    MEDDLY::domain::destroy(D);
}

// ******************************************************************
// *                 rsgen_model  protected methods                 *
// ******************************************************************


void rsgen_model::gen_explicit(const MEDDLY::dd_edge &init,
        MEDDLY::dd_edge &RS)
{
    using namespace MEDDLY;

    if (!forwd) {
        throw "can't run explicit generation backwards";
    }

#ifdef DEBUG_GENERATE
    ostream_output out(std::cerr);
#endif

    // initialize batch memory
    rangeval infty(range_special::PLUS_INFINITY, range_type::INTEGER);
    minterm_coll minterms(expl_batch, Fset);

    dd_edge unexplored(Fset);
    dd_edge explored(init);
    dd_edge batch(Fset);
    dd_edge one(Fset);
    if ('i' == set_type) {
        Fset->createConstant(1, one);
    }
    dd_edge::iterator I(explored);
    RS = init;  // add initial state(s)
    for (;;) {
        switch (set_type) {
            case '+':
                Fset->createConstant(infty, unexplored);
                break;

            case 'i':
                Fset->createConstant(-1, unexplored);
                apply(PLUS, explored, one, explored);
                // ^ add one to states to be explored.
                // -1 for unreachable becomes 0,
                // so we can iterate quickly over the set.
                break;

            default:
                Fset->createConstant(false, unexplored);
                break;
        }
        I.restart(explored);
        if (!I) break;  // done exploring!
        for (; I; ++I) {
#ifdef DEBUG_GENERATE
            out << "Exploring state: ";
            (*I).show(out);
            out << "\n";
#endif
            // what's enabled
            for (unsigned e=0; e<events.size(); e++) {
                if (!events[e]) continue;
                if (!events[e]->fire(*I, minterms.unused())) continue;
                switch (set_type) {
                    case '+':
                        minterms.unused().setValue(1+int((*I).getValue()));
                        break;

                    case 'i':
                        minterms.unused().setValue((*I).getValue());
                        // Don't add one, we already did
                        break;

                    default:
                        minterms.unused().setValue(true);
                }
#ifdef DEBUG_GENERATE
                if (events[e]->get_name()) {
                    out << "  -- (event " << events[e]->get_name() << ") --> ";
                } else {
                    out << "  -- (event " << e << ") --> ";
                }
                minterms.unused().show(out);
                out << "\n";
#endif

                // Have we seen this already in the RS
                rangeval seen;
                RS.evaluate(minterms.unused(), seen);
                switch (set_type) {
                    case '+':
                        if (!seen.isPlusInfinity()) continue;
                        break;

                    case 'i':
                        if (long(seen)>=0) continue;
                        break;

                    default:
                        if (seen) continue;
                }

                minterms.pushUnused();
                if (minterms.size() >= expl_batch) {
                    // Buffer is full; flush it
                    switch (set_type) {
                        case '+':
                            minterms.buildFunctionMin(infty, batch);
                            apply(MINIMUM, unexplored, batch, unexplored);
                            apply(MINIMUM, RS, batch, RS);
                            break;

                        case 'i':
                            minterms.buildFunctionMin(-1, batch);
                            // ^ this might be wrong
                            apply(DIST_MIN, unexplored, batch, unexplored);
                            apply(DIST_MIN, RS, batch, RS);
                            break;

                        default:
                            minterms.buildFunctionMax(false, batch);
                            apply(UNION, unexplored, batch, unexplored);
                            apply(UNION, RS, batch, RS);
                    }
                    minterms.clear();
                }
            } // for e
        } // for I

        // explored is empty.
        // Flush the buffer
        if (minterms.size()) {
            switch (set_type) {
                case '+':
                    minterms.buildFunctionMin(infty, batch);
                    apply(MINIMUM, unexplored, batch, unexplored);
                    apply(MINIMUM, RS, batch, RS);
                    break;

                case 'i':
                    minterms.buildFunctionMin(-1, batch);
                    // ^ this might be wrong
                    apply(DIST_MIN, unexplored, batch, unexplored);
                    apply(DIST_MIN, RS, batch, RS);
                    break;

                default:
                    minterms.buildFunctionMax(false, batch);
                    apply(UNION, unexplored, batch, unexplored);
                    apply(UNION, RS, batch, RS);
            }
            minterms.clear();
        }

        explored = unexplored;
    } // infinite loop
}

void rsgen_model::show_cardinality(MEDDLY::output &out,
        const MEDDLY::dd_edge &e, bool with_timing)
{
    using namespace MEDDLY;
    timer watch;
    if (approx) {
        double card;

        watch.note_time();
        apply(CARDINALITY, e, card);
        watch.note_time();

        out   << "Approximately ";
        out.put(card, 10, 6);
        out << " reachable states\n";
    }
#ifdef HAVE_LIBGMP
    else {
        mpz_t card_mpz;
        mpz_init(card_mpz);
        watch.note_time();
        apply(CARDINALITY, e, card_mpz);
        watch.note_time();

        out << "Exactly ";
        out.put_huge(card_mpz, true);
        out << " reachable states\n";

        mpz_clear(card_mpz);
    }
#endif
    if (with_timing) {
        out << "    cardinality operation took "
            << watch.get_last_seconds() << " seconds\n";
    }
    out.flush();
}

void rsgen_model::printStats(MEDDLY::output &out, const MEDDLY::forest* F)
{
    using namespace MEDDLY;
    if (F) {
        if (F->isEVPlus()) {
            out << "EV+";
        }
        if (F->isMultiTerminal()) {
            if (F->isRangeType(range_type::INTEGER)) {
                out << "MT";
            }
        }
        if (F->isForRelations()) {
            out << "MxD";
        } else {
            out << "MDD";
        }
        out << " stats:\n";
        F->reportStats(out, "\t",
            HUMAN_READABLE_MEMORY | BASIC_STATS
            | EXTRA_STATS
            | STORAGE_STATS
            // | STORAGE_DETAILED
            // | HOLE_MANAGER_STATS
            // | HOLE_MANAGER_DETAILED
        );
    }
}

void rsgen_model::write_to_file(const char* what,
        const MEDDLY::dd_edge &who, const char* outfile)
{
    if (!outfile) return;
    compressed_output writer(outfile);
    if (writer) {
        std::cout << "Writing " << what << " to " << outfile << "\n";
        D->write(writer.outstream());
        MEDDLY::mdd_writer mywriter(writer.outstream(), who.getForest());
        mywriter.writeRootEdge(who);
        mywriter.finish();
    } else {
        std::cout << "Couldn't write to file " << outfile << " for " << what << "\n";
    }
}

void rsgen_model::trad_with_checkpoints(const MEDDLY::dd_edge &relation,
        const MEDDLY::dd_edge &initial, MEDDLY::dd_edge &reachable)
{
    using namespace MEDDLY;

    //
    // Determine next checkpoint to compute
    //
    unsigned next = 0;
    for (; ; ++next) {
        std::ifstream f(chkpt_path(next));
        if (!f.good()) break;
    }

    if (0==next) {
        std::cout << "No checkpoints found\n";
        apply(COPY, initial, reachable);
        write_checkpoint(0, reachable);
    } else {
        --next;
        std::cout << "Restoring from checkpoint " << chkpt_path(next) << "\n";
        read_checkpoint(next, reachable);
    }


    do {
        ++next;
        if (verb) my_progress(next, ' ');

        //
        // Compute next states
        //
        dd_edge step(reachable);
        if (forwd) {
            apply(POST_IMAGE, reachable, relation, step);
        } else {
            apply(PRE_IMAGE, reachable, relation, step);
        }
        if (verb) my_progress(next, 'N');

        //
        // Union with reachable
        //
        switch (set_type) {
            case '+':
                apply(MINIMUM, reachable, step, step);
                break;

            case 'i':
                apply(DIST_MIN, reachable, step, step);
                break;

            default:
                apply(UNION, reachable, step, step);
        }

        if (verb) my_progress(next, ';');

        if (step == reachable) break;
        reachable = step;

        //
        // Write to checkpoint file
        //
        write_checkpoint(next, reachable);

    } while(true);

}

// ******************************************************************
// *                  rsgen_model  private methods                  *
// ******************************************************************

const char* rsgen_model::chkpt_path(unsigned iter) const
{
    if (!check_dir) {
        throw "no checkpoint directory";
    }
    static char pathname[512];
    snprintf(pathname, 512, "%s/chk.%03d.txt.xz", check_dir, iter);
    return pathname;
}

void rsgen_model::write_checkpoint(unsigned iter, const MEDDLY::dd_edge &e)
    const
{
    const char* pathname = chkpt_path(iter);
    if (verb>2) {
        std::cerr << "            writing checkpoint " << pathname << " ";
    }
    compressed_output cpout(pathname);
    if (!cpout) {
        if (verb>2) {
            std::cerr << "couldn't write\n";
        }
        return;
    }
    if (verb>2) {
        std::cerr << ".";
        std::cerr.flush();
    }

    D->write(cpout.outstream());
    MEDDLY::mdd_writer mywriter(cpout.outstream(), Fset);
    mywriter.writeRootEdge(e);
    if (verb>2) {
        std::cerr << ".";
        std::cerr.flush();
    }
    mywriter.finish();
    if (verb>2) {
        std::cerr << " done" << std::endl;
    }
}

void rsgen_model::read_checkpoint(unsigned iter, MEDDLY::dd_edge &e) const
{
    const char* pathname = chkpt_path(iter);
    std::cerr << "    restoring from checkpoint " << pathname << " ";
    compressed_input cpin(pathname);
    if (!cpin) {
        std::cerr << "couldn't read\n";
        throw "checkpoint restore error: couldn't read";
    }
    std::cerr << ".";
    std::cerr.flush();

    D->verify(cpin.instream());
    MEDDLY::mdd_reader myreader(cpin.instream(), Fset);
    std::cerr << ".";
    std::cerr.flush();
    e = myreader.getRoot(0);
    std::cerr << " done\n";
}

