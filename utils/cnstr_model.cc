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

#include "utils/cnstr_model.h"

// ******************************************************************

inline void show_comparison(MEDDLY::output &out, cnstr_model::comparison cmp)
{
    switch (cmp) {
        case cnstr_model::comparison::EQ:   return out.put("==");
        case cnstr_model::comparison::NE:   return out.put("!=");
        case cnstr_model::comparison::LT:   return out.put("< ");
        case cnstr_model::comparison::LE:   return out.put("<=");
        case cnstr_model::comparison::GT:   return out.put("> ");
        case cnstr_model::comparison::GE:   return out.put(">=");
    }
}

// ******************************************************************
// *                                                                *
// *                cnstr_model::constraint  methods                *
// *                                                                *
// ******************************************************************

cnstr_model::constraint::constraint(unsigned L) : depends(L)
{
    for (unsigned i=0; i<L; i++) {
        depends[i] = false;
    }
}

// ******************************************************************

cnstr_model::constraint::~constraint()
{
}

// ******************************************************************

void cnstr_model::constraint::show_deps(MEDDLY::output &out) const
{
    for (unsigned i=0; i<depends.size(); i++) {
        if (depends[i]) out.put('X');
        else            out.put('.');
    }
}

// ******************************************************************

bool cnstr_model::constraint::deps_less_than(const constraint* C) const
{
    if (!C) return true;
    MEDDLY_DCASSERT(depends.size() == C->depends.size());
    for (unsigned i=depends.size(); i; ) {
        --i;
        if (depends[i]) {
            if (C->depends[i]) continue;
            return false;
        }
        if (C->depends[i]) {
            return true;
        }
    }
    // equal
    return false;
}

// ******************************************************************
// *                                                                *
// *             cnstr_model::compare_two_vars  methods             *
// *                                                                *
// ******************************************************************

cnstr_model::compare_two_vars::compare_two_vars(unsigned L, unsigned _v1,
        comparison c, unsigned _v2, int d) : constraint(L)
{
    v1 = _v1;
    v2 = _v2;
    delta = d;
    cmp = c;

    add_dependency(v1);
    add_dependency(v2);
}

// ******************************************************************

void cnstr_model::compare_two_vars::build_edge(MEDDLY::dd_edge &r)
{
    using namespace MEDDLY;
    forest* mdd = r.getForest();

    if (v1 == v2) {
        //
        // Trivial case; probably shouldn't happen
        //
        if (delta) {
            mdd->createConstant(false, r);
        } else {
            mdd->createConstant(true, r);
        }
        return;
    }

    //
    // This assumes the MDD is fully reduced.
    //
    MEDDLY_DCASSERT(mdd->isFullyReduced());

    node_handle root;
    edge_value ev;

    if (v1 < v2) {
        unpacked_node* U2 = unpacked_node::newWritable(mdd, v2, FULL_ONLY);
        for (unsigned i2=0; i2<U2->getSize(); i2++) {
            unpacked_node* U1 = unpacked_node::newWritable(mdd, v1, FULL_ONLY);
            for (unsigned i1=0; i1<U1->getSize(); i1++) {
                U1->setFull(i1, mdd->handleForValue(eval(i1, i2)));
            } // for i1
            mdd->createReducedNode(U1, ev, root);
            MEDDLY_DCASSERT(ev.isVoid());
            U2->setFull(i2, root);
        } // for i2
        mdd->createReducedNode(U2, ev, root);
    } else {
        unpacked_node* U1 = unpacked_node::newWritable(mdd, v1, FULL_ONLY);
        for (unsigned i1=0; i1<U1->getSize(); i1++) {
            unpacked_node* U2 = unpacked_node::newWritable(mdd, v2, FULL_ONLY);
            for (unsigned i2=0; i2<U2->getSize(); i2++) {
                U2->setFull(i2, mdd->handleForValue(eval(i1, i2)));
            } // for i2
            mdd->createReducedNode(U2, ev, root);
            MEDDLY_DCASSERT(ev.isVoid());
            U1->setFull(i1, root);
        } // for i1
        mdd->createReducedNode(U1, ev, root);
    }
    MEDDLY_DCASSERT(ev.isVoid());
    r.set(root);
}

// ******************************************************************

void cnstr_model::compare_two_vars::show(MEDDLY::output &out) const
{
    out << "x_" << v1;
    show_comparison(out, cmp);
    out << "x_" << v2;
    if (delta<0) {
        out << delta;
    }
    if (delta>0) {
        out << "+" << delta;
    }
}

// ******************************************************************
// *                                                                *
// *                                                                *
// *                      cnstr_model  methods                      *
// *                                                                *
// *                                                                *
// ******************************************************************

// ******************************************************************
// *                   cnstr_model public methods                   *
// ******************************************************************

cnstr_model::cnstr_model()
{
    D = nullptr;
    Fset = nullptr;

    // default options

    set_solve_method('a');
    set_garbage('m');
    displaysols = '0';
#ifdef HAVE_LIBGMP
    approx = false;
#else
    approx = true;
#endif
    verb = 0;
    file_set = nullptr;
}

// ******************************************************************

cnstr_model::~cnstr_model()
{
}

// ******************************************************************

bool cnstr_model::process_arg(unsigned &i, const char** argv)
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
    // Solution heuristics
    //

    if (0==strcmp(argv[i], "-acc")) {
        set_solve_method('a');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-fold")) {
        set_solve_method('f');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-revacc")) {
        set_solve_method('r');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-small")) {
        set_solve_method('<');
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-large")) {
        set_solve_method('>');
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
    // Solutions to display
    //
    if (0==strcmp(argv[i], "-zero")) {
        displaysols='0';
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-one")) {
        displaysols='1';
        i++;
        return true;
    }
    if (0==strcmp(argv[i], "-all")) {
        displaysols='a';
        i++;
        return true;
    }

    //
    // File output
    //
    if (0==strcmp(argv[i], "-ws")) {
        i++;
        file_set = argv[i];
        i++;
        return true;
    }

    return false;
}

// ******************************************************************

void cnstr_model::show_usage(std::ostream &out, const char* arg0)
{
    show_synopsis(out, arg0, "[options]");

    show_section(out, "SOLUTION METHODS");
    show_switch(out, "-acc", "Accumulate constraints in order of level dependencies");
    show_switch(out, "-fold", "Fold constraints together in pairs");
    show_switch(out, "-revacc", "Accumulate constraints in reverse order of level dependencies");
    show_switch(out, "-small", "Add smallest constraint (in #nodes) first");
    show_switch(out, "-large", "Add largeest constraint (in #nodes) first");
    out << "\n";
    show_switch(out, "-vL", "Set verbosity level. Range is 0 (super quiet, default) to 9 (super noisy).");

    show_section(out, "GARBAGE COLLECTION");
    show_switch(out, "-ms", "mark and sweep");
    show_switch(out, "-opt", "optimistic reference counts");
    show_switch(out, "-pess", "pessimistic reference counts");

    show_section(out, "FILE OUTPUT");
    show_switch(out, "-ws file", "write the generated solution set (MDD) to file. Output filenames ending in .gz or .xz will be compressed.");

    show_section(out, "DISPLAY SOLUTIONS");
    show_switch(out, "-zero", "Do not show any solutions.");
    show_switch(out, "-one",  "Show one of the solutions (if one exists).");
    show_switch(out, "-all",  "Show all solutions.");

#ifdef HAVE_LIBGMP
    show_section(out, "CARDINALITY");
    if (approx) {
        show_switch(out, "-approx", "approximate cardinalities (default)");
        show_switch(out, "-exact", "exact cardinalities");
    } else {
        show_switch(out, "-approx", "approximate cardinalities");
        show_switch(out, "-exact", "exact cardinalities (default)");
    }
#endif
    out << "\n";
}

// ******************************************************************

