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

#include "meddly/dd_edge.h"
#include "meddly/node_marker.h"
#include "meddly/forest.h"
#include "meddly/terminal.h"

#include "meddly/io_dot.h"
#include "meddly/operators.h"

const char* edgecolor = "blue";
const char* evcolor = "purple";

inline void dot_edge_label(MEDDLY::output &s)
{
   s << " [samehead = true];";
}

inline void dot_edge_label(MEDDLY::output &s, const MEDDLY::edge_value &e)
{
    if (e.isVoid()) {
        s << " [samehead = true];";
    } else {
        s << " [label = \"";
        e.show(s, false);
        s << "\", fontcolor = \""
          << evcolor
          << "\", samehead = true];";
    }
}

// ******************************************************************
// *                                                                *
// *                 dot_maker::writing_style class                 *
// *                                                                *
// ******************************************************************

class MEDDLY::dot_maker::style : public node_marker::writing_style {
    public:
        style(const forest* F, node_marker &nm);
        virtual ~style();

        virtual void begin_level(output &s, int k);
        virtual void show_node(output &s, node_handle p);
        virtual void end_level(output &s, int k);

    private:
        std::vector <node_handle> nodes_at_this_level;
        int prev_level;
};

// ******************************************************************

MEDDLY::dot_maker::style::style(const forest* F, node_marker &)
    : node_marker::writing_style(F, SPARSE_ONLY)
{
    prev_level = 0;
}

MEDDLY::dot_maker::style::~style()
{
}

void MEDDLY::dot_maker::style::begin_level(output &s, int k)
{
    nodes_at_this_level.clear();

    const int k_rank = (k<0) ? (-k*2) -1 : k*2;
    MEDDLY_DCASSERT(k_rank > 0);

    //
    // Write the level label
    //
    s << "\n";
    s << "    L" << k_rank << " [shape=plaintext, label=\"Level: ";
    const domain* D = For->getDomain();
    const variable* v = D->getVar(unsigned(For->getVarByLevel(ABS(k))));
    if (v->hasName()) {
        s << v->getName();
    } else {
        s << ABS(k);
    }
    if (k<0) s << "'";
    s << "\"];\n";

    //
    // Chain the levels
    //
    s << "    edge[color=transparent];\n";
    s << "    L" << k_rank << " -> L" << prev_level << ";\n";
    s << "    edge [color=" << edgecolor << "];\n";
}

void MEDDLY::dot_maker::style::show_node(output &s, node_handle p)
{
    MEDDLY_DCASSERT(p>0);

    nodes_at_this_level.push_back(p);

    U->initFromNode(p);

    //
    // Draw the node
    //
    const double nodewidth = 0.25 * U->getSize();
    s << "        N" << p << " [";
    if (U->getSize()) s << "width=" << nodewidth << ", ";
    s << "label=\"";
    for (unsigned j=0; j<U->getSize(); j++) {
        if (j) s << "|";
        s << "<" << j << ">" << U->index(j);
    }
    s << "\"];\n";

    //
    // Draw the outgoing edges
    //
    for (unsigned j=0; j<U->getSize(); j++) {
        s << "            N" << p << ":" << j << " -> ";
        if (U->down(j)>0) {
            s << "N" << U->down(j);
        } else {
            s << "T" << -U->down(j);
        }
        if (U->hasEdges()) {
            dot_edge_label(s, U->edgeval(j));
        } else {
            dot_edge_label(s);
        }
        s << "\n";
    }
}

void MEDDLY::dot_maker::style::end_level(output &s, int k)
{
    const int k_rank = (k<0) ? (-k*2) -1 : k*2;
    MEDDLY_DCASSERT(k_rank > 0);

    //
    // Set these nodes to have the same rank.
    //
    s << "    { rank=same; L" << k_rank;
    for (size_t i=0; i<nodes_at_this_level.size(); i++) {
        s << " N" << nodes_at_this_level[i];
    }
    s << "; }\n";

    prev_level = k_rank;
}


// ******************************************************************
// *                                                                *
// *                       dot_maker  methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::dot_maker::dot_maker(forest* F, const char* bn)
{
    basename = bn;
    std::string fname = basename + ".dot";
    outfile.open(fname);

    For = F;
    MEDDLY_DCASSERT(For);

    nm = new node_marker(F);
    MEDDLY_DCASSERT(nm);
}

// ============================================================

MEDDLY::dot_maker::~dot_maker()
{
    delete nm;
}

// ============================================================

void MEDDLY::dot_maker::addRootEdge(const dd_edge &E)
{
    if (!nm) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
    if (!E.isAttachedTo(For)) {
        throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
    }

    roots.push_back(E);
    nm->mark(E.getNode());
}

// ============================================================

void MEDDLY::dot_maker::doneGraph()
{
    if (!nm) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }

    if (!outfile) {
        throw error(error::COULDNT_WRITE, __FILE__, __LINE__);
    }

    //
    // create output helper
    style WS(For, *nm);

    //
    // Some default options
    //
    const char* edgecolor = "blue";


    //
    // initialize the dot file
    //
    outfile << "digraph structs {\n";
    outfile << "    size=\"5,5\";\n";
    outfile << "    node [shape=record, height=0.5];\n";

    //
    // Build set of terminal nodes
    //
    std::set <node_handle> terms;
    nm->getTerminals(terms);
    for (unsigned i=0; i<roots.size(); i++) {
        if (roots[i].getNode() <= 0) {
            terms.insert(roots[i].getNode());
        }
    }

    //
    // Write the terminal nodes
    //
    outfile << "\n";
    outfile << "    L0 [shape=plaintext, label=\"Terminals\"];\n";
    outfile << "    { rank=same; L0";
    for (auto t = terms.begin(); t != terms.end(); t++)
    {
        outfile << " T" << -(*t);
    }
    outfile << "; }\n";
    for (auto t = terms.begin(); t != terms.end(); t++)
    {
        outfile << "        T" << -(*t);
        outfile << " [shape=box, style=bold, label=\"";
        terminal T(For->getTerminalType(), *t);

        switch (For->getTerminalType()) {
            case terminal_type::OMEGA:
                outfile << T.getOmega();
                break;

            case terminal_type::BOOLEAN:
                outfile << (T.getBoolean() ? "true" : "false");
                break;

            case terminal_type::INTEGER:
                outfile << T.getInteger();
                break;

            case terminal_type::REAL:
                outfile << T.getReal();
                break;

            default:
                FAIL(__FILE__, __LINE__, "Unknown terminal type");
        }

        outfile << "\"];\n";
    }

    //
    // Write all other nodes
    //
    ostream_output mout(outfile);
    nm->showByLevelsBottomUp(mout, &WS);

    //
    // Write all the root labels
    //
    mout << "\n";
    const int rootlevel = int(2 * (For->getNumVariables()+1));
    mout << "    L" << rootlevel << "[shape=plaintext, label=\"\"];\n";
    mout << "    edge [color=transparent];\n";
    mout << "    L" << rootlevel << " -> L" << rootlevel-2 << ";\n";
    mout << "    {rank=same; L" << rootlevel;
    for (unsigned i=0; i<roots.size(); i++) {
        mout << " R" << i;
    }
    mout << "; }\n";
    for (unsigned i=0; i<roots.size(); i++) {
        mout << "        R" << i << " [shape=plaintext, label=\"";
        if (roots[i].getLabel()) {
            mout << roots[i].getLabel();
        }
        mout << "\"];\n";
    }

    //
    // Root edges
    //
    mout << "    edge [color=" << edgecolor << "];\n";
    for (unsigned i=0; i<roots.size(); i++) {
        mout << "        R" << i << " -> ";
        if (roots[i].getNode() > 0) {
            mout << "N" << roots[i].getNode();
        } else {
            mout << "T" << -roots[i].getNode();
        }
        dot_edge_label(mout, roots[i].getEdgeValue());
        mout << "\n";
    }

    mout << "}\n";

    delete nm;
    nm = nullptr;
}

// ============================================================

void MEDDLY::dot_maker::runDot(const char* ext)
{
    if (nm) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }

    if (!ext) return;
    if (0==strcmp("dot", ext)) return;

    //
    // Build the command string
    //
    std::string cmd = "dot -T";
    cmd += ext;
    cmd += " -o \"";
    cmd += basename;
    cmd += ".";
    cmd += ext;
    cmd += "\" \"";
    cmd += basename;
    cmd += ".dot\"";

    //
    // Run dot on the dotfile
    //
    if (system(cmd.c_str())) {
        std::cerr << __func__ << ": Error executing DOT command:\n\t" << cmd << "\n";
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
}
