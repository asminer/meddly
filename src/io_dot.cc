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

#include "dd_edge.h"
#include "forest.h"
#include "node_marker.h"
#include "terminal.h"

#include "io_dot.h"
#include "operators.h"

#define USE_NEW_WRITE

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
        const char* edgecolor;
};

// ******************************************************************

MEDDLY::dot_maker::style::style(const forest* F, node_marker &nm)
    : node_marker::writing_style(F, SPARSE_ONLY)
{
    prev_level = 0;
    edgecolor = "blue";
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
        s << " [samehead = true];\n";
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

MEDDLY::dot_maker::dot_maker(const forest* F, const char* bn)
{
    basename = bn;
    std::string fname = basename + ".dot";
    outfile.open(fname);

    For = F;
    MEDDLY_DCASSERT(For);

    nm = new node_marker(For);
    MEDDLY_DCASSERT(nm);
}

MEDDLY::dot_maker::~dot_maker()
{
    delete nm;
}

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

#ifdef USE_NEW_WRITE

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
    outfile << "\n";
    const int rootlevel = int(2 * (For->getNumVariables()+1));
    outfile << "    L" << rootlevel << "[shape=plaintext, label=\"\"];\n";
    outfile << "    edge [color=transparent];\n";
    outfile << "    L" << rootlevel << " -> L" << rootlevel-2 << ";\n";
    outfile << "    {rank=same; L" << rootlevel;
    for (unsigned i=0; i<roots.size(); i++) {
        outfile << " R" << i;
    }
    outfile << "; }\n";
    for (unsigned i=0; i<roots.size(); i++) {
        outfile << "        R" << i << " [shape=plaintext, label=\"";
        if (roots[i].getLabel()) {
            outfile << roots[i].getLabel();
        }
        outfile << "\"];\n";
    }

    //
    // Root edges
    //
    outfile << "    edge [color=" << edgecolor << "];\n";
    for (unsigned i=0; i<roots.size(); i++) {
        outfile << "        R" << i << " -> ";
        if (roots[i].getNode() > 0) {
            outfile << "N" << roots[i].getNode();
        } else {
            outfile << "T" << -roots[i].getNode();
        }
        outfile << " [samehead = true];\n";
    }

    outfile << "}\n";
    outfile.close();

    delete nm;
    nm = nullptr;
}

// ======================================================================
#else
// ======================================================================

void MEDDLY::dot_maker::doneGraph()
{
    if (!nm) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }

    if (!outfile) {
        throw error(error::COULDNT_WRITE, __FILE__, __LINE__);
    }

    //
    // For the list of nodes at each level.
    std::vector <node_handle> nodeList;

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
    // Loop over levels, from bottom up
    //
    int prev_level = 0;
    const int TopLevel = int(For->getNumVariables());
    const domain* D = For->getDomain();
    unpacked_node* M = unpacked_node::New(For, SPARSE_ONLY);

    for (int k= For->isForRelations() ? -1 : 1; ABS(k) <= TopLevel; )
    {
        // Graph rank for this level
        const int k_rank = (k<0) ? (-k*2) -1 : k*2;
        MEDDLY_DCASSERT(k_rank > 0);

        //
        // Write the level label
        //
        outfile << "\n";
        outfile << "    L" << k_rank << " [shape=plaintext, label=\"Level: ";
        const variable* v = D->getVar(unsigned(For->getVarByLevel(ABS(k))));
        if (v->hasName()) {
            outfile << v->getName();
        } else {
            outfile << ABS(k);
        }
        if (k<0) outfile << "'";
        outfile << "\"];\n";

        //
        // Chain the levels
        //
        outfile << "    edge[color=transparent];\n";
        outfile << "    L" << k_rank << " -> L" << prev_level << ";\n";

        //
        // Grab the list of marked nodes at this level.
        //
        nodeList.clear();
        nm->getNodesAtLevel(k, nodeList);

        // Do nothing if no nodes.
        if (nodeList.size())
        {

            //
            // Set these nodes to have the same rank.
            //
            outfile << "    { rank=same; L" << k_rank;
            for (size_t i=0; i<nodeList.size(); i++) {
                outfile << " N" << nodeList[i];
            }
            outfile << "; }\n";

            //
            // Display the nodes
            //
            for (size_t i=0; i<nodeList.size(); i++) {
                M->initFromNode(nodeList[i]);
                // For->unpackNode(M, nodeList[i], SPARSE_ONLY);

                const double nodewidth = 0.25 * M->getSize();
                outfile << "        N" << nodeList[i] << " [";
                if (M->getSize()) outfile << "width=" << nodewidth << ", ";
                outfile << "label=\"";
                for (unsigned j=0; j<M->getSize(); j++) {
                    if (j) outfile << "|";
                    outfile << "<" << j << ">" << M->index(j);
                }
                outfile << "\"];\n";
            }

            //
            // Display the edges
            //
            outfile << "    edge [color=" << edgecolor << "];\n";
            for (size_t i=0; i<nodeList.size(); i++) {
                M->initFromNode(nodeList[i]);
                // For->unpackNode(M, nodeList[i], SPARSE_ONLY);

                for (unsigned j=0; j<M->getSize(); j++) {
                    outfile << "        N" << nodeList[i] << ":" << j << " -> ";
                    if (M->down(j)>0) {
                        outfile << "N" << M->down(j);
                    } else {
                        outfile << "T" << -M->down(j);
                    }
                    outfile << " [samehead = true];\n";
                }
            }

        } // if there are nodes at this level

        //
        // Get ready for next iteration
        //
        prev_level = k_rank;
        if (For->isForRelations()) {
            k = (k>=0) ? -(k+1) : -k;
        } else {
            ++k;
        }
    }

    unpacked_node::Recycle(M);

    //
    // Write all the root labels
    //
    outfile << "\n";
    const int rootlevel = int(2 * (For->getNumVariables()+1));
    outfile << "    L" << rootlevel << "[shape=plaintext, label=\"\"];\n";
    outfile << "    edge [color=transparent];\n";
    outfile << "    L" << rootlevel << " -> L" << rootlevel-2 << ";\n";
    outfile << "    {rank=same; L" << rootlevel;
    for (unsigned i=0; i<roots.size(); i++) {
        outfile << " R" << i;
    }
    outfile << "; }\n";
    for (unsigned i=0; i<roots.size(); i++) {
        outfile << "        R" << i << " [shape=plaintext, label=\"";
        if (roots[i].getLabel()) {
            outfile << roots[i].getLabel();
        }
        outfile << "\"];\n";
    }

    //
    // Root edges
    //
    outfile << "    edge [color=" << edgecolor << "];\n";
    for (unsigned i=0; i<roots.size(); i++) {
        outfile << "        R" << i << " -> ";
        if (roots[i].getNode() > 0) {
            outfile << "N" << roots[i].getNode();
        } else {
            outfile << "T" << -roots[i].getNode();
        }
        outfile << " [samehead = true];\n";
    }

    outfile << "}\n";
    outfile.close();

    delete nm;
    nm = nullptr;
}

#endif

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
