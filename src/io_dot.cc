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
#include "node_marker.h"

#include "io_dot.h"

MEDDLY::dot_maker::dot_maker(const char* bn)
{
    basename = bn;
    std::string fname = basename + ".dot";
    outfile.open(fname);

    graph_started = false;
    graph_finished = false;
}

MEDDLY::dot_maker::~dot_maker()
{
}

void MEDDLY::dot_maker::startGraph(const node_marker &nm)
{
    if (graph_started) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
    graph_started = true;

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
    const char* nodecolor = "black";


    //
    // initialize the dot file
    //
    outfile << "digraph structs {\n";
    outfile << "    size=\"5,5\";\n";
    outfile << "    node [shape=record, height=0.5];\n";

    //
    // Loop over levels, from bottom up
    //
    int prev_level = 0;
    const int TopLevel = nm.Forest()->getNumVariables();
    const domain* D = nm.Forest()->getDomain();
    unpacked_node* M = unpacked_node::New();

    for (int k= nm.Forest()->isForRelations() ? -1 : 1; ABS(k) <= TopLevel; )
    {
        // Graph rank for this level
        const int k_rank = (k<0) ? (-k*2) -1 : k*2;
        MEDDLY_DCASSERT(k_rank > 0);

        //
        // Write the level label
        //
        outfile << "\n";
        outfile << "    l" << k_rank << " [shape=plaintext, label=\"Level: ";
        const variable* v = D->getVar(nm.Forest()->getVarByLevel(ABS(k)));
        if (v->getName()) {
            outfile << v->getName();
        } else {
            outfile << ABS(k);
        }
        if (k<0) outfile << "'";
        outfile << "\"];\n";

        //
        // Chain the levels
        //
        if (prev_level) {
            outfile << "    edge[color=transparent];\n";
            outfile << "    l" << k_rank << " -> l" << prev_level << ";\n";
        }

        //
        // Grab the list of marked nodes at this level.
        //
        nodeList.clear();
        nm.getNodesAtLevel(k, nodeList);

        // Do nothing if no nodes.
        if (nodeList.size())
        {

            //
            // Set these nodes to have the same rank.
            //
            outfile << "    { rank=same; l" << k_rank;
            for (size_t i=0; i<nodeList.size(); i++) {
                outfile << " s" << nodeList[i];
            }
            outfile << ";}\n";

            //
            // Display the nodes
            //
            outfile << "    edge [color=" << nodecolor << "];\n";
            for (size_t i=0; i<nodeList.size(); i++) {
                nm.Forest()->unpackNode(M, nodeList[i], SPARSE_ONLY);

                const double nodewidth = 0.25 * M->getNNZs();
                outfile << "        s" << nodeList[i] << " [";
                if (M->getNNZs()) outfile << "width=" << nodewidth << ", ";
                outfile << "label=\"";
                for (unsigned j=0; j<M->getNNZs(); j++) {
                    if (j) outfile << "|";
                    outfile << "<" << j << ">" << M->i(j);
                    if (M->d(j)<0) outfile <<":T";    // TBD: show terminal
                }
                outfile << "\"];\n";
            }

            //
            // Display the edges
            //
            outfile << "    edge [color=" << edgecolor << "];\n";
            for (size_t i=0; i<nodeList.size(); i++) {
                nm.Forest()->unpackNode(M, nodeList[i], SPARSE_ONLY);

                for (unsigned j=0; j<M->getNNZs(); j++) {
                    if (M->d(j)<0) continue;
                    outfile << "        s" << nodeList[i] << ":"
                            << j << " -> s" << M->d(j)
                            << " [samehead = true];\n";
                }
            }

        } // if there are nodes at this level

        //
        // Get ready for next iteration
        //
        prev_level = k_rank;
        if (nm.Forest()->isForRelations()) {
            k = (k>=0) ? -(k+1) : -k;
        } else {
            ++k;
        }
    }

    unpacked_node::recycle(M);
}

void MEDDLY::dot_maker::addRootEdge(const dd_edge &E)
{
    if (!graph_started) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
    if (graph_finished) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }


    // TBD: add E to a vector, because we need to write all roots at once
    // (to have the same rank)
}

void MEDDLY::dot_maker::doneGraph()
{
    if (!graph_started) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
    if (graph_finished) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }
    graph_finished = true;

    if (!outfile) {
        throw error(error::COULDNT_WRITE, __FILE__, __LINE__);
    }

    // TBD: write all the root edges
    // and close out the graph file.
    //
    outfile << "}\n";
    outfile.close();
}

void MEDDLY::dot_maker::runDot(const char* ext)
{
    if (!graph_finished) {
        throw error(error::MISCELLANEOUS, __FILE__, __LINE__);
    }

    // TBD
}
