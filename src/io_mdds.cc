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

#include "operators.h"
#include "io_mdds.h"

#include <vector>
#include <string>

// ******************************************************************
// *                                                                *
// *                         helper methods                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    const unsigned CODELEN = 8;

    void buildCodeChars(const MEDDLY::forest* F, char* buf, char* rev)
    {
        buf[0] = 'd';
        buf[1] = 'd';
        buf[2] = '_';

        switch (F->getEdgeLabeling()) {
            case edge_labeling::MULTI_TERMINAL:
                buf[3] = 'M';
                buf[4] = 't';
                break;

            case edge_labeling::EVPLUS:
            case edge_labeling::INDEX_SET:
                buf[3] = 'E';
                buf[4] = '+';
                break;

            case edge_labeling::EVTIMES:
                buf[3] = 'E';
                buf[4] = '*';
                break;

            default:
                MEDDLY_DCASSERT(false);
        }

        if (F->isForRelations()) {
            buf[5] = 'x';
        } else {
            buf[5] = 'v';
        }

        switch (F->getRangeType()) {
            case range_type::BOOLEAN:
                buf[6] = 'b';
                break;

            case range_type::INTEGER:
                buf[6] = 'i';
                break;

            case range_type::REAL:
                buf[6] = 'r';
                break;

            default:
                MEDDLY_DCASSERT(false);
        }
        buf[7] = 0;

        rev[0] = buf[6];
        rev[1] = buf[5];
        rev[2] = buf[4];
        rev[3] = buf[3];
        rev[4] = buf[2];
        rev[5] = buf[1];
        rev[6] = buf[0];
        rev[7] = 0;
    }
}

// ******************************************************************
// *                                                                *
// *                       mdd_writer methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::mdd_writer::mdd_writer(output &s, const forest* F)
    : out(s)
{
    MEDDLY_DCASSERT(F);
    For = F;
    finished = false;
}

MEDDLY::mdd_writer::~mdd_writer()
{
    finish();
}

void MEDDLY::mdd_writer::finish()
{
    if (finished) return;
    finished = true;

    //
    // Mark all reachable nodes
    //
    node_marker M(For);
    for (unsigned i=0; i<roots.size(); i++) {
        M.mark(roots[i].getNode());
    }

    //
    // Build list of nodes to write, in bottom-up order.
    // When finished, output2handle[i] gives the handle
    // of the ith node to write.
    //
    std::vector <node_handle> output2handle;
    for (int k=1; k <= (int) For->getNumVariables(); k++) {
        if (For->isForRelations()) {
            M.getNodesAtLevel(-k, output2handle);
        }

        M.getNodesAtLevel(k, output2handle);
    } // for k
    const unsigned num_nodes = output2handle.size();

    //
    // Build the inverse list.
    // When finished, handle2output[h] gives 0 if the node with handle h
    // will not be written; otherwise it gives the order it appears in
    // the output file.
    //
    node_handle maxnode = 0;
    for (unsigned i=0; i<output2handle.size(); i++) {
        UPDATEMAX(maxnode, output2handle[i]);
    }
    std::vector <unsigned> handle2output (unsigned(1+maxnode));
    for (unsigned i=0; i<output2handle.size(); i++) {
        MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 1, output2handle[i], maxnode+1);
        handle2output[unsigned(output2handle[i])] = i+1;
    }

#ifdef DEBUG_WRITE
    std::cerr << "Got list of nodes:\n";
    for (int i=0; output2handle[i]; i++) {
        if (i) std::cerr << ", ";
        std::cerr << output2handle[i];
    }
    std::cerr << "\n";
    std::cerr << "Got inverse list:\n";
    for (int i=0; i<=maxnode; i++) {
        if (i) std::cerr << ", ";
        std::cerr << handle2output[i];
    }
    std::cerr << "\n";
#endif

    //
    // Generate code char sequence
    //
    char block[CODELEN], revblock[CODELEN];
    buildCodeChars(For, block, revblock);

    //
    // Write the nodes
    //
    out << block << " " << num_nodes << "\n";
    unpacked_node* un = unpacked_node::New(For, FULL_OR_SPARSE);
    for (unsigned i=0; i<output2handle.size(); i++) {
        out << For->getNodeLevel(output2handle[i]) << " ";
        un->initFromNode(output2handle[i]);
        // For->unpackNode(un, output2handle[i], FULL_OR_SPARSE);
        un->write(out, handle2output);
    }
    unpacked_node::Recycle(un);
    out << revblock << "\n";

    //
    // Write the actual edge pointers
    //
    out << "ptrs " << roots.size() << "\n";
    for (unsigned i=0; i<roots.size(); i++) {
        out.put('\t');
        roots[i].write(out, handle2output);
        out.put('\n');
    }
    out << "srtp\n";
    out.flush();
}

// ******************************************************************
// *                                                                *
// *                       mdd_writer methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::mdd_reader::mdd_reader(input &s, forest* F)
{
    MEDDLY_DCASSERT(F);

    rptr = 0;

    //
    // Generate code char sequence
    //
    char block[CODELEN], revblock[CODELEN];
    buildCodeChars(F, block, revblock);

    //
    // Read the entire file and store the roots.
    //

#ifdef DEBUG_READ
    try {
#endif

        s.stripWS();
        s.consumeKeyword(block);
        s.stripWS();
        unsigned num_nodes = unsigned(s.get_integer());
#ifdef DEBUG_READ_DD
        std::cerr << "Reading " << num_nodes << " nodes in " << block << " forest\n";
#endif
#ifdef DEBUG_READ
        std::cerr << "Reading " << num_nodes << " nodes in forest " << block << "\n";

#endif

        //
        // Build translation from file node# to forest node handles
        //
        std::vector <node_handle> map(1+num_nodes);

        for (unsigned node_index=1; node_index<=num_nodes; node_index++) {
            // Read this node

            //
            // read the level number
            //
            s.stripWS();
            int k = s.get_integer();
            if (!F->isValidLevel(k)) {
                throw error(error::INVALID_LEVEL, __FILE__, __LINE__);
            }

#ifdef DEBUG_READ
            std::cerr << "Reading: level " << k;
#endif

            //
            // read the node size (sparse/full)
            //
            s.stripWS();
            int rawsize = s.get_integer();
            // int n;
            unpacked_node* nb = unpacked_node::newWritable(F, k, ABS(rawsize),
                    (rawsize < 0) ? SPARSE_ONLY : FULL_ONLY);

#ifdef DEBUG_READ
            std::cerr << " rawsize " << rawsize << '\n';
#endif

            //
            // read the node
            //
            nb->read(s, map);

            //
            // Reduce the node, and update the translation
            //
            edge_value ev;
            F->createReducedNode(nb, ev, map[node_index]);

            MEDDLY_DCASSERT( !ev.isInt()    || (0==ev.getInt()) );
            MEDDLY_DCASSERT( !ev.isLong()   || (0==ev.getLong()) );
            MEDDLY_DCASSERT( !ev.isFloat()  || (0==ev.getFloat()) );
            MEDDLY_DCASSERT( !ev.isDouble() || (0==ev.getDouble()) );

#ifdef DEBUG_READ
            std::cerr << "File node " << node_index << " reduced to ";
            ostream_output myout(std::cerr);
            showNode(myout, map[node_index], SHOW_DETAILS | SHOW_INDEX);
            std::cerr << std::endl;
#endif

        } // for node_index

#ifdef DEBUG_READ
        std::cerr << "Done reading, expecting " << revblock << " keyword\n";
#endif

        //
        // match the reversed block
        //
        s.stripWS();
        s.consumeKeyword(revblock);
#ifdef DEBUG_READ
        std::cerr << "  got " << revblock << "\n";
#endif
#ifdef DEBUG_READ_DD
        std::cerr << "Finished " << block << " forest\n";
#endif

        //
        // Read the pointers
        //
        s.stripWS();
        s.consumeKeyword("ptrs");
#ifdef DEBUG_READ
        std::cerr << "Got ptrs\n";
#endif

        s.stripWS();
        unsigned num_ptrs = unsigned(s.get_integer());
#ifdef DEBUG_READ
        std::cerr << "Reading " << num_ptrs << " pointers\n";
#endif
#ifdef DEBUG_READ_DD
        std::cerr << "Reading " << num_ptrs << " pointers\n";
#endif
        for (unsigned i=0; i<num_ptrs; i++) {
            dd_edge E(F);
            E.read(s, map);
            roots.push_back(E);
        }

        s.stripWS();
        s.consumeKeyword("srtp");

#ifdef DEBUG_READ_DD
        std::cerr << "Done reading pointers\n";
#endif

        //
        // unlink map pointers
        //
        for (unsigned i=1; i<=num_nodes; i++) {
            F->unlinkNode(map[i]);
        }

#ifdef DEVELOPMENT_CODE
        F->validateIncounts(true, __FILE__, __LINE__);
#endif
#ifdef DEBUG_READ
    } // try
    catch (error& e) {
        std::cerr << "Read failed (error: " << e.getName() << ")\n";
        std::cerr << "Next few characters of file:\n";
        for (unsigned i=0; i<20; i++) {
            int c = s.get_char();
            if (EOF == c) {
                std::cerr << "EOF";
                break;
            }
            std::cerr << char(c) << " (" << c << ") ";
        }
        std::cerr << std::endl;
        throw e;
    }
#endif
}

void MEDDLY::mdd_reader::readRootEdge(dd_edge &E)
{
    if (rptr >= roots.size()) {
        throw error(error::COULDNT_READ, __FILE__, __LINE__);
    }
    E = roots[rptr];
    ++rptr;
}

