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

// #define DEBUG_READ
#define USE_NEW_WRITE

// ******************************************************************
// *                                                                *
// *                         helper methods                         *
// *                                                                *
// ******************************************************************

namespace MEDDLY {
    const unsigned CODELEN = 8;

    void buildCodeChars(const forest* F, char* buf, char* rev)
    {
        buf[0] = 'd';
        buf[1] = 'd';
        buf[2] = '_';

        switch (F->getEdgeLabeling()) {
            case edge_labeling::MULTI_TERMINAL:
                buf[3] = 'M';
                buf[4] = 't';
                break;

            case edge_labeling::INDEX_SET:
                buf[3] = 'I';
                buf[4] = '+';
                break;

            case edge_labeling::EVPLUS:
                buf[3] = 'E';
                buf[4] = '+';
                break;

            case edge_labeling::EVTIMES:
                buf[3] = 'E';
                buf[4] = '*';
                break;

            default:
                FAIL(__FILE__, __LINE__, "Unknown edge labeling");
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
                FAIL(__FILE__, __LINE__, "Unknown range type");
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

    inline void validateChar(char x, const char* legal, const char* F, unsigned L)
    {
        for (const char* p = legal; *p; ++p) {
            if (x == *p) return;
        }
#ifdef DEBUG_READ
        std::cerr << "Got char '" << x << "', expected one of: " << legal
                  << "\n\tfrom source line " << L << "\n";
#endif
        throw error(error::INVALID_FILE, F, L);
    }

    forest* code2forest(const char* code, domain* D)
    {
        // sanity checks
        validateChar(code[0], "d", __FILE__, __LINE__);
        validateChar(code[1], "d", __FILE__, __LINE__);
        validateChar(code[2], "_", __FILE__, __LINE__);

        //
        // determine edge labeling
        //
        edge_labeling EL = edge_labeling::MULTI_TERMINAL;
        validateChar(code[3], "MIE", __FILE__, __LINE__);
        switch (code[3]) {
            case 'M':
                validateChar(code[4], "t", __FILE__, __LINE__);
                EL = edge_labeling::MULTI_TERMINAL;
                break;

            case 'I':
                validateChar(code[4], "+", __FILE__, __LINE__);
                EL = edge_labeling::INDEX_SET;
                break;

            case 'E':
                validateChar(code[4], "+*", __FILE__, __LINE__);
                if ('+' == code[4]) {
                    EL = edge_labeling::EVPLUS;
                } else {
                    EL = edge_labeling::EVTIMES;
                }
                break;

            default:
                FAIL(__FILE__, __LINE__, "Unknown edge label");
        }

        //
        // determine set/relation
        //
        validateChar(code[5], "xv", __FILE__, __LINE__);
        set_or_rel SR = ('x' == code[5]) ? RELATION : SET;

        //
        // determine range type
        //
        validateChar(code[6], "bir", __FILE__, __LINE__);
        range_type RT = range_type::BOOLEAN;
        switch (code[6]) {
            case 'b':
                RT = range_type::BOOLEAN;
                break;

            case 'i':
                RT = range_type::INTEGER;
                break;

            case 'r':
                RT = range_type::REAL;
                break;

            default:
                FAIL(__FILE__, __LINE__, "Unknown range type");
        }

        return forest::create(D, SR, RT, EL);
    }
}

// ******************************************************************
// *                                                                *
// *                mdd_writer::writing_style  class                *
// *                                                                *
// ******************************************************************

class MEDDLY::mdd_writer::style : public node_marker::writing_style {
    public:
        style(const forest* F, node_marker &nm);
        virtual ~style();

        virtual void begin_level(output &s, int k);
        virtual void show_node(output &s, node_handle p);
        virtual void end_level(output &s, int k);

        inline size_t numWritten() const {
            return written;
        }
        const std::vector <size_t>& getH2O() const {
            return handle2output;
        }

    private:
        std::vector <size_t> handle2output;
        size_t written;
};

// ******************************************************************

MEDDLY::mdd_writer::style::style(const forest* F, node_marker &nm)
    : node_marker::writing_style(F, FULL_OR_SPARSE)
{
    handle2output.resize(1+nm.lastMarked(), 0);
    written = 0;
}

MEDDLY::mdd_writer::style::~style()
{
}

void MEDDLY::mdd_writer::style::begin_level(output &s, int k)
{
}

void MEDDLY::mdd_writer::style::show_node(output &s, node_handle p)
{
    MEDDLY_DCASSERT(p>0);
    MEDDLY_DCASSERT(0 == handle2output.at(p));
    handle2output[p] = ++written;

    s << For->getNodeLevel(p) << " ";
    U->initFromNode(p);
    U->write(s, handle2output);
}

void MEDDLY::mdd_writer::style::end_level(output &s, int k)
{
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

#ifdef USE_NEW_WRITE

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
    // Generate code char sequence
    //
    char block[CODELEN], revblock[CODELEN];
    buildCodeChars(For, block, revblock);

    //
    // Write the nodes
    //
    const size_t num_nodes = M.countMarked();
    out << block << " " << num_nodes << "\n";
    out.flush();
    style WS(For, M);
    M.showByLevelsBottomUp(out, &WS);
    MEDDLY_DCASSERT(num_nodes == WS.numWritten());
    out << revblock << "\n";

    //
    // Write the actual edge pointers
    //
    out << "ptrs " << roots.size() << "\n";
    for (unsigned i=0; i<roots.size(); i++) {
        out.put('\t');
        roots[i].write(out, WS.getH2O());
        out.put('\n');
    }
    out << "srtp\n";
    out.flush();
}

// ======================================================================
#else
// ======================================================================

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
    const size_t num_nodes = output2handle.size();

    //
    // Build the inverse list.
    // When finished, handle2output[h] gives 0 if the node with handle h
    // will not be written; otherwise it gives the order it appears in
    // the output file.
    //
    node_handle maxnode = 0;
    for (size_t i=0; i<output2handle.size(); i++) {
        UPDATEMAX(maxnode, output2handle[i]);
    }
    std::vector <size_t> handle2output (1+maxnode);
    for (size_t i=0; i<output2handle.size(); i++) {
        MEDDLY_CHECK_RANGE(1, output2handle[i], maxnode+1);
        handle2output[output2handle[i]] = i+1;
    }

#ifdef DEBUG_WRITE
    std::cerr << "Got list of nodes:\n";
    for (size_t i=0; output2handle[i]; i++) {
        if (i) std::cerr << ", ";
        std::cerr << output2handle[i];
    }
    std::cerr << "\n";
    std::cerr << "Got inverse list:\n";
    for (node_handle i=0; i<=maxnode; i++) {
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
    for (size_t i=0; i<output2handle.size(); i++) {
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

#endif

// ******************************************************************
// *                                                                *
// *                       mdd_reader methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::mdd_reader::mdd_reader(input &s, forest* F, output* prog)
{
    MEDDLY_DCASSERT(F);
    For = F;
    rptr = 0;
    file_nodes = 0;

    //
    // Generate code char sequence
    //
    char block[CODELEN], revblock[CODELEN];
    buildCodeChars(F, block, revblock);

    //
    // Read the forest code sequence
    //

#ifdef DEBUG_READ
    try {
#endif
        s.stripWS();
        s.consumeKeyword(block);

#ifdef DEBUG_READ
    } // try
    catch (error& e) {
        handleError(s, e);
    }
#endif

    readAfterForest(s, prog);
}

MEDDLY::mdd_reader::mdd_reader(input &s, domain* D, output* prog)
{
    rptr = 0;
    file_nodes = 0;

    //
    // Read code char sequence
    //
    char block[CODELEN];
#ifdef DEBUG_READ
    try {
#endif
        s.stripWS();

        block[CODELEN-1] = 0;
        for (unsigned i=0; i<CODELEN; i++) {
            int c = s.get_char();
            if ( (EOF == c) || (' ' == c) || ('\t' == c) || ('\n' == c) || ('\r' == c) )
            {
                block[i] = 0;
                break;
            }
            block[i] = c;
        }
#ifdef DEBUG_READ
        std::cerr << "Read code block: '" << block << "'\n";
#endif

        For = code2forest(block, D);

#ifdef DEBUG_READ
    } // try
    catch (error &e) {
        handleError(s, e);
    }
#endif

    readAfterForest(s, prog);
}

void MEDDLY::mdd_reader::readRootEdge(dd_edge &E)
{
    if (rptr >= roots.size()) {
        throw error(error::COULDNT_READ, __FILE__, __LINE__);
    }
    E = roots[rptr];
    ++rptr;
}

void MEDDLY::mdd_reader::readAfterForest(input &s, output* prog)
{
    MEDDLY_DCASSERT(For);

    //
    // Generate code char sequence
    //
    char block[CODELEN], revblock[CODELEN];
    buildCodeChars(For, block, revblock);

    //
    // Read the entire file and store the roots.
    //

#ifdef DEBUG_READ
    try {
#endif
        //
        // The MDD type has already been consumed
        // Pick up with the number of nodes
        //

        s.stripWS();
        file_nodes = size_t(s.get_integer());
#ifdef DEBUG_READ_DD
        std::cerr << "Reading " << file_nodes << " nodes in " << block << " forest\n";
#endif
#ifdef DEBUG_READ
        std::cerr << "Reading " << file_nodes << " nodes in forest " << block << "\n";

#endif
        if (prog) {
            (*prog) << "Reading " << file_nodes << " nodes in forest " << block << "\n";
        }

        //
        // Build translation from file node# to forest node handles
        //
        std::vector <node_handle> map(1+file_nodes);

        for (size_t node_index=1; node_index<=file_nodes; node_index++) {
            // Read this node

            //
            // read the level number
            //
            s.stripWS();
            int k = s.get_integer();
            if (!For->isValidLevel(k)) {
                throw error(error::INVALID_LEVEL, __FILE__, __LINE__);
            }

#ifdef DEBUG_READ
            std::cerr << "Reading #" << node_index << ": level " << k;
#endif

            //
            // read the node size (sparse/full)
            //
            s.stripWS();
            int rawsize = s.get_integer();
            // int n;
            unpacked_node* nb = unpacked_node::newWritable(For, k, ABS(rawsize),
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
            For->createReducedNode(nb, ev, map[node_index]);

            MEDDLY_DCASSERT( !ev.isInt()    || (0==int(ev)) );
            MEDDLY_DCASSERT( !ev.isLong()   || (0==long(ev)) );
            MEDDLY_DCASSERT( !ev.isFloat()  || (0==float(ev)) );
            MEDDLY_DCASSERT( !ev.isDouble() || (0==double(ev)) );

#ifdef DEBUG_READ
            std::cerr << "File node " << node_index << " reduced to ";
            ostream_output myout(std::cerr);
            For->showNode(myout, map[node_index], SHOW_DETAILS | SHOW_INDEX);
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
        if (prog) {
            (*prog) << "Finished " << block << " forest\n";
        }

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
        if (prog) {
            (*prog) << "Reading " << num_ptrs << " pointers\n";
        }
        for (unsigned i=0; i<num_ptrs; i++) {
            dd_edge E(For);
            E.read(s, map);
            roots.push_back(E);
        }

        s.stripWS();
        s.consumeKeyword("srtp");

#ifdef DEBUG_READ_DD
        std::cerr << "Done reading pointers\n";
#endif
        if (prog) {
            (*prog) << "Done reading pointers\n";
        }

        //
        // unlink map pointers
        //
        for (size_t i=1; i<=file_nodes; i++) {
            For->unlinkNode(map[i]);
        }

#ifdef DEVELOPMENT_CODE
        For->validateIncounts(true, __FILE__, __LINE__, "mdd reader");
#endif
#ifdef DEBUG_READ
    } // try
    catch (error& e) {
        handleError(s, e);
    }
#endif
}

void MEDDLY::mdd_reader::handleError(input &s, error &e) const
{
#ifdef DEBUG_READ
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
#endif
}
