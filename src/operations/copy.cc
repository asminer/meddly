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

#include "../defines.h"
#include "copy.h"

#include "../oper_unary.h"
#include "../ct_vector.h"
#include "../forest_levels.h"
#include "../forest_edgerules.h"

namespace MEDDLY {
    class copy_inforest;
    class copy_MT;
    class copy_EV_fast;

    class COPY_factory;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                      copy_inforest  class                      *
// *                                                                *
// ******************************************************************

/// Copy within a forest.
/// Provided for uniformity and because users mignt not
/// know any better.
///
class MEDDLY::copy_inforest : public unary_operation {
    public:
        copy_inforest(forest* f);
        virtual ~copy_inforest();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);
};

// ******************************************************************

MEDDLY::copy_inforest::copy_inforest(forest* f)
    : unary_operation(f, f)
{
}

MEDDLY::copy_inforest::~copy_inforest()
{
}

void MEDDLY::copy_inforest::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
    cv = av;
    cp = resF->linkNode(ap);
}


// ******************************************************************
// *                                                                *
// *                         copy_MT  class                         *
// *                                                                *
// ******************************************************************

/// New, hopefully faster, copy operation for multi-terminal DDs.
/// The target forest can be multi-terminal or edge-valued :)
///
class MEDDLY::copy_MT : public unary_operation {
    public:
        copy_MT(forest* arg, forest* res);
        virtual ~copy_MT();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

    protected:

        /// Implement compute(), recursively.
        void _compute(int L, unsigned in,
                node_handle A, edge_value &cv, node_handle &cp);

        void traceout(const edge_value &v, node_handle p)
        {
#ifdef TRACE
            if (v.isVoid()) {
                out << p;
            } else {
                out << "<";
                v.show(out);
                out << ", " << p << ">";
            }
#endif
        }

    private:
        ct_entry_type* ct;

#ifdef TRACE
        ostream_output out;
#endif
};

// ******************************************************************

MEDDLY::copy_MT::copy_MT(forest* arg, forest* res)
    : unary_operation(arg, res)
#ifdef TRACE
      , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    if (!arg->isMultiTerminal()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    ct = new ct_entry_type("copy_mt");
    ct->setFixed(arg);
    if (res->isMultiTerminal()) {
        ct->setResult(res);
    } else {
        ct->setResult(res->getEdgeType(), res);
    }
    ct->doneBuilding();
}

MEDDLY::copy_MT::~copy_MT()
{
    ct->markForDestroy();
}

void MEDDLY::copy_MT::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    _compute(L, in, ap, cv, cp);
}

void MEDDLY::copy_MT::_compute(int L, unsigned in,
        node_handle A, edge_value& cv, node_handle &cp)
{
    //
    // Terminal case
    //
    if (argF->isTerminalNode(A)) {
        //
        // Translate A to result forest terminal
        //
        bool abool;
        int aint;
        float afloat;
        if (resF->isMultiTerminal()) {
            //
            // Result is MT
            //
            cv.set();
            switch (resF->getTerminalType()) {
                case terminal_type::BOOLEAN:
                    argF->getValueFromHandle(A, abool);
                    cp = resF->handleForValue(abool);
                    break;

                case terminal_type::INTEGER:
                    argF->getValueFromHandle(A, aint);
                    cp = resF->handleForValue(aint);
                    break;

                case terminal_type::REAL:
                    argF->getValueFromHandle(A, afloat);
                    cp = resF->handleForValue(afloat);
                    break;

                default:
                    FAIL(__FILE__, __LINE__, "Unknown terminal type");
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            } // switch
        } else {
            //
            // Translate A to result forest edge value
            //
            switch (resF->getEdgeType()) {
                case edge_type::INT:
                case edge_type::LONG:
                    argF->getValueFromHandle(A, aint);
                    resF->getEdgeForValue(aint, cv, cp);
                    break;

                case edge_type::FLOAT:
                case edge_type::DOUBLE:
                    argF->getValueFromHandle(A, afloat);
                    resF->getEdgeForValue(afloat, cv, cp);
                    break;

                default:
                    FAIL(__FILE__, __LINE__, "Unknown edge type");
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            } // switch
        } // MT vs EV
        //
        // Add nodes up to level L.
        //
        if (argF->isIdentityReduced()) {
            cp = resF->makeIdentitiesTo(cp, 0, L, in);
        } else {
            cp = resF->makeRedundantsTo(cp, 0, L);
        }
        return;
    } // A is terminal

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
#ifdef TRACE
    out << "copy_MT::_compute(" << A << ")\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(ct->getKeySize());
    ct_vector res(ct->getResultSize());
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        //
        // compute table 'hit'
        //
        if (resF->isMultiTerminal()) {
            cv.set();
            cp = resF->linkNode(res[0].getN());
        } else {
            res[0].get(cv);
            cp = resF->linkNode(res[1].getN());
        }
#ifdef TRACE
        out << "  CT hit ";
        traceout(cv, cp);
        out << "\n";
        out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
        //
        // done: compute table 'hit'
        //
    } else {
        //
        // compute table 'miss'; do computation
        //

        //
        // Initialize unpacked nodes
        //
        unpacked_node* Au = unpacked_node::newFromNode(argF, A, FULL_ONLY);
        unpacked_node* Cu = unpacked_node::newWritable(resF, Alevel, FULL_ONLY);
        MEDDLY_DCASSERT(Au->getSize() == Cu->getSize());
#ifdef TRACE
        out << "A: ";
        Au->show(out, true);
        out.indent_more();
        out.put('\n');
#endif
        const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

        //
        // Recurse over child edges
        //
        for (unsigned i=0; i<Cu->getSize(); i++) {
            edge_value v;
            node_handle d;
            _compute(Cnextlevel, i, Au->down(i), v, d);
            if (resF->isMultiTerminal()) {
                MEDDLY_DCASSERT(v.isVoid());
                Cu->setFull(i, d);
            } else {
                Cu->setFull(i, v, d);
            }
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "copy_MT::_compute(" << A << ") done\n";
        out << "A: ";
        Au->show(out, true);
        out << "\nC: ";
        Cu->show(out, true);
        out << "\n";
#endif

        //
        // Reduce
        //
        resF->createReducedNode(Cu, cv, cp);
#ifdef TRACE
        out << "reduced to ";
        traceout(cv, cp);
        out << ": ";
        resF->showNode(out, cp, SHOW_DETAILS);
        out << "\n";
#endif

        //
        // Add to CT
        //
        if (resF->isMultiTerminal()) {
            MEDDLY_DCASSERT(cv.isVoid());
            res[0].setN(cp);
        } else {
            res[0].set(cv);
            res[1].setN(cp);
        }
        ct->addCT(key, res);

        //
        // Cleanup
        //
        unpacked_node::Recycle(Au);
        // Cu is recycled when we reduce it :)

        //
        // done: compute table 'miss'
        //
    }

    //
    // Adjust result for singletons, added identities/redundants.
    // Do this for both CT hits and misses.
    //
    const int Clevel = resF->getNodeLevel(cp);
    if (Clevel == L) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        cp = resF->redirectSingleton(in, cp);
        return;
    }

    //
    // Add nodes but only from Alevel;
    // if Clevel is below Alevel it means nodes were eliminated
    // in the result forest.
    //
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, Alevel, L, in);
    } else {
        cp = resF->makeRedundantsTo(cp, Alevel, L);
    }
}

// ******************************************************************
// *                                                                *
// *                       copy_EV_fast class                       *
// *                                                                *
// ******************************************************************

/// Fast copy operation for EV: no 'push down' required.
/// Requires the same operation and no information loss.
///
class MEDDLY::copy_EV_fast : public unary_operation {
    public:
        copy_EV_fast(forest* arg, forest* res);
        virtual ~copy_EV_fast();

        virtual void compute(int L, unsigned in,
                const edge_value &av, node_handle ap,
                edge_value &cv, node_handle &cp);

    protected:

        // Recursive implementation
        void _compute(int L, unsigned in,
                node_handle ap, node_handle &cp);

        void traceout(const edge_value &v, node_handle p)
        {
#ifdef TRACE
            if (v.isVoid()) {
                out << p;
            } else {
                out << "<";
                v.show(out);
                out << ", " << p << ">";
            }
#endif
        }

    private:
        ct_entry_type* ct;

#ifdef TRACE
        ostream_output out;
#endif
};

// ******************************************************************

MEDDLY::copy_EV_fast::copy_EV_fast(forest* arg, forest* res)
    : unary_operation(arg, res)
#ifdef TRACE
      , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    if (arg->isMultiTerminal() || res->isMultiTerminal()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    ct = new ct_entry_type("copy_ev_fast");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::copy_EV_fast::~copy_EV_fast()
{
    ct->markForDestroy();
}

void MEDDLY::copy_EV_fast::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    MEDDLY_DCASSERT(!av.isVoid());
    _compute(L, in, ap, cp);
    cv = edge_value(resF->getEdgeType(), av);
}

void MEDDLY::copy_EV_fast::_compute(int L, unsigned in,
        node_handle A, node_handle &cp)
{
    //
    // Terminal case
    //
    if (argF->isTerminalNode(A)) {
        //
        // TBD: will we ever need to translate terminal nodes?
        //
        cp = A;
        //
        // Add nodes up to level L.
        //
        if (argF->isIdentityReduced()) {
            cp = resF->makeIdentitiesTo(cp, 0, L, in);
        } else {
            cp = resF->makeRedundantsTo(cp, 0, L);
        }
        return;
    } // A is terminal

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
#ifdef TRACE
    out << "copy_EV_fast::_compute(" << A << ")\n";
#endif


    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        //
        // compute table 'hit'
        //
        cp = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "  CT hit " << cp << "\n";
        out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
        //
        // done: compute table 'hit'
        //
    } else {
        //
        // compute table 'miss'; do computation
        //

        //
        // Initialize unpacked nodes
        //
        unpacked_node* Au = unpacked_node::newFromNode(argF, A, SPARSE_ONLY);
        unpacked_node* Cu = unpacked_node::newWritable(resF, Alevel, Au->getSize(), SPARSE_ONLY);
#ifdef TRACE
        out << "A: ";
        Au->show(out, true);
        out.indent_more();
        out.put('\n');
#endif
        const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

        //
        // Recurse over child edges
        //
        for (unsigned z=0; z<Cu->getSize(); z++) {
            node_handle d;
            _compute(Cnextlevel, Au->index(z), Au->down(z), d);
            edge_value v(resF->getEdgeType(), Au->edgeval(z));
            Cu->setSparse(z, Au->index(z), v, d);
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "copy_EV_fast::_compute(" << A << ") done\n";
        out << "A: ";
        Au->show(out, true);
        out << "\nC: ";
        Cu->show(out, true);
        out << "\n";
#endif

        //
        // Reduce
        //
        edge_value cv;
        resF->createReducedNode(Cu, cv, cp);
#ifdef TRACE
        out << "reduced to ";
        traceout(cv, cp);
        out << ": ";
        resF->showNode(out, cp, SHOW_DETAILS);
        out << "\n";
#endif
        MEDDLY_DCASSERT(cv == resF->getNoOpEdge());

        //
        // Add to CT
        //
        res[0].setN(cp);
        ct->addCT(key, res);

        //
        // Cleanup
        //
        unpacked_node::Recycle(Au);
        // Cu is recycled when we reduce it :)

        //
        // done: compute table 'miss'
        //
    }

    //
    // Adjust result for singletons, added identities/redundants.
    // Do this for both CT hits and misses.
    //
    const int Clevel = resF->getNodeLevel(cp);
    if (Clevel == L) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        cp = resF->redirectSingleton(in, cp);
        return;
    }

    //
    // Add nodes but only from Alevel;
    // if Clevel is below Alevel it means nodes were eliminated
    // in the result forest.
    //
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, Alevel, L, in);
    } else {
        cp = resF->makeRedundantsTo(cp, Alevel, L);
    }
}

// ******************************************************************
// *                                                                *
// *                         copy_EV  class                         *
// *                                                                *
// ******************************************************************

/// New, hopefully faster, copy operation for edge-valued DDs.
/// The target forest can be multi-terminal or edge-valued :)
/// This is the 'slow' version where we push values down.
///
namespace MEDDLY {

    template <class EdgeOp>
    class copy_EV : public unary_operation {
        public:
            copy_EV(forest* arg, forest* res);
            virtual ~copy_EV();

            virtual void compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    edge_value &cv, node_handle &cp);

        protected:

            // Recursive implementation
            void _compute(int L, unsigned in,
                    const edge_value &av, node_handle ap,
                    edge_value &cv, node_handle &cp);

            void traceout(const edge_value &v, node_handle p)
            {
#ifdef TRACE
                if (v.isVoid()) {
                    out << p;
                } else {
                    out << "<";
                    v.show(out);
                    out << ", " << p << ">";
                }
#endif
            }

        private:
            ct_entry_type* ct;

#ifdef TRACE
            ostream_output out;
#endif
    };

}; // namespace MEDDLY

// ******************************************************************

template <class EdgeOp>
MEDDLY::copy_EV<EdgeOp>::copy_EV(forest* arg, forest* res)
    : unary_operation(arg, res)
#ifdef TRACE
      , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__);
    if (arg->isMultiTerminal()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    ct = new ct_entry_type("copy_ev");
    ct->setFixed(arg->getEdgeType(), arg);
    if (res->isMultiTerminal()) {
        ct->setResult(res);
    } else {
        ct->setResult(res->getEdgeType(), res);
    }
    ct->doneBuilding();
}

template <class EdgeOp>
MEDDLY::copy_EV<EdgeOp>::~copy_EV()
{
    ct->markForDestroy();
}

template <class EdgeOp>
void MEDDLY::copy_EV<EdgeOp>::compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    _compute(L, in, av, ap, cv, cp);
}


template <class EdgeOp>
void MEDDLY::copy_EV<EdgeOp>::_compute(int L, unsigned in,
        const edge_value &av, node_handle ap, edge_value& cv, node_handle &cp)
{
    //
    // Terminal case
    //
    if (argF->isTerminalNode(ap)) {
        //
        // TBD: Check which terminal we reached
        //

        // if (OMEGA_INFINITY == ap) then what???

        if (resF->isMultiTerminal()) {
            //
            // Result is MT, copy edge value into terminal
            //
            bool abool;
            int aint;
            float afloat;
            cv.set();
            switch (resF->getTerminalType()) {
                case terminal_type::BOOLEAN:
                    av.copyInto(abool);
                    cp = resF->handleForValue(abool);
                    break;

                case terminal_type::INTEGER:
                    av.copyInto(aint);
                    cp = resF->handleForValue(aint);
                    break;

                case terminal_type::REAL:
                    av.copyInto(afloat);
                    cp = resF->handleForValue(afloat);
                    break;

                default:
                    FAIL(__FILE__, __LINE__, "Unknown terminal type");
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            } // switch
        } else {
            //
            // Result is EV, copy edge value over
            //
            int aint;
            long along;
            float afloat;
            double adouble;
            switch (resF->getEdgeType()) {
                case edge_type::INT:
                    av.copyInto(aint);
                    resF->getEdgeForValue(aint, cv, cp);
                    break;

                case edge_type::LONG:
                    av.copyInto(along);
                    resF->getEdgeForValue(along, cv, cp);
                    break;

                case edge_type::FLOAT:
                    av.copyInto(afloat);
                    resF->getEdgeForValue(afloat, cv, cp);
                    break;

                case edge_type::DOUBLE:
                    av.copyInto(adouble);
                    resF->getEdgeForValue(adouble, cv, cp);
                    break;

                default:
                    FAIL(__FILE__, __LINE__, "Unknown edge type");
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            } // switch
        }
        //
        // Stack nodes above terminal node as needed
        //
        if (argF->isIdentityReduced()) {
            cp = resF->makeIdentitiesTo(cp, 0, L, in);
        } else {
            cp = resF->makeRedundantsTo(cp, 0, L);
        }
        return;
    } // ap is terminal

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(ap);
#ifdef TRACE
    out << "copy_EV::_compute(";
    traceout(av, ap);
    out << ")\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(2);
    ct_vector res( resF->isMultiTerminal() ? 1 : 2 );
    key[0].set(av);
    key[1].setN(ap);
    if (ct->findCT(key, res)) {
        //
        // compute table 'hit'
        //
        if (resF->isMultiTerminal()) {
            cv.set();
            cp = resF->linkNode(res[0].getN());
        } else {
            res[0].get(cv);
            cp = resF->linkNode(res[1].getN());
        }
#ifdef TRACE
        out << "  CT hit ";
        traceout(cv, cp);
        out << "\n";
        out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
        //
        // done: compute table 'hit'
        //
    } else {
        //
        // compute table 'miss'; do computation
        //

        //
        // Initialize unpacked nodes
        //
        unpacked_node* Au = unpacked_node::newFromNode(argF, ap, FULL_ONLY);
        unpacked_node* Cu = unpacked_node::newWritable(resF, Alevel, FULL_ONLY);
        MEDDLY_DCASSERT(Au->getSize() == Cu->getSize());
#ifdef TRACE
        out << "A: ";
        Au->show(out, true);
        out.indent_more();
        out.put('\n');
#endif
        const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

        //
        // Recurse over child edges
        //
        for (unsigned i=0; i<Cu->getSize(); i++) {
            edge_value v;
            node_handle d;
            _compute(Cnextlevel, i, EdgeOp::applyOp(av, Au->edgeval(i)),
                    Au->down(i), v, d);
            if (resF->isMultiTerminal()) {
                MEDDLY_DCASSERT(v.isVoid());
                Cu->setFull(i, d);
            } else {
                Cu->setFull(i, v, d);
            }
        }

#ifdef TRACE
        out.indent_less();
        out.put('\n');
        out << "copy_EV::_compute(";
        traceout(av, ap);
        out << ") done\n";
        out << "A: ";
        Au->show(out, true);
        out << "\nC: ";
        Cu->show(out, true);
        out << "\n";
#endif

        //
        // Reduce
        //
        resF->createReducedNode(Cu, cv, cp);
#ifdef TRACE
        out << "reduced to ";
        traceout(cv, cp);
        out << ": ";
        resF->showNode(out, cp, SHOW_DETAILS);
        out << "\n";
#endif

        //
        // Add to CT
        //
        if (resF->isMultiTerminal()) {
            MEDDLY_DCASSERT(cv.isVoid());
            res[0].setN(cp);
        } else {
            res[0].set(cv);
            res[1].setN(cp);
        }
        ct->addCT(key, res);

        //
        // Cleanup
        //
        unpacked_node::Recycle(Au);
        // Cu is recycled when we reduce it :)

        //
        // done: compute table 'miss'
        //
    }

    //
    // Adjust result for singletons, added identities/redundants.
    // Do this for both CT hits and misses.
    //
    const int Clevel = resF->getNodeLevel(cp);
    if (Clevel == L) {
        //
        // We don't need to add identities or redundants;
        // just check if we need to avoid a singleton.
        //
        cp = resF->redirectSingleton(in, cp);
        return;
    }

    //
    // Add nodes but only from Alevel;
    // if Clevel is below Alevel it means nodes were eliminated
    // in the result forest.
    //
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, Alevel, L, in);
    } else {
        cp = resF->makeRedundantsTo(cp, Alevel, L);
    }
}

// ******************************************************************
// *                                                                *
// *                       COPY_factory class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::COPY_factory : public unary_factory {
    public:
        virtual void setup();
        virtual unary_operation* build(forest* arg, forest* res);
};

void MEDDLY::COPY_factory::setup()
{
    _setup("COPY", "Copy functions, within the same forest or across forests. The argument and result forests must have the same domain, and must both be sets or both be relations. The function ranges may be different.");
}

MEDDLY::unary_operation* MEDDLY::COPY_factory::build(forest* arg, forest* res)
{
    if (!arg || !res) {
        return nullptr;
    }

    unary_operation* uop =  find(arg, res);
    if (uop) {
        return uop;
    }

    if (arg->isForRelations() != res->isForRelations()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    //
    // Super fast case: same forests
    //
    if (arg == res) {
        return add(new copy_inforest(res));
    }

    //
    // Source is MT?
    //
    if (arg->isMultiTerminal())
    {
        return add(new copy_MT(arg, res));
    }


    //
    // "fast" EV to EV copies (same operation, no info loss)
    //
    if ( ((arg->isEVPlus() || arg->isIndexSet()) && res->isEVPlus())
            || (arg->isEVTimes() && res->isEVTimes()) )
    {
        switch (arg->getRangeType()) {
            case range_type::INTEGER:
                return add( new copy_EV_fast(arg, res) );

            case range_type::REAL:
                if (res->getRangeType() == range_type::REAL) {
                    return add( new copy_EV_fast(arg, res) );
                }

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);

        }; // switch
    }

    //
    // "slow" (push-down) EV source
    //
    if (arg->isEVPlus() || arg->isIndexSet()) {
        switch (arg->getRangeType()) {
            case range_type::INTEGER:
                return add(
                    new copy_EV< EdgeOp_plus<long> >(arg, res)
                );

            case range_type::REAL:
                return add(
                    new copy_EV< EdgeOp_plus<float> >(arg, res)
                );

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        };
    }

    if (arg->isEVTimes()) {
        switch (arg->getRangeType()) {
            case range_type::INTEGER:
                return add(
                    new copy_EV< EdgeOp_times<long> >(arg, res)
                );

            case range_type::REAL:
                return add(
                    new copy_EV< EdgeOp_times<float> >(arg, res)
                );

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        };
    }


    //
    // Catch all for any other cases
    //
    return nullptr;
}


// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_factory& MEDDLY::COPY()
{
    static COPY_factory F;
    return F;
}

