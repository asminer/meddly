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

#include "../ct_entry_key.h"        // remove when we can
#include "../ct_entry_result.h"     // remove when we can
#include "../compute_table.h"       // remove when we can

#include "../oper_unary.h"
#include "../ct_vector.h"

namespace MEDDLY {
    class copy_inforest;
    class copy_MT;
    class copy_EV_fast;

    unary_list COPY_cache;
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

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
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

void MEDDLY::copy_inforest::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
    cv = av;
    ap = resF->linkNode(cp);
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

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp);

    protected:

        /*
           Recursive copy.

           This will correctly build a copy at the same level as
           node A, but in the result forest.
           It is the caller's responsibility to add any nodes
           above this one, or to check if it is a singleton node
           (if the target forest is identity reduced).

           @param   A   Source node

           @param   cv  on output: edge value for copy
           @param   cp  on output: target node for copy
        */
        void _compute(node_handle A, edge_value &cv, node_handle &cp);

        void traceout(const edge_value &v, node_handle p) const
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
    if (res->isForRelations()) {
        checkAllRelations(__FILE__, __LINE__, RELATION);
    } else {
        checkAllRelations(__FILE__, __LINE__, SET);
    }
    if (!arg->isMultiTerminal()) {
        throw error(error::TYPE_MISMATCH);
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

void MEDDLY::copy_MT::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    _compute(ap, cv, cp);

    unsigned aplevel = argF->getNodeLevel(ap);
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, aplevel, L, -1);
    } else {
        cp = resF->makeRedundantsTo(cp, aplevel, L);
    }
}

void MEDDLY::copy_MT::_compute(node_handle A, edge_value& cv, node_handle &cp)
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
                    return;

                case terminal_type::INTEGER:
                    argF->getValueFromHandle(A, aint);
                    cp = resF->handleForValue(aint);
                    return;

                case terminal_type::REAL:
                    argF->getValueFromHandle(A, afloat);
                    cp = resF->handleForValue(afloat);
                    return;

                default:
                    MEDDLY_DCASSERT(false);
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            } // switch
        } else {
            //
            // Result is EV
            //
            switch (resF->getEdgeType()) {
                case edge_type::INT:
                case edge_type::LONG:
                    argF->getValueFromHandle(A, aint);
                    resF->getEdgeForValue(aint, cv, cp);
                    return;

                case edge_type::FLOAT:
                case edge_type::DOUBLE:
                    argF->getValueFromHandle(A, afloat);
                    resF->getEdgeForValue(afloat, cv, cp);
                    return;

                default:
                    MEDDLY_DCASSERT(false);
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            } // switch
        } // MT vs EV
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
    ct_vector key(1);
    ct_vector res( resF->isMultiTerminal() ? 1 : 2 );
    key[0].setN(A);
    if (ct->findCT(key, res)) {
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
        return;
    }

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = argF->newUnpacked(A, FULL_ONLY);
    unpacked_node* Cu = unpacked_node::newFull(resF, Alevel, Au->getSize());
#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
#endif

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

    for (unsigned i=0; i<Cu->getSize(); i++) {
        int Audlevel = argF->getNodeLevel(Au->down(i));
        edge_value v;
        node_handle d;
        _compute(Au->down(i), v, d);
        node_handle dc;
        if (argF->isIdentityReduced()) {
            dc = resF->makeIdentitiesTo(d, Audlevel, Cnextlevel, i);
        } else {
            dc = resF->makeRedundantsTo(d, Audlevel, Cnextlevel);
        }
#ifdef TRACE
        if (dc != d) {
            out << "built chain from " << d << " to " << dc << "\n";
        }
#endif
        d = resF->redirectSingleton(
#ifdef DEVELOPMENT_CODE
                Cu->getLevel(),
#endif
                i, dc
        );
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
    << ": ";
    resF->showNode(out, cp, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
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

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp);

    protected:

        /*
           Recursive copy.

           This will correctly build a copy at the same level as
           node A, but in the result forest.
           It is the caller's responsibility to add any nodes
           above this one, or to check if it is a singleton node
           (if the target forest is identity reduced).

           @param   A   Source node

           @param   cv  on output: edge value for copy
           @param   cp  on output: target node for copy
        */
        void _compute(node_handle A, node_handle &cp);

        void traceout(const edge_value &v, node_handle p) const
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
    if (res->isForRelations()) {
        checkAllRelations(__FILE__, __LINE__, RELATION);
    } else {
        checkAllRelations(__FILE__, __LINE__, SET);
    }
    if (arg->isMultiTerminal() || res->isMultiTerminal()) {
        throw error(error::TYPE_MISMATCH);
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

void MEDDLY::copy_EV_fast::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    MEDDLY_DCASSERT(!av.isVoid());
    _compute(ap, cp);

    unsigned aplevel = argF->getNodeLevel(ap);
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, aplevel, L, -1);
    } else {
        cp = resF->makeRedundantsTo(cp, aplevel, L);
    }

    cv = edge_value(resF->getEdgeType(), av);
}

void MEDDLY::copy_EV_fast::_compute(node_handle A, node_handle &cp)
{
    //
    // Terminal case
    //
    if (argF->isTerminalNode(A)) {
        //
        // TBD: will we ever need to translate terminal nodes?
        //
        cp = A;
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
        cp = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "  CT hit ";
        traceout(cv, cp);
        out << "\n";
        out << "  at level " << resF->getNodeLevel(cp) << "\n";
#endif
        return;
    }

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = argF->newUnpacked(A, SPARSE_ONLY);
    unpacked_node* Cu = unpacked_node::newSparse(resF, Alevel, Au->getSize());
#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
#endif

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

    for (unsigned z=0; z<Cu->getSize(); z++) {
        int Audlevel = argF->getNodeLevel(Au->down(z));
        node_handle d;
        _compute(Au->down(z), d);
        node_handle dc;
        if (argF->isIdentityReduced()) {
            dc = resF->makeIdentitiesTo(d, Audlevel, Cnextlevel, Au->index(z));
        } else {
            dc = resF->makeRedundantsTo(d, Audlevel, Cnextlevel);
        }
#ifdef TRACE
        if (dc != d) {
            out << "built chain from " << d << " to " << dc << "\n";
        }
#endif
        d = resF->redirectSingleton(
#ifdef DEVELOPMENT_CODE
                Cu->getLevel(),
#endif
                Au->index(z), dc
        );
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
    << ": ";
    resF->showNode(out, cp, SHOW_DETAILS);
    out << "\n";
#endif
    MEDDLY_DCASSERT(cv == resF->getNoOpEdge());

    //
    // Save result in CT
    //
    res[0].setN(cp);
    ct->addCT(key, res);

    //
    // Cleanup
    //
    unpacked_node::Recycle(Au);
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

            virtual void compute(const edge_value &av, node_handle ap,
                    int L,
                    edge_value &cv, node_handle &cp);

        protected:

            /*
               Recursive copy.

               This will correctly build a copy at the same level as
               node A, but in the result forest.
               It is the caller's responsibility to add any nodes
               above this one, or to check if it is a singleton node
               (if the target forest is identity reduced).

               @param   av  source edge value
               @param   ap  source node

               @param   cv  on output: edge value for copy
               @param   cp  on output: target node for copy
               */
            void _compute(edge_value av, node_handle ap,
                    edge_value &cv, node_handle &cp);

            void traceout(const edge_value &v, node_handle p) const
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
    if (res->isForRelations()) {
        checkAllRelations(__FILE__, __LINE__, RELATION);
    } else {
        checkAllRelations(__FILE__, __LINE__, SET);
    }
    if (arg->isMultiTerminal()) {
        throw error(error::TYPE_MISMATCH);
    }

    ct = new ct_entry_type("copy_mt");
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
void MEDDLY::copy_EV<EdgeOp>::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    _compute(av, ap, cv, cp);

    unsigned aplevel = argF->getNodeLevel(ap);
    if (argF->isIdentityReduced()) {
        cp = resF->makeIdentitiesTo(cp, aplevel, L, -1);
    } else {
        cp = resF->makeRedundantsTo(cp, aplevel, L);
    }
}

template <class EdgeOp>
void MEDDLY::copy_EV<EdgeOp>::_compute(edge_value av, node_handle ap,
        edge_value& cv, node_handle &cp)
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
                    return;

                case terminal_type::INTEGER:
                    av.copyInto(aint);
                    cp = resF->handleForValue(aint);
                    return;

                case terminal_type::REAL:
                    av.copyInto(afloat);
                    cp = resF->handleForValue(afloat);
                    return;

                default:
                    MEDDLY_DCASSERT(false);
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
                    return;

                case edge_type::LONG:
                    av.copyInto(along);
                    resF->getEdgeForValue(along, cv, cp);
                    return;

                case edge_type::FLOAT:
                    av.copyInto(afloat);
                    resF->getEdgeForValue(afloat, cv, cp);
                    return;

                case edge_type::DOUBLE:
                    av.copyInto(adouble);
                    resF->getEdgeForValue(adouble, cv, cp);
                    return;

                default:
                    MEDDLY_DCASSERT(false);
                    throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
            } // switch
        }

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
        return;
    }


    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = argF->newUnpacked(ap, FULL_ONLY);
    unpacked_node* Cu = unpacked_node::newFull(resF, Alevel, Au->getSize());
#ifdef TRACE
    out << "A: ";
    Au->show(out, true);
#endif

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    const int Cnextlevel = resF->isForRelations()
            ? MXD_levels::downLevel(Alevel)
            : MDD_levels::downLevel(Alevel);

    for (unsigned i=0; i<Cu->getSize(); i++) {
        int Audlevel = argF->getNodeLevel(Au->down(i));
        edge_value v;
        node_handle d;
        _compute(EdgeOp::accumulate(av, Au->edgeval(i)), Au->down(i), v, d);
        node_handle dc;
        if (argF->isIdentityReduced()) {
            dc = resF->makeIdentitiesTo(d, Audlevel, Cnextlevel, i);
        } else {
            dc = resF->makeRedundantsTo(d, Audlevel, Cnextlevel);
        }
#ifdef TRACE
        if (dc != d) {
            out << "built chain from " << d << " to " << dc << "\n";
        }
#endif
        d = resF->redirectSingleton(
#ifdef DEVELOPMENT_CODE
                Cu->getLevel(),
#endif
                i, dc
        );
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
    << ": ";
    resF->showNode(out, cp, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
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
}



// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_operation* MEDDLY::COPY(forest* arg, forest* res)
{
    if (!arg || !res) {
        return nullptr;
    }

    unary_operation* uop =  COPY_cache.find(arg, res);
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
        return COPY_cache.add(new copy_inforest(res));
    }

    //
    // Source is MT?
    //
    if (arg->isMultiTerminal())
    {
        return COPY_cache.add(new copy_MT(arg, res));
    }


    //
    // "fast" EV to EV copies (same operation, no info loss)
    //
    if ( ((arg->isEVPlus() || arg->isIndexSet()) && res->isEVPlus())
            || (arg->isEVTimes() && res->isEVTimes()) )
    {
        switch (arg->getRangeType()) {
            case range_type::INTEGER:
                return COPY_cache.add( new copy_EV_fast(arg, res) );

            case range_type::REAL:
                if (res->getRangeType() == range_type::REAL) {
                    return COPY_cache.add( new copy_EV_fast(arg, res) );
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
                return COPY_cache.add(
                    new copy_EV< EdgeOp_plus<long> >(arg, res)
                );

            case range_type::REAL:
                return COPY_cache.add(
                    new copy_EV< EdgeOp_plus<float> >(arg, res)
                );

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        };
    }

    if (arg->isEVTimes()) {
        switch (arg->getRangeType()) {
            case range_type::INTEGER:
                return COPY_cache.add(
                    new copy_EV< EdgeOp_times<long> >(arg, res)
                );

            case range_type::REAL:
                return COPY_cache.add(
                    new copy_EV< EdgeOp_times<float> >(arg, res)
                );

            default:
                throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
        };
    }


    //
    // Catch all for any other cases
    //
    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

void MEDDLY::COPY_init()
{
    COPY_cache.reset("Copy");
}

void MEDDLY::COPY_done()
{
    MEDDLY_DCASSERT(COPY_cache.isEmpty());
}

