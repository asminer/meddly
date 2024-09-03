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
#include "complement.h"

#include "../ct_vector.h"
#include "../compute_table.h"
#include "../oper_unary.h"

namespace MEDDLY {
    class compl_mt;

    class compl_mdd;
    class compl_mxd;

    unary_list COMPL_cache;
};

// #define TRACE

#ifdef TRACE
#include "../operators.h"
#endif

// ******************************************************************
// *                                                                *
// *                         compl_mt class                         *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mt : public unary_operation {
    public:
        compl_mt(forest* arg, forest* res);
        virtual ~compl_mt();

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        /// Build complement of A, put result in C.
        /// The level of C will be the level of A, or lower.
        void _compute(node_handle A, node_handle &C);

        /*
            Build a complemented identity pattern:
                [ p 1 1 ... 1 ]
                [ 1 p 1 ... 1 ]
                [ 1 1 p ... 1 ]
                [ :         : ]
                [ 1 1 1 ... p ]


                @param  p   Diagonal pointer
                @param  K   Start the pattern above this level
                @param  L   Stop the pattern at this level
        */
        node_handle _identity_complement(node_handle p, int K, int L);

        inline node_handle identity_complement(node_handle p, int K, int L)
        {
            MEDDLY_DCASSERT(L>=0);
            MEDDLY_DCASSERT(K>=0);
            if (0==L) return p;
            if (K==L) return p;
            MEDDLY_DCASSERT(L>=K);
            return _identity_complement(p, K, L);
        }

    private:
        ct_entry_type* ct;

#ifdef TRACE
        ostream_output out;
#endif
};

// ******************************************************************

MEDDLY::compl_mt::compl_mt(forest* arg, forest* res)
    : unary_operation(arg, res)
#ifdef TRACE
      , out(std::cout)
#endif
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("mdd_complement");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::compl_mt::~compl_mt()
{
    ct->markForDestroy();
}

void MEDDLY::compl_mt::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
#ifdef TRACE
    out.indentation(0);
#endif
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    _compute(ap, cp);

    int aplevel = argF->getNodeLevel(ap);
    if (argF->isIdentityReduced()) {
        cp = identity_complement(cp, aplevel, L);
    } else {
        cp = resF->makeRedundantsTo(cp, aplevel, L);
    }
}

void MEDDLY::compl_mt::_compute(node_handle A, node_handle &C)
{
    //
    // Terminal cases.
    //
    if (argF->isTerminalNode(A)) {
        bool ta;
        argF->getValueFromHandle(A, ta);
        C = resF->handleForValue(!ta);
        return;
    }

    //
    // Determine level information
    //
#ifdef TRACE
    out << "compl_mt::_compute(" << A << ")\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        C = resF->linkNode(res[0].getN());
#ifdef TRACE
        out << "CT hit " << res[0].getN() << "\n";
        out << "  at level " << resF->getNodeLevel(C) << "\n";
#endif
        return;
    }

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = argF->newUnpacked(A, FULL_ONLY);
    const int Alevel = argF->getNodeLevel(A);
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
        node_handle d;
        _compute(Au->down(i), d);
        node_handle dc;
        if (argF->isIdentityReduced()) {
            dc = identity_complement(d, Audlevel, Cnextlevel);
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
        Cu->setFull(i, d);
    } // for i

#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "compl_mt::_compute(" << A << ") done\n";
    out << "A: ";
    Au->show(out, true);
    out << "\nC: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce
    //
    edge_value v;
    resF->createReducedNode(Cu, v, C);
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);

    //
    // Cleanup
    //
    unpacked_node::Recycle(Au);
}

MEDDLY::node_handle MEDDLY::compl_mt::_identity_complement(node_handle p,
    int K, int L)
{
    //
    // TBD: need to check for singleton edges at the bottom,
    // but only if K<0
    //
    // look at forest::makeIdentitiesTo
    //
    // maybe add this to forest? to build generic patterns
    //  [ p q q q ]
    //  [ q p q q ]
    //  [ q q p q ]
    //  [ q q q p ]
    //
    MEDDLY_DCASSERT(L>0);
    MEDDLY_DCASSERT(K>=0);
    unpacked_node* Uun;
    unpacked_node* Upr;
    edge_value ev;

    terminal ONE(true);
    node_handle chain_to_one = resF->makeRedundantsTo(ONE.getHandle(), 0, K);

    for (K++; K<=L; K++) {
        Uun = unpacked_node::newFull(resF, K, resF->getLevelSize(K));

        for (unsigned i=0; i<Uun->getSize(); i++) {
            Upr = unpacked_node::newFull(resF, -K, Uun->getSize());
            for (unsigned j=0; j<Uun->getSize(); j++) {
                Upr->setFull(j, resF->linkNode(chain_to_one));
            }
            Upr->setFull(i, (i ? resF->linkNode(p) : p));
            node_handle h;
            resF->createReducedNode(Upr, ev, h);
            MEDDLY_DCASSERT(ev.isVoid());
            Uun->setFull(i, h);
        }

        resF->createReducedNode(Uun, ev, p);
        MEDDLY_DCASSERT(ev.isVoid());

        chain_to_one = resF->makeRedundantsTo(chain_to_one, K-1, K);
    } // for k
    resF->unlinkNode(chain_to_one);
    return p;
}


// ******************************************************************
// *                                                                *
// *                        compl_mdd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mdd : public unary_operation {
    public:
        compl_mdd(forest* arg, forest* res);
        virtual ~compl_mdd();

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, int L);

    private:
        ct_entry_type* ct;
};

// ******************************************************************

MEDDLY::compl_mdd::compl_mdd(forest* arg, forest* res)
    : unary_operation(arg, res)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, SET);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("mdd_complement");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::compl_mdd::~compl_mdd()
{
    ct->markForDestroy();
}

void MEDDLY::compl_mdd::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    cp = _compute(ap, L);
}

MEDDLY::node_handle MEDDLY::compl_mdd::_compute(node_handle A, int L)
{
    //
    // Terminal cases
    //
    if (argF->isTerminalNode(A)) {
        bool ta;
        argF->getValueFromHandle(A, ta);
        return resF->makeRedundantsTo(resF->handleForValue(!ta), 0, L);
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        return resF->makeRedundantsTo(resF->linkNode(res[0].getN()), Alevel, L);
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = argF->newUnpacked(A, FULL_ONLY);
    unpacked_node* Cu = unpacked_node::newFull(resF, Alevel, Au->getSize());

    //
    // Build result node
    //
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i, _compute(Au->down(i), Cu->getLevel()-1));
    }

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);

    return resF->makeRedundantsTo(C, Alevel, L);
}

// ******************************************************************
// *                                                                *
// *                        compl_mxd  class                        *
// *                                                                *
// ******************************************************************

class MEDDLY::compl_mxd : public unary_operation {
    public:
        compl_mxd(forest* arg, forest* res);
        virtual ~compl_mxd();

        virtual void compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp);

    protected:
        node_handle _compute(node_handle A, int L);

        node_handle _compute_primed(int in, node_handle A, int L);

        /*
            Build a complemented identity pattern:
                [ p 1 1 ... 1 ]
                [ 1 p 1 ... 1 ]
                [ 1 1 p ... 1 ]
                [ :         : ]
                [ 1 1 1 ... p ]


                @param  p   Diagonal pointer
                @param  K   Start the pattern above this level
                @param  L   Stop the pattern at this level
        */
        node_handle _identity_complement(node_handle p, int K, int L);

    protected:
        inline unpacked_node* patternNode(forest* F, int L, int in,
                node_handle A, node_storage_flags fs)
        {
            MEDDLY_DCASSERT(L<0);
            if (0==A || F->isFullyReduced()) {
                return unpacked_node::newRedundant(F, L, A, fs);
            }
            if (F->isIdentityReduced()) {
                return unpacked_node::newIdentity(F, L, in, A, fs);
            }
            std::cout << "L: " << L << "\n";
            std::cout << "A: " << A << "\n";
            std::cout << "A.level: " << F->getNodeLevel(A) << "\n";
            MEDDLY_DCASSERT(false);
        }

    protected:
        inline node_handle identity_complement(node_handle p, int K, int L)
        {
            MEDDLY_DCASSERT(L>=0);
            MEDDLY_DCASSERT(K>=0);
            if (0==L) return p;
            if (K==L) return p;
            MEDDLY_DCASSERT(L>=K);
            return _identity_complement(p, K, L);
        }

    private:
        ct_entry_type* ct;

        ostream_output out;
};

// ******************************************************************

MEDDLY::compl_mxd::compl_mxd(forest* arg, forest* res)
    : unary_operation(arg, res), out(std::cout)
{
    checkDomains(__FILE__, __LINE__);
    checkAllRelations(__FILE__, __LINE__, RELATION);
    checkAllRanges(__FILE__, __LINE__, range_type::BOOLEAN);
    checkAllLabelings(__FILE__, __LINE__, edge_labeling::MULTI_TERMINAL);

    ct = new ct_entry_type("mxd_complement");
    ct->setFixed(arg);
    ct->setResult(res);
    ct->doneBuilding();
}

MEDDLY::compl_mxd::~compl_mxd()
{
    ct->markForDestroy();
}

void MEDDLY::compl_mxd::compute(const edge_value &av, node_handle ap,
                int L,
                edge_value &cv, node_handle &cp)
{
    MEDDLY_DCASSERT(av.isVoid());
    cv.set();
    out.indentation(0);
    cp = _compute(ap, L);
}

MEDDLY::node_handle MEDDLY::compl_mxd::_compute(node_handle A, int L)
{
    //
    // Terminal cases
    //
    if (0==A) {
        // Terminal zero;
        // that's a matrix of all zeroes.
        return resF->makeRedundantsTo(resF->handleForValue(true), 0, L);
    }

    if (A < 0) {
        // Terminal one; that's either a matrix of all ones
        // or an identity matrix.

        if (0==L || argF->isFullyReduced()) {
            return resF->makeRedundantsTo(resF->handleForValue(false), 0, L);
        }

        MEDDLY_DCASSERT(argF->isIdentityReduced());
        return identity_complement(0, 0, L);
    }

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
    const int Clevel = ABS(Alevel);

#ifdef TRACE
    out << "compl_mxd::_compute(" << A << ", " << L << ")\n";
    out << A << " level " << Alevel << "\n";
    out << "result level " << Clevel << " before chain\n";
#endif

    //
    // Check compute table
    //
    ct_vector key(1);
    ct_vector res(1);
    key[0].setN(A);
    if (ct->findCT(key, res)) {
        node_handle C = resF->linkNode(res[0].getN());
        if (argF->isIdentityReduced()) {
            C = identity_complement(C, Clevel, L);
        } else {
            C = resF->makeRedundantsTo(C, Clevel, L);
        }
#ifdef TRACE
        out << "CT hit " << res[0].getN() << "\n";
        out << "after chain " << C << "\n";
#endif
        return C;
    }

    //
    // Do computation
    //

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   unpacked_node::newRedundant(argF, Clevel, A, FULL_ONLY)
        :   argF->newUnpacked(A, FULL_ONLY);
    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i,
            _compute_primed(int(i), Au->down(i), MXD_levels::downLevel(Clevel))
        );
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "compl_mxd::_compute(" << A << ", " << L
              << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C);
    MEDDLY_DCASSERT(dummy.isVoid());
#ifdef TRACE
    out << "reduced to " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif

    //
    // Save result in CT
    //
    res[0].setN(C);
    ct->addCT(key, res);

    if (argF->isIdentityReduced()) {
        C = identity_complement(C, Clevel, L);
    } else {
        C = resF->makeRedundantsTo(C, Clevel, L);
    }

#ifdef TRACE
    out << "chain to level " << L << " = " << C << ": ";
    resF->showNode(out, C, SHOW_DETAILS);
    out << "\n";
#endif
    return C;
}

MEDDLY::node_handle MEDDLY::compl_mxd::_compute_primed(int in,
            node_handle A, const int Clevel)
{
    MEDDLY_DCASSERT(Clevel<0);

    //
    // Determine level information
    //
    const int Alevel = argF->getNodeLevel(A);
#ifdef TRACE
    out << "compl_mxd::_compute_primed(" << in << ", " << A << ", " << Clevel << ")\n";
    out << A << " level " << Alevel << "\n";
#endif

    //
    // Initialize unpacked nodes
    //
    unpacked_node* Au = (Alevel != Clevel)
        ?   patternNode(argF, Clevel, in, A, FULL_ONLY)
        :   argF->newUnpacked(A, FULL_ONLY);

    unpacked_node* Cu = unpacked_node::newFull(resF, Clevel, Au->getSize());

    //
    // Build result node
    //
#ifdef TRACE
    out.indent_more();
    out.put('\n');
#endif
    for (unsigned i=0; i<Cu->getSize(); i++) {
        Cu->setFull(i, _compute(Au->down(i), MXD_levels::downLevel(Clevel)));
    }
#ifdef TRACE
    out.indent_less();
    out.put('\n');
    out << "compl_mxd::_compute_primed(" << in << ", " << A << ", " << Clevel << ") done\n";
    out << "  A: ";
    Au->show(out, true);
    out << "\n  C: ";
    Cu->show(out, true);
    out << "\n";
#endif

    //
    // Reduce / cleanup
    //
    unpacked_node::Recycle(Au);
    edge_value dummy;
    node_handle C;
    resF->createReducedNode(Cu, dummy, C, in);
    MEDDLY_DCASSERT(dummy.isVoid());

#ifdef TRACE
    out << "compl_mxd::_compute_primed(" << in << ", " << A << ", "
              << Clevel << ") = " << C << "\n  ";
    resF->showNode(out, C, SHOW_DETAILS);
    out.put('\n');
#endif
    return C;
}

MEDDLY::node_handle MEDDLY::compl_mxd::_identity_complement(node_handle p,
    int K, int L)
{
    MEDDLY_DCASSERT(L>0);
    MEDDLY_DCASSERT(K>=0);
    unpacked_node* Uun;
    unpacked_node* Upr;
    edge_value ev;

    terminal ONE(true);
    node_handle chain_to_one = resF->makeRedundantsTo(ONE.getHandle(), 0, K);

    for (K++; K<=L; K++) {
        Uun = unpacked_node::newFull(resF, K, resF->getLevelSize(K));

        for (unsigned i=0; i<Uun->getSize(); i++) {
            Upr = unpacked_node::newFull(resF, -K, Uun->getSize());
            for (unsigned j=0; j<Uun->getSize(); j++) {
                Upr->setFull(j, resF->linkNode(chain_to_one));
            }
            Upr->setFull(i, (i ? resF->linkNode(p) : p));
            node_handle h;
            resF->createReducedNode(Upr, ev, h);
            MEDDLY_DCASSERT(ev.isVoid());
            Uun->setFull(i, h);
        }

        resF->createReducedNode(Uun, ev, p);
        MEDDLY_DCASSERT(ev.isVoid());

        chain_to_one = resF->makeRedundantsTo(chain_to_one, K-1, K);
    } // for k
    resF->unlinkNode(chain_to_one);
    return p;
}

// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::unary_operation* MEDDLY::COMPLEMENT(forest* arg, forest* res)
{
    if (!arg || !res) return nullptr;

    unary_operation* uop = COMPL_cache.find(arg, res);
    if (uop) {
        return uop;
    }

    if (arg->isForRelations()) {
        return COMPL_cache.add(new compl_mxd(arg, res));
    } else {
        return COMPL_cache.add(new compl_mdd(arg, res));
    }
}

void MEDDLY::COMPLEMENT_init()
{
    COMPL_cache.reset("Complement");
}

void MEDDLY::COMPLEMENT_done()
{
    MEDDLY_DCASSERT( COMPL_cache.isEmpty() );
}

