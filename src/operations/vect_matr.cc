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
#include "vect_matr.h"

#include "../forest.h"
#include "../forest_levels.h"
#include "../oper_numer.h"

namespace MEDDLY {
    class base_evplus_mt;

    class VM_evplus_mt;
    class MV_evplus_mt;

    inline bool isEvPlusStyle(const forest* f) {
        return f->isEVPlus() || f->isIndexSet();
    }
};

// ******************************************************************
// *                                                                *
// *                      base_evplus_mt class                      *
// *                                                                *
// ******************************************************************

class MEDDLY::base_evplus_mt : public numerical_operation {
    public:
        base_evplus_mt(const char* name, const dd_edge &x_ind,
                const dd_edge& A, const dd_edge &y_ind);

        virtual ~base_evplus_mt();

        virtual void compute(double* y, const double* x);

        virtual void compute_r(int ht, double* y, node_handle y_ind,
                const double* x, node_handle x_ind, node_handle A) = 0;

    protected:
        const forest* fx;
        const forest* fA;
        const forest* fy;
        node_handle x_root;
        node_handle A_root;
        node_handle y_root;
        int L;

        inline virtual bool checkForestCompatibility() const
        {
            auto o1 = fx->variableOrder();
            auto o2 = fA->variableOrder();
            auto o3 = fy->variableOrder();
            return o1->is_compatible_with(*o2) && o1->is_compatible_with(*o3);
        }
};

// ******************************************************************

MEDDLY::base_evplus_mt::base_evplus_mt(const char* name,
  const dd_edge &x_ind, const dd_edge& A, const dd_edge &y_ind)
 : numerical_operation(name)
{
    fx = x_ind.getForest();
    fA = A.getForest();
    fy = y_ind.getForest();
    MEDDLY_DCASSERT(fx);
    MEDDLY_DCASSERT(fA);
    MEDDLY_DCASSERT(fy);
    // everyone must use the same domain
    if  (       (fx->getDomain() != fy->getDomain())
            ||  (fx->getDomain() != fA->getDomain())  )
    {
        throw error(error::DOMAIN_MISMATCH, __FILE__, __LINE__);
    }

    // Check edge types
    if (
            (fy->getRangeType() != range_type::INTEGER)
            || (fy->isForRelations())
            || (fx->getRangeType() != range_type::INTEGER)
            || (fx->isForRelations())
            || (fA->getRangeType() != range_type::REAL)
            || (!fA->isForRelations())
       )
    {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    // A can't be fully reduced.
    if (fA->isFullyReduced()) {
        throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    }

    // For now, fy and fx must be Indexed sets or EVPLUS forests.
    if ( !isEvPlusStyle(fy) || !isEvPlusStyle(fx) ) {
        throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
    }

    x_root = x_ind.getNode();
    A_root = A.getNode();
    y_root = y_ind.getNode();
    L = fx->getMaxLevelIndex();
}

MEDDLY::base_evplus_mt::~base_evplus_mt()
{
}

void MEDDLY::base_evplus_mt::compute(double* y, const double* x)
{
    if (!checkForestCompatibility()) {
        throw error(error::INVALID_OPERATION, __FILE__, __LINE__);
    }
    compute_r(L, y, y_root, x, x_root, A_root);
}

// ******************************************************************
// *                                                                *
// *                       VM_evplus_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::VM_evplus_mt : public base_evplus_mt {
    public:
        VM_evplus_mt(const char* name, const dd_edge &x_ind,
                const dd_edge& A, const dd_edge &y_ind);

        virtual void compute_r(int k, double* y, node_handle y_ind,
                const double* x, node_handle x_ind, node_handle A);

        void comp_pr(int k, double* y, node_handle y_ind, const double* x,
                node_handle x_ind, unsigned ain, node_handle A);

};

MEDDLY::VM_evplus_mt::VM_evplus_mt(const char* name,
        const dd_edge &x_ind, const dd_edge& A, const dd_edge &y_ind)
    : base_evplus_mt(name, x_ind, A, y_ind)
{
}

void MEDDLY::VM_evplus_mt::compute_r(int k, double* y, node_handle y_ind,
    const double* x, node_handle x_ind, node_handle a)
{
    // Handles the unprimed levels of a
    if (0==k) {
        y[0] += x[0] * fA->getRealFromHandle(a);
        return;
    }

    // It should be impossible for an indexing function to skip levels, right?
    MEDDLY_DCASSERT(fx->getNodeLevel(x_ind) == k);
    MEDDLY_DCASSERT(fy->getNodeLevel(y_ind) == k);
    int aLevel = fA->getNodeLevel(a);

    //
    // A is identity matrix times a constant; exploit that if we can
    //
    if (0==aLevel && (x_ind == y_ind)) {
        if (fx == fy && fx->isIndexSet()) {
            // yes we can
            float v = fA->getRealFromHandle(a);
            for (long i = fx->getIndexSetCardinality(x_ind)-1; i>=0; i--) {
                y[i] += x[i] * v;
            }
            return;
        }
    }

    //
    // Check if A is an identity node
    //
    if (ABS(aLevel) < k) {
        // Init sparse readers
        unpacked_node* xR = unpacked_node::newFromNode(fx, x_ind, SPARSE_ONLY);
        unpacked_node* yR = unpacked_node::newFromNode(fy, y_ind, SPARSE_ONLY);

        unsigned xp = 0;
        unsigned yp = 0;
        for (;;) {
            if (xR->index(xp) < yR->index(yp)) {
                xp++;
                if (xp >= xR->getSize()) break;
                continue;
            }
            if (xR->index(xp) > yR->index(yp)) {
                yp++;
                if (yp >= yR->getSize()) break;
                continue;
            }
            // match, need to recurse
            compute_r(k-1, y + yR->edgeval(yp).getLong(), yR->down(yp),
                        x + xR->edgeval(xp).getLong(), xR->down(xp), a);
            xp++;
            if (xp >= xR->getSize()) break;
            yp++;
            if (yp >= yR->getSize()) break;
        } // for (;;)

        // Cleanup
        unpacked_node::Recycle(yR);
        unpacked_node::Recycle(xR);

        // Done
        return;
    }

    //
    // A is not an identity node.
    //

    // Init sparse readers
    unpacked_node* aR = unpacked_node::New(fA, SPARSE_ONLY);
    if (aLevel == k) {
        aR->initFromNode(a);
    } else {
        aR->initRedundant(k, a);
    }

    unpacked_node* xR = unpacked_node::newFromNode(fx, x_ind, SPARSE_ONLY);

    unsigned xp = 0;
    unsigned ap = 0;
    for (;;) {
        if (aR->index(ap) < xR->index(xp)) {
            ap++;
            if (ap >= aR->getSize()) break;
            continue;
        }
        if (aR->index(ap) > xR->index(xp)) {
            xp++;
            if (xp >= xR->getSize()) break;
            continue;
        }
        // match, need to recurse
        comp_pr(k, y, y_ind, x + xR->edgeval(xp).getLong(),
                    xR->down(xp), aR->index(ap), aR->down(ap));
        ap++;
        if (ap >= aR->getSize()) break;
        xp++;
        if (xp >= xR->getSize()) break;
    } // for (;;)

    // Cleanup
    unpacked_node::Recycle(xR);
    unpacked_node::Recycle(aR);
}

void MEDDLY::VM_evplus_mt::comp_pr(int k, double* y, node_handle y_ind,
        const double* x, node_handle x_ind, unsigned ain, node_handle a)
{
    // Handles the primed levels of A
    if (0==k) {
        y[0] += x[0] * fA->getRealFromHandle(a);
        return;
    }

    // Init sparse readers
    unpacked_node* aR = unpacked_node::New(fA, SPARSE_ONLY);
    if (fA->getNodeLevel(a) == -k) {
        aR->initFromNode(a);
    } else {
        aR->initIdentity(k, ain, a);
    }

    unpacked_node* yR = unpacked_node::newFromNode(fy, y_ind, SPARSE_ONLY);


    unsigned yp = 0;
    unsigned ap = 0;
    for (;;) {
        if (aR->index(ap) < yR->index(yp)) {
            ap++;
            if (ap >= aR->getSize()) break;
            continue;
        }
        if (aR->index(ap) > yR->index(yp)) {
            yp++;
            if (yp >= yR->getSize()) break;
            continue;
        }
        // match, need to recurse
        compute_r(k-1, y + yR->edgeval(yp).getLong(), yR->down(yp),
                    x, x_ind, aR->down(ap));
        ap++;
        if (ap >= aR->getSize()) break;
        yp++;
        if (yp >= yR->getSize()) break;
    } // for (;;)

    // Cleanup
    unpacked_node::Recycle(yR);
    unpacked_node::Recycle(aR);
}


// ******************************************************************
// *                                                                *
// *                       MV_evplus_mt class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::MV_evplus_mt : public base_evplus_mt {
    public:
        MV_evplus_mt(const char* name, const dd_edge &x_ind,
                const dd_edge& A, const dd_edge &y_ind);

        virtual void compute_r(int k, double* y, node_handle y_ind,
                const double* x, node_handle x_ind, node_handle A);

        void comp_pr(int k, double* y, node_handle y_ind, const double* x,
                node_handle x_ind, unsigned ain, node_handle A);

};

MEDDLY::MV_evplus_mt::MV_evplus_mt(const char* name,
        const dd_edge &x_ind, const dd_edge& A, const dd_edge &y_ind)
    : base_evplus_mt(name, x_ind, A, y_ind)
{
}

void MEDDLY::MV_evplus_mt::compute_r(int k, double* y, node_handle y_ind,
        const double* x, node_handle x_ind, node_handle a)
{
    // Handles the unprimed levels of a
    if (0==k) {
        y[0] += x[0] * fA->getRealFromHandle(a);
        return;
    }

    // It should be impossible for an indexing function to skip levels, right?
    MEDDLY_DCASSERT(fx->getNodeLevel(x_ind) == k);
    MEDDLY_DCASSERT(fy->getNodeLevel(y_ind) == k);
    int aLevel = fA->getNodeLevel(a);

    //
    // A is identity matrix times a constant; exploit that if we can
    //
    if (0==aLevel && (x_ind == y_ind)) {
        if (fx == fy && fx->isIndexSet()) {
            // yes we can
            float v = fA->getRealFromHandle(a);
            for (long i = fy->getIndexSetCardinality(y_ind)-1; i>=0; i--) {
                y[i] += x[i] * v;
            }
            return;
        }
    }

    //
    // Check if a is an identity node
    //
    if (ABS(aLevel) < k) {
        // Init sparse readers
        unpacked_node* xR = unpacked_node::newFromNode(fx, x_ind, SPARSE_ONLY);
        unpacked_node* yR = unpacked_node::newFromNode(fy, y_ind, SPARSE_ONLY);

        unsigned xp = 0;
        unsigned yp = 0;
        for (;;) {
            if (xR->index(xp) < yR->index(yp)) {
                xp++;
                if (xp >= xR->getSize()) break;
                continue;
            }
            if (xR->index(xp) > yR->index(yp)) {
                yp++;
                if (yp >= yR->getSize()) break;
                continue;
            }
            // match, need to recurse
            compute_r(k-1, y + yR->edgeval(yp).getLong(), yR->down(yp),
                        x + xR->edgeval(xp).getLong(), xR->down(xp), a);
            xp++;
            if (xp >= xR->getSize()) break;
            yp++;
            if (yp >= yR->getSize()) break;
        } // for (;;)

        // Cleanup
        unpacked_node::Recycle(yR);
        unpacked_node::Recycle(xR);

        // Done
        return;
    }

    //
    // A is not an identity node.
    //

    // Init sparse readers
    unpacked_node* aR = unpacked_node::New(fA, SPARSE_ONLY);
    if (aLevel == k) {
        aR->initFromNode(a);
    } else {
        aR->initRedundant(k, a);
    }

    unpacked_node* yR = unpacked_node::newFromNode(fy, y_ind, SPARSE_ONLY);


    unsigned yp = 0;
    unsigned ap = 0;
    for (;;) {
        if (aR->index(ap) < yR->index(yp)) {
            ap++;
            if (ap >= aR->getSize()) break;
            continue;
        }
        if (aR->index(ap) > yR->index(yp)) {
            yp++;
            if (yp >= yR->getSize()) break;
            continue;
        }
        // match, need to recurse
        comp_pr(k, y + yR->edgeval(yp).getLong(), yR->down(yp),
                    x, x_ind, aR->index(ap), aR->down(ap));
        ap++;
        if (ap >= aR->getSize()) break;
        yp++;
        if (yp >= yR->getSize()) break;
    } // for (;;)

    // Cleanup
    unpacked_node::Recycle(yR);
    unpacked_node::Recycle(aR);
}

void MEDDLY::MV_evplus_mt::comp_pr(int k, double* y, node_handle y_ind,
        const double* x, node_handle x_ind, unsigned ain, node_handle a)
{
    // Handles the primed levels of A
    if (0==k) {
        y[0] += x[0] * fA->getRealFromHandle(a);
        return;
    }

    // Init sparse readers
    unpacked_node* aR = unpacked_node::New(fA, SPARSE_ONLY);
    if (fA->getNodeLevel(a) == -k) {
        aR->initFromNode(a);
    } else {
        aR->initIdentity(k, ain, a);
    }

    unpacked_node* xR = unpacked_node::newFromNode(fx, x_ind, SPARSE_ONLY);


    unsigned xp = 0;
    unsigned ap = 0;
    for (;;) {
        if (aR->index(ap) < xR->index(xp)) {
            ap++;
            if (ap >= aR->getSize()) break;
            continue;
        }
        if (aR->index(ap) > xR->index(xp)) {
            xp++;
            if (xp >= xR->getSize()) break;
            continue;
        }
        // match, need to recurse
        compute_r(k-1, y, y_ind,
                    x + xR->edgeval(xp).getLong(), xR->down(xp), aR->down(ap));
        ap++;
        if (ap >= aR->getSize()) break;
        xp++;
        if (xp >= xR->getSize()) break;
    } // for (;;)

    // Cleanup
    unpacked_node::Recycle(xR);
    unpacked_node::Recycle(aR);
}



// ******************************************************************
// *                                                                *
// *                           Front  end                           *
// *                                                                *
// ******************************************************************

MEDDLY::numerical_operation* MEDDLY::EXPLVECT_MATR_MULT(const dd_edge &xind,
        const dd_edge &A, const dd_edge &yind)
{
    const forest* fA = A.getForest();
    switch (fA->getEdgeLabeling()) {
        case edge_labeling::MULTI_TERMINAL:
            return new VM_evplus_mt("VectMatrMult", xind, A, yind);

        case edge_labeling::EVTIMES:
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    };
}

MEDDLY::numerical_operation* MEDDLY::MATR_EXPLVECT_MULT(const dd_edge &xind,
        const dd_edge &A, const dd_edge &yind)
{
    const forest* fA = A.getForest();
    switch (fA->getEdgeLabeling()) {
        case edge_labeling::MULTI_TERMINAL:
            return new MV_evplus_mt("MatrVectMult", xind, A, yind);

        case edge_labeling::EVTIMES:
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);

        default:
            throw error(error::TYPE_MISMATCH, __FILE__, __LINE__);
    };
}
