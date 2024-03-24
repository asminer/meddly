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

#include "oper.h"

#include "ct_entry_result.h"
#include "compute_table.h"
#include "ct_initializer.h"

// ******************************************************************
// *                       operation  statics                       *
// ******************************************************************

MEDDLY::compute_table* MEDDLY::operation::Monolithic_CT;
MEDDLY::operation** MEDDLY::operation::op_list;
unsigned* MEDDLY::operation::op_holes;
unsigned MEDDLY::operation::list_size;
unsigned MEDDLY::operation::list_alloc;
unsigned MEDDLY::operation::free_list;

// ******************************************************************
// *                                                                *
// *                       operation  methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::operation::operation(opname* n, unsigned et_slots)
{
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Creating operation %p\n", this);
    fflush(stdout);
#endif
    theOpName = n;
    num_etids = et_slots;

    is_marked_for_deletion = false;
    next = nullptr;

    registerOperation(*this);

    //
    // Delay CT initialization!
    // The derived class hasn't set up the entry types yet!
    //
    CT = 0;

    //
    // Set up slots to save our entry_types.
    //
    if (et_slots) {
        etype = new ct_entry_type* [et_slots];
        for (unsigned i=0; i<et_slots; i++) {
            etype[i] = 0;
        }
    } else {
        etype = 0;
    }

    //
    // Allocate CTresults
    //
    if (et_slots) {
        CTresult = new ct_entry_result [et_slots];
    } else {
        CTresult = 0;
    }

    //
    // Allocate our slots
    //
    compute_table::registerOp(this, et_slots);
}



MEDDLY::operation::~operation()
{
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Deleting operation %p %s\n", this, getName());
    fflush(stdout);
#endif

    if (CT) {
        for (unsigned i=0; i<num_etids; i++) {
            if (CT[i] != Monolithic_CT) {
                delete CT[i];
            }
        }
        delete[] CT;
        CT = nullptr;
    }
    // Don't delete the entries in etype, they're owned by compute_table.
    delete[] etype;
    delete[] CTresult;
    compute_table::unregisterOp(this, num_etids);

    unregisterOperation(*this);
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Deleted operation %p %s\n", this, getName());
    fflush(stdout);
#endif
}


void MEDDLY::operation::removeStaleComputeTableEntries()
{
    bool has_monolithic = false;
    if (CT) {
        for (unsigned i=0; i<num_etids; i++) {
            if (0==CT[i]) continue;
            if (CT[i]->isOperationTable()) {
                CT[i]->removeStales();
            } else {
                has_monolithic = true;
            }
        }
    }
    if (has_monolithic) {
        Monolithic_CT->removeStales();
    }
}

void MEDDLY::operation::removeAllComputeTableEntries()
{
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Removing entries for operation %p %s\n", this, getName());
    fflush(stdout);
#endif
    if (is_marked_for_deletion) return;
    is_marked_for_deletion = true;
    for (unsigned i=0; i<num_etids; i++) {
        etype[i]->markForDeletion();
    }
    removeStaleComputeTableEntries();
    for (unsigned i=0; i<num_etids; i++) {
        etype[i]->unmarkForDeletion();
    }
    is_marked_for_deletion = false;
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Removed entries for operation %p %s\n", this, getName());
    fflush(stdout);
#endif
}

void MEDDLY::operation::removeStalesFromMonolithic()
{
    // DON'T make this inlined;
    // we don't want it in the header file because
    // it requires compute_table.h
    if (Monolithic_CT) Monolithic_CT->removeStales();
}

void MEDDLY::operation::removeAllFromMonolithic()
{
    // DON'T make this inlined;
    // we don't want it in the header file because
    // it requires compute_table.h
    if (Monolithic_CT) Monolithic_CT->removeAll();
}

void MEDDLY::operation::countAllNodeEntries(const forest* f, size_t* counts)
{
    if (Monolithic_CT) {
        Monolithic_CT->countNodeEntries(f, counts);
    }
    for (unsigned i=0; i<list_size; i++) {
        if (op_list[i]) {
            op_list[i]->countCTEntries(f, counts);
        }
    }
}

void MEDDLY::operation::countCTEntries(const forest* f, size_t* counts)
    const
{
    if (CT) {
        for (unsigned i=0; i<num_etids; i++) {
            if (0==CT[i]) continue;
            if (CT[i]->isOperationTable()) {
                CT[i]->countNodeEntries(f, counts);
            }
        }
    }
}


void MEDDLY::operation::showMonolithicComputeTable(output &s, int verbLevel)
{
    // DON'T make this inlined;
    // we don't want it in the header file because
    // it requires compute_table.h
    if (Monolithic_CT) Monolithic_CT->show(s, verbLevel);
}


void MEDDLY::operation::showAllComputeTables(output &s, int verbLevel)
{
    if (Monolithic_CT) {
        Monolithic_CT->show(s, verbLevel);
        return;
    }
    for (unsigned i=0; i<list_size; i++) {
        if (op_list[i]) {
            op_list[i]->showComputeTable(s, verbLevel);
        }
    }
}

void MEDDLY::operation::purgeAllMarked()
{
    removeStalesFromMonolithic();
    for (unsigned i=0; i<list_size; i++) {
        if (!op_list[i]) continue;
        if (op_list[i]->isMarkedForDeletion()) {
            destroyOperation(op_list[i]);
        }
    }
}

void MEDDLY::operation::showComputeTable(output &s, int verbLevel) const
{
    bool has_monolithic = false;
    if (CT) {
        for (unsigned i=0; i<num_etids; i++) {
            if (0==CT[i]) continue;
            if (CT[i]->isOperationTable()) {
                CT[i]->show(s, verbLevel);
            } else {
                has_monolithic = true;
            }
        }
    }
    if (has_monolithic) {
        Monolithic_CT->show(s, verbLevel);
    }
}


//
// Protected
//

void MEDDLY::operation::markForDeletion()
{
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Marking operation %p %s for deletion\n", this, getName());
    fflush(stdout);
#endif
    if (is_marked_for_deletion) return;
    is_marked_for_deletion = true;
    for (unsigned i=0; i<num_etids; i++) {
        etype[i]->markForDeletion();
    }
    if (CT) {
        for (unsigned i=0; i<num_etids; i++) {
            if (CT[i] && CT[i]->isOperationTable()) CT[i]->removeStales();
        }
    }
}

void MEDDLY::operation::registerInForest(MEDDLY::forest* f)
{
    if (f) f->registerOperation(this);
}

void MEDDLY::operation::unregisterInForest(MEDDLY::forest* f)
{
    if (f) f->unregisterOperation(this);
}

void MEDDLY::operation::registerEntryType(unsigned slot, ct_entry_type* et)
{
    MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, slot, num_etids);
    MEDDLY_DCASSERT(etype);
    MEDDLY_DCASSERT(0==etype[slot]);
    etype[slot] = et;
    compute_table::registerEntryType(first_etid + slot, et);
}

void MEDDLY::operation::buildCTs()
{
    if (0==num_etids) return;

    CT = new compute_table* [num_etids];

    if (Monolithic_CT) {
        for (unsigned i=0; i<num_etids; i++) {
            CT[i] = Monolithic_CT;
        }
    } else {
        for (unsigned i=0; i<num_etids; i++) {
            CT[i] = ct_initializer::createForOp(this, i);
        }
    }

    //
    // Initialize CTresults
    //
    for (unsigned i=0; i<num_etids; i++) {
        CTresult[i].initialize(etype[i]);
    }

    //
    // Most operations use only one slot
    //
    CT0 = CT[0];
}


//
// Private
//

void MEDDLY::operation::initializeStatics()
{
  op_list = nullptr;
  op_holes = nullptr;
  list_size = 0;
  list_alloc = 0;
  free_list = 0;
  Monolithic_CT = 0;
}

void MEDDLY::operation::destroyAllOps()
{
  for (unsigned i=0; i<list_size; i++) delete op_list[i];
  free(op_list);
  free(op_holes);
  initializeStatics();
}

void MEDDLY::operation::registerOperation(operation &o)
{
    //
    // assign an index to this operation
    //
    if (free_list) {
        o.oplist_index = free_list;
        free_list = op_holes[free_list];
    } else {
        if (list_size >= list_alloc) {
            unsigned nla = list_alloc + 256;
            op_list = (operation**) realloc(op_list, nla * sizeof(void*));
            op_holes = (unsigned*) realloc(op_holes, nla * sizeof(unsigned));
            if (0==op_list || 0==op_holes) {
                throw error(error::INSUFFICIENT_MEMORY, __FILE__, __LINE__);
            }
            list_alloc = nla;
            for (unsigned i=list_size; i<list_alloc; i++) {
                op_list[i] = 0;
                op_holes[i] = 0;
            }
        }
        if (0==list_size) {
            // Never use slot 0
            list_size++;
        }
        o.oplist_index = list_size;
        list_size++;
    }
    op_list[o.oplist_index] = &o;
}

void MEDDLY::operation::unregisterOperation(operation &o)
{
    if (o.oplist_index) {
        MEDDLY_DCASSERT(op_list[o.oplist_index] == &o);
        op_list[o.oplist_index] = nullptr;
        op_holes[o.oplist_index] = free_list;
        free_list = o.oplist_index;
    }
}

// ******************************************************************
// *                                                                *
// *                      front-end  functions                      *
// *                                                                *
// ******************************************************************

void MEDDLY::destroyOperation(MEDDLY::operation* &op)
{
    if (!op) return;
    opname* parent = op->getParent();
    if (parent) parent->removeOperationFromCache(op);
    if (!op->isMarkedForDeletion()) {
        op->markForDeletion();
        operation::removeStalesFromMonolithic();
    }
    delete op;
    op = nullptr;
}

