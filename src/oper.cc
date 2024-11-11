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

std::vector <MEDDLY::operation*> MEDDLY::operation::op_list;
std::vector <unsigned> MEDDLY::operation::free_list;

// ******************************************************************
// *                                                                *
// *                       operation  methods                       *
// *                                                                *
// ******************************************************************

MEDDLY::operation::operation(const char* n, unsigned et_slots)
{
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Creating operation %p\n", this);
    fflush(stdout);
#endif
    name = n;

    registerOperation(*this);

    //
    // Delay CT initialization!
    // The derived class hasn't set up the entry types yet!
    //
#ifdef ALLOW_DEPRECATED_0_17_6
    CT = 0;
    num_etids = et_slots;

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
#endif

    //
    // Initialize list of forests
    //
    FList.clear();
}

MEDDLY::operation::operation()
{
    name = nullptr;
    registerOperation(*this);
    FList.clear();

#ifdef ALLOW_DEPRECATED_0_17_6
    CT = nullptr;
    etype = nullptr;
    CTresult = nullptr;
    num_etids = 0;
#endif
}


MEDDLY::operation::~operation()
{
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Deleting operation %p %s\n", this, getName());
    fflush(stdout);
#endif

#ifdef ALLOW_DEPRECATED_0_17_6
    if (CT) {
        delete[] CT;
        CT = nullptr;
    }

    //
    // Eventually operation destructors
    // will destroy these
    //
    for (unsigned i=0; i<num_etids; i++) {
        if (etype[i]) {
            etype[i]->markForDestroy();
        }
    }

    // Don't delete the entries in etype, they're owned by compute_table.
    delete[] etype;
    delete[] CTresult;
#endif

    unregisterOperation(*this);
#ifdef DEBUG_CLEANUP
    fprintf(stdout, "Deleted operation %p %s\n", this, getName());
    fflush(stdout);
#endif
}


void MEDDLY::operation::destroyAllWithForest(const forest* f)
{
    if (!f) return;
    for (unsigned i=0; i<op_list.size(); i++) {
        if (!op_list[i]) continue;
        operation* op = op_list[i];
        //
        // Check if this operation uses f
        //
        bool uses_f = false;
        for (unsigned j=0; j<op->FList.size(); j++) {
            if (f->FID() != op->FList[j]) continue;
            uses_f = true;
            break;
        }
        if (!uses_f) continue;

        delete op;
        MEDDLY_DCASSERT(nullptr == op_list[i]);
    }
}


//
// Protected
//

void MEDDLY::operation::registerInForest(MEDDLY::forest* f)
{
    if (!f) return;
#ifdef FOREST_OPN_REGISTRY
    if (f) f->registerOperation(this);
#endif
    //
    // See if FList already contains this FID
    //
    for (unsigned i=0; i<FList.size(); i++) {
        if (f->FID() == FList[i]) return;
    }

    //
    // Nope, add it
    //
    FList.push_back(f->FID());
}

void MEDDLY::operation::unregisterInForest(MEDDLY::forest* f)
{
#ifdef FOREST_OPN_REGISTRY
    if (f) f->unregisterOperation(this);
#endif
    //
    // It is safe to NOT update FList
    //
}

#ifdef ALLOW_DEPRECATED_0_17_6
void MEDDLY::operation::registerEntryType(unsigned slot, ct_entry_type* et)
{
    MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, slot, num_etids);
    MEDDLY_DCASSERT(etype);
    MEDDLY_DCASSERT(0==etype[slot]);
    etype[slot] = et;
}

void MEDDLY::operation::buildCTs()
{
    if (0==num_etids) return;

    //
    // TBD: for now; eventually remove!
    CT = new compute_table* [num_etids];

    for (unsigned i=0; i<num_etids; i++) {
        CT[i] = etype[i]->getCT();
    }
    //
    // Most operations use only one slot
    //
    CT0 = CT[0];

    //
    // Initialize CTresults
    //
    for (unsigned i=0; i<num_etids; i++) {
        CTresult[i].initialize(etype[i]);
    }
}
#endif


//
// Private
//

void MEDDLY::operation::initializeStatics()
{
    //
    // Global operation registry
    //
    op_list.clear();
    op_list.push_back(nullptr);
    free_list.clear();
}

void MEDDLY::operation::destroyAllOps()
{
    unsigned i = op_list.size();
    while (i) {
        --i;
#ifdef DEVELOPMENT_CODE
        delete op_list.at(i);
#else
        delete op_list[i];
#endif
        MEDDLY_DCASSERT(nullptr == op_list[i]);
    }
    op_list.clear();
    free_list.clear();
}

void MEDDLY::operation::registerOperation(operation &o)
{
    if (free_list.size()) {
        unsigned u = free_list.back();
        free_list.pop_back();
        MEDDLY_DCASSERT(nullptr == op_list.at(u));
        o.oplist_index = u;
        op_list[u] = &o;
    } else {
        o.oplist_index = op_list.size();
        op_list.push_back(&o);
    }
}

void MEDDLY::operation::unregisterOperation(operation &o)
{
    if (!o.oplist_index) return;
    MEDDLY_DCASSERT(op_list[o.oplist_index] == &o);
    op_list[o.oplist_index] = nullptr;
    if (o.oplist_index == op_list.size()-1) {
        // last operator, shrink
        op_list.pop_back();
    } else {
        // Not last operator, add to free list
        free_list.push_back(o.oplist_index);
    }
}

