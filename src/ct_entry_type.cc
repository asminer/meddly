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

#include "ct_entry_type.h"
#include "ct_entry_key.h"
#include "ct_initializer.h"
#include "error.h"
#include "forest.h"

#include "io.h"
#include "operators.h"

// #define DEBUG_ENTRY_TYPE

// ******************************************************************
// *                                                                *
// *                       ct_itemtype methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::ct_itemtype::ct_itemtype(char c)
{
    switch (c) {
        case 'N':   type = ct_typeID::NODE;         break;
        case 'I':   type = ct_typeID::INTEGER;      break;
        case 'L':   type = ct_typeID::LONG;         break;
        case 'F':   type = ct_typeID::FLOAT;        break;
        case 'D':   type = ct_typeID::DOUBLE;       break;
        case 'G':   type = ct_typeID::GENERIC;      break;
        default :   type = ct_typeID::ERROR;
    }
    typeUpdate(nullptr);
}

char MEDDLY::ct_itemtype::getTypeChar() const
{
    switch (type) {
        case ct_typeID::NODE    : return 'N';
        case ct_typeID::INTEGER : return 'I';
        case ct_typeID::LONG    : return 'L';
        case ct_typeID::FLOAT   : return 'F';
        case ct_typeID::DOUBLE  : return 'D';
        case ct_typeID::GENERIC : return 'G';
        default                 : return '?';
    }
}

void MEDDLY::ct_itemtype::show(output &s) const
{
    s.put(getTypeChar());
    if (ct_typeID::NODE != type) return;
    s.put("(f ");
#ifdef USE_FID
    s.put(nodeFID);
#else
    s.put(nodeF ? nodeF->FID() : 0);
#endif
    s.put(')');
}

void MEDDLY::ct_itemtype::typeUpdate(forest* f)
{
#ifdef USE_FID
    nodeFID = f ? f->FID() : 0;
#else
    nodeF = f;
#endif
    switch (type) {
        case ct_typeID::NODE    :
                                    twoslots = sizeof(node_handle)>sizeof(int);
                                    should_hash = true;
                                    return;

        case ct_typeID::INTEGER :
                                    twoslots = false;
                                    should_hash = true;
                                    return;

        case ct_typeID::LONG    :
                                    twoslots = true;
                                    should_hash = true;
                                    return;

        case ct_typeID::FLOAT   :
                                    twoslots = sizeof(float) > sizeof(int);
                                    should_hash = false;
                                    return;

        case ct_typeID::DOUBLE  :
                                    twoslots = sizeof(double) > sizeof(int);
                                    should_hash = false;
                                    return;

        case ct_typeID::GENERIC :
                                    twoslots = true;
                                    should_hash = false;
                                    return;

        default                 :
                                    twoslots = false;
                                    should_hash = false;
    }
}

// **********************************************************************
// *                                                                    *
// *                         ct_object  methods                         *
// *                                                                    *
// **********************************************************************

MEDDLY::ct_object::ct_object()
{
}

MEDDLY::ct_object::~ct_object()
{
}

void MEDDLY::ct_object::show(output &s) const
{
    s.put_hex((unsigned long) this);
    s.put(' ');
    s.put('G');
}

// **********************************************************************
// *                                                                    *
// *                       ct_entry_type  methods                       *
// *                                                                    *
// **********************************************************************

MEDDLY::ct_entry_type::~ct_entry_type()
{
    if (CT) {
        if (CT->isOperationTable()) {
            delete CT;
        }
        // otherwise CT is the monolithic table; don't delete it.
    }
    unregisterEntry(this);
}



#ifdef ALLOW_DEPRECATED_0_17_6
MEDDLY::ct_entry_type::ct_entry_type(const char* _name, const char* pattern)
{
    CT = nullptr;
    registerEntry(this);
    name = _name;
    is_marked_for_deletion = false;
    updatable_result = false;

    bool saw_dot = false;
    bool saw_colon = false;

    unsigned dot_slot = 0;
    unsigned colon_slot = 0;

    //
    // First scan: find '.' and ':'
    //

    unsigned length;
    for (length=0; pattern[length]; length++) {
        if ('.' == pattern[length]) {
            if (saw_dot || saw_colon) {
                throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
            }
            saw_dot = true;
            dot_slot = length;
            continue;
        }
        if (':' == pattern[length]) {
            if (saw_colon) {
                throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
            }
            saw_colon = true;
            colon_slot = length;
            continue;
        }
    }
    if (!saw_colon) {
        throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    }
    if (!saw_dot) {
        dot_slot = colon_slot;
    }

    //
    // Build fixed, repeating, and result arrays
    //
    for (unsigned i=0; i<dot_slot; i++) {
        key_fixed.push_back(pattern[i]);
    }

    for (unsigned i=1+dot_slot; i<colon_slot; i++) {
        key_repeating.push_back(pattern[i]);
    }

    for (unsigned i=1+colon_slot; i<length; i++) {
        result.push_back(pattern[i]);
    }

    //
    // "Done" building
    //
    buildCT();


#ifdef DEBUG_ENTRY_TYPE
    ostream_output sout(std::cout);
    sout << "Built #" << etID << " entry type " << name << " with pattern '"
         << pattern << "'\n";
    sout << "        Key start: \"";
    for (unsigned i=0; i<key_fixed.size(); i++) {
        key_fixed[i].show(sout);
    }
    sout << "\"  (" << fixed_intslots << " intslots)\n";
    sout << "        Key repeat: \"";
    for (unsigned i=0; i<key_repeating.size(); i++) {
        key_repeating[i].show(sout);
    }
    sout << "\"  (" << repeating_intslots << " intslots)\n";
    sout << "        Result: \"";
    for (unsigned i=0; i<result.size(); i++) {
        result[i].show(sout);
    }
    sout << "\"  (" << result_intslots << " intslots)\n";
#endif
}

void MEDDLY::ct_entry_type::setForestForSlot(unsigned i, forest* f)
{
    if (i<key_fixed.size()) {
        key_fixed[i].setForest(f);
        return;
    }
    i -= key_fixed.size();

    if (key_repeating.size()) {
        // adjust for the .
        i--;
        if (i<key_repeating.size()) {
            key_repeating[i].setForest(f);
            return;
        }
    }

    i -= key_repeating.size();
    // adjust for :
    i--;

    if (i < result.size()) {
        result[i].setForest(f);
        return;
    }

    // i is too large
    throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
}

#endif // allow_deprecated_0_17_6

MEDDLY::ct_entry_type::ct_entry_type(const char* _name)
{
    CT = nullptr;
    registerEntry(this);
    name = _name;
    is_marked_for_deletion = false;

    updatable_result = false;

    // fixed_bytes = 0;
    fixed_intslots = 0;

    // repeating_bytes = 0;
    repeating_intslots = 0;

    // result_bytes = 0;
    result_intslots = 0;
}


bool MEDDLY::ct_entry_type::hasForest(const forest* f) const
{
    unsigned i;
    for (i=0; i<key_fixed.size(); i++) {
        if (key_fixed[i].hasForest(f)) return true;
    }
    for (i=0; i<key_repeating.size(); i++) {
        if (key_repeating[i].hasForest(f)) return true;
    }
    for (i=0; i<result.size(); i++) {
        if (result[i].hasForest(f)) return true;
    }
    return false;
}


void MEDDLY::ct_entry_type::clearForestCTBits(std::vector <bool> &skipF) const
{
    unsigned i;
    for (i=0; i<key_fixed.size(); i++) {
        forest* f = key_fixed[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < skipF.size());
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = true;
    }
    for (i=0; i<key_repeating.size(); i++) {
        forest* f = key_repeating[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < skipF.size());
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = true;
    }
    for (i=0; i<result.size(); i++) {
        forest* f = result[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < skipF.size());
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = true;
    }
}

void MEDDLY::ct_entry_type::clearAllForestCTBits(std::vector <bool> &skipF)
{
    for (unsigned i=1; i<all_entries.size(); i++) {
        if (all_entries[i]) {
            all_entries[i]->clearForestCTBits(skipF);
        }
    }
}

void MEDDLY::ct_entry_type::sweepForestCTBits(std::vector <bool> &whichF) const
{
    unsigned i;
    for (i=0; i<key_fixed.size(); i++) {
        forest* f = key_fixed[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < whichF.size());
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
    for (i=0; i<key_repeating.size(); i++) {
        forest* f = key_repeating[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < whichF.size());
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
    for (i=0; i<result.size(); i++) {
        forest* f = result[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < whichF.size());
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
}

void MEDDLY::ct_entry_type::sweepAllForestCTBits(std::vector <bool> &whichF)
{
    for (unsigned i=1; i<all_entries.size(); i++) {
        if (all_entries[i]) {
            all_entries[i]->sweepForestCTBits(whichF);
        }
    }
}


void MEDDLY::ct_entry_type::show(output &s) const
{
    s.put(name);
    s.put(' ');
    s.put('[');
    unsigned i;
    for (i=0; i<key_fixed.size(); i++) {
        key_fixed[i].show(s);
    }
    s.put('.');
    for (i=0; i<key_repeating.size(); i++) {
        key_repeating[i].show(s);
    }
    s.put(':');
    for (i=0; i<result.size(); i++) {
        result[i].show(s);
    }
    s.put(']');
    s.put('\n');
}

void MEDDLY::ct_entry_type::invalidateForest(const forest* f)
{
#ifndef USE_FID
    for (unsigned i=0; i<key_fixed.size(); i++) {
        key_fixed[i].invalidateForest(f);
    }
    for (unsigned i=0; i<key_repeating.size(); i++) {
        key_repeating[i].invalidateForest(f);
    }
    for (unsigned i=0; i<result.size(); i++) {
        result[i].invalidateForest(f);
    }
#endif
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Entry registry static methods and members
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::vector <MEDDLY::ct_entry_type*> MEDDLY::ct_entry_type::all_entries;

void MEDDLY::ct_entry_type::initStatics()
{
    all_entries.clear();
    all_entries.push_back(nullptr);
    // ^ reserve index 0 for 'no entry'
}

void MEDDLY::ct_entry_type::doneStatics()
{
    all_entries.clear();
}

void MEDDLY::ct_entry_type::invalidateAllWithForest(const forest* f)
{
#ifndef USE_FID
    for (unsigned i=0; i<all_entries.size(); i++) {
        if (all_entries[i]) {
            all_entries[i]->invalidateForest(f);
        }
    }
#endif
}

void MEDDLY::ct_entry_type::removeAllCTEntriesWithForest(const forest* f)
{
    for (unsigned i=0; i<all_entries.size(); i++) {
        if (!all_entries[i]) continue;
        if (!all_entries[i]->CT->isOperationTable()) continue;
        if (!all_entries[i]->hasForest(f)) continue;

        //
        // Still here?
        // We have an operation cache that uses forest f.
        // Clear it.
        all_entries[i]->CT->removeAll();
    }
}

/*
void MEDDLY::ct_entry_type::countNodeEntries(const forest* f, std::vector <unsigned long> &counts)
{
    compute_table::countMonolithicNodeEntries(f, counts);
    for (unsigned i=0; i<all_entries.size(); i++) {
        if (!all_entries[i]) continue;
        if (!all_entries[i]->CT->isOperationTable()) continue;
        if (!all_entries[i]->hasForest(f)) continue;

        all_entries[i]->CT->countNodeEntries(f, counts);
    }
}

void MEDDLY::ct_entry_type::showAllComputeTables(output &s, int verbLevel)
{
    if (compute_table::showMonolithicComputeTable(s, verbLevel)) return;
    for (unsigned i=0; i<all_entries.size(); i++) {
        if (!all_entries[i]) continue;
        if (!all_entries[i]->CT->isOperationTable()) continue;
        all_entries[i]->CT->show(s, verbLevel);
    }
}
*/

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Private helpers
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

bool MEDDLY::ct_entry_type::keyIsOurs(const ct_entry_key *k) const
{
    MEDDLY_DCASSERT(k);
    return k->getET() == this;
}

void MEDDLY::ct_entry_type::buildCT()
{
    MEDDLY_DCASSERT(!CT);

    //
    // Preprocess fixed portion of key
    //
    fixed_intslots = 0;
    fixed_entirely_hashable = true;
    for (unsigned i=0; i<key_fixed.size(); i++) {
        fixed_intslots += key_fixed[i].intslots();
        if (!key_fixed[i].shouldBeHashed()) {
            fixed_entirely_hashable = false;
        }
    }

    //
    // Preprocess repeating portion of key
    //
    repeating_intslots = 0;
    repeating_entirely_hashable = true;
    for (unsigned i=0; i<key_repeating.size(); i++) {
        repeating_intslots += key_repeating[i].intslots();
        if (!key_repeating[i].shouldBeHashed()) {
            repeating_entirely_hashable = false;
        }
    }

    //
    // Preprocess result
    //
    result_intslots = 0;
    for (unsigned i=0; i<result.size(); i++) {
        result_intslots += result[i].intslots();
    }

    //
    // Build CT or use monolithic as appropriate
    //
    if (compute_table::Monolithic()) {
        CT = compute_table::Monolithic();
    } else {
        CT = ct_initializer::createForOp(this);
    }
}

