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
    nodeFID = 0;
    switch (c) {
        case 'N':   type = ct_typeID::NODE;         break;
        case 'I':   type = ct_typeID::INTEGER;      break;
        case 'L':   type = ct_typeID::LONG;         break;
        case 'F':   type = ct_typeID::FLOAT;        break;
        case 'D':   type = ct_typeID::DOUBLE;       break;
        case 'G':   type = ct_typeID::GENERIC;      break;
        default :   type = ct_typeID::ERROR;
    }
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
    s.put(nodeFID);
    s.put(')');
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

#ifdef ALLOW_DEPRECATED_0_17_6
MEDDLY::ct_entry_type::ct_entry_type(const char* _name, const char* pattern)
{
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
    countFixed();

    for (unsigned i=1+dot_slot; i<colon_slot; i++) {
        key_repeating.push_back(pattern[i]);
    }
    countRepeating();

    for (unsigned i=1+colon_slot; i<length; i++) {
        result.push_back(pattern[i]);
    }
    countResult();


#ifdef DEBUG_ENTRY_TYPE
    ostream_output sout(std::cout);
    sout << "Built entry type " << name << " with pattern '"
         << pattern << "'\n";
    sout << "Key start: \"";
    for (unsigned i=0; i<key_fixed.size(); i++) {
        key_fixed[i].show(sout);
    }
    sout << "\"  (" << fixed_bytes << " bytes)\n";
    sout << "Key repeat: \"";
    for (unsigned i=0; i<key_repeating.size(); i++) {
        key_repeating[i].show(sout);
    }
    sout << "\"  (" << repeating_bytes << " bytes)\n";
    sout << "Result: \"";
    for (unsigned i=0; i<result.size(); i++) {
        result[i].show(sout);
    }
    sout << "\"  (" << result_bytes << " bytes)\n";
#endif
}

#endif // allow_deprecated_0_17_6

MEDDLY::ct_entry_type::ct_entry_type(const char* _name)
{
    name = _name;
    is_marked_for_deletion = false;

    updatable_result = false;

    fixed_bytes = 0;
    fixed_intslots = 0;

    repeating_bytes = 0;
    repeating_intslots = 0;

    result_bytes = 0;
    result_intslots = 0;
}

MEDDLY::ct_entry_type::~ct_entry_type()
{
}

#ifdef ALLOW_DEPRECATED_0_17_6

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



void MEDDLY::ct_entry_type::clearForestCTBits(bool* skipF, unsigned N) const
{
    unsigned i;
    for (i=0; i<key_fixed.size(); i++) {
        forest* f = key_fixed[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = 1;
    }
    for (i=0; i<key_repeating.size(); i++) {
        forest* f = key_repeating[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = 1;
    }
    for (i=0; i<result.size(); i++) {
        forest* f = result[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = 1;
    }
}

void MEDDLY::ct_entry_type::sweepForestCTBits(bool* whichF, unsigned N) const
{
    unsigned i;
    for (i=0; i<key_fixed.size(); i++) {
        forest* f = key_fixed[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
    for (i=0; i<key_repeating.size(); i++) {
        forest* f = key_repeating[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
    for (i=0; i<result.size(); i++) {
        forest* f = result[i].rawForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
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

void MEDDLY::ct_entry_type::countFixed()
{
    fixed_bytes = 0;
    fixed_intslots = 0;
    for (unsigned i=0; i<key_fixed.size(); i++) {
        fixed_bytes += key_fixed[i].bytes();
        fixed_intslots += key_fixed[i].intslots();
    }
}

void MEDDLY::ct_entry_type::countRepeating()
{
    repeating_bytes = 0;
    repeating_intslots = 0;
    for (unsigned i=0; i<key_repeating.size(); i++) {
        repeating_bytes += key_repeating[i].bytes();
        repeating_intslots += key_repeating[i].intslots();
    }
}

void MEDDLY::ct_entry_type::countResult()
{
    result_bytes = 0;
    result_intslots = 0;
    for (unsigned i=0; i<result.size(); i++) {
        result_bytes += result[i].bytes();
        result_intslots += result[i].intslots();
    }
}

