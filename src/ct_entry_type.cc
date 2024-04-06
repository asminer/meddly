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

// #define DEBUG_ENTRY_TYPE

// ******************************************************************
// *                                                                *
// *                       ct_itemtype methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::ct_itemtype::ct_itemtype(char c)
{
    nodeFor = nullptr;
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

// **********************************************************************
// *                                                                    *
// *                          Helper functions                          *
// *                                                                    *
// **********************************************************************

inline unsigned bytes4typeID(MEDDLY::ct_typeID t)
{
    switch (t) {
        case MEDDLY::ct_typeID::NODE    : return sizeof(MEDDLY::node_handle);
        case MEDDLY::ct_typeID::INTEGER : return sizeof(int);
        case MEDDLY::ct_typeID::LONG    : return sizeof(long);
        case MEDDLY::ct_typeID::FLOAT   : return sizeof(float);
        case MEDDLY::ct_typeID::DOUBLE  : return sizeof(double);
        case MEDDLY::ct_typeID::GENERIC : return sizeof(MEDDLY::ct_object*);
        default                         : return 0;
    }
}

inline char typeID2char(MEDDLY::ct_typeID t)
{
    switch (t) {
        case MEDDLY::ct_typeID::NODE    : return 'N';
        case MEDDLY::ct_typeID::INTEGER : return 'I';
        case MEDDLY::ct_typeID::LONG    : return 'L';
        case MEDDLY::ct_typeID::FLOAT   : return 'F';
        case MEDDLY::ct_typeID::DOUBLE  : return 'D';
        case MEDDLY::ct_typeID::GENERIC : return 'G';
        default                         : return '?';
    }
}

inline MEDDLY::ct_typeID char2typeID(char c)
{
    switch (c) {
        case 'N':   return MEDDLY::ct_typeID::NODE;
        case 'I':   return MEDDLY::ct_typeID::INTEGER;
        case 'L':   return MEDDLY::ct_typeID::LONG;
        case 'F':   return MEDDLY::ct_typeID::FLOAT;
        case 'D':   return MEDDLY::ct_typeID::DOUBLE;
        case 'G':   return MEDDLY::ct_typeID::GENERIC;
        default :   return MEDDLY::ct_typeID::ERROR;
    }
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

    //
    // Determine lengths for each portion
    //

    if (!saw_colon) {
        throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    }

    unsigned len_ks_type = 0;
    len_kr_type = 0;
    len_r_type  = 0;

    if (saw_dot) {
        // "012.456:89"
        len_ks_type = dot_slot;
        MEDDLY_DCASSERT(colon_slot > dot_slot);
        len_kr_type = (colon_slot - dot_slot) - 1;
    } else {
        // "01234:67"
        len_ks_type = colon_slot;
        len_kr_type = 0;
    }
    len_r_type = (length - colon_slot) - 1;

    if (
        (saw_dot && 0==len_kr_type)             // "foo.:bar" is bad
        ||
        (len_ks_type + len_kr_type == 0)        // no key?
        ||
        (len_r_type == 0)                       // no result?
      )
    {
        throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
    }

    //
    // Build fixed, starting portion of key
    //
    fixed_bytes = 0;
    if (len_ks_type) {
        key_fixed.resize(len_ks_type);
        for (unsigned i=0; i<len_ks_type; i++) {
            key_fixed[i] = ct_itemtype(pattern[i]);
            fixed_bytes += key_fixed[i].bytes();
        }
    }

    /*
    if (len_ks_type) {
        ks_type = new ct_typeID[len_ks_type];
        ks_forest = new forest*[len_ks_type];
        for (unsigned i=0; i<len_ks_type; i++) {
            ks_type[i] = char2typeID(pattern[i]);
            fixed_bytes += bytes4typeID(ks_type[i]);
            ks_forest[i] = 0;
        }
    } else {
        // This is possible if the pattern begins with .
        ks_type = 0;
        ks_forest = 0;
    }
    */

    //
    // Build repeating portion of key
    //
    kr_bytes = 0;
    if (len_kr_type) {
        kr_type = new ct_typeID[len_kr_type];
        kr_forest = new forest*[len_kr_type];
        for (unsigned i=0; i<len_kr_type; i++) {
            kr_type[i] = char2typeID(pattern[i + dot_slot + 1]);
            kr_bytes += bytes4typeID(kr_type[i]);
            kr_forest[i] = 0;
        }
    } else {
        kr_type = 0;
        kr_forest = 0;
    }

    //
    // Build result
    //
    MEDDLY_DCASSERT(len_r_type);
    r_bytes = 0;
    r_type = new ct_typeID[len_r_type];
    r_forest = new forest*[len_r_type];
    for (unsigned i=0; i<len_r_type; i++) {
        r_type[i] = char2typeID(pattern[i + colon_slot + 1]);
        r_bytes += bytes4typeID(r_type[i]);
        r_forest[i] = 0;
    }

#ifdef DEBUG_ENTRY_TYPE
    printf("Built entry type %s with pattern '%s'\n", name, pattern);
    printf("Key start: \"");
    for (unsigned i=0; i<key_fixed.size(); i++) {
        fputc(key_fixed[i].getTypeChar(), stdout);
    }
    printf("\"  (%u bytes)\n", fixed_bytes);
    printf("Key repeat: \"");
    for (unsigned i=0; i<len_kr_type; i++) {
        fputc(typeID2char(kr_type[i]), stdout);
    }
    printf("\"  (%u bytes)\n", kr_bytes);
    printf("Result: \"");
    for (unsigned i=0; i<len_r_type; i++) {
        fputc(typeID2char(r_type[i]), stdout);
    }
    printf("\"  (%u bytes)\n", r_bytes);
#endif
}

#endif // allow_deprecated_0_17_6

MEDDLY::ct_entry_type::ct_entry_type(const char* _name)
{
    name = _name;
    is_marked_for_deletion = false;

    updatable_result = false;

    throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
}

MEDDLY::ct_entry_type::~ct_entry_type()
{
    delete[] kr_type;
    delete[] kr_forest;
    delete[] r_type;
    delete[] r_forest;
}

#ifdef ALLOW_DEPRECATED_0_17_6

void MEDDLY::ct_entry_type::setForestForSlot(unsigned i, forest* f)
{
    if (i<key_fixed.size()) {
        key_fixed[i].setForest(f);
        /*
        if (ct_typeID::NODE != ks_type[i]) {
            throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
        }
        ks_forest[i] = f;
        */
        return;
    }
    i -= key_fixed.size();

    if (len_kr_type) {
        // adjust for the .
        i--;
        if (i<len_kr_type) {
            if (ct_typeID::NODE != kr_type[i]) {
                throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
            }
            kr_forest[i] = f;
            return;
        }
    }

    i -= len_kr_type;
    // adjust for :
    i--;

    if (i < len_r_type) {
        if (ct_typeID::NODE != r_type[i]) {
            throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
        }
        r_forest[i] = f;
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
        forest* f = key_fixed[i].getForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = 1;
    }
    for (i=0; i<len_kr_type; i++) {
        forest* f = kr_forest[i];
        if (0==f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (skipF[ f->FID() ]) continue;
        f->clearAllCacheBits();
        skipF[ f->FID() ] = 1;
    }
    for (i=0; i<len_r_type; i++) {
        forest* f = r_forest[i];
        if (0==f) continue;
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
        forest* f = key_fixed[i].getForest();
        if (!f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
    for (i=0; i<len_kr_type; i++) {
        forest* f = kr_forest[i];
        if (0==f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
    for (i=0; i<len_r_type; i++) {
        forest* f = r_forest[i];
        if (0==f) continue;
        MEDDLY_DCASSERT(f->FID() < N);
        if (whichF[ f->FID() ]) {
            f->sweepAllCacheBits();
            whichF[ f->FID() ] = 0;
        }
    }
}

