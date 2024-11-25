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

#include "minterms.h"

#include "io.h"
#include "domain.h"

// #define DEBUG_COLLECT_FIRST
// #include "operators.h"

// ******************************************************************
// *                                                                *
// *                        minterm  methods                        *
// *                                                                *
// ******************************************************************


MEDDLY::minterm::minterm(const domain* D, set_or_rel sr,
        node_storage_flags fs) : termval(true)
{
    _D = D;
    if (D) {
        num_vars = D->getNumVariables();
    } else {
        num_vars = 0;
    }
    for_relations = sr;
    sparse = (fs == SPARSE_ONLY);
    if (!sparse) {
        _from.resize(1+num_vars, DONT_CARE);
        if (for_relations) {
            _to.resize(1+num_vars, DONT_CARE);
        }
    }
}

bool MEDDLY::minterm::hasSparseVars(const std::vector <unsigned> &vars) const
{
    if (isFull()) return false;
    if (_index.size() != vars.size()) return false;
    for (unsigned i=0; i<vars.size(); i++) {
        if (_index[i] != vars[i]) return false;
    }
    return true;
}

void MEDDLY::minterm::show(output &s) const
{
    if (isFull()) {
        //
        // Full storage
        //
        s.put("[ top");
        for (unsigned i=num_vars; i; i--) {
            s.put(", ");
            if (_from[i] < 0)   s.put('x');
            else                s.put(_from[i]);
            if (isForSets()) continue;
            s.put("->");
            if (DONT_CARE == _to[i])        s.put('x');
            else if (DONT_CHANGE == _to[i]) s.put('i');
            else                            s.put(_to[i]);
        }
        s.put(", bot]");
    } else {
        //
        // Sparse storage
        //
        s.put('(');
        for (unsigned z=0; z<size(); z++) {
            if (z) s.put(", ");
            s.put(whichVar(z));
            s.put(':');
            if (_from[z] < 0)   s.put('x');
            else                s.put(_from[z]);
            if (isForSets()) continue;
            s.put("->");
            if (DONT_CARE == _to[z])        s.put('x');
            else if (DONT_CHANGE == _to[z]) s.put('i');
            else                            s.put(_to[z]);
        }
        s.put(')');
    }
}

// ******************************************************************
// *                                                                *
// *                      minterm_coll methods                      *
// *                                                                *
// ******************************************************************

MEDDLY::minterm_coll::minterm_coll(const domain* D, set_or_rel sr)
{
    _D = D;
    num_vars = D ? D->getNumVariables() : 0;
    for_relations = sr;
    sparse = false;
    freelist = nullptr;
}

MEDDLY::minterm_coll::minterm_coll(const domain* D, set_or_rel sr,
        const std::vector<unsigned> &varlist)
{
    _D = D;
    num_vars = D ? D->getNumVariables() : 0;
    for_relations = sr;
    sparse = true;
    freelist = nullptr;

    _varlist.resize(varlist.size());
    for (unsigned i=0; i<varlist.size(); i++) {
        _varlist[i] = varlist[i];
    }
}

MEDDLY::minterm_coll::~minterm_coll()
{
    clear();
    while (freelist) {
        minterm* m = freelist;
        freelist = freelist->next;
        delete m;
    }
}

void MEDDLY::minterm_coll::clear()
{
    while (_mtlist.size()) {
        recycle( _mtlist.back() );
        _mtlist.back() = nullptr;
        _mtlist.pop_back();
    }
}

void MEDDLY::minterm_coll::show(output &s,
        const char* pre, const char* post) const
{
    for (unsigned i=0; i<_mtlist.size(); i++) {
        if (pre) {
            s.put(pre);
        } else {
            s.put("  ");
            s.put(i, 4);
            s.put("  ");
        }
        _mtlist[i]->show(s);
        s.put(post);
    }
}

bool MEDDLY::minterm_coll::_check(minterm* m) const
{
    MEDDLY_DCASSERT(m->getDomain() == _D);
    MEDDLY_DCASSERT(m->isForRelations() == for_relations);
    if (sparse) {
        MEDDLY_DCASSERT(!m->isFull());
        MEDDLY_DCASSERT(m->hasSparseVars(_varlist));
    } else {
        MEDDLY_DCASSERT(m->isFull());
    }
    return true;
}


int MEDDLY::minterm_coll::_collect_first(int L, unsigned low,
        unsigned hi, unsigned& lend)
{
    int val = at(low)->var(L);
    bool val_dc;
    if ((DONT_CARE == val) && isForRelations() && (L>0)) {
        val_dc = (DONT_CHANGE == at(low)->var(-L));
    } else {
        val_dc = false;
    }
    lend = low+1;
    unsigned hptr = hi;

    for (;;)
    {
        //
        // Find smallest element different from val
        //
        for (;;) {
            const int vl = at(lend)->var(L);
            if (vl != val) break;
            if ((DONT_CARE == vl) && isForRelations() && (L>0)) {
                const bool vl_dc = (DONT_CHANGE == at(lend)->var(-L));
                if (val_dc != vl_dc) break;
            }
            ++lend;
            if (lend >= hptr) {
#ifdef DEBUG_COLLECT_FIRST
                FILE_output out(stdout);
                out << "Collected first on range [" << low << ", "
                    << hi << "):\n";
                for (unsigned i=low; i<hi; i++) {
                    if (i==lend) out << "    ======================\n";
                    out << "    ";
                    out.put(long(i), 2);
                    out << "  ";
                    at(i)->show(out);
                    out << "\n";
                }
                out << "val " << val << "\n";
#endif
                return val;
            }
        }

        //
        // Find high element to swap with
        //
        for (;;) {
            const int vh = at(hptr-1)->var(L);
            if (vh == val) {
                if ((DONT_CARE == vh) && isForRelations() && (L>0)) {
                    const bool vh_dc = (DONT_CHANGE == at(hptr-1)->var(-L));
                    if (val_dc == vh_dc) break;
                } else {
                    break;
                }
            }
            --hptr;
            if (lend >= hptr) {
#ifdef DEBUG_COLLECT_FIRST
                FILE_output out(stdout);
                out << "Collected first on range [" << low << ", "
                    << hi << "):\n";
                for (unsigned i=low; i<hi; i++) {
                    if (i==lend) out << "    ======================\n";
                    out << "    ";
                    out.put(long(i), 2);
                    out << "  ";
                    at(i)->show(out);
                    out << "\n";
                }
                out << "val " << val << "\n";
#endif
                return val;
            }
        }

        //
        // Swap elements
        //
        minterm* tmp = _mtlist[lend];
        _mtlist[lend] = _mtlist[hptr-1];
        _mtlist[hptr-1] = tmp;
    }

}

