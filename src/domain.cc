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

#include "defines.h"
#include "domain.h"
#include "io.h"

#include "initializer.h"

#include "unique_table.h"
#include "operators.h"

#include "oper.h"

// #define DUMP_ON_FOREST_DESTROY

// ******************************************************************
// *                         domain statics                         *
// ******************************************************************

MEDDLY::domain* MEDDLY::domain::domain_list;

// ******************************************************************
// *                                                                *
// *                         domain methods                         *
// *                                                                *
// ******************************************************************

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Public static
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MEDDLY::domain* MEDDLY::domain::create(variable** vars, unsigned N)
{
    if (!initializer_list::libraryIsRunning()) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
#ifdef ALLOW_DEPRECATED_0_17_2
    return new expert_domain(vars, N);
#else
    return new domain(vars, N);
#endif
}

MEDDLY::domain* MEDDLY::domain::createBottomUp(const int* bounds, unsigned N)
{
    if (!initializer_list::libraryIsRunning()) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
#ifdef ALLOW_DEPRECATED_0_17_2
    domain* d = new expert_domain(nullptr, 0);
#else
    domain* d = new domain(nullptr, 0);
#endif
    d->createVariablesBottomUp(bounds, N);
    return d;
}


void MEDDLY::domain::destroy(MEDDLY::domain* &d)
{
    if (!d) return;
    if (!initializer_list::libraryIsRunning()) {
        throw error(error::UNINITIALIZED, __FILE__, __LINE__);
    }
    d->markForDeletion();
    // operation::purgeAllMarked();

    //
    // Remove d from the domain_list
    //
    domain* dp = d->prev;
    domain* dn = d->next;
    if (dp) {
        MEDDLY_DCASSERT(dp->next == d);
        MEDDLY_DCASSERT(domain_list != d);
        dp->next = dn;
    } else {
        MEDDLY_DCASSERT(domain_list == d);
        domain_list = dn;
    }
    if (dn) {
        MEDDLY_DCASSERT(dn->prev == d);
        dn->prev = dp;
    }
    delete d;
    d = nullptr;
}


void MEDDLY::domain::testMarkAllDomains(bool mark)
{
#ifdef DEVELOPMENT_CODE
    domain* p = nullptr;
#endif
    for (domain* d = domain_list; d; d=d->next) {
        MEDDLY_DCASSERT(d->prev == p);
        d->is_marked_for_deletion = mark;
#ifdef DEVELOPMENT_CODE
        p = d;
#endif
    }
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// I/O methods
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::domain::write(output &s) const
{
    s << "dom\n" << nVars << "\n";
    for (unsigned i=nVars; i; i--) {
        s.put(long(vars[i]->getBound(false)));
        s.put(' ');
        MEDDLY_DCASSERT(vars[i]->getBound(false) == vars[i]->getBound(true));
    }
    s << "\nmod\n";
}

void MEDDLY::domain::read(input &s)
{
    // domain must be empty -- no variables defined so far
    if (!forestReg.empty() || nVars) {
        throw error(error::DOMAIN_NOT_EMPTY, __FILE__, __LINE__);
    }

    s.stripWS();
    s.consumeKeyword("dom");
    s.stripWS();
    nVars = unsigned(s.get_integer());
    if (nVars) {
        vars = new variable*[1+nVars];
        vars[0] = nullptr;
    } else {
        vars = nullptr;
    }
    for (unsigned i=nVars; i; i--) {
        s.stripWS();
        long bound = s.get_integer();
        vars[nVars-i+1] = new variable(bound, "");
        vars[i]->addToList(this);
    }
    s.stripWS();
    s.consumeKeyword("mod");
}

void MEDDLY::domain::showInfo(output &strm)
{
    // list variables handles, their bounds and heights.
    strm << "Domain info:\n";
    strm << "  #variables: " << nVars << "\n";
    strm << "  Variables listed in height-order (ascending):\n";
    strm << "    height\t\tname\t\textensible\t\tbound\t\tprime-bound\n";
    for (unsigned i = 1; i < nVars + 1; ++i) {
        strm    << "    " << i << "\t\t" << vars[i]->getName().c_str()
                << "\t\t" << (vars[i]->isExtensible()? "yes": "no")
                << "\t\t" << vars[i]->getBound(false)
                << "\t\t" << vars[i]->getBound(true) << "\n";
  }
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Domain building
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::domain::createVariablesBottomUp(const int* bounds, unsigned N)
{
    // domain must be empty -- no variables defined so far
    if (hasForests() || nVars != 0) {
        throw error(error::DOMAIN_NOT_EMPTY, __FILE__, __LINE__);
    }

    vars = new variable*[1+N];
    nVars = N;

    vars[0] = nullptr;
    for (unsigned i=1; i<=N; i++) {
        vars[i] = new variable(bounds[i-1], "");
        vars[i]->addToList(this);
    }

    // Create the default variable order
    var_orders.clear();
    int* defaultOrder = new int[N + 1];
    defaultOrder[0] = 0;
    for (unsigned i = 1; i <= N; i++) {
        defaultOrder[i] = int(i);
    }
    default_var_order = std::make_shared<variable_order>(defaultOrder, N);
    delete[] defaultOrder;
    var_orders.push_back(default_var_order);
}


void MEDDLY::domain::createVariablesTopDown(const int* bounds, unsigned N)
{
    // domain must be empty -- no variables defined so far
    if (hasForests() || nVars != 0) {
        throw error(error::DOMAIN_NOT_EMPTY, __FILE__, __LINE__);
    }

    vars = new variable*[1+N];
    nVars = N;

    vars[0] = nullptr;
    for (unsigned i=N; i; i--) {
        vars[N-i+1] = MEDDLY::createVariable(bounds[i], "");
        vars[i]->addToList(this);
    }

    // Create the default variable order
    var_orders.clear();
    int* defaultOrder = new int[N + 1];
    defaultOrder[0] = 0;
    for (unsigned i = N; i >= 1; i--) {
        defaultOrder[N - i + 1] = int(i);
    }
    default_var_order = std::make_shared<variable_order>(defaultOrder, N);
    delete[] defaultOrder;
    var_orders.push_back(default_var_order);
}


unsigned MEDDLY::domain::findLevelOfVariable(const variable *v) const
{
    // TBD: more efficient implementation based on binary search?
    unsigned i;
    for (i=nVars; i; i--) {
        if (vars[i] == v) break;
    }
    return i;
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Variable ordering/reordering
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

std::shared_ptr<const MEDDLY::variable_order>
MEDDLY::domain::makeVariableOrder(const int* order)
{
    cleanVariableOrders();

    for (const auto& p : var_orders) {
        if (p->is_compatible_with(order)) {
            return p;
        }
    }

    std::shared_ptr<const variable_order>
        p = std::make_shared<variable_order>(order, getNumVariables());

    var_orders.push_back(p);
    return p;
}

std::shared_ptr<const MEDDLY::variable_order>
MEDDLY::domain::makeVariableOrder(const variable_order& order)
{
    cleanVariableOrders();

    for (const auto& p : var_orders) {
        if (p->is_compatible_with(order)) {
            return p;
        }
    }

    std::shared_ptr<const variable_order>
        p = std::make_shared<variable_order>(order);

    var_orders.push_back(p);
    return p;
}

void MEDDLY::domain::cleanVariableOrders()
{
    // var_orders[0] is reserved
    size_t i = 1;
    while (i < var_orders.size()) {
        if (var_orders[i].use_count() == 1) {
            var_orders[i] = var_orders.back();
            var_orders.pop_back();
        }
        else {
            i++;
        }
    }
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Forest registry.
// Not inlined to hide forest details from our header file.
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::domain::registerForest(forest* f)
{
    MEDDLY_DCASSERT(f);
#ifdef DEBUG_CLEANUP
    std::cerr << "In domain " << this << ": registering forest " << f->FID() << "\n";
#endif
    forestReg.insert(f->FID());
}

void MEDDLY::domain::unregisterForest(forest* f)
{
    // Don't bother with the registry if we're marked for deletion.
    // Also, this prevents us from trying to change the container
    // while we're iterating through it in our destructor.
    if (is_marked_for_deletion) return;

    MEDDLY_DCASSERT(f);
#ifdef DEBUG_CLEANUP
    std::cerr << "In domain " << this << ": unregistering forest " << f->FID() << "\n";
#endif
    forestReg.erase(f->FID());
}

void MEDDLY::domain::markForDeletion()
{
    if (is_marked_for_deletion) return;
#ifdef DEBUG_CLEANUP
    std::cerr << "Marking domain " << this << " for deletion\n";
#endif
    is_marked_for_deletion = true;


    for (auto it = forestReg.begin(); it != forestReg.end(); ++it) {
        forest* f = forest::getForestWithID(*it);
        if (f) {
            f->markForDeletion();
        }
    }
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// domain constructor/destructor
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++



MEDDLY::domain::domain(variable** v, unsigned N)
{
    vars = v;
    nVars = N;
    for (unsigned i=1; i<=N; i++) {
        vars[i]->addToList(this);
    }
    is_marked_for_deletion = false;

    //
    // Create the default variable order
    //
    int* defaultOrder = new int[N + 1];
    for (unsigned i = 0; i < N + 1; i++) {
        defaultOrder[i] = int(i);
    }
    default_var_order = std::make_shared<variable_order>(defaultOrder, N);
    delete[] defaultOrder;
    var_orders.push_back(default_var_order);

    //
    // Add myself to the master list of domains
    //
    if (domain_list) {
        domain_list->prev = this;
    }
    prev = nullptr;
    next = domain_list;
    domain_list = this;
}

MEDDLY::domain::~domain()
{
    //
    // Delete all forests using this domain
    //
    MEDDLY_DCASSERT(is_marked_for_deletion);
    for (auto it=forestReg.begin(); it!=forestReg.end(); ++it) {
        forest* f = forest::getForestWithID(*it);
        if (!f) continue;
#ifdef DUMP_ON_FOREST_DESTROY
        MEDDLY_DCASSERT(ef);
        std::cerr << "Destroying forest #" << *it << "\n";
        ostream_output s(std::cerr);
        f->dump(s, SHOW_DETAILS);
#endif
        delete f;
    }

    //
    // Delete my variables
    //
    for (unsigned i=1; i<=nVars; i++) {
        vars[i]->removeFromList(this);
    }
    delete[] vars;

    //
    // DON'T remove myself from the master list
    // That's done in domain::destroy()
    //
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Domain registry
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MEDDLY::domain::initDomList()
{
    domain_list = nullptr;
}

void MEDDLY::domain::markDomList()
{
    for (domain* d = domain_list; d; d=d->next) {
        d->markForDeletion();
    }
}

void MEDDLY::domain::deleteDomList()
{
    while (domain_list) {
        domain* dn = domain_list->next;
        delete domain_list;
        domain_list = dn;
    }
    MEDDLY_DCASSERT(!domain_list);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// OLD: Forest building
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef ALLOW_DEPRECATED_0_17_2

MEDDLY::forest* MEDDLY::domain::createForest(bool rel, range_type t,
    edge_labeling e, const policies &p, int* level_reduction_rule, int tv)
{
    return forest::create(this, rel, t, e, p, level_reduction_rule, tv);
}

MEDDLY::forest*
MEDDLY::domain
::createForest(bool rel, range_type t, edge_labeling e)
{
    return forest::create(this, rel, t, e);
}

#endif
