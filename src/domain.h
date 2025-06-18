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

#ifndef MEDDLY_DOMAIN_H
#define MEDDLY_DOMAIN_H

#include "variable.h"
#include "policies.h"
#include <vector>
#include <memory>
#include <set> // for the forest registry

namespace MEDDLY {
    class initializer_list;
    class variable_order;

    class domain;

    class input;
    class output;

    class forest;
    class relforest;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                          domain class                          *
// *                                                                *
// *                                                                *
// ******************************************************************


/** Domain class.
    A domain is an ordered collection of variables,
    along with a rich set of operations for adding and removing variables.
    A variable may be shared in more than one domain.

    When a domain is destroyed, all of its forests are destroyed.
*/
class MEDDLY::domain {
    public:
        /** Front-end function to create a domain with the given variables.
                @param  vars    List of variables, in order.
                                vars[i] gives the variable at level i.
                                Note that vars[0] should be 0.
                @param  N       Number of variables.
                                vars[N] refers to the top-most variable.

                @return A new domain.
        */
        static domain* create(variable** vars=nullptr, unsigned N=0);

        /** Front-end function to create a domain with given variable bounds.
            Equivalent to creating an empty domain and then building the
            domain bottom up.

                @param  bounds  variable bounds.  bounds[i] gives the bound
                                for the variable at level i+1.  If bound<=0,
                                the variable is marked as extensible,
                                with initial bound as abs(bound).
                                Note: an extensible variable has a
                                range [1 .. +infinity].

                @param  N       Number of variables.

                @return A new domain.
        */
        static domain* createBottomUp(const int* bounds, unsigned N);

        /**
            Destroy a domain.  Will destroy all forests associated
            with the domain, and remove this domain from the domain_list.
        */
        static void destroy(domain* &d);

        /**
            For domain testing only.
            Mark/unmark all domains.
        */
        static void testMarkAllDomains(bool mark);

    private:
        /// Number of variables.
        unsigned nVars;

        /// Array of pointers to variables.
        /// Dimension is 1+nVars and element 0 is always null,
        /// so that vars[level] gives the variable for a given level.
        variable** vars;

    public:
        //
        // Inlined Getters for variables and variable information.
        //

        /// Get the number of variables in this domain.
        inline unsigned getNumVariables() const {
            return nVars;
        }

        /// @return The variable at level \a lev. Const version.
        inline const variable* getVar(unsigned lev) const {
            MEDDLY_CHECK_RANGE(1U, lev, 1+nVars);
            return vars[lev];
        }

        /// @return The variable at level \a lev.
        inline variable* getVar(unsigned lev) {
            MEDDLY_CHECK_RANGE(1U, lev, 1+nVars);
            return vars[lev];
        }

        /** Get the specified bound of a variable.
            @param  lev     Level number, should be 1 for bottom-most
                            and getNumVariables() for top-most.
            @param  prime   If prime is true, get the bound for
                            the primed variable.
            @return         The bound set for variable at level \a lev.
        */
        inline int getVariableBound(unsigned lev, bool prime = false) const {
            const variable* v = getVar(lev);
            MEDDLY_DCASSERT(v);
            return v->getBound(prime);
        }

        /**
            Enlarge the possible values for a variable.
            This could modify all nodes in all forests, depending on the
            choice of reduction rule.
                @param  lev     Level number, should be 1 for bottom-most
                                and getNumVariables() for top-most.

                @param  prime   If prime is true, enlarge the bound for
                                the primed variable only, otherwise both
                                the primed and unprimed are enlarged.

                @param  b       New bound, if less than the current bound
                                an error code is returned.  If bound<=0,
                                the variable is marked as extensible,
                                with initial bound as abs(bound).
                                Note: an extensible variable has
                                a range [1 .. +infinity].
        */
        inline void enlargeVariableBound(unsigned lev, bool prime, int b) {
            variable* v = getVar(lev);
            MEDDLY_DCASSERT(v);
            v->enlargeBound(prime, b);
        }


        /**
            Shrink the possible values for a variable.
            This could modify all nodes in all forests, depending on the
            choice of reduction rule.

                @param  lev     Level number, should be 1 for bottom-most
                                and getNumVariables() for top-most.

                @param  b       New bound.

                @param  force   If \a b is too small, and information will
                                be lost, proceed anyway if \a force is true,
                                otherwise throw an error.
        */
        inline void shrinkVariableBound(unsigned lev, int b, bool force) {
            variable* v = getVar(lev);
            MEDDLY_DCASSERT(v);
            v->shrinkBound(b, force);
        }


    public:
        //
        // File I/O to save and recover a domain
        //

        /**
            Write the domain to a file in a format that can be read back later.
                @param  s   Stream to write to

                @throws     COULDNT_WRITE, if writing failed
        */
        void write(output &s) const;

        /**
            Initialize the domain from data in a file.
            Allows reconstruction of a domain that we saved using \a write().
            The domain should be empty.
                @param  s   Stream to read from

                @throws     INVALID_FILE, if the file does not match
                            what we expect
        */
        void read(input &s);

        /**
            Display lots of information about the domain.
            This is primarily for aid in debugging.
                @param  strm    Stream to write to.
        */
        void showInfo(output &strm);

    public:
        //
        // Domain building
        //

        /**
            Create all variables at once, from the bottom up.
            Requires the domain to be "empty" (containing no variables or
            forests).

                @param  bounds  variable bounds.  bounds[i] gives the bound
                                for the variable at level i+1.  If bound<=0,
                                the variable is marked as extensible,
                                with initial bound as abs(bound).
                                Note: an extensible variable has a
                                range [1 .. +infinity].

                @param  N       Number of variables.
        */
        void createVariablesBottomUp(const int* bounds, unsigned N);

        /**
            Create all variables at once, from the top down.
            Requires the domain to be "empty" (containing no variables or
            forests).

                @param  bounds  Current variable bounds.  bounds[0] gives
                                the bound for the top-most variable, and
                                bounds[N-1] gives the bound for the bottom-most
                                variable.  If bound<=0, the variable is marked
                                as extensible, with initial bound as abs(bound).
                                Note: an extensible variable has a
                                range [1 .. +infinity].

                @param  N       Number of variables.
        */
        void createVariablesTopDown(const int* bounds, unsigned N);

        /** Find the level of a given variable.
                @param  v   Variable to search for.

                @return 0,  if the variable was not found;
                        i,  with getVar(i) == v, otherwise.
        */
        unsigned findLevelOfVariable(const variable* v) const;

    private:
        //
        // Variable ordering/reordering members
        //

        // var_orders[0] is reserved to store the default variable order
        std::vector< std::shared_ptr<const variable_order> > var_orders;
        std::shared_ptr<const variable_order> default_var_order;

    public:
        //
        // Variable ordering/reordering methods
        //

        std::shared_ptr<const variable_order> makeVariableOrder(const int* order);
        std::shared_ptr<const variable_order> makeVariableOrder(const variable_order& order);
        inline std::shared_ptr<const variable_order> makeDefaultVariableOrder() {
            return default_var_order;
        }
        void cleanVariableOrders();

    private:
        //
        // Forest registry members
        //

        /// Registry of all forests with this domain.
        /// Stored as a set of FIDs.
        std::set<unsigned> forestReg;

        /// Registry of all relation forests with this domain.
        /// Stored as a set of FIDs.
        std::set<unsigned> relforestReg;

        /// Are we about to delete the domain.
        bool is_marked_for_deletion;


    private:
        //
        // Forest registry methods
        //

        /// Does the registry contain any forests?
        inline bool hasForests() const {
            return !forestReg.empty();
        }

        /// Add f to the forest registry.
        void registerForest(forest* f);

        /// Remove f from the forest registry.
        void unregisterForest(forest* f);

        friend class forest; // calls the above two methods

        /// Add f to the relation forest registry.
        void registerRelforest(relforest* f);

        /// Remove f from the relation forest registry.
        void unregisterRelforest(relforest* f);

        friend class relforest; // calls the above two methods

        /// Mark this domain, and all forests in our registry, for deletion
        void markForDeletion();

    public:

        /** Are we marked for deletion?
            Not sure if this is useful anywhere; right now it is used
            only for testing the domain list
        */
        inline bool isMarkedForDeletion() const {
            return is_marked_for_deletion;
        }


    protected:
        /// Constructor.
        domain(variable** v, unsigned N);

        /// Destructor.
        virtual ~domain();


    private:
        //
        // Members and methods for the registry of all domains.
        //

        /// domain before us in the registry
        domain* prev;
        /// domain after us in the registry
        domain* next;

        /// Global domain registry
        static domain* domain_list;

        /// Initialize the domain registry (on library initialization)
        static void initDomList();
        /// Mark all domains in the registry
        static void markDomList();
        /// Delete all domains in the registry (on library cleanup)
        static void deleteDomList();

        /*
            Initializer_list will call initDomList() and deleteDomList().
        */
        friend class initializer_list;

    public:

    // --------------------------------------------------------------------
    //
    // Ideas for future expansion.
    //

#if 0
        /**
            Insert a new variable.
                @param  lev     Level to insert above; use 0 for a
                                new bottom-most level.
                @param  v       Variable to insert.
        */
        inline void insertVariableAboveLevel(unsigned lev, variable* v)
        {
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }

        /**
            Add a new variable with bound 1.
            Can be used when the domain already has forests, in which case
            all forests are modified as appropriate.
                @param  below   Placement information: the new variable will
                                appear immediately above the level \a below.
        */
        inline void createVariable(int below)
        {
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }


        /** Remove a variable at the specified level.
            An error is thrown if the variable size is not 1.
            Use shrinkVariableBound() to make the bound 1.
            All forests are modified as appropriate.
                @param  lev   Level number.
        */
        inline void removeVariableAtLevel(unsigned lev)
        {
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }


        /** Swap the locations of variables in forests.
            I.e., changes the variable ordering of all forests with this domain.
                @param  lev1    Level of first variable.
                @param  lev2    Level of second variable.
        */
        inline void swapOrderOfVariables(unsigned lev1, unsigned lev2)
        {
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }

        /** Find the actual bound of a variable.
                @param  lev     Level number.
                @return         The smallest shrinkable bound before
                                information loss for variable at level \a lev.
        */
        inline unsigned findVariableBound(unsigned vh) const
        {
            throw error(error::NOT_IMPLEMENTED, __FILE__, __LINE__);
        }
#endif

};


#endif
