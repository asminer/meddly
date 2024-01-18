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

#ifndef MEDDLY_VARIABLE_H
#define MEDDLY_VARIABLE_H

#include <vector>
#include <string>

namespace MEDDLY {
    class variable;
    class domain;
};

// ******************************************************************
// *                                                                *
// *                         variable class                         *
// *                                                                *
// ******************************************************************

/** Variable class.
    Abstract base class.
    A variable consists of an optional name, and a bound.
    A single variable object is used to describe both
    the primed and unprimed versions of the variable.

    Note1: variables are automatically deleted when
    removed from their last domain.

    Additional features are provided in the expert interface.
*/
class MEDDLY::variable {
    public:
        variable(int bound, std::string name);
        virtual ~variable();
        inline void setName(std::string newname) {
            name = newname;
        }
        inline bool hasName() const {
            return name.length()>0;
        }
        inline const std::string& getName() const {
            return name;
        }
        inline int getBound(bool primed) const {
            return primed ? pr_bound : un_bound;
        }
        inline bool isExtensible() const {
            return is_extensible;
        }

        /// Update our list of domains: add \a d.
        inline void addToList(domain* d) {
            domlist.push_back(d);
        }

        /// Update our list of domains: remove \a d.
        void removeFromList(const domain* d);

        /** Enlarge the possible values for a variable.
            This could modify all nodes in all forests, depending on the
            choice of reduction rule.
            @param  prime   If prime is true, enlarge the bound for
                            the primed variable only, otherwise both
                            the primed and unprimed are enlarged.
            @param  b       New bound, if less than the current bound
                            an error is thrown. If bound<=0, the variable
                            is marked as extensible, with initial bound as
                            abs(bound).  Note: an extensible variable has
                            a range [0 .. +infinity].
        */
        void enlargeBound(bool prime, int b);

        /** Shrink the possible values for a variable.
            This could modify all nodes in all forests, depending on the
            choice of reduction rule.
            @param  b       New bound, if more than the current bound
                            an error is thrown.
            @param  force   If \a b is too small, and information will be lost,
                            proceed anyway if \a force is true, otherwise
                            return an error code.
        */
        void shrinkBound(int b, bool force);

    private:
        // List of all domains this variable belongs to
        std::vector <domain*> domlist;

        std::string name;
        int un_bound;
        int pr_bound;
        bool is_extensible;
};


namespace MEDDLY {
    /** Front-end function to create a variable.
        Provided for backward compatability; one can instead
        simply use the constructor.
            @param  bound   The initial bound for the variable.
                            If bound<=0, the variable is marked as extensible,
                            with initial bound as abs(bound).
                            Note: an extensible variable has
                            range [0 .. +infinity].

            @param  name    Variable name (used only in display / debugging), or 0.
            @return A new variable, or 0 on error.
    */
    inline variable* createVariable(int bound, std::string name)
    {
        return new variable(bound, name);
    }
};


#endif
