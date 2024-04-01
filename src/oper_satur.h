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

#ifndef MEDDLY_OPER_SATUR_H
#define MEDDLY_OPER_SATUR_H

#include "oper.h"
#include "forest.h"

namespace MEDDLY {
    class ct_object;
    class saturation_operation;
    class forest;
};


// ******************************************************************
// *                                                                *
// *                   saturation_operation class                   *
// *                                                                *
// ******************************************************************

/** Mechanism to apply a saturation operation in a specific forest.
    Basically, a unary operation but without needing to 'cache'
    operations by forests, because the relation is fixed.
    Specific operations will be derived from this class.
*/
class MEDDLY::saturation_operation : public operation {
    public:
        saturation_operation(const char* name, unsigned et_slots,
                forest* arg, forest* res);

    protected:
        virtual ~saturation_operation();

    protected:
        // Fairly standard checks; call these in the operator's constructor.
        //
        /// Make sure all three domains are the same.
        inline void checkDomains(const char* file, unsigned line) const
        {
            if  ( (argF->getDomain() != resF->getDomain()) )
            {
                throw error(error::DOMAIN_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments set/relation status matches
        inline void checkAllRelations(const char* file, unsigned line,
                set_or_rel a) const
        {
            if  (
                    (argF->isForRelations() != a)  ||
                    (resF->isForRelations() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure all arguments match the edge labeling rule
        inline void checkAllLabelings(const char* file, unsigned line,
                edge_labeling a) const
        {
            if  (
                    (argF->getEdgeLabeling() != a)  ||
                    (resF->getEdgeLabeling() != a)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }
        /// Make sure the arguments match the range type
        inline void checkAllRanges(const char* file, unsigned line,
                range_type rt) const
        {
            if  (
                    (argF->getRangeType() != rt)  ||
                    (resF->getRangeType() != rt)
                )
            {
                throw error(error::TYPE_MISMATCH, file, line);
            }
        }

    protected:
        virtual bool checkForestCompatibility() const;

    public:
        virtual void compute(const dd_edge &arg, dd_edge &res) = 0;

        virtual bool isReachable(const dd_edge& a, const dd_edge& constraint);

    protected:
        forest* argF;
        forest* resF;
};


#endif // #include guard
