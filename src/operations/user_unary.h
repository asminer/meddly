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

#ifndef MEDDLY_USER_UNARY_H
#define MEDDLY_USER_UNARY_H

#include "../oper_unary.h"

#define TRACE_USER

namespace MEDDLY {
    class user_unary_operation : public unary_operation {
        public:
            user_unary_operation(forest* arg, forest* res,
                    user_defined_unary F);

            virtual ~user_unary_operation();

        private:
            user_defined_unary F;
            ct_entry_type* ct;
#ifdef TRACE_USER
            ostream_output out;
            unsigned top_count;
#endif
    };
};

#endif

