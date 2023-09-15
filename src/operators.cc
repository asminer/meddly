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

#include "operators.h"
#include "ops_builtin.h"
#include "opname.h"
#include "error.h"
#include "forest.h"

namespace MEDDLY {

    dd_edge operator+(const dd_edge &e1, const dd_edge &e2)
    {
        if (! e1.sameForest(e2)) {
            throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
        }

        dd_edge res(e1.getForest());

        if (res.getForest()->isRangeType(range_type::BOOLEAN)) {
            apply(UNION, e1, e2, res);
        } else {
            apply(PLUS, e1, e2, res);
        }
        return res;
    }

    dd_edge operator-(const dd_edge &e1, const dd_edge &e2)
    {
        if (! e1.sameForest(e2)) {
            throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
        }

        dd_edge res(e1.getForest());

        if (res.getForest()->isRangeType(range_type::BOOLEAN)) {
            apply(DIFFERENCE, e1, e2, res);
        } else {
            apply(MINUS, e1, e2, res);
        }
        return res;
    }

    dd_edge operator*(const dd_edge &e1, const dd_edge &e2)
    {
        if (! e1.sameForest(e2)) {
            throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
        }

        dd_edge res(e1.getForest());

        if (res.getForest()->isRangeType(range_type::BOOLEAN)) {
            apply(INTERSECTION, e1, e2, res);
        } else {
            apply(MULTIPLY, e1, e2, res);
        }
        return res;
    }

    dd_edge operator/(const dd_edge &e1, const dd_edge &e2)
    {
        if (! e1.sameForest(e2)) {
            throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
        }

        dd_edge res(e1.getForest());
        apply(DIVIDE, e1, e2, res);
        return res;
    }

    dd_edge operator&(const dd_edge &e1, const dd_edge &e2)
    {
        if (! e1.sameForest(e2)) {
            throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
        }

        dd_edge res(e1.getForest());
        apply(INTERSECTION, e1, e2, res);
        return res;
    }

    dd_edge operator|(const dd_edge &e1, const dd_edge &e2)
    {
        if (! e1.sameForest(e2)) {
            throw error(error::FOREST_MISMATCH, __FILE__, __LINE__);
        }

        dd_edge res(e1.getForest());
        apply(UNION, e1, e2, res);
        return res;
    }

    dd_edge operator!(const dd_edge &e)
    {
        dd_edge res(e.getForest());
        apply(COMPLEMENT, e, res);
        return res;
    }


    dd_edge operator+=(dd_edge &res, const dd_edge &upd)
    {
        if (res.getForest()->isRangeType(range_type::BOOLEAN)) {
            apply(UNION, res, upd, res);
        } else {
            apply(PLUS, res, upd, res);
        }
        return res;
    }

    dd_edge operator-=(dd_edge &res, const dd_edge &upd)
    {
        if (res.getForest()->isRangeType(range_type::BOOLEAN)) {
            apply(DIFFERENCE, res, upd, res);
        } else {
            apply(MINUS, res, upd, res);
        }
        return res;
    }

    dd_edge operator*=(dd_edge &res, const dd_edge &upd)
    {
        if (res.getForest()->isRangeType(range_type::BOOLEAN)) {
            apply(INTERSECTION, res, upd, res);
        } else {
            apply(MULTIPLY, res, upd, res);
        }
        return res;
    }

    dd_edge operator/=(dd_edge &res, const dd_edge &upd)
    {
        apply(DIVIDE, res, upd, res);
        return res;
    }

    dd_edge operator&=(dd_edge &res, const dd_edge &upd)
    {
        apply(INTERSECTION, res, upd, res);
        return res;
    }

    dd_edge operator|=(dd_edge &res, const dd_edge &upd)
    {
        apply(UNION, res, upd, res);
        return res;
    }

};
