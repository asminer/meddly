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

#ifndef MEDDLY_FOREST_LEVELS_H
#define MEDDLY_FOREST_LEVELS_H

#include "defines.h"

namespace MEDDLY {
    class MDD_levels;
    class MXD_levels;
}

// ******************************************************************
// *                                                                *
// *                        MDD_levels class                        *
// *                                                                *
// ******************************************************************

/**
    Small class for MDD level behavior,
    used for template parameters in operations.
 */
class MEDDLY::MDD_levels {
    public:

        /// Go "down a level"
        static inline int downLevel(int k) {
            return k-1;
        }

        ///  Go "up a level"
        static inline int upLevel(int k) {
            return k+1;
        }

        /// Determine the top of two levels
        static inline int topLevel(int k1, int k2) {
            return MAX(k1, k2);
        }
};


// ******************************************************************
// *                                                                *
// *                        MXD_levels class                        *
// *                                                                *
// ******************************************************************

/**
    Small class for MXD level behavior,
    used for template parameters in operations.
 */
class MEDDLY::MXD_levels {
    public:

        /// Go "down a level"
        static inline int downLevel(int k) {
            return (k>0) ? (-k) : (-k-1);
        }

        ///  Go "up a level"
        static inline int upLevel(int k) {
            return (k<0) ? (-k) : (-k-1);
        }

        /// Determine the top of two levels
        static inline int topLevel(int k1, int k2) {
            if (ABS(k1) == ABS(k2)) {
                return MAX(k1, k2);
            }
            return (ABS(k1) > ABS(k2)) ? k1 : k2;
        }

        /// Determine the top unprimed level of two levels.
        /// This is ABS(topLevel(k1, k2)) but computed more efficiently.
        static inline int topUnprimed(int k1, int k2) {
            return MAX(ABS(k1), ABS(k2));
        }
};

#endif // include guard

