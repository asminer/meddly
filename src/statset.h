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

#ifndef MEDDLY_STATSET_H
#define MEDDLY_STATSET_H

#include "defines.h"

namespace MEDDLY {
    struct statset;
};


/// Collection of various stats for performance measurement
struct MEDDLY::statset {
    /// Number of times we scanned for reachable nodes
    long reachable_scans;
    /// Number of times a dead node was resurrected.
    long reclaimed_nodes;
    /// Number of times the forest storage array was compacted
    long num_compactions;
    /// Number of times the garbage collector ran.
    long garbage_collections;

#ifdef TRACK_UNREACHABLE_NODES
    /// Current number of unreachable (disconnected) nodes
    long unreachable_nodes;
#endif

    /// Current number of connected nodes
    long active_nodes;

    /// Peak number of active nodes
    long peak_active;

    /// Constructor. Zeroes out everything.
    statset();

    // nice helpers

    inline void incActive(long b) {
        active_nodes += b;
        if (active_nodes > peak_active) {
            peak_active = active_nodes;
        }
        MEDDLY_DCASSERT(active_nodes >= 0);
    }

    inline void decActive(long b) {
        active_nodes -= b;
        MEDDLY_DCASSERT(active_nodes >= 0);
    }
};

#endif
