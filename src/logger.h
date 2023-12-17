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

#ifndef LOGGER_H
#define LOGGER_H

namespace MEDDLY {
    class logger;
    class forest;
};

/**
    Abstract base class for logging of forest stats.
    The idea is to allow an external tool to read this information,
    and display stats in real time.  Potentially much more detailed
    that the stats collected above.
*/
class MEDDLY::logger {
        bool nfix;
        bool node_counts;
        bool time_stamps;

        long startsec;
        long startusec;
    public:
        logger();
        virtual ~logger();

    /*
        Settings.
        Cannot be changed once we attach to a forest.
    */
    public:
        inline bool recordingNodeCounts() const {
            return node_counts;
        }
        inline void recordNodeCounts() {
            if (nfix) node_counts = true;
        }
        inline void ignoreNodeCounts() {
            if (nfix) node_counts = false;
        }


        inline bool recordingTimeStamps() const {
            return time_stamps;
        }
        inline void recordTimeStamps() {
            if (nfix) time_stamps = true;
        }
        inline void ignoreTimeStamps() {
            if (nfix) time_stamps = false;
        }

    /*
        Hooks, used in various places.
        Must be overridden in derived classes.
    */
    public:
        /**
            Insert a comment string.
            May be ignored depending on the file format.

              @param  comment   String to insert.
        */
        virtual void addComment(const char* comment) = 0;

        /**
            Start a new phase of computation.
            May be ignored depending on the file format.

              @param  f         Forest this applies to.
              @param  comment   Info to display about this phase.
        */
        virtual void newPhase(const forest* f, const char* comment) = 0;

        /**
            Called once, when the logger is attached to a forest.
            Must call method fixLogger().

              @param  f       Forest info to log
              @param  name    Forest name; better for displaying if we
                              have multiple forests.
        */
        virtual void logForestInfo(const forest* f, const char* name) = 0;

        /**
            Change node counts by specified amounts.
        */
        virtual void addToActiveNodeCount(const forest* f, int level, long delta) = 0;

    /*
        Helpers for derived classes
    */
    protected:
        /* Use this for generating timestamps. */
        void currentTime(long &sec, long &usec);

        /* Call this inside logForestInfo() */
        inline void fixLogger() {
            nfix = false;
        }

};

#endif
