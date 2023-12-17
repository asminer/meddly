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

#include "logger.h"

// ******************************************************************
// *                                                                *
// *                         logger methods                         *
// *                                                                *
// ******************************************************************

MEDDLY::logger::logger()
{
    nfix = true;
    /*
        All other settings to false;
        the default logger does nothing!
    */
    node_counts = false;
    time_stamps = false;

    /* Do this now, regardless */
#ifdef HAVE_SYS_TIME_H
    static struct timeval curr_time;
    static struct timezone tz;
    gettimeofday(&curr_time, &tz);
    startsec = curr_time.tv_sec;
    startusec = curr_time.tv_usec;
#else
    startsec = 0;
    startusec = 0;
#endif
}

MEDDLY::logger::~logger()
{
}

void MEDDLY::logger::currentTime(long &sec, long &usec)
{
#ifdef HAVE_SYS_TIME_H
    static struct timeval curr_time;
    static struct timezone tz;
    gettimeofday(&curr_time, &tz);
    sec = curr_time.tv_sec - startsec;
    usec = curr_time.tv_usec - startusec;
    if (usec<0) {
        sec--;
        usec += 1000000;
    }
#else
    sec = 0;
    usec = 0;
#endif
}

