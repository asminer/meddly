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



/*
 * \brief A timer class.
 *
 * Timer is started when a timer object is initialized.
 * note_time() is used to note the current time.
 * get_last_interval() is the time between the last two instances when the
 * time was noted (in microseconds)
 *
 * */

#ifndef TIMERS_H
#define TIMERS_H

#include "../config.h"

#define SHOW_WHICH_TIMER

#ifdef SHOW_WHICH_TIMER
#include <iostream>
#endif

#if HAVE_GETRUSAGE

/*
    Timer based on getrusage,
    preferred version.
*/

#include <sys/resource.h>

class timer
{
        struct rusage usage;
	    struct timeval prev_time;
	    long last_interval;

    public:
	    timer()
	    {
#ifdef SHOW_WHICH_TIMER
            std::cerr << "getrusage timer\n";
#endif
            getrusage(RUSAGE_SELF, &usage);
		    prev_time = usage.ru_utime;
	    }

	    inline void note_time()
	    {
            getrusage(RUSAGE_SELF, &usage);
            update(usage.ru_utime);
	    }

	    inline long get_last_interval()
	    {
		    return last_interval;
	    }

        inline double get_last_seconds()
        {
            return last_interval / 1000000.0;
        }

    private:
        inline void update(const struct timeval &curr)
        {
		    last_interval = (curr.tv_sec - prev_time.tv_sec) * 1000000;
		    last_interval += curr.tv_usec - prev_time.tv_usec;
		    prev_time = curr;
        }
};

#else
#if HAVE_GETTIMEOFDAY

/*
    Timer based on gettimeofday.
*/

#include <sys/time.h>
#include <time.h>

class timer
{
	    struct timezone time_zone;
	    struct timeval curr_time;
	    struct timeval prev_time;
	    long last_interval;

    public:
	    timer()
	    {
#ifdef SHOW_WHICH_TIMER
            std::cerr << "gettimeofday timer\n";
#endif
		    gettimeofday(&curr_time, &time_zone);
		    prev_time = curr_time;
	    }

	    inline void note_time()
	    {
		    gettimeofday(&curr_time, &time_zone);
            update(curr_time);
	    }

	    inline long get_last_interval()
	    {
		    return last_interval;
	    }

        inline double get_last_seconds()
        {
            return last_interval / 1000000.0;
        }

private:
        inline void update(const struct timeval &curr)
        {
		    last_interval = (curr.tv_sec - prev_time.tv_sec) * 1000000;
		    last_interval += curr.tv_usec - prev_time.tv_usec;
		    prev_time = curr;
        }
};


#else

/*
    Can't make a timer.
    Fail safe version.
*/

class timer
{
    public:
	    timer()
	    {
#ifdef SHOW_WHICH_TIMER
            std::cerr << "failsafe (do nothing) timer\n";
#endif
	    }

	    inline void note_time()
	    {
	    }

	    inline long get_last_interval()
	    {
		    return 0;
	    }

        inline double get_last_seconds()
        {
            return 0.0;
        }
};

#endif
#endif
#endif

