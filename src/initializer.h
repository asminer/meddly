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

#ifndef MEDDLY_INITIALIZER_H
#define MEDDLY_INITIALIZER_H

namespace MEDDLY {
    class initializer_list;

    /**
        Build list of initializers for Meddly.
        Custom-built initialization lists will usually include this list.
            @param    prev  Initializers to execute before the default list;
                            can be null.

            @return   List of initializers.

        TBD: maybe move this into "library.h" and "library.cc"?
    */
    initializer_list* defaultInitializerList(initializer_list* prev);


    /** Initialize the library with custom settings.
        Should be called before using any other functions.
            @param  L   List of initializers.  Will execute the "setup()"
                        methods in order now, and the "cleanup()" methods
                        in reverse order on library cleanup.
    */
    void initialize(initializer_list* L);

};

// ******************************************************************
// *                                                                *
// *                     initializer_list class                     *
// *                                                                *
// ******************************************************************

/** Mechanism for initializing and/or cleaning up library structures.
    Any user additions to the library should utilize this class.
    Derive a class from this one, provide the \a setup and \a cleanup
    methods.
    Implementation in meddly.cc
*/
class MEDDLY::initializer_list {
    public:
        /**
            Constructor.
            Takes the initializer(s) to run before this one.
            Cleanup runs in the reverse order.
        */
        initializer_list(initializer_list* previous);
        virtual ~initializer_list();

        /**
            Run all setup methods for the list of initializers,
            "previous first".
        */
        void setupAll();

        /**
            Run all cleanup methods for the list of initializers,
            "previous last".
        */
        void cleanupAll();

    protected:
        virtual void setup() = 0;
        virtual void cleanup() = 0;

    private:
        initializer_list* previous;
};

#endif
