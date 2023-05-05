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

    // ******************************************************************
    // *                  library management functions                  *
    // ******************************************************************

    /**
        Build list of initializers for Meddly.
        Custom-built initialization lists will usually include this list.
            @param    prev  Initializers to execute before the default list;
                            can be null.

            @return   List of initializers.
    */
    initializer_list* defaultInitializerList(initializer_list* prev);


    /** Initialize the library with custom settings.
        Should be called before using any other functions.
            @param  L   List of initializers.  Will execute the "setup()"
                        methods in order now, and the "cleanup()" methods
                        in reverse order on library cleanup.
    */
    void initialize(initializer_list* L);

    /** Initialize the library with default settings.
        See meddly_expert.h for functions to initialize
        the library with non-default settings.
        Should be called before using any other functions.
    */
    void initialize();

    /** Clean up the library.
        Can be called to free memory used by the library;
        after it is called, the library may be initialized again.
    */
    void cleanup();

    /** Get the information about the library.
        @param  what  Determines the type of information to obtain.
        @return A human-readable information string.
                The string depends on parameter \a what, as follows.
                0: Title string, e.g., "MDD library version 0.0.0"
                1: Copyright string
                2: License string
                3: Reference url
                4: String with library features
                5: Release date
                Anything else: null string.
    */
    const char* getLibraryInfo(int what = 0);

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
            Initialize the library from the given list.
        */
        static void initializeLibrary(initializer_list* L);

        /// Did we initialize the library?
        static inline bool libraryIsRunning() { return isRunning; }

        /**
            Clean up the library using the same list
            passed to initializeLibrary().
        */
        static void cleanupLibrary();

    protected:
        virtual void setup() = 0;
        virtual void cleanup() = 0;

    private:
        initializer_list* previous;
        static initializer_list* meddlyInitializers;
        static bool isRunning;
};

#endif
