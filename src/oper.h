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

#ifndef MEDDLY_OPER_H
#define MEDDLY_OPER_H

#include "defines.h"
#include "io.h"
#include "error.h"

#include <vector>

namespace MEDDLY {
    class operation;
    class compute_table;
    class forest;

    class initializer_list;
    class ct_initializer;
    class ct_entry_type;
    class ct_entry_result;

    /// Argument and result types for apply operations.
    enum class opnd_type {
        FOREST      = 0,
        BOOLEAN     = 1,
        INTEGER     = 2,
        REAL        = 3,
        HUGEINT     = 4
    };

    /// Progress callback functions, used by operation class.
    typedef void (*operation_progress)(unsigned itno, char status);
};

// ******************************************************************
// *                                                                *
// *                        operation  class                        *
// *                                                                *
// ******************************************************************

/** Generic operation.
    Operations are tied to specific forests.
    Necessary for compute table entries.
*/
class MEDDLY::operation {
        friend class initializer_list;
    public:
        struct option {
                char    type;
                char    opt_c;
                long    opt_i;
                double  opt_f;
            public:
                option() {
                    type=' ';
                }
                option(char c) {
                    type='c';
                    opt_c = c;
                }
                option(long i) {
                    type='i';
                    opt_i = i;
                }
                option(double f) {
                    type='f';
                    opt_f = f;
                }
        };
    public:
        /** Constructor (OLD)
                @param  n           Operation name, for debugging
                @param  et_slots    Number of different compute table entry
                                    types used by this operation.
                                    Derived class constructors must register
                                    exactly this many entry types.
        */
        operation(const char* n, unsigned et_slots=0);

        /// New constructor.
        operation();


        /** Destructor.
            Safe to call directly now :)
        */
        virtual ~operation();

#ifdef ALLOW_DEPRECATED_0_17_6
        /// Safely destroy the given operation.
        inline static void destroy(operation* op) {
            delete op;
        }
#endif

        /// Safely destroy all operations associated with the given forest.
        /// Called in forest destructor.
        static void destroyAllWithForest(const forest* f);

        /// Get the name of this operation; for display / debugging
        inline const char* getName() const { return name; }


        /// Set the function to call for progress updates.
        /// Only used by some operations, like traditional
        /// reachability set construction.
        inline void setProgressNotifier(operation_progress p) {
            np_func = p;
        }

        /// Does the operation use the progress notifier?
        inline bool usesProgressNotifier() const {
            return uses_progress;
        }

        /// How many user-adjustable options are there for this
        inline unsigned getNumOptions() const {
            return num_options;
        }

        /** Get the option type.
                @param  i   Option number, between 0 and #options-1.
                @return One of the following.
                    ' ':    There is no option.
                    'c':    The option type is 'char'
                    'i':    The option type is 'long'
                    'f':    The option type is 'double'
         */
        inline char getOptionType(unsigned i) const {
            if (i >= num_options) return ' ';
            return options[i];
        }

        /** Set an option.
                @param  i   Option number, between 0 and #options-1.
                @param  x   Value
                @return true on success, false otherwise.
                @throws INVALID_OPTION error on type mismatch.
        */
        inline void setOption(unsigned i, option x) {
            if (getOptionType(i) != x.type) {
                throw error(error::INVALID_OPTION, __FILE__, __LINE__);
            }
            _setOption(i, x);
        }

        /** Get an actual option. Override in derived classes;
            the default behavior returns an empty option.
        */
        virtual option getOption(unsigned i) const;

        /// Display the operation name and options.
        virtual void showNameAndOptions(output &s) const;

    protected:
        /// Set the name after constructing.
        inline void setName(const char* n) {
            if (!n) return;
            MEDDLY_DCASSERT(!name);
            name = n;
        }

        /// Set the option types. Default is nullptr.
        inline void setOptionTypes(const char* optypes) {
            options = optypes;
            num_options = 0;
            while (options[num_options]) {
                ++num_options;
            }
        }

        /** Set an actual option. Override in derived classes;
            the default behavior throws an error.
        */
        virtual bool _setOption(unsigned i, option x);

        /// Indicate that the progress notifier might be used.
        /// (Default: it's never used.)
        inline void usesProgress(bool up = true) {
            uses_progress = up;
        }

        /// Call the progress function, if there is one.
        /// Convention:
        ///     @param  it      Iteration number if that makes sense;
        ///                     otherwise, operation dependent.
        ///
        ///     @param  status  Where are we in the iteration:
        ///                         ' ': just starting
        ///                         ';': just finishing
        ///                     all others are operation dependent.
        ///
        /// Most operations will never call this.
        ///
        inline void notifyProgress(unsigned it, char status) const
        {
            if (np_func) {
                np_func(it, status);
            }
        }

        void registerInForest(forest* f);
        void unregisterInForest(forest* f);

#ifdef ALLOW_DEPRECATED_0_17_6
        void registerEntryType(unsigned slot, ct_entry_type* et);
        void buildCTs();
#endif


    // ------------------------------------------------------------
    public:  // public  methods for the operation registry
    // ------------------------------------------------------------
        /// Unique ID > 0 for the operation.
        inline unsigned getID() const { return oplist_index; }

        /// Find the operation with the given ID.
        inline static operation* getOpWithID(unsigned i) {
#ifdef DEVELOPMENT_CODE
            return op_list.at(i);
#else
            return op_list[i];
#endif
        }

        /// How many operations are active?
        inline static unsigned getOpListSize() {
            return op_list.size();
        }


    // ------------------------------------------------------------
    private: // private methods for the operation registry
    // ------------------------------------------------------------
        /// Called during library initialization.
        static void initializeStatics();

        /// Called during library cleanup.
        static void destroyAllOps();

        static void registerOperation(operation &o);
        static void unregisterOperation(operation &o);

    // ------------------------------------------------------------
    private: // private members for the operation registry
    // ------------------------------------------------------------
        unsigned oplist_index;
        static std::vector <operation*> op_list;
        static std::vector <unsigned>   free_list;

    //
    // Members
    //

    protected:
#ifdef ALLOW_DEPRECATED_0_17_6
        /// Compute table to use (for entry type 0), if any.
        compute_table* CT0;
        /// Array of compute tables, one per entry type.
        compute_table** CT;

        //
        // TBD: remove all of these
        //

        /** Array of entry types.
            Owned by the compute_table class; we have
            these pointers for convenience.
        */
        ct_entry_type** etype;
        /** Array of entry results.
            Use these during computation.
            We only ever need one result per entry type.
        */
        ct_entry_result* CTresult;

        /**
            Number of entry_types needed by this operation.
            This gives the dimension of arrays CT and etype.
        */
        unsigned num_etids;
#endif

    private:
        const char* name;
        const char* options;
        unsigned num_options;
        /// List of forest IDs associated with this operation.
        std::vector <unsigned> FList;

        operation_progress np_func;

        bool uses_progress;
};

#endif // #include guard
