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

#ifndef MEDDLY_CT_ENTRY_TYPE_H
#define MEDDLY_CT_ENTRY_TYPE_H

#include "defines.h"

#include "old_meddly.h"

namespace MEDDLY {
    class ct_object;
    class ct_entry_type;
    class compute_table;

    class expert_forest;

    enum class ct_typeID {
        ERROR = 0,
        NODE = 1,
        INTEGER = 2,
        LONG = 3,
        FLOAT = 4,
        DOUBLE = 5,
        GENERIC = 6 // ct_object
    };
    inline int intOf(ct_typeID t) {
        return static_cast<typename std::underlying_type<ct_typeID>::type>(t);
    }

    union ct_entry_item {
        int I;
        unsigned int U;
        long L;
        unsigned long UL;
        node_handle N;
        float F;
        double D;
        ct_object* G;
    };

};

// ******************************************************************
// *                                                                *
// *                         ct_object class                        *
// *                                                                *
// ******************************************************************

/** Generic objects in compute tables.
    Used for things other than dd_edges and simple types.
*/
class MEDDLY::ct_object {
    public:
        ct_object();
        virtual ~ct_object();
        virtual opnd_type getType() = 0;
};


// ******************************************************************
// *                                                                *
// *                      ct_entry_type  class                      *
// *                                                                *
// ******************************************************************

/**
  Type information about entries.
  Usually there is one type of entry for each operation,
  but there could be more than one type.

  These are built by operations and then registered
  with the compute table.
*/
class MEDDLY::ct_entry_type {
    public:
        /** Constructor.
              @param  name    Name of the entry type; used only for displaying
                              CT entries (usually while debugging).

              @param  pattern Pattern for an entry.  The following characters
                              are supported in this string:
                                'N': node (in a forest)
                                'I': int
                                'L': long
                                'F': float
                                'D': double
                                'G': pointer to a ct_object
                                ':': separates key portion from result portion;
                                     must appear exactly once
                                '.': for repeating entries; can appear at most
                                     once.  Everything between '.' and ':' can
                                     repeat zero or more times.

              @throws INVALID_ARGUMENT if pattern is illegal
        */
        ct_entry_type(const char* name, const char* pattern);
        ~ct_entry_type();

        /**
            Set the forest for 'N' items in the pattern.
              @param  i   Slot.  Character i in the pattern must be 'N'.
              @param  f   Forest.
        */
        void setForestForSlot(unsigned i, expert_forest* f);

        /** Clear CT bits for any forests this entry type uses.
              @param  skipF   If skipF[i] is true, then we do nothing
                              for forests with ID i.  We set this to
                              true after clearing forest with ID i to
                              prevent clearing the bits twice.

              @param  N       Size of in_use array, for sanity checks.
        */
        void clearForestCTBits(bool* skipF, unsigned N) const;

        /** Notify forests that we're done marking CT bits.
            The forests can choose to start the sweep phase if they like.
              @param  whichF  If whichF[i] is true, then we notify the
                              forest with ID i, and set whichF[i] to false.
                              This prevents notifying a forest twice.

              @param  N       Size of in_use array, for sanity checks.
        */
        void sweepForestCTBits(bool* skipF, unsigned N) const;


        //
        // Inlines; mostly used by compute tables
        //

        /// Unique ID for use in compute tables
        inline unsigned getID() const { return etID; }

        /**
            Results might be overwritten.
            Indicate that in these entries, the result portion of the
            entry might be updated.
            The CT will make storage decisions based on this.
            If this is never called, then we assume results cannot
            be updated; this is true for most operations.
        */
        inline void mightUpdateResults() { updatable_result = true; }

        /**
            Is the result portion updatable?
        */
        inline bool isResultUpdatable() const { return updatable_result; }

        /**
            Name; for debugging/printing
        */
        inline const char* getName() const { return name; }


        /**
            Does this entry type allow repetitions in the key?
            I.e., was there a '.' in the pattern?
        */
        inline bool isRepeating() const { return len_kr_type; }


        /**
            Get the number of items in the key.
              @param  reps  Number of repetitions.
                            If this is not a repeating type,
                            then this is ignored.

              @return Total number of slots in the key.
        */
        inline unsigned getKeySize(unsigned reps) const {
            return len_ks_type + (reps * len_kr_type);
        }

        /**
            Get the number of bytes in the key.
              @param  reps  Number of repetitions.
                            If this is not a repeating type,
                            then this is ignored.

              @return Total number of bytes required for the key.
        */
        inline unsigned getKeyBytes(unsigned reps) const {
            return ks_bytes + (reps * kr_bytes);
        }


        /**
            Get the type for item i in the key.
            Automatically handles repetitions.
              @param  i   Slot number, between 0 and getKeySize().

              @param  t   On output, the type for item i.
              @param  f   If t is 'N', the forest for item i.
                          Otherwise, null.
        */
        inline void getKeyType(unsigned i, ct_typeID &t, expert_forest* &f)
        const {
            if (i<len_ks_type) {
                MEDDLY_DCASSERT(ks_type);
                MEDDLY_DCASSERT(ks_forest);

                t = ks_type[i];
                f = ks_forest[i];
            } else {
                MEDDLY_DCASSERT(len_kr_type);
                MEDDLY_DCASSERT(kr_type);
                MEDDLY_DCASSERT(kr_forest);

                i -= len_ks_type;
                i %= len_kr_type;
                t = kr_type[i];
                f = kr_forest[i];
            }
        }


        /**
            Get the type for item i in the key.
            Automatically handles repetitions.
              @param  i   Slot number, between 0 and getKeySize().
        */
        inline ct_typeID getKeyType(unsigned i) const {
            if (i<len_ks_type) {
                MEDDLY_DCASSERT(ks_type);
                return ks_type[i];
            } else {
                MEDDLY_DCASSERT(len_kr_type);
                MEDDLY_DCASSERT(kr_type);
                i -= len_ks_type;
                i %= len_kr_type;
                return kr_type[i];
            }
        }


        /**
            Get the forest for item i in the key.
            Automatically handles repetitions.
              @param  i   Slot number, between 0 and getKeySize().
              @return     Forest for that slot, or 0 if the type
                          is not 'N'.
        */
        inline expert_forest* getKeyForest(unsigned i) const {
            if (i<len_ks_type) {
                MEDDLY_DCASSERT(ks_forest);
                return ks_forest[i];
            } else {
                MEDDLY_DCASSERT(len_kr_type);
                MEDDLY_DCASSERT(kr_forest);
                i -= len_ks_type;
                i %= len_kr_type;
                return kr_forest[i];
            }
        }

        /**
            Get the number of items in the result
        */
        inline unsigned getResultSize() const { return len_r_type; }

        /**
            Get the number of bytes in the result
        */
        inline unsigned getResultBytes() const { return r_bytes; }


        /**
            Get the type for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().

              @param  t   On output, the type for item i.
              @param  f   If t is 'N', the forest for item i.
                          Otherwise, null.
        */
        inline void getResultType(unsigned i, ct_typeID &t, expert_forest* &f)
        const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, len_r_type);
            MEDDLY_DCASSERT(r_type);
            MEDDLY_DCASSERT(r_forest);
            t = r_type[i];
            f = r_forest[i];
        }


        /**
            Get the type for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().
        */
        inline ct_typeID getResultType(unsigned i) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, len_r_type);
            MEDDLY_DCASSERT(r_type);
            return r_type[i];
        }


        /**
            Get the forest for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().
              @return     Forest for that slot, or 0 if the type
                          is not 'N'.
        */
        inline expert_forest* getResultForest(unsigned i) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0, i, len_r_type);
            MEDDLY_DCASSERT(r_forest);
            return r_forest[i];
        }


        /// Mark for deletion
        inline void markForDeletion() { is_marked_for_deletion = true; }


        /// Unmark for deletion
        inline void unmarkForDeletion() { is_marked_for_deletion = false; }


        /// Should we remove all CT entries of this type?
        inline bool isMarkedForDeletion() const {
            return is_marked_for_deletion;
        }

    private:
        /// Unique ID, set by compute table
        unsigned etID;

        /// Name; for displaying CT entries
        const char* name;

        /// Starting portion of key pattern.
        ct_typeID* ks_type;

        /// Forests in starting portion of key.
        expert_forest** ks_forest;

        /// Length of ks_type and ks_forest arrays.
        unsigned len_ks_type;

        /// Total bytes in the starting portion of the key.
        unsigned ks_bytes;


        /// Repeating portion of key pattern (or null for no repeats).
        ct_typeID* kr_type;

        /// Forests in repeating portion of key (or null).
        expert_forest** kr_forest;

        /// Length of kr_type and kr_forest arrays (zero if no repeats).
        unsigned len_kr_type;

        /// Total bytes in the repeating portion of the key.
        unsigned kr_bytes;


        /// Result pattern
        ct_typeID* r_type;

        /// Forests in result
        expert_forest** r_forest;

        /// Length of r_type and r_forest arrays.
        unsigned len_r_type;

        /// Total bytes in the result.
        unsigned r_bytes;

        /// Can the result be changed later
        bool updatable_result;

        /// For deleting all entries of this type
        bool is_marked_for_deletion;

        friend class compute_table;
};

#endif
