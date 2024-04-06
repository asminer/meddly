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
#include "edge_value.h"
#include "oper.h"
#include <type_traits>

// #define OLD_TYPE_IFACE

namespace MEDDLY {
    class ct_object;
    class ct_entry_type;
    class compute_table;

    class forest;
    class output;

    enum class ct_typeID {
        ERROR   = 0,
        NODE    = 1,
        INTEGER = 2,
        LONG    = 3,
        FLOAT   = 4,
        DOUBLE  = 5,
        GENERIC = 6 // ct_object
    };

    // TBD should this be unsigned?
    inline int intOf(ct_typeID t) {
        return static_cast<typename std::underlying_type<ct_typeID>::type>(t);
    }

    class ct_itemtype;

    // TBD move this because it is entry related
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
// *                        ct_itemtype class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::ct_itemtype {
    public:
        ct_itemtype() {
            type = ct_typeID::ERROR;
            nodeFor = nullptr;
        }
        /** Set type based on the following codes:
                'N': node (in a forest)
                'I': int
                'L': long
                'F': float
                'D': double
                'G': pointer to a ct_object
        */
        ct_itemtype(char code);

        ct_itemtype(forest* f) {
            type = ct_typeID::NODE;
            nodeFor = f;
        }
        ct_itemtype(ct_typeID t) {
            MEDDLY_DCASSERT(t != ct_typeID::NODE);
            type = t;
            nodeFor = nullptr;
        }
        ct_itemtype(edge_type et) {
            nodeFor = nullptr;
            switch (et) {
                case edge_type::INT:
                    type = ct_typeID::INTEGER;
                    break;

                case edge_type::LONG:
                    type = ct_typeID::LONG;
                    break;

                case edge_type::FLOAT:
                    type = ct_typeID::FLOAT;
                    break;

                case edge_type::DOUBLE:
                    type = ct_typeID::DOUBLE;
                    break;

                default:
                    type = ct_typeID::ERROR;
            }
        }

        /// Get the type of this item.
        inline ct_typeID getType() const {
            return type;
        }
        /// Check if the item has the specified type.
        inline bool hasType(ct_typeID t) const {
            return t == type;
        }
        /// Get the type as a character;
        /// for the old-style 'pattern' interface.
        char getTypeChar() const;

        /// Make sure this is a node type item, and return its forest.
        inline forest* getForest() const {
            MEDDLY_DCASSERT(ct_typeID::NODE == type);
            return nodeFor;
        }
        /// Return the raw forest (no checks)
        inline forest* rawForest() const {
            return nodeFor;
        }
        /// Check if this item is associated with forest f.
        inline bool hasForest(const forest* f) const {
            return f == nodeFor;
        }
        /// Check if this item has an associated forest.
        inline bool hasForest() const {
            return nodeFor;
        }
        /// Number of bytes required to store this item in a CT
        inline unsigned bytes() const {
            static const unsigned sizes[] = {
                0,                      // ERROR   = 0,
                sizeof(node_handle),    // NODE    = 1,
                sizeof(int),            // INTEGER = 2,
                sizeof(long),           // LONG    = 3,
                sizeof(float),          // FLOAT   = 4,
                sizeof(double),         // DOUBLE  = 5,
                sizeof(ct_object*)      // GENERIC = 6 // ct_object
            };
            return sizes[ getTypeInt() ];
        }
        /// Number of integer slots needed to store this item in a CT
        inline unsigned intSlots() const {
            static const unsigned sizes[] = {
                0,                                  // ERROR   = 0,
                sizeof(node_handle) / sizeof(int),  // NODE    = 1,
                sizeof(int)         / sizeof(int),  // INTEGER = 2,
                sizeof(long)        / sizeof(int),  // LONG    = 3,
                sizeof(float)       / sizeof(int),  // FLOAT   = 4,
                sizeof(double)      / sizeof(int),  // DOUBLE  = 5,
                sizeof(ct_object*)  / sizeof(int)   // GENERIC = 6 // ct_object
            };
            return sizes[ getTypeInt() ];
        }

        /// Set the forest; should be called when building an operation.
        inline void setForest(forest* f) {
            if (type != ct_typeID::NODE) {
                throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
            }
            MEDDLY_DCASSERT(!nodeFor);
            nodeFor = f;
        }

    public:
        void show(output &s) const;

    protected:
        inline unsigned getTypeInt() const {
            unsigned u =
                static_cast <typename std::underlying_type<ct_typeID>::type>
                    (type);
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, u, 7u);
            return u;
        }

    private:
        ct_typeID   type;
        forest      *nodeFor;
};

// ******************************************************************
// *                                                                *
// *                         ct_object class                        *
// *                                                                *
// ******************************************************************

// TBD move this
//
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
        friend class compute_table;
    public:
#ifdef ALLOW_DEPRECATED_0_17_6
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

        /**
            Set the forest for 'N' items in the pattern.
              @param  i   Slot.  Character i in the pattern must be 'N'.
              @param  f   Forest.
        */
        void setForestForSlot(unsigned i, forest* f);
#endif
        /**
            New constructor.

              @param  name    Name of the entry type; used only for displaying
                              CT entries (usually while debugging).

            Use this constructor if you are going to specify the
            key and result types using methods below:

                set_fkey()      for setting the fixed portion of the key,
                                if any.
                set_rkey()      for setting the repeating portion of the key,
                                if any.
                set_result()    for setting the result portion of the key.
        */
        ct_entry_type(const char* name);
        ~ct_entry_type();

        inline void set_fkey(ct_itemtype a) { }
        inline void set_fkey(ct_itemtype a, ct_itemtype b) { }
        inline void append_fkey(ct_itemtype a) { }

        inline void set_rkey(ct_itemtype a) { }
        inline void set_rkey(ct_itemtype a, ct_itemtype b) { }
        inline void append_rkey(ct_itemtype a) { }

        inline void set_result(ct_itemtype a) { }
        inline void set_result(ct_itemtype a, ct_itemtype b) { }
        inline void append_result(ct_itemtype a) { }

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
        inline bool isRepeating() const { return key_repeating.size(); }


        /**
            Get the number of items in the key.
              @param  reps  Number of repetitions.
                            If this is not a repeating type,
                            then this is ignored.

              @return Total number of slots in the key.
        */
        inline unsigned getKeySize(unsigned reps) const {
            return key_fixed.size() + (reps * key_repeating.size());
        }

        /**
            Get the number of bytes in the key.
              @param  reps  Number of repetitions.
                            If this is not a repeating type,
                            then this is ignored.

              @return Total number of bytes required for the key.
        */
        inline unsigned getKeyBytes(unsigned reps) const {
            return fixed_bytes + (reps * repeating_bytes);
        }


#ifdef OLD_TYPE_IFACE
        /**
            Get the type for item i in the key.
            Automatically handles repetitions.
              @param  i   Slot number, between 0 and getKeySize().

              @param  t   On output, the type for item i.
              @param  f   If t is 'N', the forest for item i.
                          Otherwise, null.
        */
        inline void getKeyType(unsigned i, ct_typeID &t, forest* &f)
        const {
            if (i<key_fixed.size()) {
                t = key_fixed[i].getType();
                f = key_fixed[i].rawForest();
            } else {
                MEDDLY_DCASSERT(key_repeating.size());

                i -= key_fixed.size();
                i %= key_repeating.size();
                t = key_repeating[i].getType();
                f = key_repeating[i].rawForest();
            }
        }


        /**
            Get the type for item i in the key.
            Automatically handles repetitions.
              @param  i   Slot number, between 0 and getKeySize().
        */
        inline ct_typeID getKeyType(unsigned i) const {
            if (i<key_fixed.size()) {
                return key_fixed[i].getType();
            } else {
                MEDDLY_DCASSERT(key_repeating.size());
                i -= key_fixed.size();
                i %= key_repeating.size();
                return key_repeating[i].getType();
            }
        }


        /**
            Get the forest for item i in the key.
            Automatically handles repetitions.
              @param  i   Slot number, between 0 and getKeySize().
              @return     Forest for that slot, or 0 if the type
                          is not 'N'.
        */
        inline forest* getKeyForest(unsigned i) const {
            if (i<key_fixed.size()) {
                return key_fixed[i].rawForest();
            } else {
                MEDDLY_DCASSERT(key_repeating.size());
                i -= key_fixed.size();
                i %= key_repeating.size();
                return key_repeating[i].rawForest();
            }
        }
#else
        /**
            Get the type for item i in the key.
            Automatically handles repetitions.
              @param  i   Slot number, between 0 and getKeySize().
              @return     Type info for slot i
        */
        inline const ct_itemtype& getKeyType(unsigned i) const {
            if (i<key_fixed.size()) {
                return key_fixed[i];
            } else {
                MEDDLY_DCASSERT(key_repeating.size());
                i -= key_fixed.size();
                i %= key_repeating.size();
                return key_repeating[i];
            }
        }

#endif

        /**
            Get the number of items in the result
        */
        inline unsigned getResultSize() const { return result.size(); }

        /**
            Get the number of bytes in the result
        */
        inline unsigned getResultBytes() const { return result_bytes; }

#ifdef OLD_TYPE_IFACE
        /**
            Get the type for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().

              @param  t   On output, the type for item i.
              @param  f   If t is 'N', the forest for item i.
                          Otherwise, null.
        */
        inline void getResultType(unsigned i, ct_typeID &t, forest* &f)
        const
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i, unsigned(result.size()));
            t = result[i].getType();
            f = result[i].rawForest();
        }


        /**
            Get the type for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().
        */
        inline ct_typeID getResultType(unsigned i) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i, unsigned(result.size()));
            return result[i].getType();
        }


        /**
            Get the forest for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().
              @return     Forest for that slot, or 0 if the type
                          is not 'N'.
        */
        inline forest* getResultForest(unsigned i) const {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__, 0u, i, unsigned(result.size()));
            return result[i].rawForest();
        }
#else
        /**
            Get the type for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().
        */
        inline const ct_itemtype& getResultType(unsigned i) const
        {
            MEDDLY::CHECK_RANGE(__FILE__, __LINE__,
                    0u, i, unsigned(result.size()));
            return result[i];
        }
#endif

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

        /// Fixed, initial portion of the key.
        std::vector <ct_itemtype> key_fixed;

        /// Total bytes in the starting portion of the key.
        unsigned fixed_bytes;

        /// Repeating portion of the key.
        std::vector <ct_itemtype> key_repeating;

        /// Total bytes in the repeating portion of the key.
        unsigned repeating_bytes;

        /// Result pattern
        std::vector <ct_itemtype> result;

        /// Total bytes in the result.
        unsigned result_bytes;

        /// Can the result be changed later
        bool updatable_result;

        /// For deleting all entries of this type
        bool is_marked_for_deletion;
};

#endif
