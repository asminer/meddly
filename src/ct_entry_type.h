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
#include "compute_table.h"
#include "forest.h"

#define DELETE_ON_ZERO

namespace MEDDLY {
    class ct_object;
    class ct_entry_type;
#ifdef ALLOW_DEPRECATED_0_17_6
    class ct_entry_key;
    class ct_entry_result;
#endif
    class ct_vector;
    // class compute_table;

    class output;
    class relforest;

    enum class ct_typeID {
        ERROR   = 0,
        NODE    = 1,
        INTEGER = 2,
        LONG    = 3,
        FLOAT   = 4,
        DOUBLE  = 5,
        GENERIC = 6 // ct_object
    };

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

        unsigned raw[2];
    };

};

// #define USE_FID

// ******************************************************************
// *                                                                *
// *                        ct_itemtype class                       *
// *                                                                *
// ******************************************************************

class MEDDLY::ct_itemtype {
    public:
        /// Default constructor
        ct_itemtype() {
            type = ct_typeID::ERROR;
            typeUpdate(nullptr);
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

        /// Set the type to 'node' and assign the forest.
        ct_itemtype(forest* f) {
            type = ct_typeID::NODE;
            typeUpdate(f);
        }
        /// Set the type from the actual type enum.
        ct_itemtype(ct_typeID t) {
            ASSERT(__FILE__, __LINE__, t != ct_typeID::NODE);
            type = t;
            typeUpdate(nullptr);
        }
        /// Set the type from an edge value type
        ct_itemtype(edge_type et) {
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
            typeUpdate(nullptr);
        }

        /// Get the type of this item.
        inline ct_typeID getType() const {
            return type;
        }
        /// Check if the item has the specified type.
        inline bool hasType(ct_typeID t) const {
            return t == type;
        }
        /// Are we a node type (most common to check for)
        inline bool hasNodeType() const {
            return ct_typeID::NODE == type;
        }

        /// Get the type as a character;
        /// for the old-style 'pattern' interface.
        char getTypeChar() const;

        /// Return the raw forest (no checks)
        inline forest* rawForest() const {
#ifdef USE_FID
            return forest::getForestWithID(nodeFID);
#else
            return nodeF;
#endif
        }
        /// Check if this item is associated with forest f.
        inline bool hasForest(const forest* f) const {
            ASSERT(__FILE__, __LINE__, f);
#ifdef USE_FID
            return f->FID() == nodeFID;
#else
            return f == nodeF;
#endif
        }
        /// Number of integer slots needed to store this item in a CT
        inline unsigned intslots() const {
            return 1+twoslots;
        }
        /// Are two integer slots required
        inline bool requiresTwoSlots() const {
            return twoslots;
        }
        /// Should this type be hashed
        inline bool shouldBeHashed() const {
            return should_hash;
        }

        /// Set the forest; should be called when building an operation.
        inline void setForest(forest* f) {
            if (type != ct_typeID::NODE) {
                throw error(error::INVALID_ARGUMENT, __FILE__, __LINE__);
            }
#ifdef USE_FID
            ASSERT(__FILE__, __LINE__, !nodeFID);
            nodeFID = f ? f->FID() : 0;
#else
            ASSERT(__FILE__, __LINE__, !nodeF);
            nodeF = f;
#endif
        }

        /// Invalidate our forest if it is f.
        inline void invalidateForest(const forest* f) {
#ifndef USE_FID
            if (nodeF == f) {
                nodeF = nullptr;
            }
#endif
        }

        /// Notify forest that the given node is being
        /// added to a CT entry
        inline void cacheNode(node_handle n) const
        {
            ASSERT(__FILE__, __LINE__, ct_typeID::NODE == type);
#ifdef USE_FID
            forest* nodeF = forest::getForestWithID(nodeFID);
#endif
            ASSERT(__FILE__, __LINE__, nodeF);
            nodeF->cacheNode(n);
        }

        /// Notify forest that the given node is being
        /// removed as a CT entry
        inline void uncacheNode(node_handle n) const
        {
            ASSERT(__FILE__, __LINE__, ct_typeID::NODE == type);
#ifdef USE_FID
            forest* nodeF = forest::getForestWithID(nodeFID);
#endif
            if (nodeF) {
                nodeF->uncacheNode(n);
            }
        }

        /// Check if a CT entry is dead
        inline bool isDeadEntry(node_handle n) const
        {
            ASSERT(__FILE__, __LINE__, ct_typeID::NODE == type);
#ifdef USE_FID
            forest* nodeF = forest::getForestWithID(nodeFID);
#endif
            return nodeF ? nodeF->isDeadEntry(n) : true;
        }

        /// Check if a CT entry is stale.
        /// If not, and if mark is true, we'll mark it.
        inline bool isStaleEntry(node_handle n, bool mark) const
        {
            ASSERT(__FILE__, __LINE__, ct_typeID::NODE == type);
#ifdef USE_FID
            forest* nodeF = forest::getForestWithID(nodeFID);
#endif
            if (nodeF) {
                if (mark) {
                    if (nodeF->isStaleEntry(n)) return true;
                    nodeF->setCacheBit(n);
                    return false;
                } else {
                    return nodeF->isStaleEntry(n);
                }
            }
            return true;
        }

    public:
        /// Display; used for debugging
        void show(output &s) const;

    protected:
        void typeUpdate(forest* f);

    private:
        ct_typeID   type;
#ifdef USE_FID
        unsigned    nodeFID;
#else
        forest*     nodeF;
#endif

        bool        twoslots;
        bool        should_hash;
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
        ~ct_entry_type();

    public:
        // ***************************************************************
        //
        // Methods for building the entry type.
        //
        // ***************************************************************

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

                setFixed()      for setting the fixed portion of the key,
                                if any.
                setRepeat()     for setting the repeating portion of the key,
                                if any.
                setResult()     for setting the result portion of the key.

                doneBuilding()  call once when we're done building.
        */
        ct_entry_type(const char* name=nullptr);

        inline void setFixed(ct_itemtype a) {
            ASSERT(__FILE__, __LINE__, !CT);
            key_fixed.resize(1);
            key_fixed[0] = a;
        }
        inline void setFixed(ct_itemtype a, ct_itemtype b) {
            ASSERT(__FILE__, __LINE__, !CT);
            key_fixed.resize(2);
            key_fixed[0] = a;
            key_fixed[1] = b;
        }
        inline void setFixed(ct_itemtype a, ct_itemtype b, ct_itemtype c) {
            ASSERT(__FILE__, __LINE__, !CT);
            key_fixed.resize(3);
            key_fixed[0] = a;
            key_fixed[1] = b;
            key_fixed[2] = c;
        }
        inline void setFixed(ct_itemtype a, ct_itemtype b,
                ct_itemtype c, ct_itemtype d)
        {
            ASSERT(__FILE__, __LINE__, !CT);
            key_fixed.resize(4);
            key_fixed[0] = a;
            key_fixed[1] = b;
            key_fixed[2] = c;
            key_fixed[3] = d;
        }
        inline void appendFixed(ct_itemtype a) {
            ASSERT(__FILE__, __LINE__, !CT);
            key_fixed.push_back(a);
        }

        inline void setRepeat(ct_itemtype a) {
            ASSERT(__FILE__, __LINE__, !CT);
            key_repeating.resize(1);
            key_repeating[0] = a;
        }
        inline void setRepeat(ct_itemtype a, ct_itemtype b) {
            ASSERT(__FILE__, __LINE__, !CT);
            key_repeating.resize(2);
            key_repeating[0] = a;
            key_repeating[1] = b;
        }
        inline void appendRepeat(ct_itemtype a) {
            ASSERT(__FILE__, __LINE__, !CT);
            key_repeating.push_back(a);
        }

        inline void setResult(ct_itemtype a) {
            ASSERT(__FILE__, __LINE__, !CT);
            result.resize(1);
            result[0] = a;
        }
        inline void setResult(ct_itemtype a, ct_itemtype b) {
            ASSERT(__FILE__, __LINE__, !CT);
            result.resize(2);
            result[0] = a;
            result[1] = b;
        }
        inline void appendResult(ct_itemtype a) {
            ASSERT(__FILE__, __LINE__, !CT);
            result.push_back(a);
        }

        /// Call this (at least) once after setting the key and result items.
        inline void doneBuilding() {
            if (!CT) buildCT();
        }


    public:
        // ***************************************************************
        //
        // Convenience methods for use in operations
        //
        // ***************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6
        /**
            Find an entry in the compute table of this type.
                @param  key     Key to search for.
                @param  res     Where to store the result, if any.
         */
        inline void findCT(ct_entry_key* key, ct_entry_result &res) {
            ASSERT(__FILE__, __LINE__, CT);
            ASSERT(__FILE__, __LINE__, keyIsOurs(key));

            CT->find(key, res);
        }

        /**
            Add a compute table entry of this type.
                @param  key   Key portion of the entry.  Will be recycled.
                @param  res   Result portion of the entry.
         */
        inline void addCT(ct_entry_key* key, const ct_entry_result &res) {
            ASSERT(__FILE__, __LINE__, CT);
            ASSERT(__FILE__, __LINE__, keyIsOurs(key));

            CT->addEntry(key, res);
        }

        /**
            Done with a key.
            Call this after an unsuccessful findCT(),
            if we're not going to follow it with an addCT().
        */
        inline void noaddCT(ct_entry_key* key) {
            ASSERT(__FILE__, __LINE__, CT);
            ASSERT(__FILE__, __LINE__, keyIsOurs(key));

            CT->doneKey(key);
            CT->recycle(key);
        }

        inline compute_table* getCT() const {
            return CT;
        }

#endif
        //
        // New interface
        //

        /** Find an entry in the compute table of this type.
                @param  key   Key to search for.
                @param  res   Where to store the result, if found.
                @return true, if a result was found; false otherwise.
        */
        inline bool findCT(ct_vector &key, ct_vector &res)
        {
            ASSERT(__FILE__, __LINE__, CT);

            return CT->find(*this, key, res);
        }

        /**
            Add an entry (key plus result) to the compute table of this type.
                @param  key   Key portion of the entry.
                @param  res   Result portion of the entry.
        */
        inline void addCT(ct_vector &key, const ct_vector &res)
        {
            ASSERT(__FILE__, __LINE__, CT);

            CT->addEntry(*this, key, res);
        }

        /**
            Done with a key.
            Call this after an unsuccessful findCT(),
            if we're not going to follow it with an addCT().
        */
        inline void noaddCT(ct_vector &key) {
            ASSERT(__FILE__, __LINE__, CT);
            CT->doneKey(key);
        }

    public:
        // ***************************************************************
        //
        // Methods primarily for compute tables
        //
        // ***************************************************************

        /**
            Return true if the key or result of this entry
            contains nodes in forest f.
        */
        bool hasForest(const forest* f) const;

        /** Clear CT bits for any forests this entry type uses.
              @param  skipF   If skipF[i] is true, then we do nothing
                              for forests with ID i.  We set this to
                              true after clearing forest with ID i to
                              prevent clearing the bits twice.
        */
        void clearForestCTBits(std::vector <bool> &skipF) const;

        /** Calls clearForestCTBits() for all registered entries.
            @param  skipF   Passed to clearForestCTBits().
        */
        static void clearAllForestCTBits(std::vector <bool> &skipF);

        /** Notify forests that we're done marking CT bits.
            The forests can choose to start the sweep phase if they like.
              @param  whichF  If whichF[i] is true, then we notify the
                              forest with ID i, and set whichF[i] to false.
                              This prevents notifying a forest twice.
        */
        void sweepForestCTBits(std::vector <bool> &whichF) const;

        /** Calls sweepForestCTBits() for all registered entries.
            @param  whichF   Passed to sweepForestCTBits().
        */
        static void sweepAllForestCTBits(std::vector <bool> &whichF);

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
            Set the name; allows empty constructors.
        */
        inline void setName(const char* n) { name = n; }

        /**
            Does this entry type allow repetitions in the key?
            I.e., was there a '.' in the pattern?
        */
        inline bool isRepeating() const { return key_repeating.size(); }

        /**
            Can we have a key of the given size?
        */
        inline bool canHaveKeySize(unsigned sz) const {
            if (sz < key_fixed.size()) return false;
            const unsigned rsize = sz - key_fixed.size();
            if (rsize) {
                // there's some leftover key after the fixed portion
                // make sure it can be divided into the repeating portion.
                if (key_repeating.size()) {
                    return 0 == (rsize % key_repeating.size());
                } else {
                    // no repeating portion
                    return false;
                }
            }
            // no leftover key after the fixed portion
            return true;
        }

        /**
            Get number of repetitions for a key of the given size.
            Assumes canHaveKeySize() is true.
        */
        inline unsigned repeatsForKeySize(unsigned sz) const {
            if (sz <= key_fixed.size()) return 0;
            const unsigned rsize = sz - key_fixed.size();
            return rsize / key_repeating.size();
        }

        /**
            Get the number of items in the key.
              @param  reps  Number of repetitions.
                            If this is not a repeating type,
                            then this is ignored.

              @return Total number of slots in the key.
        */
        inline unsigned getKeySize(unsigned reps=0) const {
            return key_fixed.size() + (reps * key_repeating.size());
        }

        /**
            Get the number of integer slots required for the key.
              @param  reps  Number of repetitions.
                            If this is not a repeating type,
                            then this is ignored.

              @return Total number of integer slots required for the key.
        */
        inline unsigned getKeyIntslots(unsigned reps) const {
            return fixed_intslots + (reps * repeating_intslots);
        }

        /**
            Is the entire key hashable?
                @param  reps    Number of repetitions.
                                If this is not a repeating type,
                                then this is ignored.
                                (In practice, it only matters if
                                there are zero repetitions,
                                or more than that.)

                @return true, if every key item is hashable
                              (method shouldBeHashed() returns true).
        */
        inline bool isEntireKeyHashable(unsigned reps) const {
            return fixed_entirely_hashable
                    &&
                    ( (0==reps) || repeating_entirely_hashable );
        }

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
                ASSERT(__FILE__, __LINE__, key_repeating.size());
                i -= key_fixed.size();
                i %= key_repeating.size();
                return key_repeating[i];
            }
        }

        /**
            Get the number of items in the result
        */
        inline unsigned getResultSize() const { return result.size(); }

        /**
            Get the number of integer slots in the result
        */
        inline unsigned getResultIntslots() const { return result_intslots; }

        /**
            Get the type for item i in the result.
              @param  i   Slot number, between 0 and getResultSize().
        */
        inline const ct_itemtype& getResultType(unsigned i) const
        {
            CHECK_RANGE(__FILE__, __LINE__,
                    0u, i, unsigned(result.size()));
            return result[i];
        }

        /// Mark for deletion
        inline void markForDeletion() { is_marked_for_deletion = true; }

        /// Mark for destroy
        inline void markForDestroy() {
            markForDeletion();
            destroyWhenEmpty = true;
#ifdef DELETE_ON_ZERO
            if (0==numEntries) {
                delete this;
            }
#endif
        }

        /// Unmark for deletion
        inline void unmarkForDeletion() { is_marked_for_deletion = false; }


        /// Should we remove all CT entries of this type?
        inline bool isMarkedForDeletion() const {
            return is_marked_for_deletion;
        }

        /// Display info, for debugging
        void show(output &s) const;

        /// Update all of our items that forest f has been deleted
        void invalidateForest(const forest* f);

        /// Get current number of entries
        inline unsigned long getNumEntries() const {
            return numEntries;
        }

    public:
        // ***************************************************************
        //
        // Registry of all entries: interface
        //
        // ***************************************************************

        static void initStatics();
        static void doneStatics();

        static inline const ct_entry_type* getEntryType(unsigned etid)
        {
#ifdef DEVELOPMENT_CODE
            return all_entries.at(etid);
#else
            return all_entries[etid];
#endif
        }

        /**
            Notify all items in all entries that forest f is invalid.
            Should be called automatically when a forest is destroyed.
         */
        static void invalidateAllWithForest(const forest* f);


        /**
            Immediately remove ALL entries involving forest f.
            Required for things like variable reordering.
        */
        static void removeAllCTEntriesWithForest(const forest* f);

        /// Display info for all registered entries
        static void showAll(output &s);

        /// CT is adding another entry of type etid
        static inline void incEntries(unsigned etid) {
            ASSERT(__FILE__, __LINE__, all_entries.at(etid));
            all_entries[etid]->numEntries++;
        }

        /// CT is removing an entry of this type
        static inline void decEntries(unsigned etid) {
            ASSERT(__FILE__, __LINE__, all_entries.at(etid));
            ASSERT(__FILE__, __LINE__, all_entries[etid]->numEntries);
            all_entries[etid]->numEntries--;
#ifdef DELETE_ON_ZERO
            if (all_entries[etid]->destroyWhenEmpty) {
                if (0==all_entries[etid]->numEntries) {
                    delete all_entries[etid];
                    ASSERT(__FILE__, __LINE__, nullptr == all_entries[etid]);
                }
            }
#endif
        }

        /**
            Check if ANY registered entry type uses a ct_object.
            If so, we must clean the monolithic compute table.
         */
        static bool mightHaveCTObjects();


    private:
        // ***************************************************************
        //
        // Registry of all entries: back end
        //
        // ***************************************************************


        /// Unique ID, for the life of the library.
        /// Guaranteed not to be 0.
        unsigned etID;

        /// Global registry of all entries
        static std::vector <ct_entry_type*> all_entries;

        static inline void registerEntry(ct_entry_type* et) {
            if (et) {
                et->etID = all_entries.size();
                all_entries.push_back(et);
#ifdef DEBUG_CLEANUP
                std::cout << "registered ct entry " << et << ", #" << et->etID << "\n";
#endif
            }
        }

        static inline void unregisterEntry(ct_entry_type* et) {
            if (et) {
#ifdef DEVELOPMENT_CODE
                all_entries.at(et->etID) = nullptr;
#else
                all_entries[et->etID] = nullptr;
#endif
#ifdef DEBUG_CLEANUP
                std::cout << "unregistered ct entry " << et << ", #" << et->etID << "\n";
#endif
            }
        }


    private:
        // ***************************************************************
        //
        // Helper methods and ugly details
        //
        // ***************************************************************

#ifdef ALLOW_DEPRECATED_0_17_6
        /// Make sure k has us as entry type.
        bool keyIsOurs(const ct_entry_key *k) const;
#endif

        /// Return true if the entry has some valid forest
        bool hasValidForest() const;

        /// Return true if an entry might have a ct_object
        bool hasCTObject() const;

        /// Set counts; build CT.
        void buildCT();

        /// Compute table to use for this entry
        compute_table* CT;

        /// Do we 'own' the compute table?
        /// I.e., is this an operation-specific (not monolithic) CT?
        bool CT_is_ours;

        /// Name; for displaying CT entries
        const char* name;

        /// How many entries of this type exist in some CT?
        unsigned long numEntries;

        /// Fixed, initial portion of the key.
        std::vector <ct_itemtype> key_fixed;

        /// Number of integer slots required for the fixed key
        unsigned fixed_intslots;

        /// Are all items in the fixed key, hashable?
        bool fixed_entirely_hashable;

        /// Repeating portion of the key.
        std::vector <ct_itemtype> key_repeating;

        /// Number of integer slots required for the repeating key
        unsigned repeating_intslots;

        /// Are all items in the repeating key, hashable?
        bool repeating_entirely_hashable;

        /// Result pattern
        std::vector <ct_itemtype> result;

        /// Number of integer slots required for the result
        unsigned result_intslots;

        /// Can the result be changed later
        bool updatable_result;

        /// For deleting all entries of this type
        bool is_marked_for_deletion;

        /// Delete this as soon as there are no entries
        bool destroyWhenEmpty;
};

#endif
