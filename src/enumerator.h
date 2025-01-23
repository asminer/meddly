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

#ifndef MEDDLY_ENUMERATOR_H
#define MEDDLY_ENUMERATOR_H

#ifdef ALLOW_DEPRECATED_0_17_7

namespace MEDDLY {
    class dd_edge;
    class enumerator;
    class forest;
    class unpacked_node;
};

// ******************************************************************
// *                                                                *
// *                                                                *
// *                        enumerator class                        *
// *                                                                *
// *                                                                *
// ******************************************************************

/** Class for enumerating values encoded by a dd-edge.
    Effectively, a single class encapsulating
    several possible iterators.
*/
class MEDDLY::enumerator {
  public:
    class iterator {
      public:
        iterator(const forest* F);
        virtual ~iterator();

        inline bool build_error() const {
            return 0==F;
        }

        virtual bool start(const dd_edge &e);
        virtual bool start(const dd_edge &e, const int* m);
        virtual bool next() = 0;

        /**
            Return the highest level changed during the last increment.
        */
        inline int levelChanged() const {
            return level_change;
        }


        /** Get the current variable assignments.
            For variable i, use index i for the
            unprimed variable, and index -i for the primed variable.
        */
        inline const int* getAssignments() const {
            return index;
        }

        /** Get primed assignments.
            It is much faster to use getAssigments()
            and look at the negative indexes;
            however, this works.
        */
        const int* getPrimedAssignments();

        /// For integer-ranged edges, get the current non-zero value.
        virtual void getValue(int& edgeValue) const;

        /// For integer-ranged edges, get the current non-zero value.
        virtual void getValue(long& edgeValue) const;

        /// For real-ranged edges, get the current non-zero value.
        virtual void getValue(float& edgeValue) const;



      protected:
        // Current parent forest.
        const forest* F;
        // Path, as list of unpacked nodes
        unpacked_node*    rawpath;
        unpacked_node*    path;   // rawpath, shifted so we can use path[-k]
        // Path nnz pointers
        int*      rawnzp;
        int*      nzp;   // rawnzp, shifted so we can use nzp[-k]
        // Path indexes
        int*      rawindex;
        int*      index;  // rawindex, shifted so we can use index[-k]
        //
        int       minLevel; // 1 or -#vars, depending.
        int       maxLevel; // #vars
        //
        int       level_change;

      private:
        // Used only by getPrimedAssignments.
        int*      prindex;
    };

  public:
    /// Enumerator type.
    enum type {
      EMPTY,
      FULL,
      ROW_FIXED,
      COL_FIXED
    };

  public:
    /// Empty constructor.
    enumerator();

    /// Proper constructor.
    enumerator(type t, const forest* F);

    /** What is usually wanted constructor.
        Equivalent to enumerator(FULL, e.getForest())
        followed by start(e).
    */
    enumerator(const dd_edge &e);

    /// Destructor.
    ~enumerator();

    /// Re-initialize
    void init(type t, const forest* F);

    inline operator bool() const {
        return is_valid;
    }

    inline void operator++() {
        is_valid &= I->next();
    }



    /**
        Start iterating through edge e.
    */
    void start(const dd_edge &e);

    /** Start iterating through edge e.
        The unprimed variables will be fixed to
        the given values.
          @param  e         Edge to iterate.
                            Must be a relation.
          @param  minterm   Array of dimension 1+vars in e.
                            minterm[k] gives the fixed variable
                            assignment for (unprimed) variable k.
    */
    void startFixedRow(const dd_edge &e, const int* minterm);

    /** Start iterating through edge e.
        The primed variables will be fixed to
        the given values.
          @param  e         Edge to iterate.
                            Must be a relation.
          @param  minterm   Array of dimension 1+vars in e.
                            minterm[k] gives the fixed variable
                            assignment for (unprimed) variable k.
    */
    void startFixedColumn(const dd_edge &e, const int* minterm);


    /** Get the current variable assignments.
        For variable i, use index i for the
        unprimed variable, and index -i for the primed variable.
    */
    inline const int* getAssignments() const {
        if (I && is_valid) return I->getAssignments(); else return 0;
    }

    /// Get the current primed variable assignments.
    inline const int* getPrimedAssignments() const {
        if (I && is_valid) return I->getPrimedAssignments(); else return 0;
    }

    inline void getValue(int &v) const {
        if (I && is_valid) I->getValue(v);
    }
    inline void getValue(long &v) const {
        if (I && is_valid) I->getValue(v);
    }
    inline void getValue(float &v) const {
        if (I && is_valid) I->getValue(v);
    }

    inline int levelChanged() const {
        if (I) return I->levelChanged();
        return 0;
    }

    inline type getType() const {
        return T;
    }

  private:
    iterator* I;
    bool is_valid;
    type T;
};

#endif

#endif
