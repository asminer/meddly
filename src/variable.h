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

#ifndef MEDDLY_VARIABLE_H
#define MEDDLY_VARIABLE_H

namespace MEDDLY {
    class variable;
    class expert_variable;
    class domain;
};

//
// TBD: merge variable and expert_variable classes
//

// ******************************************************************
// *                                                                *
// *                         variable class                         *
// *                                                                *
// ******************************************************************

/** Variable class.
    Abstract base class.
    A variable consists of an optional name, and a bound.
    A single variable object is used to describe both
    the primed and unprimed versions of the variable.

    Note1: variables are automatically deleted when
    removed from their last domain.

    Additional features are provided in the expert interface.
*/
class MEDDLY::variable {
  protected:
    variable(int bound, char* name);
    virtual ~variable();
  public:
    int getBound(bool primed) const;
    const char* getName() const;
    void setName(char* newname);
    bool isExtensible() const;
  protected:
    int un_bound;
    int pr_bound;
    bool is_extensible;
  private:
    char* name;
};

// ******************************************************************
// *                    inlined variable methods                    *
// ******************************************************************

inline int MEDDLY::variable::getBound(bool primed) const {
  return primed ? pr_bound : un_bound;
}
inline const char* MEDDLY::variable::getName() const { return name; }
inline bool MEDDLY::variable::isExtensible() const {
  return is_extensible;
}

// ******************************************************************
// *                                                                *
// *                     expert_variable  class                     *
// *                                                                *
// ******************************************************************

class MEDDLY::expert_variable : public variable {
  public:
    expert_variable(int b, char* n);

    /// Update our list of domains: add \a d.
    void addToList(domain* d);
    /// Update our list of domains: remove \a d.
    void removeFromList(const domain* d);

    /** Enlarge the possible values for a variable.
      This could modify all nodes in all forests, depending on the
      choice of reduction rule.
      @param  prime   If prime is true, enlarge the bound for
                      the primed variable only, otherwise both
                      the primed and unprimed are enlarged.
      @param  b       New bound, if less than the current bound
                      an error is thrown.
                      If bound<=0, the variable is marked as extensible,
                      with initial bound as abs(bound).
                      Note: an extensible variable has a range [1 .. +infinity].
    */
    void enlargeBound(bool prime, int b);

    /** Shrink the possible values for a variable.
      This could modify all nodes in all forests, depending on the
      choice of reduction rule.
      @param  b       New bound, if more than the current bound
                      an error is thrown.
      @param  force   If \a b is too small, and information will be lost,
                      proceed anyway if \a force is true, otherwise
                      return an error code.
    */
    void shrinkBound(int b, bool force);

  private:
    domain** domlist;
    int dl_alloc;
    int dl_used;

    virtual ~expert_variable();
};

// ******************************************************************
// *                inlined expert_variable  methods                *
// ******************************************************************

#endif
